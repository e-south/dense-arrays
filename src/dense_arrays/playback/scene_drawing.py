"""Compose graph, duplex, legends, and visible evidence within a figure.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING

from .duplex_drawing import draw_duplex
from .duplex_frames import DuplexFrames, duplex_frame_for_axis
from .graph.geometry import (
    KMER_FONT_FAMILY,
    KMER_FONT_SIZE_PT,
    KMER_FONT_WEIGHT,
    matplotlib_layout_spec,
)
from .graph.layout import build_graph_scene
from .graph.projection import project_explanation_graph
from .graph.routing import quadratic_arc_length, route_graph_scene
from .graph_drawing import draw_graph
from .presentation import PlaybackDocument, resolve_distance_brackets, resolve_evidence
from .theme import RESTING_COLOR, RESTING_TEXT_COLOR, blend_color
from .timeline import placement_progress

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

_PAPER = "#ffffff"
_GRAPH_TEXT = "#1F2423"
_EVIDENCE_LINE_WIDTH = 155
_MAX_EVIDENCE_FRACTION = 0.25
_MAX_NATIVE_DISTANCE_BRACKETS = 2


def _graph_font_size_for_cap_height(
    cap_height_px: float | None, figure_dpi: float
) -> float:
    if cap_height_px is None:
        return KMER_FONT_SIZE_PT
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    properties = FontProperties(
        family=KMER_FONT_FAMILY,
        size=1.0,
        weight=KMER_FONT_WEIGHT,
    )
    unit_cap_height_pt = float(
        TextPath((0, 0), "ACGT", prop=properties).get_extents().height
    )
    if unit_cap_height_pt <= 0:
        return KMER_FONT_SIZE_PT
    return cap_height_px * 72.0 / (float(figure_dpi) * unit_cap_height_pt)


def draw_document(
    document: PlaybackDocument,
    *,
    transition_index: int,
    progress: float,
    figure: Figure,
    duplex_frames: DuplexFrames | None = None,
) -> None:
    """Compose a single scene state with its visible evidence and caller labels."""
    if not document.plan.steps:
        msg = "playback document requires at least one step"
        raise ValueError(msg)
    transition_index = max(0, min(transition_index, len(document.plan.steps)))
    progress = max(0.0, min(progress, 1.0))
    figure.clear()
    figure.set_facecolor(_PAPER)
    graph_axis, duplex_axis, legend_axis = document_axes(figure, document)
    duplex_frame = None
    displayed_duplex_cap_height = None
    if duplex_frames is not None:
        duplex_frame, displayed_duplex_cap_height = duplex_frame_for_axis(
            duplex_frames,
            document,
            transition_index,
            progress,
            duplex_axis,
        )
    if graph_axis is not None:
        draw_graph(
            graph_axis,
            document,
            transition_index,
            progress,
            kmer_font_size_pt=_graph_font_size_for_cap_height(
                displayed_duplex_cap_height,
                figure.dpi,
            ),
        )
    if duplex_frames is None:
        draw_duplex(duplex_axis, document, transition_index, progress)
    else:
        duplex_axis.imshow(
            duplex_frame,
            interpolation="lanczos",
            resample=True,
        )
        duplex_axis.axis("off")
    if legend_axis is not None:
        _draw_legend(legend_axis, document, transition_index, progress)
    if duplex_frames is None or not duplex_frames.renders_distance_brackets:
        draw_distance_brackets(duplex_axis, document, transition_index, progress)
    lines = evidence_lines(document)
    line_height = min(
        11 / (72 * figure.get_figheight()),
        (_MAX_EVIDENCE_FRACTION - 0.05) / len(lines),
    )
    for index, (line, color) in enumerate(reversed(lines)):
        figure.text(
            0.02,
            0.025 + index * line_height,
            line,
            color=color,
            ha="left",
            va="bottom",
            fontsize=min(8, line_height * figure.get_figheight() * 72 * 8 / 11),
        )


def evidence_lines(document: PlaybackDocument) -> tuple[tuple[str, str], ...]:
    """Bound canvas evidence to three lines, disclosing metadata for full details."""
    evidence = resolve_evidence(document)
    lines = [(evidence.qualification, "#59635f")]
    if evidence.constraints:
        text = "; ".join(evidence.constraints)
        if len(text) > _EVIDENCE_LINE_WIDTH:
            text = (
                f"FAILED {len(evidence.constraints)} distance constraints; "
                "full actual/required results in metadata"
            )
        lines.append((text, "#9d2525"))
    if evidence.notices:
        text = "; ".join(evidence.notices)
        if len(text) > _EVIDENCE_LINE_WIDTH:
            excerpt = textwrap.shorten(text, width=100, placeholder="…")
            text = f"{excerpt} ({len(evidence.notices)} notices; full text in metadata)"
        lines.append((text, "#59635f"))
    return tuple(lines)


def draw_distance_brackets(
    axis: Axes, document: PlaybackDocument, transition_index: int, progress: float
) -> None:
    """Draw declared distances below a native or producer-rendered duplex."""
    brackets = resolve_distance_brackets(document)
    if not brackets:
        if document.presentation.show_distance_bracket == "always":
            axis.text(
                0.5,
                0.015,
                "No declared distance constraints",
                transform=axis.transAxes,
                ha="center",
                fontsize=8,
                color=blend_color(
                    RESTING_TEXT_COLOR,
                    "#59635f",
                    placement_progress(0, transition_index, progress),
                ),
            )
        return
    if len(brackets) > _MAX_NATIVE_DISTANCE_BRACKETS or any(
        len(bracket.label) > _EVIDENCE_LINE_WIDTH for bracket in brackets
    ):
        axis.text(
            0.5,
            0.04,
            f"{len(brackets)} declared distances; full results in metadata",
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontsize=7,
            color=blend_color(
                RESTING_TEXT_COLOR,
                "#9d2525" if any(bracket.failed for bracket in brackets) else "#59635f",
                placement_progress(
                    len(document.plan.steps) - 1, transition_index, progress
                ),
            ),
            clip_on=True,
        )
        return
    length = len(document.plan.realized_sequence)
    for index, bracket in enumerate(brackets):
        y = 0.04 + index * 0.10
        x1, x2 = (
            0.08 + bracket.start / length * 0.84,
            0.08 + bracket.end / length * 0.84,
        )
        result = document.plan.constraint_results[index]
        emphasis = min(
            placement_progress(step_index, transition_index, progress)
            for step_index, step in enumerate(document.plan.steps)
            if step.placement_id
            in (result.upstream_placement_id, result.downstream_placement_id)
        )
        color = blend_color(
            RESTING_COLOR, "#9d2525" if bracket.failed else "#59635f", emphasis
        )
        axis.plot(
            (x1, x1, x2, x2),
            (y + 0.025, y, y, y + 0.025),
            transform=axis.transAxes,
            color=color,
            linewidth=1,
            gid="distance-bracket",
        )
        axis.text(
            (x1 + x2) / 2,
            y + 0.03,
            bracket.label,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontsize=7,
            color=color,
        )


def document_axes(
    figure: Figure, document: PlaybackDocument
) -> tuple[Axes | None, Axes, Axes | None]:
    """Allocate axes while reserving physical space for scene evidence."""
    header_height = 0.015
    footer_height = min(
        _MAX_EVIDENCE_FRACTION,
        (0.06 + len(evidence_lines(document)) * 11 / 72) / figure.get_figheight(),
    )
    has_legend = bool(document.presentation.legend_entries)
    graph_fraction = float(document.presentation.graph_fraction)
    show_graph = document.presentation.graph_detail != "none"
    if not show_graph and has_legend:
        grid = figure.add_gridspec(
            2,
            1,
            height_ratios=(0.82, 0.18),
            left=0.004,
            right=0.997,
            top=1.0 - header_height,
            bottom=footer_height,
            hspace=0.0,
        )
        return None, figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[1, 0])
    if not show_graph:
        grid = figure.add_gridspec(
            1, 1, left=0.02, right=0.98, top=1.0 - header_height, bottom=footer_height
        )
        return None, figure.add_subplot(grid[0, 0]), None
    if has_legend:
        grid = figure.add_gridspec(
            2,
            2,
            width_ratios=(graph_fraction, 1.0 - graph_fraction),
            height_ratios=(0.82, 0.18),
            left=0.004,
            right=0.997,
            top=1.0 - header_height,
            bottom=footer_height,
            wspace=0.001,
            hspace=0.0,
        )
        graph_axis = figure.add_subplot(grid[:, 0])
        duplex_axis = figure.add_subplot(grid[0, 1])
        legend_axis = figure.add_subplot(grid[1, 1])
        return graph_axis, duplex_axis, legend_axis
    grid = figure.add_gridspec(
        1,
        2,
        width_ratios=(graph_fraction, 1.0 - graph_fraction),
        left=0.004,
        right=0.997,
        top=1.0 - header_height,
        bottom=footer_height,
        wspace=0.001,
    )
    graph_axis = figure.add_subplot(grid[0, 0])
    duplex_axis = figure.add_subplot(grid[0, 1])
    legend_axis = None
    return graph_axis, duplex_axis, legend_axis


def _draw_legend(
    axis: Axes, document: PlaybackDocument, transition_index: int, progress: float
) -> None:
    entries = document.presentation.legend_entries
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.0, 1.0)
    axis.axis("off")
    if not entries:
        return
    group_span = min(0.70, 0.175 * len(entries))
    segment = group_span / len(entries)
    group_start = (1.0 - group_span) / 2.0
    for index, entry in enumerate(entries):
        emphasis = placement_progress(0, transition_index, progress)
        center = group_start + segment * (index + 0.5)
        axis.scatter(
            (center - 0.070,),
            (0.5,),
            transform=axis.transAxes,
            marker="s",
            s=112,
            facecolor=blend_color(RESTING_COLOR, entry.color, emphasis),
            edgecolor="none",
            linewidth=0.0,
        )
        axis.text(
            center - 0.047,
            0.5,
            entry.label,
            transform=axis.transAxes,
            ha="left",
            va="center",
            color=blend_color(RESTING_TEXT_COLOR, _GRAPH_TEXT, emphasis),
            fontsize=11.5,
            family=KMER_FONT_FAMILY,
            fontweight="normal",
        )


def transition_frame_counts(
    document: PlaybackDocument,
    figure: Figure,
    *,
    fps: int,
    seconds_per_step: float,
) -> tuple[int, ...]:
    """Allocate transition frames by visible traversal arc length."""
    figure.clear()
    graph_axis, _duplex_axis, _legend_axis = document_axes(figure, document)
    if graph_axis is None or not resolve_evidence(document).animate_chain:
        count = max(1, round(fps * seconds_per_step))
        figure.clear()
        return (count,) * (len(document.plan.steps) + 1)
    semantic_graph = project_explanation_graph(document.plan)
    layout_spec = matplotlib_layout_spec(graph_axis, semantic_graph)
    scene = build_graph_scene(document.plan, layout_spec=layout_spec)
    routes = route_graph_scene(scene)
    lengths = tuple(
        max(quadratic_arc_length(routed.curve), 1.0) for routed in routes.traversal
    )
    expected = len(document.plan.steps) + 1
    if len(lengths) != expected:
        msg = f"playback requires {expected} traversal edges, found {len(lengths)}"
        raise ValueError(msg)
    mean_length = sum(lengths) / len(lengths)
    counts = tuple(
        max(1, round(fps * seconds_per_step * length / mean_length))
        for length in lengths
    )
    figure.clear()
    return counts
