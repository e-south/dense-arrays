"""Short-wide publication renderers for synchronized path and duplex playback."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import Any

from .graph.routing import (
    EDGE_LABEL_FONT_SIZE_PT,
    quadratic_arc_length,
    quadratic_arc_t,
    quadratic_segment,
)
from .graph_layout import (
    END_NODE_ID,
    KMER_FONT_FAMILY,
    KMER_FONT_SIZE_PT,
    KMER_FONT_WEIGHT,
    START_NODE_ID,
    build_graph_scene,
    matplotlib_layout_spec,
    project_explanation_graph,
    quadratic_point,
    route_graph_scene,
    step_color,
)
from .html import PlaybackDocument
from .theme import (
    constraint_relation_color,
    step_text_color,
)

DuplexFrameRenderer = Callable[[PlaybackDocument, int], Any]
_PAPER = "#ffffff"
_DUPLEX_RASTER_OVERSAMPLE = 2.0
_INK = "#4b5563"
_LINE = "#9cc9c1"
_TRAVERSED = "#50635f"
_ACTIVE = "#167a70"
_GRAPH_TEXT = "#1F2423"
_TERMINAL_FONT_SIZE_PT = 12.4
_COMPLEMENTS = str.maketrans("ATCGRYSWKMBDHVN", "TAGCYRSWMKVHDBN")


def _require_documents(documents: tuple[PlaybackDocument, ...]) -> None:
    if not documents:
        raise ValueError("at least one playback document is required")


def _smoothstep(progress: float) -> float:
    progress = max(0.0, min(1.0, progress))
    return progress * progress * (3.0 - (2.0 * progress))


def _matplotlib_path(curve):
    from matplotlib.path import Path

    return Path(
        (
            curve.visible_start,
            curve.visible_control or curve.control,
            curve.visible_end,
        ),
        (Path.MOVETO, Path.CURVE3, Path.CURVE3),
    )


def _matplotlib_segment_path(points):
    from matplotlib.path import Path

    return Path(points, (Path.MOVETO, Path.CURVE3, Path.CURVE3))


def _draw_graph(
    axis,
    document: PlaybackDocument,
    transition_index: int,
    progress: float,
    *,
    kmer_font_size_pt: float = KMER_FONT_SIZE_PT,
) -> None:
    from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

    steps = document.plan.steps
    semantic_graph = project_explanation_graph(document.plan)
    layout_spec = matplotlib_layout_spec(
        axis,
        semantic_graph,
        kmer_font_size_pt=kmer_font_size_pt,
    )
    scene = build_graph_scene(document.plan, layout_spec=layout_spec)
    routes = route_graph_scene(scene)
    axis.set_xlim(0, layout_spec.viewport.width_pt)
    axis.set_ylim(0, layout_spec.viewport.height_pt)
    axis.set_aspect("equal", adjustable="box")
    axis.axis("off")

    visible_context = tuple(
        routed
        for routed in routes.context
        if document.presentation.graph_detail == "full"
        and routed.edge.relation_kind != "declared_constraint"
    )
    for routed in visible_context:
        curve = routed.curve
        declared_constraint = routed.edge.relation_kind == "declared_constraint"
        axis.add_patch(
            FancyArrowPatch(
                path=_matplotlib_path(curve),
                arrowstyle="-|>",
                mutation_scale=5.4,
                linewidth=1.7 if declared_constraint else 0.85,
                color=(
                    constraint_relation_color(document.presentation.color_profile)
                    if declared_constraint
                    else _LINE
                ),
                alpha=0.82 if declared_constraint else 0.34,
                zorder=0,
            )
        )

    active_geometry = None
    for routed in routes.traversal:
        edge = routed.edge
        edge_index = edge.traversal_index
        if edge_index is None:
            raise ValueError("traversal edge is missing its timeline index")
        completed = edge_index < transition_index or (
            edge_index == transition_index and progress >= 0.999
        )
        active = edge_index == transition_index and progress < 0.999
        color = _TRAVERSED if completed else _LINE
        width = 2.4 if completed else 1.45
        curve = routed.curve
        axis.add_patch(
            FancyArrowPatch(
                path=_matplotlib_path(curve),
                arrowstyle="-|>",
                mutation_scale=7.8,
                linewidth=width,
                color=color,
                alpha=1.0 if completed or active else 0.72,
                zorder=1,
            )
        )
        active_t = None
        if active:
            active_t = quadratic_arc_t(curve, progress)
            if active_t > curve.visible_t_start:
                partial_end = min(active_t, curve.visible_t_end)
                if partial_end > curve.visible_t_start:
                    axis.add_patch(
                        FancyArrowPatch(
                            path=_matplotlib_segment_path(
                                quadratic_segment(
                                    curve,
                                    curve.visible_t_start,
                                    partial_end,
                                )
                            ),
                            arrowstyle="-",
                            linewidth=2.4,
                            color=_ACTIVE,
                            alpha=1.0,
                            zorder=1.2,
                        )
                    )
        if (
            document.presentation.show_edge_costs
            and edge.added_bases is not None
            and routed.label_position is not None
        ):
            axis.text(
                routed.label_position[0],
                routed.label_position[1],
                str(edge.added_bases),
                ha="center",
                va="center",
                color=_GRAPH_TEXT,
                fontsize=routed.label_font_size or EDGE_LABEL_FONT_SIZE_PT,
                family=KMER_FONT_FAMILY,
                bbox={
                    "boxstyle": "round,pad=0.10",
                    "facecolor": _PAPER,
                    "edgecolor": _LINE,
                    "linewidth": 0.7,
                    "alpha": 0.97,
                },
                zorder=2.0,
            )
        if active:
            active_geometry = (curve, active_t)

    for index, step in enumerate(steps):
        geometry = scene.geometry(step.placement_id)
        x, y = scene.position(step.placement_id)
        color = step_color(step, index, document.presentation.color_profile)
        axis.add_patch(
            FancyBboxPatch(
                (x - geometry.width_pt / 2, y - geometry.height_pt / 2),
                geometry.width_pt,
                geometry.height_pt,
                boxstyle="round,pad=0.0,rounding_size=1.5",
                facecolor=color,
                edgecolor="none",
                linewidth=0.0,
                alpha=1.0,
                zorder=3,
            )
        )
        axis.text(
            x,
            y,
            step.placement_sequence,
            ha="center",
            va="center",
            color=step_text_color(
                step,
                index,
                document.presentation.color_profile,
            ),
            fontsize=kmer_font_size_pt,
            family=KMER_FONT_FAMILY,
            fontweight=KMER_FONT_WEIGHT,
            zorder=4,
        )

    for label, node_id in (("Start", START_NODE_ID), ("End", END_NODE_ID)):
        x, y = scene.position(node_id)
        geometry = scene.geometry(node_id)
        radius = geometry.width_pt / 2.0
        axis.add_patch(
            Circle(
                (x, y),
                radius=radius,
                facecolor=_PAPER,
                edgecolor=_INK,
                linewidth=1.25,
                zorder=4,
            )
        )
        axis.text(
            x,
            y + radius + 3.5,
            label,
            ha="center",
            va="bottom",
            color=_GRAPH_TEXT,
            fontsize=_TERMINAL_FONT_SIZE_PT,
            family=KMER_FONT_FAMILY,
            zorder=4,
        )

    ordering_status = getattr(
        document.plan.ordering_status, "value", str(document.plan.ordering_status)
    )
    if active_geometry is not None and ordering_status != "layout_only":
        curve, active_t = active_geometry
        point = quadratic_point(
            curve.motion_start,
            curve.control,
            curve.motion_end,
            active_t,
        )
        axis.scatter(
            (point[0],),
            (point[1],),
            s=26,
            facecolor=_ACTIVE,
            edgecolor="none",
            linewidth=0.0,
            zorder=1.5,
        )
    if document.presentation.show_authority_notice:
        axis.text(
            0.5,
            0.015,
            "Realized order",
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            color="#6B7280",
            fontsize=10.5,
            family=KMER_FONT_FAMILY,
            fontweight="normal",
            zorder=6,
        )


def _draw_fallback_duplex(axis, document: PlaybackDocument, step_index: int) -> None:
    from matplotlib.patches import FancyBboxPatch

    plan = document.plan
    sequence = plan.realized_sequence
    current = plan.steps[step_index]
    reveal_end = max(step.end for step in plan.steps[: step_index + 1])
    complement = sequence.translate(_COMPLEMENTS)
    length = len(sequence)
    axis.set_xlim(-4, length + 3)
    axis.set_ylim(-2.2, 2.2)
    axis.axis("off")
    for index, step in enumerate(plan.steps[: step_index + 1]):
        y = (
            1.05 + (index % 2) * 0.48
            if step.orientation != "rev"
            else -1.30 - (index % 2) * 0.48
        )
        color = step_color(step, index, document.presentation.color_profile)
        axis.add_patch(
            FancyBboxPatch(
                (step.start, y),
                step.end - step.start,
                0.38,
                boxstyle="round,pad=0.01,rounding_size=0.08",
                facecolor=color,
                edgecolor=_ACTIVE if index == step_index else color,
                linewidth=2.6 if index == step_index else 1.2,
            )
        )
        axis.text(
            (step.start + step.end) / 2,
            y + 0.19,
            step.placement_sequence,
            ha="center",
            va="center",
            color="white",
            fontsize=5.8,
            family=KMER_FONT_FAMILY,
        )
    font_size = max(5.2, min(10.5, 830 / max(1, length)))
    for index in range(reveal_end):
        color = _ACTIVE if current.start <= index < current.end else _INK
        axis.text(
            index + 0.5,
            0.30,
            sequence[index],
            ha="center",
            va="center",
            color=color,
            fontsize=font_size,
            family=KMER_FONT_FAMILY,
        )
        axis.text(
            index + 0.5,
            -0.30,
            complement[index],
            ha="center",
            va="center",
            color=color,
            fontsize=font_size,
            family=KMER_FONT_FAMILY,
        )
    axis.text(-2.0, 0.30, "5'", ha="center", va="center", color=_INK, fontsize=11)
    axis.text(-2.0, -0.30, "3'", ha="center", va="center", color=_INK, fontsize=11)
    axis.text(length + 0.7, 0.30, "3'", ha="left", va="center", color=_INK, fontsize=11)
    axis.text(
        length + 0.7, -0.30, "5'", ha="left", va="center", color=_INK, fontsize=11
    )


def _duplex_frame_for_axis(
    renderer: DuplexFrameRenderer,
    document: PlaybackDocument,
    transition_index: int,
    progress: float,
    axis,
):
    """Prepare a duplex frame at its final display resolution."""
    import numpy as np
    from PIL import Image

    frame = _duplex_transition_frame(
        renderer,
        document,
        transition_index,
        progress,
    )
    source_height, source_width = frame.shape[:2]
    axis_bounds = axis.get_window_extent()
    scale = min(
        float(axis_bounds.width) / max(source_width, 1),
        float(axis_bounds.height) / max(source_height, 1),
        1.0,
    )
    raster_scale = min(1.0, scale * _DUPLEX_RASTER_OVERSAMPLE)
    target_width = max(1, round(source_width * raster_scale))
    target_height = max(1, round(source_height * raster_scale))
    native_cap_height = _duplex_native_cap_height_px(renderer)
    displayed_cap_height = (
        native_cap_height * scale if native_cap_height is not None else None
    )
    if target_width == source_width and target_height == source_height:
        return frame, displayed_cap_height
    resized = Image.fromarray(frame).resize(
        (target_width, target_height),
        resample=Image.Resampling.LANCZOS,
    )
    return np.asarray(resized), displayed_cap_height


def _duplex_native_cap_height_px(renderer: DuplexFrameRenderer) -> float | None:
    owner = getattr(renderer, "__self__", None)
    value = getattr(owner, "native_nucleotide_cap_height_px", None)
    if value is None:
        return None
    numeric = float(value)
    return numeric if numeric > 0 else None


def _preferred_figure_height_inches(renderer: DuplexFrameRenderer) -> float:
    owner = getattr(renderer, "__self__", None)
    value = getattr(owner, "preferred_figure_height_inches", 2.4)
    numeric = float(value)
    if not 0.0 < numeric < float("inf"):
        raise ValueError("preferred_figure_height_inches must be finite and positive")
    return numeric


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


def _draw_document(
    document: PlaybackDocument,
    *,
    transition_index: int,
    progress: float,
    figure,
    duplex_frame_renderer: DuplexFrameRenderer | None = None,
) -> None:
    if not document.plan.steps:
        raise ValueError("playback document requires at least one step")
    transition_index = max(0, min(transition_index, len(document.plan.steps)))
    progress = max(0.0, min(progress, 1.0))
    figure.clear()
    figure.set_facecolor(_PAPER)
    graph_axis, duplex_axis, legend_axis = _document_axes(figure, document)
    duplex_frame = None
    displayed_duplex_cap_height = None
    if duplex_frame_renderer is not None:
        duplex_frame, displayed_duplex_cap_height = _duplex_frame_for_axis(
            duplex_frame_renderer,
            document,
            transition_index,
            progress,
            duplex_axis,
        )
    if graph_axis is not None:
        _draw_graph(
            graph_axis,
            document,
            transition_index,
            progress,
            kmer_font_size_pt=_graph_font_size_for_cap_height(
                displayed_duplex_cap_height,
                figure.dpi,
            ),
        )
    if duplex_frame_renderer is None:
        step_index = min(transition_index, len(document.plan.steps) - 1)
        _draw_fallback_duplex(duplex_axis, document, step_index)
    else:
        duplex_axis.imshow(
            duplex_frame,
            interpolation="lanczos",
            resample=True,
        )
        duplex_axis.axis("off")
    if legend_axis is not None:
        _draw_legend(legend_axis, document)


def _document_axes(figure, document: PlaybackDocument):
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
            top=0.985,
            bottom=0.015,
            hspace=0.0,
        )
        return None, figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[1, 0])
    if not show_graph:
        return None, figure.add_subplot(1, 1, 1), None
    if has_legend:
        grid = figure.add_gridspec(
            2,
            2,
            width_ratios=(graph_fraction, 1.0 - graph_fraction),
            height_ratios=(0.82, 0.18),
            left=0.004,
            right=0.997,
            top=0.985,
            bottom=0.015,
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
        top=0.985,
        bottom=0.015,
        wspace=0.001,
    )
    graph_axis = figure.add_subplot(grid[0, 0])
    duplex_axis = figure.add_subplot(grid[0, 1])
    legend_axis = None
    return graph_axis, duplex_axis, legend_axis


def _draw_legend(axis, document: PlaybackDocument) -> None:
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
        center = group_start + segment * (index + 0.5)
        axis.scatter(
            (center - 0.070,),
            (0.5,),
            transform=axis.transAxes,
            marker="s",
            s=112,
            facecolor=entry.color,
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
            color=_GRAPH_TEXT,
            fontsize=11.5,
            family=KMER_FONT_FAMILY,
            fontweight="normal",
        )


def _duplex_transition_frame(
    renderer: DuplexFrameRenderer,
    document: PlaybackDocument,
    transition_index: int,
    progress: float,
):
    import numpy as np

    final_index = len(document.plan.steps) - 1
    if transition_index > final_index:
        return renderer(document, final_index)
    current = np.asarray(renderer(document, transition_index), dtype=np.float32)
    previous = (
        np.full_like(current, 255.0)
        if transition_index == 0
        else np.asarray(renderer(document, transition_index - 1), dtype=np.float32)
    )
    if previous.shape != current.shape:
        return current.astype(np.uint8)
    if progress <= 0.0:
        return previous.astype(np.uint8)
    if progress >= 0.999:
        return current.astype(np.uint8)
    difference = np.max(np.abs(current[..., :3] - previous[..., :3]), axis=2)
    content_mask = (difference > 1.5).astype(np.float32)
    if not np.any(content_mask):
        return current.astype(np.uint8)
    step = document.plan.steps[transition_index]
    orientation = getattr(step.orientation, "value", step.orientation)
    direction = 1 if str(orientation) == "rev" else -1
    settle_start = 0.42
    if progress <= settle_start:
        return previous.astype(np.uint8)
    local_progress = min(1.0, (progress - settle_start) / (1.0 - settle_start))
    overshoot = 1.4
    shifted = local_progress - 1.0
    settled = 1.0 + (overshoot + 1.0) * shifted**3 + overshoot * shifted**2
    offset = round(direction * 28.0 * (1.0 - settled))
    opacity = _smoothstep(min(1.0, local_progress / 0.55))
    output = previous.copy()
    height = current.shape[0]
    source_start = max(0, -offset)
    source_end = min(height, height - offset)
    destination_start = source_start + offset
    destination_end = source_end + offset
    if source_start >= source_end:
        return output.astype(np.uint8)
    incoming = current[source_start:source_end]
    alpha = content_mask[source_start:source_end, :, None] * opacity
    destination = output[destination_start:destination_end]
    output[destination_start:destination_end] = (destination * (1.0 - alpha)) + (
        incoming * alpha
    )
    return np.clip(output, 0, 255).astype(np.uint8)


def _transition_frame_counts(
    document: PlaybackDocument,
    figure,
    *,
    fps: int,
    seconds_per_step: float,
) -> tuple[int, ...]:
    figure.clear()
    graph_axis, _duplex_axis, _legend_axis = _document_axes(figure, document)
    if graph_axis is None:
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
        raise ValueError(
            f"playback requires {expected} traversal edges, found {len(lengths)}"
        )
    mean_length = sum(lengths) / len(lengths)
    counts = tuple(
        max(1, round(fps * seconds_per_step * length / mean_length))
        for length in lengths
    )
    figure.clear()
    return counts


def render_collection_poster_png(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    dpi: int = 180,
    duplex_frame_renderer: DuplexFrameRenderer | None = None,
) -> Path:
    _require_documents(documents)
    if dpi <= 0:
        raise ValueError("dpi must be positive")
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(
        figsize=(16, _preferred_figure_height_inches(duplex_frame_renderer)),
        facecolor=_PAPER,
    )
    document = documents[0]
    _draw_document(
        document,
        transition_index=len(document.plan.steps),
        progress=1.0,
        figure=figure,
        duplex_frame_renderer=duplex_frame_renderer,
    )
    figure.savefig(output_path, dpi=dpi, facecolor=_PAPER)
    plt.close(figure)
    return output_path


def render_collection_mp4(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    fps: int = 30,
    seconds_per_step: float = 0.70,
    hold_seconds: float = 0.75,
    lead_seconds: float = 0.25,
    scene_transition_seconds: float = 0.0,
    duplex_frame_renderer: DuplexFrameRenderer | None = None,
) -> Path:
    _require_documents(documents)
    if (
        fps <= 0
        or seconds_per_step <= 0
        or hold_seconds < 0
        or lead_seconds < 0
        or scene_transition_seconds < 0
    ):
        raise ValueError(
            "fps and seconds_per_step must be positive; hold_seconds, lead_seconds, "
            "and scene_transition_seconds must be non-negative"
        )
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(
        figsize=(16, _preferred_figure_height_inches(duplex_frame_renderer)),
        facecolor=_PAPER,
    )
    writer = FFMpegWriter(
        fps=fps,
        codec="h264",
        metadata={"title": "Dense-array solution playback"},
        extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
    )
    hold_frames = max(1, round(fps * hold_seconds))
    lead_frames = max(1, round(fps * lead_seconds))
    transition_frames = round(fps * scene_transition_seconds)
    with writer.saving(figure, str(output_path), dpi=150):
        for document_index, document in enumerate(documents):
            transition_frame_counts = _transition_frame_counts(
                document,
                figure,
                fps=fps,
                seconds_per_step=seconds_per_step,
            )
            _draw_document(
                document,
                transition_index=0,
                progress=0.0,
                figure=figure,
                duplex_frame_renderer=duplex_frame_renderer,
            )
            if document_index > 0 and transition_frames:
                for frame_index in range(transition_frames):
                    _draw_document(
                        document,
                        transition_index=0,
                        progress=0.0,
                        figure=figure,
                        duplex_frame_renderer=duplex_frame_renderer,
                    )
                    figure.add_artist(
                        __import__(
                            "matplotlib.patches", fromlist=["Rectangle"]
                        ).Rectangle(
                            (0.0, 0.0),
                            1.0,
                            1.0,
                            transform=figure.transFigure,
                            facecolor=_PAPER,
                            edgecolor="none",
                            alpha=1.0 - ((frame_index + 1) / transition_frames),
                            zorder=1000,
                        )
                    )
                    writer.grab_frame(facecolor=_PAPER)
            for _ in range(lead_frames):
                writer.grab_frame(facecolor=_PAPER)
            for transition_index, transition_frames in enumerate(
                transition_frame_counts
            ):
                for frame_index in range(transition_frames):
                    progress = (frame_index + 1) / transition_frames
                    _draw_document(
                        document,
                        transition_index=transition_index,
                        progress=progress,
                        figure=figure,
                        duplex_frame_renderer=duplex_frame_renderer,
                    )
                    writer.grab_frame(facecolor=_PAPER)
            _draw_document(
                document,
                transition_index=len(document.plan.steps),
                progress=1.0,
                figure=figure,
                duplex_frame_renderer=duplex_frame_renderer,
            )
            for _ in range(hold_frames):
                writer.grab_frame(facecolor=_PAPER)
            if document_index < len(documents) - 1 and transition_frames:
                for frame_index in range(transition_frames):
                    _draw_document(
                        document,
                        transition_index=len(document.plan.steps),
                        progress=1.0,
                        figure=figure,
                        duplex_frame_renderer=duplex_frame_renderer,
                    )
                    figure.add_artist(
                        __import__(
                            "matplotlib.patches", fromlist=["Rectangle"]
                        ).Rectangle(
                            (0.0, 0.0),
                            1.0,
                            1.0,
                            transform=figure.transFigure,
                            facecolor=_PAPER,
                            edgecolor="none",
                            alpha=(frame_index + 1) / transition_frames,
                            zorder=1000,
                        )
                    )
                    writer.grab_frame(facecolor=_PAPER)
    plt.close(figure)
    return output_path


def render_collection_gif(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    fps: int = 15,
    seconds_per_step: float = 0.70,
    hold_seconds: float = 0.70,
    lead_seconds: float = 0.25,
    scene_transition_seconds: float = 0.0,
    duplex_frame_renderer: DuplexFrameRenderer | None = None,
) -> Path:
    _require_documents(documents)
    if (
        fps <= 0
        or seconds_per_step <= 0
        or hold_seconds < 0
        or lead_seconds < 0
        or scene_transition_seconds < 0
    ):
        raise ValueError(
            "fps and seconds_per_step must be positive; hold_seconds, lead_seconds, and scene_transition_seconds must be non-negative"
        )
    import matplotlib.pyplot as plt
    from matplotlib.animation import PillowWriter

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure = plt.figure(
        figsize=(16, _preferred_figure_height_inches(duplex_frame_renderer)),
        facecolor=_PAPER,
    )
    writer = PillowWriter(fps=fps)
    hold_frames = max(1, round(fps * hold_seconds))
    lead_frames = max(1, round(fps * lead_seconds))
    with writer.saving(figure, str(output_path), dpi=100):
        for document in documents:
            transition_frame_counts = _transition_frame_counts(
                document,
                figure,
                fps=fps,
                seconds_per_step=seconds_per_step,
            )
            _draw_document(
                document,
                transition_index=0,
                progress=0.0,
                figure=figure,
                duplex_frame_renderer=duplex_frame_renderer,
            )
            for _ in range(lead_frames):
                writer.grab_frame(facecolor=_PAPER)
            for transition_index, transition_frames in enumerate(
                transition_frame_counts
            ):
                for frame_index in range(transition_frames):
                    progress = (frame_index + 1) / transition_frames
                    _draw_document(
                        document,
                        transition_index=transition_index,
                        progress=progress,
                        figure=figure,
                        duplex_frame_renderer=duplex_frame_renderer,
                    )
                    writer.grab_frame(facecolor=_PAPER)
            _draw_document(
                document,
                transition_index=len(document.plan.steps),
                progress=1.0,
                figure=figure,
                duplex_frame_renderer=duplex_frame_renderer,
            )
            for _ in range(hold_frames):
                writer.grab_frame(facecolor=_PAPER)
    plt.close(figure)
    return output_path
