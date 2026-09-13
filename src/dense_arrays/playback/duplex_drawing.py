"""Draw sequence, complement, placement bars, and caller labels.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .duplex_geometry import FEATURE_HEIGHT, draw_nucleotide, fit_duplex_grid
from .theme import RESTING_COLOR, RESTING_TEXT_COLOR, blend_color
from .timeline import complement_sequence, placement_progress

if TYPE_CHECKING:
    from collections.abc import Sequence

    from matplotlib.axes import Axes

    from .models import PlaybackStep
    from .presentation import PlaybackDocument

_INK = "#4b5563"
_TRACK_PITCH = 1.1
_CAPTION_TRACK_PITCH = 1.35


def _placement_tracks(
    steps: Sequence[PlaybackStep], *, with_labels: bool
) -> tuple[float, ...]:
    """Allocate non-overlapping lanes once from the complete placement plan."""
    lane_ends: dict[bool, list[int]] = {False: [], True: []}
    pitch = _CAPTION_TRACK_PITCH if with_labels else _TRACK_PITCH
    positions = []
    for step in steps:
        reverse = step.orientation == "rev"
        ends = lane_ends[reverse]
        lane = next(
            (index for index, end in enumerate(ends) if end <= step.start), len(ends)
        )
        if lane == len(ends):
            ends.append(step.end)
        else:
            ends[lane] = step.end
        positions.append(-1.45 - lane * pitch if reverse else 1.05 + lane * pitch)
    return tuple(positions)


def draw_duplex(
    axis: Axes,
    document: PlaybackDocument,
    step_index: int,
    progress: float = 1.0,
    *,
    bottom_padding_pt: float = 0,
) -> float:
    """Draw the fixed duplex and return its nucleotide cap height in pixels."""
    from matplotlib.patches import FancyBboxPatch

    plan = document.plan
    sequence = plan.realized_sequence
    complement = complement_sequence(sequence)
    length = len(sequence)
    tracks = _placement_tracks(plan.steps, with_labels=bool(document.label_overrides))
    geometry = fit_duplex_grid(
        axis,
        sequence,
        (min(-2.2, min(tracks) - 0.75), max(2.2, max(tracks) + 1.25)),
        bottom_padding_pt=bottom_padding_pt,
    )
    axis.axis("off")
    for index, step in enumerate(plan.steps):
        y = tracks[index]
        emphasis = placement_progress(index, step_index, progress)
        color = blend_color(RESTING_COLOR, document.step_color(index), emphasis)
        axis.add_patch(
            FancyBboxPatch(
                (step.start, y),
                step.end - step.start,
                FEATURE_HEIGHT,
                boxstyle="round,pad=0.01,rounding_size=0.08",
                facecolor=color,
                edgecolor="none",
                linewidth=0,
            )
        )
        feature_sequence = (
            complement_sequence(step.placement_sequence)
            if step.orientation == "rev"
            else step.placement_sequence
        )
        for offset, base in enumerate(feature_sequence):
            draw_nucleotide(
                axis,
                (step.start + offset + 0.5, y + FEATURE_HEIGHT / 2),
                base,
                blend_color(RESTING_TEXT_COLOR, "#FFFFFF", emphasis),
                geometry,
            )
        if step.placement_id in document.label_overrides:
            reverse = step.orientation == "rev"
            axis.text(
                (step.start + step.end) / 2,
                y - 0.1 if reverse else y + FEATURE_HEIGHT + 0.1,
                document.label_overrides[step.placement_id],
                ha="center",
                va="top" if reverse else "bottom",
                color=blend_color(RESTING_TEXT_COLOR, _INK, emphasis),
                fontsize=geometry.label_font_size_pt,
                family=geometry.font.get_family(),
            )
    coordinate_steps = {
        coordinate: index
        for index, step in enumerate(plan.steps)
        for span in step.added_spans
        for coordinate in range(span.start, span.end)
    }
    for index in range(length):
        placement = coordinate_steps.get(index)
        emphasis = (
            0
            if placement is None
            else placement_progress(placement, step_index, progress)
        )
        color = blend_color(RESTING_COLOR, _INK, emphasis)
        draw_nucleotide(axis, (index + 0.5, 0.30), sequence[index], color, geometry)
        draw_nucleotide(axis, (index + 0.5, -0.30), complement[index], color, geometry)
    terminal_color = blend_color(
        RESTING_TEXT_COLOR, _INK, placement_progress(0, step_index, progress)
    )
    axis.text(
        -2.0, 0.30, "5'", ha="center", va="center", color=terminal_color, fontsize=11
    )
    axis.text(
        -2.0, -0.30, "3'", ha="center", va="center", color=terminal_color, fontsize=11
    )
    axis.text(
        length + 0.7,
        0.30,
        "3'",
        ha="left",
        va="center",
        color=terminal_color,
        fontsize=11,
    )
    axis.text(
        length + 0.7,
        -0.30,
        "5'",
        ha="left",
        va="center",
        color=terminal_color,
        fontsize=11,
    )
    return geometry.cap_height_px
