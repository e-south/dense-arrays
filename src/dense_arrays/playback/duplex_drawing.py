"""Draw sequence, complement, placement bars, and caller labels.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .graph.geometry import KMER_FONT_FAMILY
from .timeline import complement_sequence, current_added_indices, revealed_indices

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from .presentation import PlaybackDocument

_ACTIVE = "#167a70"
_INK = "#4b5563"


def draw_duplex(axis: Axes, document: PlaybackDocument, step_index: int) -> None:
    """Draw the nucleotide reveal mask and complete placement bars."""
    from matplotlib.patches import FancyBboxPatch

    plan = document.plan
    sequence = plan.realized_sequence
    revealed = revealed_indices(plan.steps, step_index)
    current_added = frozenset(current_added_indices(plan.steps, step_index))
    complement = complement_sequence(sequence)
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
        color = document.step_color(index)
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
        axis.text(
            (step.start + step.end) / 2,
            y + 0.48,
            document.step_label(index),
            ha="center",
            va="bottom",
            color=_INK,
            fontsize=7,
        )
    font_size = max(5.2, min(10.5, 830 / max(1, length)))
    for index in sorted(revealed):
        color = _ACTIVE if index in current_added else _INK
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
