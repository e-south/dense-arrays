"""Shared reveal and complement semantics for playback renderers.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

    from .models import PlaybackStep

_IUPAC_COMPLEMENTS = str.maketrans("ATCGRYSWKMBDHVN", "TAGCYRSWMKVHDBN")


def placement_progress(index: int, transition_index: int, progress: float) -> float:
    """Return smooth emphasis from neutral context through a placed step."""
    if index < transition_index:
        return 1.0
    if index > transition_index:
        return 0.0
    progress = max(0.0, min(progress, 1.0))
    return progress * progress * (3.0 - 2.0 * progress)


def revealed_indices(steps: Sequence[PlaybackStep], step_index: int) -> tuple[int, ...]:
    """Return every coordinate revealed through ``step_index``."""
    return tuple(
        sorted(
            {
                index
                for step in steps[: step_index + 1]
                for span in step.added_spans
                for index in range(span.start, span.end)
            }
        )
    )


def current_added_indices(
    steps: Sequence[PlaybackStep], step_index: int
) -> tuple[int, ...]:
    """Return only coordinates first introduced by the current step."""
    return tuple(
        index
        for span in steps[step_index].added_spans
        for index in range(span.start, span.end)
    )


def complement_sequence(sequence: str) -> str:
    """Return the coordinate-aligned IUPAC DNA complement."""
    return sequence.translate(_IUPAC_COMPLEMENTS)
