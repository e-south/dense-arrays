"""Shared reveal and complement semantics for playback renderers."""

from __future__ import annotations

from collections.abc import Sequence

from .models import PlaybackStep

_IUPAC_COMPLEMENTS = str.maketrans("ATCGRYSWKMBDHVN", "TAGCYRSWMKVHDBN")


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
