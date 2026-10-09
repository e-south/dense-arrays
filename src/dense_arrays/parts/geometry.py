"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/geometry.py

Project a source motif core through its part's realized placement.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.realized import Orientation

if TYPE_CHECKING:
    from dense_arrays.parts.models import Part
    from dense_arrays.realized import Placement


def core_interval(
    part: Part, placement: Placement
) -> tuple[int | None, int | None, str | None]:
    """Compose the core's source strand with the part's placement orientation."""
    start, end, strand = part.core_start, part.core_end, part.core_orientation
    if start is not None:
        if placement.orientation == Orientation.REVERSE:
            start, end = len(part.sequence) - end, len(part.sequence) - start
            strand = "reverse" if strand == "forward" else "forward"
        start, end = placement.start + start, placement.start + end
    return start, end, strand
