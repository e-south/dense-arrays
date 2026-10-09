"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/geometry.py

Validate supplied-part joins and project final-sequence core coordinates.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.realized import Orientation, RealizedArray
from dense_arrays.sequence import reverse_complement

if TYPE_CHECKING:
    from collections.abc import Mapping

    from dense_arrays.parts import Part


def validate_array(array: RealizedArray, parts: Mapping[str, Part]) -> None:
    """Check each occurrence against its bound source part and explicit strand."""
    if not isinstance(array, RealizedArray):
        msg = "collection arrays must be RealizedArray records"
        raise TypeError(msg)
    if array.coordinate_space != "realized_sequence":
        msg = "collections require realized_sequence coordinates"
        raise ValueError(msg)
    for placement in array.placements:
        if placement.feature_id not in parts:
            msg = f"array {array.source_id!r}: unknown part {placement.feature_id!r}"
            raise ValueError(msg)
        part = parts[placement.feature_id]
        if placement.orientation == Orientation.UNSPECIFIED:
            msg = (
                "bound part placements require explicit forward or reverse orientation"
            )
            raise ValueError(msg)
        oriented = (
            reverse_complement(part.sequence)
            if placement.orientation == Orientation.REVERSE
            else part.sequence
        )
        if placement.sequence != oriented:
            msg = f"array {array.source_id!r}: placement differs from its oriented part"
            raise ValueError(msg)
