"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/projections.py

Scalar sequence and placement records for interoperable library handoffs.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING

from dense_arrays.realized import Orientation

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping

    from dense_arrays.artifacts.records import Design
    from dense_arrays.parts.models import Part


@dataclass(frozen=True)
class SequenceRecord:
    """One final sequence, joined by design identity rather than DNA equality."""

    design_ref: str
    sequence_id: str
    plan_id: str
    sequence: str
    length: int
    gc_fraction: float

    def to_dict(self) -> dict[str, object]:
        """Encode scalar fields without losing source or sequence identity."""
        return {"schema": "dense_arrays.sequence_record.v1", **asdict(self)}


@dataclass(frozen=True)
class PlacementRecord:
    """One realized placement with zero-based half-open final coordinates."""

    design_ref: str
    sequence_id: str
    plan_id: str
    placement_id: str
    collection_id: str
    part_id: str
    group: str | None
    sequence: str
    orientation: str
    start: int
    end: int
    core_start: int | None
    core_end: int | None
    core_orientation: str | None

    def to_dict(self) -> dict[str, object]:
        """Keep unknown core annotations null, including in tabular exports."""
        return {"schema": "dense_arrays.placement_record.v1", **asdict(self)}


def project(
    design: Design, view: str, parts: Mapping[str, Part], collection_id: str
) -> Iterator[Design | SequenceRecord | PlacementRecord]:
    """Project stored evidence; do not solve, sample or infer new annotations."""
    if view == "designs":
        yield design
    elif view == "sequences":
        sequence = design.realized.sequence
        yield SequenceRecord(
            design.reference,
            design.sequence_id,
            design.plan_id,
            sequence,
            len(sequence),
            (sequence.count("G") + sequence.count("C")) / len(sequence),
        )
    else:
        for placement in design.realized.placements:
            part = parts[placement.feature_id]
            start, end, orientation = (
                part.core_start,
                part.core_end,
                part.core_orientation,
            )
            if start is not None:
                if placement.orientation == Orientation.REVERSE:
                    start, end = len(part.sequence) - end, len(part.sequence) - start
                    orientation = "reverse" if orientation == "forward" else "forward"
                start, end = placement.start + start, placement.start + end
            direction = {
                Orientation.FORWARD: "forward",
                Orientation.REVERSE: "reverse",
                Orientation.UNSPECIFIED: "unspecified",
            }[placement.orientation]
            yield PlacementRecord(
                design.reference,
                design.sequence_id,
                design.plan_id,
                placement.placement_id,
                collection_id,
                part.part_id,
                part.group,
                placement.sequence,
                direction,
                placement.start,
                placement.end,
                start,
                end,
                orientation,
            )
