"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/projections.py

Sequence and placement tables for supplied-array collections.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING

from dense_arrays.parts.geometry import core_interval
from dense_arrays.parts.serialization import part_to_dict
from dense_arrays.realized import Orientation

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping

    from dense_arrays.parts import Part

    from .models import ArrayRecord


@dataclass(frozen=True)
class ArraySequence:
    """One sequence with its preserved source ID and collection namespace."""

    collection_id: str
    array_id: str
    sequence_id: str
    sequence: str
    length: int
    gc_fraction: float

    def to_dict(self) -> dict[str, object]:
        """Encode a scalar table row without construction-run fields."""
        return {"schema": "dense_arrays.array_sequence.v1", **asdict(self)}


@dataclass(frozen=True)
class ArrayPlacement:
    """A supplied occurrence and its core in final-sequence coordinates."""

    collection_id: str
    array_id: str
    sequence_id: str
    placement_id: str
    part_id: str
    kind: str
    group: str | None
    sequence: str
    orientation: str
    start: int
    end: int
    core_start: int | None
    core_end: int | None
    core_orientation: str | None

    def to_dict(self) -> dict[str, object]:
        """Preserve unknown core annotations as null rather than inferred spans."""
        return {"schema": "dense_arrays.array_placement.v1", **asdict(self)}


@dataclass(frozen=True)
class CollectionPart:
    """One entry in the complete bound catalog, including unused supplied parts."""

    collection_id: str
    part: Part

    def to_dict(self) -> dict[str, object]:
        """Keep original source/core annotations and caller metadata intact."""
        return {
            "schema": "dense_arrays.collection_part.v1",
            "collection_id": self.collection_id,
            **part_to_dict(self.part),
        }


def project(
    record: ArrayRecord, view: str, parts: Mapping[str, Part]
) -> Iterator[ArrayRecord | ArraySequence | ArrayPlacement]:
    """Project geometry only; no scoring, generation or runtime inference occurs."""
    if view == "arrays":
        yield record
    elif view == "sequences":
        dna = record.realized.sequence
        yield ArraySequence(
            record.collection_id,
            record.array_id,
            record.sequence_id,
            dna,
            len(dna),
            (dna.count("G") + dna.count("C")) / len(dna),
        )
    else:
        for placement in record.realized.placements:
            part = parts[placement.feature_id]
            start, end, strand = core_interval(part, placement)
            yield ArrayPlacement(
                record.collection_id,
                record.array_id,
                record.sequence_id,
                placement.placement_id,
                part.part_id,
                placement.kind.value,
                part.group,
                placement.sequence,
                "reverse"
                if placement.orientation == Orientation.REVERSE
                else "forward",
                placement.start,
                placement.end,
                start,
                end,
                strand,
            )
