"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/models.py

Supplied-array identities and evidence boundaries, independent of execution.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    mutable_json,
    records,
    required_text,
    semantic_digest,
)
from dense_arrays.artifacts.reading import ReadCost, ReadLimits
from dense_arrays.parts import Part
from dense_arrays.playback.serialization import realized_array_to_dict

if TYPE_CHECKING:
    from dense_arrays.realized import RealizedArray

SCHEMA = "dense_arrays.array_collection.v1"
BOUNDARY = "supplied_sequences_parts_and_placements"
MANIFEST = "collection.json"
DATABASE = "arrays.sqlite3"


@dataclass(frozen=True)
class ArrayCollection:
    """Parts and an iterable of supplied arrays, consumed once during publication."""

    parts: tuple[Part, ...]
    arrays: Iterable[RealizedArray]
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Freeze the part catalog while leaving the array stream incremental."""
        parts = records(self.parts, Part, field_name="parts")
        if len({part.part_id for part in parts}) != len(parts):
            msg = "collection part IDs must be unique; use explicit source namespaces"
            raise ValueError(msg)
        if not isinstance(self.arrays, Iterable) or isinstance(
            self.arrays, (str, bytes)
        ):
            msg = "arrays must be an iterable of RealizedArray records"
            raise TypeError(msg)
        object.__setattr__(self, "parts", parts)
        object.__setattr__(self, "provenance", immutable_json_mapping(self.provenance))


@dataclass(frozen=True)
class ArrayFilter:
    """Match array identities and arrays containing any specified part or group."""

    array_ids: tuple[str, ...] = ()
    part_ids: tuple[str, ...] = ()
    groups: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Reject duplicate or blank selectors; preserve supplied group labels."""
        for name in ("array_ids", "part_ids", "groups"):
            value = getattr(self, name)
            if not isinstance(value, (tuple, list)):
                msg = f"{name} must be a list or tuple"
                raise TypeError(msg)
            for item in value:
                required_text(item, field_name=name)
            if len(set(value)) != len(value):
                msg = f"{name} must not repeat identities"
                raise ValueError(msg)
            object.__setattr__(self, name, tuple(value))

    def to_dict(self) -> dict[str, object]:
        """Bind cursors and exports to the exact predicate."""
        return {
            name: list(getattr(self, name))
            for name in ("array_ids", "part_ids", "groups")
        }


def sequence_identity(sequence: str) -> str:
    """Use the same versioned exact-DNA identity as generated designs."""
    return semantic_digest({"schema": "dense_arrays.sequence.v1", "sequence": sequence})


@dataclass(frozen=True)
class ArrayRecord:
    """One supplied realization scoped to its immutable collection."""

    collection_id: str
    realized: RealizedArray

    @property
    def array_id(self) -> str:
        """Preserve the caller's array identity independently of sequence equality."""
        return self.realized.source_id

    @property
    def sequence_id(self) -> str:
        """Identify exact sequence equality without merging annotations."""
        return sequence_identity(self.realized.sequence)

    def to_dict(self) -> dict[str, object]:
        """Serialize geometry without inserting plan, solver or attempt fields."""
        return {
            "schema": "dense_arrays.array_record.v1",
            "collection_id": self.collection_id,
            "sequence_id": self.sequence_id,
            "realized": realized_array_to_dict(self.realized),
        }


@dataclass(frozen=True)
class CollectionSummary:
    """Committed inventory; verification is explicit and scoped to supplied data."""

    manifest: Mapping[str, object]
    verified: bool = False
    read_limits: ReadLimits = field(default_factory=ReadLimits)

    def __post_init__(self) -> None:
        """Detach stored metadata from mutable dictionaries."""
        object.__setattr__(self, "manifest", immutable_json_mapping(self.manifest))

    @property
    def collection_id(self) -> str:
        """Return the content-bound collection identity."""
        return self.manifest["collection_id"]

    @property
    def arrays(self) -> int:
        """Count supplied realizations, including separately identified equal DNA."""
        return self.manifest["arrays"]

    @property
    def placements(self) -> int:
        """Count annotation occurrences, including overlaps and repeated parts."""
        return self.manifest["placements"]

    @property
    def parts(self) -> int:
        """Count all supplied parts, including unused parts."""
        return self.manifest["parts"]

    @property
    def cost(self) -> ReadCost:
        """Describe the metadata-only summary read."""
        return ReadCost(
            self.collection_id, 0, "manifest", "summary", 1, self.read_limits
        )

    @property
    def verification_cost(self) -> ReadCost:
        """Expose the complete catalog/array scan and database byte bound."""
        return ReadCost(
            self.collection_id,
            0,
            "scan",
            "verification",
            self.parts + self.arrays,
            self.read_limits,
            bytes_estimate=self.manifest["database"]["bytes"],
        )

    def to_dict(self) -> dict[str, object]:
        """Keep original source provenance separate from collection export software."""
        return {**mutable_json(self.manifest), "verified": self.verified}
