"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/libraries.py

Frozen accepted libraries and explicit sequence-exclusion policies.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    integer,
    object_fields,
    records,
    required_text,
    semantic_digest,
)
from dense_arrays.planning.lineage import ParentRun

if TYPE_CHECKING:
    from pathlib import Path

_REFERENCE_COMPONENTS = 3


@dataclass(frozen=True)
class ExcludedDesign:
    """One exact final-sequence exclusion with its original design evidence."""

    design_ref: str
    sequence_id: str
    record_digest: str
    cell_id: str = "default"

    def __post_init__(self) -> None:
        """Validate the recorded cell identity and original design reference."""
        required_text(self.design_ref, field_name="design_ref")
        digest(self.sequence_id, field_name="sequence_id")
        digest(self.record_digest, field_name="record_digest")
        if (
            not self.cell_id
            or len(self.design_ref.split("/")) != _REFERENCE_COMPONENTS
            or self.design_ref.split("/")[1] != self.cell_id
        ):
            msg = "excluded design requires a full reference matching its cell"
            raise ValueError(msg)


@dataclass(frozen=True)
class ParentLibrary:
    """One cell of a verified parent snapshot, with its inherited exclusions."""

    run_id: str
    plan_id: str
    revision: int
    state: str
    target: int
    accepted: int
    accepted_digest: str
    exclusions: tuple[ExcludedDesign, ...]
    cell_id: str = "default"

    def __post_init__(self) -> None:
        """Preserve separate parent attainment and ancestor exclusion identities."""
        required_text(self.run_id, field_name="parent.run_id")
        digest(self.plan_id, field_name="parent.plan_id")
        digest(self.accepted_digest, field_name="parent.accepted_digest")
        integer(self.revision, field_name="parent.revision", minimum=0)
        integer(self.target, field_name="parent.target", minimum=0)
        integer(self.accepted, field_name="parent.accepted", minimum=0)
        required_text(self.cell_id, field_name="parent.cell_id")
        if (
            self.state not in {"completed", "stopped", "failed", "inactive"}
            or self.accepted > self.target
            or (self.state in {"completed", "inactive"})
            != (self.accepted == self.target)
            or (self.state == "inactive" and self.target != 0)
        ):
            msg = "parent must have consistent terminal attainment"
            raise ValueError(msg)
        exclusions = records(self.exclusions, ExcludedDesign, field_name="exclusions")
        if any(e.cell_id != self.cell_id for e in exclusions):
            msg = "parent exclusions must belong to the declared cell"
            raise ValueError(msg)
        if len({e.sequence_id for e in exclusions}) != len(exclusions) or len(
            {e.design_ref for e in exclusions}
        ) != len(exclusions):
            msg = "parent exclusions must have unique sequence and design identities"
            raise ValueError(msg)
        object.__setattr__(self, "exclusions", exclusions)
        direct = tuple(
            e for e in exclusions if e.design_ref.startswith(self.run_id + "/")
        )
        if (
            len(direct) != self.accepted
            or accepted_digest(direct) != self.accepted_digest
        ):
            msg = "parent accepted digest or count does not match its exclusions"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Encode all ancestor exclusions so later execution needs no parent path."""
        value = {
            f.name: getattr(self, f.name)
            for f in fields(self)
            if f.name != "exclusions"
        }
        value["exclusions"] = [asdict(e) for e in self.exclusions]
        if self.cell_id == "default":
            value.pop("cell_id")
        return {"schema": "dense_arrays.parent_library.v1", **value}

    @classmethod
    def from_dict(cls, value: object) -> ParentLibrary:
        """Reject unsupported lineage or changed exclusion evidence."""
        data = object_fields(
            value,
            {
                "schema",
                "run_id",
                "plan_id",
                "revision",
                "state",
                "target",
                "accepted",
                "accepted_digest",
                "exclusions",
                "cell_id",
            },
            "parent library",
        )
        if data.pop("schema", None) != "dense_arrays.parent_library.v1":
            msg = "unsupported parent-library schema"
            raise ValueError(msg)
        data.setdefault("cell_id", "default")
        if set(data) != {f.name for f in fields(cls)}:
            msg = "parent library requires every canonical field"
            raise ValueError(msg)
        if not isinstance(data["exclusions"], list):
            msg = "parent exclusions must be an array"
            raise TypeError(msg)
        data["exclusions"] = tuple(
            ExcludedDesign(
                **object_fields(
                    e,
                    {"design_ref", "sequence_id", "record_digest", "cell_id"},
                    "exclusion",
                )
            )
            for e in data["exclusions"]
        )
        return cls(**data)


def accepted_digest(exclusions: tuple[ExcludedDesign, ...]) -> str:
    """Bind the ordered accepted records, not only their sequence equality class."""
    return semantic_digest(
        {
            "schema": "dense_arrays.accepted_library.v1",
            "designs": [asdict(e) for e in exclusions],
        }
    )


@dataclass(frozen=True)
class LibraryExclusion:
    """Exclude a declared accepted library under an explicit cell mapping."""

    source: ParentRun | ParentLibrary
    cell_mapping: Mapping[str, str]
    uniqueness: str

    def __post_init__(self) -> None:
        """Validate the supported scope before opening any source artifact."""
        if not isinstance(self.source, (ParentRun, ParentLibrary)):
            msg = "exclusion source must be ParentRun or a frozen ParentLibrary"
            raise TypeError(msg)
        cell = (
            self.source.cell_id if isinstance(self.source, ParentLibrary) else "default"
        )
        if not isinstance(self.cell_mapping, Mapping) or dict(self.cell_mapping) != {
            cell: cell
        }:
            msg = "exclusion cell_mapping must preserve the source cell identity"
            raise ValueError(msg)
        if self.uniqueness != "exact_sequence_per_cell.v1":
            msg = "exclusion uniqueness must be exact_sequence_per_cell.v1"
            raise ValueError(msg)
        object.__setattr__(
            self, "cell_mapping", MappingProxyType(dict(self.cell_mapping))
        )

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Encode an unresolved source locator or a complete frozen snapshot."""
        source = self.source
        return {
            "source": source.to_dict()
            if isinstance(source, ParentLibrary)
            else {
                "run": str(source.run)
                if base is None
                else os.path.relpath(source.run.absolute(), base)
            },
            "cell_mapping": dict(self.cell_mapping),
            "uniqueness": self.uniqueness,
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> LibraryExclusion:
        """Read explicit source, mapping and uniqueness without guessing policies."""
        data = object_fields(value, {"source", "cell_mapping", "uniqueness"}, "exclude")
        source = data.pop("source", None)
        if isinstance(source, Mapping) and source.get("schema") is not None:
            source = ParentLibrary.from_dict(source)
        else:
            source = ParentRun(**object_fields(source, {"run"}, "exclude.source"))
            if base is not None and not source.run.is_absolute():
                source = ParentRun(base / source.run)
        return cls(source=source, **data)


def exclusion_count(value: object) -> int:
    """Count wire identities before constructing a frozen exclusion collection."""
    data = object_fields(value, {"source", "cell_mapping", "uniqueness"}, "exclude")
    source = data.get("source")
    if not isinstance(source, Mapping) or not isinstance(
        source.get("exclusions"), list
    ):
        msg = "resolved exclusion requires a frozen library with an exclusions array"
        raise TypeError(msg)
    return len(source["exclusions"])
