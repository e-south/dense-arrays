"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/provenance.py

Immutable import evidence retained with normalized parts and resolved plans.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    immutable_json_mapping,
    integer,
    mutable_json,
    object_fields,
    records,
    required_text,
)
from dense_arrays.parts.filters import PartFilter
from dense_arrays.parts.models import Normalization

if TYPE_CHECKING:
    from collections.abc import Mapping

IMPORT_SCHEMA = "dense_arrays.import.v1"


@dataclass(frozen=True)
class Transformation:
    """One explicit normalization with its original source location and values."""

    row: int
    field: str
    column: str
    before: str
    after: str

    def __post_init__(self) -> None:
        """Require an actual change and an explicit one-based logical row."""
        integer(self.row, field_name="transformation.row", minimum=1)
        for name in ("field", "column", "before", "after"):
            required_text(getattr(self, name), field_name=f"transformation.{name}")
        if self.field != "sequence" or self.before == self.after:
            msg = "transformations require a changed sequence field"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Return exact source evidence without normalizing it again."""
        return {
            name: getattr(self, name)
            for name in ("row", "field", "column", "before", "after")
        }


@dataclass(frozen=True, repr=False)
class ImportReport:
    """Bounded display with complete explicit source transformations on demand."""

    kind: str
    rows: int
    columns: Mapping[str, str] = field(default_factory=dict)
    ignored_columns: tuple[str, ...] = ()
    transformations: tuple[Transformation, ...] = ()
    normalization: Normalization = field(default_factory=Normalization)
    id_policy: str = "provided"
    collection_id: str | None = None
    selection: PartFilter | None = None

    def __post_init__(self) -> None:
        """Freeze evidence and reject incomplete or contradictory import contracts."""
        if self.kind not in {"inline", "table", "pool"} or self.id_policy not in {
            "provided",
            "row",
        }:
            msg = "unsupported import kind or ID policy"
            raise ValueError(msg)
        integer(self.rows, field_name="import.rows", minimum=0)
        if self.kind == "pool":
            digest(self.collection_id, field_name="collection_id")
        elif self.collection_id is not None or self.selection is not None:
            msg = "collection identity and selection require a pool source"
            raise ValueError(msg)
        if self.selection is not None and not isinstance(self.selection, PartFilter):
            msg = "import selection must be PartFilter"
            raise TypeError(msg)
        if not isinstance(self.normalization, Normalization):
            msg = "import normalization must be Normalization"
            raise TypeError(msg)
        object.__setattr__(self, "columns", immutable_json_mapping(self.columns))
        for name, column in self.columns.items():
            required_text(name, field_name="import.field")
            required_text(column, field_name="import.column")
        object.__setattr__(
            self,
            "ignored_columns",
            records(self.ignored_columns, str, field_name="ignored_columns"),
        )
        object.__setattr__(
            self,
            "transformations",
            records(self.transformations, Transformation, field_name="transformations"),
        )
        if any(
            t.row > self.rows or self.columns.get(t.field) != t.column
            for t in self.transformations
        ):
            msg = "transformation location does not match the import report"
            raise ValueError(msg)

    def __repr__(self) -> str:
        """Display counts without exposing all source values or transformations."""
        return (
            f"ImportReport({self.kind}, {self.rows} rows, "
            f"{len(self.transformations)} transformations)"
        )

    def to_dict(self) -> dict[str, object]:
        """Serialize complete normalization evidence under its own schema."""
        return {
            "schema": IMPORT_SCHEMA,
            "kind": self.kind,
            "rows": self.rows,
            "columns": mutable_json(self.columns),
            "ignored_columns": list(self.ignored_columns),
            "transformations": [t.to_dict() for t in self.transformations],
            "normalization": {
                "uppercase": self.normalization.uppercase,
                "trim_outer_whitespace": self.normalization.trim_outer_whitespace,
            },
            "id_policy": self.id_policy,
            "collection_id": self.collection_id,
            "selection": None if self.selection is None else self.selection.to_dict(),
        }

    @classmethod
    def from_dict(cls, value: object) -> ImportReport:
        """Read a complete versioned report without applying request defaults."""
        keys = {
            "schema",
            "kind",
            "rows",
            "columns",
            "ignored_columns",
            "transformations",
            "normalization",
            "id_policy",
            "collection_id",
            "selection",
        }
        data = object_fields(value, keys, "import_report")
        if set(data) != keys or data.pop("schema") != IMPORT_SCHEMA:
            msg = f"unsupported or incomplete import schema; supported: {IMPORT_SCHEMA}"
            raise ValueError(msg)
        if not isinstance(data["transformations"], list):
            msg = "transformations must be an array"
            raise TypeError(msg)
        data["transformations"] = tuple(
            Transformation(
                **object_fields(
                    t, {"row", "field", "column", "before", "after"}, "transformation"
                )
            )
            for t in data["transformations"]
        )
        data["normalization"] = Normalization(**data["normalization"])
        if data["selection"] is not None:
            data["selection"] = PartFilter.from_dict(data["selection"])
        report = cls(**data)
        if report.to_dict() != value:
            msg = "import report requires every canonical field"
            raise ValueError(msg)
        return report
