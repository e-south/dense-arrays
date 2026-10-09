"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/models.py

Immutable identities and source declarations for eligible parts.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    integer,
    required_text,
)
from dense_arrays.problem import motif_library
from dense_arrays.sequence import reverse_complement

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence


@dataclass(frozen=True)
class Part:
    """One eligible supplied occurrence; equal sequences retain distinct IDs."""

    part_id: str
    sequence: str
    group: str | None = None
    source: str | None = None
    core_start: int | None = None
    core_end: int | None = None
    core_orientation: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate sequence/core coordinates and detach mutable caller metadata."""
        required_text(self.part_id, field_name="part_id")
        motif_library([self.sequence])
        for name in ("group", "source"):
            value = getattr(self, name)
            if value is not None:
                required_text(value, field_name=name)
        self._validate_core()
        object.__setattr__(self, "metadata", immutable_json_mapping(self.metadata))

    def _validate_core(self) -> None:
        core = (self.core_start, self.core_end, self.core_orientation)
        if all(value is None for value in core):
            return
        if any(value is None for value in core):
            msg = "core annotations require start, end and orientation"
            raise ValueError(msg)

        start = integer(self.core_start, field_name="core_start", minimum=0)
        end = integer(self.core_end, field_name="core_end", minimum=1)
        if not start < end <= len(self.sequence):
            msg = "core interval must fit the supplied sequence"
            raise ValueError(msg)
        if self.core_orientation not in {"forward", "reverse"}:
            msg = "core_orientation must be forward or reverse"
            raise ValueError(msg)

    @property
    def core_sequence(self) -> str | None:
        """Return the annotated core in its declared orientation, or unknown."""
        if self.core_start is None:
            return None
        core = self.sequence[self.core_start : self.core_end]
        return reverse_complement(core) if self.core_orientation == "reverse" else core


@dataclass(frozen=True)
class PartSelector:
    """Select supplied occurrences by explicit IDs or caller-defined groups."""

    part_ids: tuple[str, ...] = ()
    groups: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Freeze selectors and reject ambiguous or empty predicates."""
        for name in ("part_ids", "groups"):
            values = getattr(self, name)
            if not isinstance(values, (tuple, list)):
                msg = f"{name} must be a list or tuple of labels"
                raise TypeError(msg)
            for value in values:
                required_text(value, field_name=name)
            if len(values) != len(set(values)):
                msg = f"{name} must not repeat labels"
                raise ValueError(msg)
            object.__setattr__(self, name, tuple(values))
        if bool(self.part_ids) == bool(self.groups):
            msg = "select either part_ids or groups, with at least one label"
            raise ValueError(msg)

    def indices(self, parts: Sequence[Part]) -> tuple[int, ...]:
        """Resolve the selector and fail for every unknown label."""
        field_name = "part_id" if self.part_ids else "group"
        labels = set(self.part_ids or self.groups)
        available = {getattr(part, field_name) for part in parts}
        if missing := labels - available:
            msg = f"unknown {field_name} values: {sorted(missing)}"
            raise ValueError(msg)
        return tuple(
            i for i, part in enumerate(parts) if getattr(part, field_name) in labels
        )


@dataclass(frozen=True)
class Normalization:
    """Opt-in DNA transformations; internal whitespace is never normalized."""

    uppercase: bool = False
    trim_outer_whitespace: bool = False

    def __post_init__(self) -> None:
        """Reject implicit truth-value coercion."""
        if not all(
            isinstance(v, bool) for v in (self.uppercase, self.trim_outer_whitespace)
        ):
            msg = "normalization flags must be booleans"
            raise TypeError(msg)


@dataclass(frozen=True)
class PartTable:
    """A declared table format, column mapping and normalization policy."""

    table: str | Path
    format: str
    columns: Mapping[str, str] = field(default_factory=dict)
    metadata_columns: Mapping[str, str] = field(default_factory=dict)
    id_policy: str = "provided"
    normalization: Normalization = field(default_factory=Normalization)
    sheet: str | None = None

    def __post_init__(self) -> None:
        """Freeze mappings and reject unsupported import policy."""
        object.__setattr__(self, "table", Path(self.table))
        if self.format not in {"csv", "tsv", "parquet", "xlsx"}:
            msg = "format must be csv, tsv, parquet or xlsx"
            raise ValueError(msg)
        if self.sheet is not None:
            required_text(self.sheet, field_name="sheet")
            if self.format != "xlsx":
                msg = "sheet is only supported for xlsx inputs"
                raise ValueError(msg)
        if self.id_policy not in {"provided", "row"}:
            msg = "id_policy must be provided or row"
            raise ValueError(msg)
        if not isinstance(self.normalization, Normalization):
            msg = "normalization must be Normalization"
            raise TypeError(msg)
        fields = set(Part.__dataclass_fields__) - {"metadata"}
        if unknown := set(self.columns) - fields:
            msg = f"unknown columns fields: {sorted(unknown)}"
            raise ValueError(msg)
        for name in ("columns", "metadata_columns"):
            mapping = immutable_json_mapping(getattr(self, name))
            for key, value in mapping.items():
                required_text(key, field_name=name)
                required_text(value, field_name=f"{name}.{key}")
            object.__setattr__(self, name, mapping)
        if set(self.metadata_columns) & fields:
            msg = "metadata_columns cannot override part fields"
            raise ValueError(msg)
