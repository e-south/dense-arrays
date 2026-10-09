"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/filters.py

Typed predicates over supplied part identity and available numeric evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import dataclass, field
from numbers import Real
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import object_fields, records, required_text

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from dense_arrays.parts.models import Part


@dataclass(frozen=True)
class Range:
    """Inclusive numeric bounds; missing observations never match."""

    min: float | None = None
    max: float | None = None

    def __post_init__(self) -> None:
        """Require finite nonboolean bounds without coercing values."""
        if self.min is None and self.max is None:
            msg = "range requires min or max"
            raise ValueError(msg)
        for value in (self.min, self.max):
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(value)
            ):
                msg = "range bounds must be finite numbers"
                raise ValueError(msg)
        if self.min is not None and self.max is not None and self.min > self.max:
            msg = "range minimum exceeds maximum"
            raise ValueError(msg)

    def matches(self, value: float | None) -> bool:
        """Apply inclusive bounds to available evidence only."""
        return (
            value is not None
            and (self.min is None or value >= self.min)
            and (self.max is None or value <= self.max)
        )


@dataclass(frozen=True)
class PartFilter:
    """OR within IDs/groups, AND across fields and metric ranges."""

    part_ids: tuple[str, ...] = ()
    groups: tuple[str, ...] = ()
    metrics: Mapping[str, Range] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Freeze predicates and reject duplicate or untyped selectors."""
        for name in ("part_ids", "groups"):
            values = records(getattr(self, name), str, field_name=name)
            for value in values:
                required_text(value, field_name=name)
            if len(set(values)) != len(values):
                msg = f"{name} must not repeat labels"
                raise ValueError(msg)
            object.__setattr__(self, name, values)
        values = dict(self.metrics)
        for key, value in values.items():
            required_text(key, field_name="metric")
            if not isinstance(value, Range):
                msg = "metric predicates must be Range values"
                raise TypeError(msg)
        object.__setattr__(self, "metrics", MappingProxyType(values))

    def validate(self, parts: Sequence[Part]) -> None:
        """Reject unknown IDs/groups and unavailable metrics before evaluating rows."""
        self.validate_available({p.part_id for p in parts}, {p.group for p in parts})

    def validate_available(self, part_ids: set[str], groups: set[str | None]) -> None:
        """Validate against explicit identities or equivalent indexed lookups."""
        for name, labels, available in (
            ("part IDs", self.part_ids, part_ids),
            ("groups", self.groups, groups),
        ):
            if missing := set(labels) - available:
                msg = f"unknown {name}: {sorted(missing)}"
                raise ValueError(msg)
        if missing := set(self.metrics) - {"length"}:
            msg = f"unavailable part metrics: {sorted(missing)}; supported: length"
            raise ValueError(msg)

    def matches(self, part: Part) -> bool:
        """Match one part after validating the filter against its collection."""
        return (
            (not self.part_ids or part.part_id in self.part_ids)
            and (not self.groups or part.group in self.groups)
            and all(
                interval.matches(len(part.sequence))
                for interval in self.metrics.values()
            )
        )

    def to_dict(self) -> dict[str, object]:
        """Serialize explicit predicates under the part-filter schema."""
        return {
            "schema": "dense_arrays.part-filter.v1",
            "part_ids": list(self.part_ids),
            "groups": list(self.groups),
            "metrics": {
                name: {"min": r.min, "max": r.max} for name, r in self.metrics.items()
            },
        }

    @classmethod
    def from_dict(cls, value: object, *, declared: bool = True) -> PartFilter:
        """Parse a declared file filter or an embedded request predicate."""
        data = object_fields(
            value, {"schema", "part_ids", "groups", "metrics"}, "part filter"
        )
        if declared and data.pop("schema", None) != "dense_arrays.part-filter.v1":
            msg = "unsupported part-filter schema"
            raise ValueError(msg)
        if not declared and "schema" in data:
            msg = "embedded part filters do not have a schema field"
            raise ValueError(msg)
        data["metrics"] = {
            name: Range(**object_fields(r, {"min", "max"}, "range"))
            for name, r in data.get("metrics", {}).items()
        }
        return cls(**data)
