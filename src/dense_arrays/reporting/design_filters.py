"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/design_filters.py

Typed selection of accepted designs by identity, placements and composition.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import object_fields, records, required_text
from dense_arrays.parts.filters import Range
from dense_arrays.reporting.metrics import design_metrics

if TYPE_CHECKING:
    from collections.abc import Mapping

    from dense_arrays.artifacts.records import Design
    from dense_arrays.parts.models import Part

METRICS = {"length", "gc_fraction", "placement_count", "packing_density"}


@dataclass(frozen=True)
class DesignFilter:
    """OR within identity fields, AND across fields and inclusive metric ranges."""

    design_ids: tuple[str, ...] = ()
    cells: tuple[str, ...] = ()
    part_ids: tuple[str, ...] = ()
    groups: tuple[str, ...] = ()
    metrics: Mapping[str, Range] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Freeze explicitly typed selectors without interpreting expressions."""
        for name in ("design_ids", "cells", "part_ids", "groups"):
            values = records(getattr(self, name), str, field_name=name)
            for value in values:
                required_text(value, field_name=name)
            if len(set(values)) != len(values):
                msg = f"{name} must not repeat labels"
                raise ValueError(msg)
            object.__setattr__(self, name, values)
        metrics = dict(self.metrics)
        if unknown := set(metrics) - METRICS:
            msg = f"unavailable design metrics: {sorted(unknown)}"
            raise ValueError(msg)
        if any(not isinstance(value, Range) for value in metrics.values()):
            msg = "metric predicates must be Range values"
            raise TypeError(msg)
        object.__setattr__(self, "metrics", MappingProxyType(metrics))

    @property
    def identities(self) -> int:
        """Number of explicit identity selectors retained by this predicate."""
        return sum(
            len(getattr(self, name))
            for name in ("design_ids", "cells", "part_ids", "groups")
        )

    def matches(
        self, design: Design, parts: Mapping[str, Part], collection_id: str = ""
    ) -> bool:
        """Use selected placements, never incidental sequence substring matches."""
        if self.design_ids and not {design.design_id, design.reference}.intersection(
            self.design_ids
        ):
            return False
        if self.cells and not {
            design.cell_id,
            f"{design.run_id}/{design.cell_id}",
        }.intersection(self.cells):
            return False
        selected = {p.feature_id for p in design.realized.placements}
        aliases = selected | {f"{collection_id}/{p}" for p in selected}
        if self.part_ids and not aliases.intersection(self.part_ids):
            return False
        if self.groups and not {parts[p].group for p in selected}.intersection(
            self.groups
        ):
            return False
        if not self.metrics:
            return True
        values = design_metrics(design.realized)
        values["packing_density"] = values["density"]
        return all(
            interval.matches(values[name]) for name, interval in self.metrics.items()
        )

    def to_dict(self) -> dict[str, object]:
        """Encode the same declarative predicate used by both interfaces."""
        return {
            "schema": "dense_arrays.design-filter.v1",
            **{
                name: list(getattr(self, name))
                for name in ("design_ids", "cells", "part_ids", "groups")
            },
            "metrics": {
                name: {"min": r.min, "max": r.max} for name, r in self.metrics.items()
            },
        }

    @classmethod
    def from_dict(cls, value: object) -> DesignFilter:
        """Reject unknown fields, versions and untyped metric bounds."""
        data = object_fields(
            value,
            {"schema", "design_ids", "cells", "part_ids", "groups", "metrics"},
            "design filter",
        )
        if data.pop("schema", None) != "dense_arrays.design-filter.v1":
            msg = "unsupported design-filter schema"
            raise ValueError(msg)
        data["metrics"] = {
            name: Range(**object_fields(r, {"min", "max"}, "range"))
            for name, r in data.get("metrics", {}).items()
        }
        return cls(**data)
