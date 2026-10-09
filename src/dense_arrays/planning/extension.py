"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/extension.py

Explicit additional targets and immutable parent-library exclusions.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, object_fields, required_text
from dense_arrays.planning.lineage import ParentRun
from dense_arrays.planning.models import Limits

if TYPE_CHECKING:
    from pathlib import Path

EXTENSION_SCHEMA = "dense_arrays.extension.v1"


@dataclass(frozen=True)
class ExtensionSpec:
    """Generate additional designs under unchanged rules and new explicit effort."""

    parent: ParentRun
    additional: int | Mapping[str, int]
    limits: Limits
    seed: int

    def __post_init__(self) -> None:
        """Require typed lineage, a positive additional target, effort and seed."""
        if not isinstance(self.parent, ParentRun) or not isinstance(
            self.limits, Limits
        ):
            msg = "extension requires ParentRun and Limits"
            raise TypeError(msg)
        if isinstance(self.additional, Mapping):
            for cell, count in self.additional.items():
                required_text(cell, field_name="additional cell")
                integer(count, field_name=f"additional.{cell}", minimum=0)
            if not self.additional or not sum(self.additional.values()):
                msg = "additional cell targets must include a positive count"
                raise ValueError(msg)
            object.__setattr__(
                self, "additional", MappingProxyType(dict(self.additional))
            )
        else:
            integer(self.additional, field_name="additional", minimum=1)
        integer(self.seed, field_name="seed", minimum=0)

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Export explicit extension intent without resolving its parent run."""
        return {
            "schema": EXTENSION_SCHEMA,
            "parent": {
                "run": str(self.parent.run)
                if base is None
                else os.path.relpath(self.parent.run.absolute(), base)
            },
            "additional": dict(self.additional)
            if isinstance(self.additional, Mapping)
            else self.additional,
            "limits": self.limits.to_dict(),
            "seed": self.seed,
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> ExtensionSpec:
        """Read a strict extension request with no implicit seed or effort reset."""
        data = object_fields(
            value, {"schema", "parent", "additional", "limits", "seed"}, "extension"
        )
        if data.pop("schema", None) != EXTENSION_SCHEMA:
            msg = "unsupported extension schema"
            raise ValueError(msg)
        if missing := {"parent", "additional", "limits", "seed"} - set(data):
            msg = f"extension requires explicit fields: {sorted(missing)}"
            raise ValueError(msg)
        parent = ParentRun(**object_fields(data.pop("parent"), {"run"}, "parent"))
        if base is not None and not parent.run.is_absolute():
            parent = ParentRun(base / parent.run)
        limits = Limits(
            **object_fields(
                data.pop("limits"),
                {"attempts", "active_seconds", "solver_seconds", "model_pairs"},
                "limits",
            )
        )
        return cls(parent=parent, limits=limits, **data)
