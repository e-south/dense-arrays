"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/matrices/allocation.py

Deterministic targets with explicit inactive cells and no redistribution.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

from dense_arrays._record_validation import integer, required_text


@dataclass(frozen=True)
class Allocation:
    """Choose one per-cell count, complete named counts, or an allocated total."""

    per_cell: int | None = None
    counts: Mapping[str, int] | None = None
    total: int | None = None
    policy: str | None = None
    zero_cells: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Reject competing declarations and freeze named count records."""
        if sum(v is not None for v in (self.per_cell, self.counts, self.total)) != 1:
            msg = "allocation requires exactly one of per_cell, counts or total"
            raise ValueError(msg)
        for name in ("per_cell", "total"):
            if getattr(self, name) is not None:
                integer(getattr(self, name), field_name=f"allocation.{name}", minimum=0)
        if self.counts is not None:
            object.__setattr__(self, "counts", _counts(self.counts))
        if not isinstance(self.zero_cells, (tuple, list)) or any(
            not isinstance(v, str) or not v for v in self.zero_cells
        ):
            msg = "allocation.zero_cells must contain cell identities"
            raise TypeError(msg)
        if len(set(self.zero_cells)) != len(self.zero_cells):
            msg = "allocation.zero_cells repeats a cell"
            raise ValueError(msg)
        object.__setattr__(self, "zero_cells", tuple(self.zero_cells))
        if self.total is not None:
            if self.policy != "balanced":
                msg = "total allocation requires policy='balanced'"
                raise ValueError(msg)
        elif self.policy is not None or self.zero_cells:
            msg = "policy and zero_cells apply only to total allocation"
            raise ValueError(msg)

    @property
    def policy_version(self) -> str:
        """Name the explicit allocation rule used by the resolved plan."""
        return "ordered_balanced.v1" if self.total is not None else "explicit_counts.v1"

    def resolve(self, cells: tuple[str, ...]) -> tuple[int, ...]:
        """Allocate over the declared expansion order, retaining every zero target."""
        if not cells or len(set(cells)) != len(cells):
            msg = "allocation requires nonempty unique cell identities"
            raise ValueError(msg)
        if self.per_cell is not None:
            return (self.per_cell,) * len(cells)
        if self.counts is not None:
            if set(self.counts) != set(cells):
                msg = "allocation.counts must name every expanded cell exactly once"
                raise ValueError(msg)
            return tuple(self.counts[cell] for cell in cells)
        if set(self.zero_cells) - set(cells):
            msg = "allocation.zero_cells names unknown cells"
            raise ValueError(msg)
        active = len(cells) - len(self.zero_cells)
        if self.total < active or (active == 0 and self.total != 0):
            msg = (
                "total cannot cover active cells; declare zero-target cells explicitly"
            )
            raise ValueError(msg)
        quotient, remainder = divmod(self.total, active) if active else (0, 0)
        result = []
        for cell in cells:
            if cell in self.zero_cells:
                result.append(0)
            else:
                result.append(quotient + (remainder > 0))
                remainder = max(0, remainder - 1)
        return tuple(result)

    def to_dict(self) -> dict[str, object]:
        """Encode all choices without applying allocation to an unknown domain."""
        return {
            "per_cell": self.per_cell,
            "counts": None if self.counts is None else dict(self.counts),
            "total": self.total,
            "policy": self.policy,
            "zero_cells": list(self.zero_cells),
        }


def _counts(value: object) -> Mapping[str, int]:
    if not isinstance(value, Mapping) or not value:
        msg = "allocation.counts must be a nonempty cell mapping"
        raise TypeError(msg)
    for name, count in value.items():
        required_text(name, field_name="allocation cell")
        integer(count, field_name=f"allocation.counts.{name}", minimum=0)
    return MappingProxyType(dict(value))
