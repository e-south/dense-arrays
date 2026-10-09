"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/selections/requests.py

Explicit library allocations with versioned, declarative wire contracts.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

from dense_arrays._record_validation import integer, object_fields, required_text
from dense_arrays.reporting.design_filters import DesignFilter

SELECTION_SCHEMA = "dense_arrays.library-selection.v1"


@dataclass(frozen=True)
class Take:
    """Choose a total count or explicit cell quotas, without replacement."""

    count: int | None = None
    per_cell: Mapping[str, int] | None = None
    policy: str = "first"
    seed: int | None = None
    shortfall: str = "error"

    def __post_init__(self) -> None:
        """Reject ambiguous allocations and ignored randomness options."""
        if (self.count is None) == (self.per_cell is None):
            msg = "take requires exactly one of count or per_cell"
            raise ValueError(msg)
        if self.count is not None:
            integer(self.count, field_name="count", minimum=0)
        if self.per_cell is not None:
            if not isinstance(self.per_cell, Mapping):
                msg = (
                    "per_cell must be an explicit mapping of cell references to counts"
                )
                raise TypeError(msg)
            for cell, count in self.per_cell.items():
                required_text(cell, field_name="cell reference")
                integer(count, field_name=f"per_cell[{cell}]", minimum=0)
            object.__setattr__(self, "per_cell", MappingProxyType(dict(self.per_cell)))
        if self.policy not in {"first", "random"}:
            msg = "take policy must be first or random"
            raise ValueError(msg)
        if self.policy == "random":
            integer(self.seed, field_name="random selection seed", minimum=0)
        elif self.seed is not None:
            msg = "seed applies only to random selection"
            raise ValueError(msg)
        if self.shortfall not in {"error", "allow_partial"}:
            msg = "shortfall must be error or allow_partial"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Keep zero allocations explicit; unlisted cells have no quota."""
        return {
            **(
                {"count": self.count}
                if self.count is not None
                else {"per_cell": dict(self.per_cell)}
            ),
            "policy": self.policy,
            **({"seed": self.seed} if self.seed is not None else {}),
            "shortfall": self.shortfall,
        }

    @classmethod
    def from_dict(cls, value: object) -> Take:
        """Read a strict allocation without coercing counts or policy labels."""
        return cls(
            **object_fields(
                value, {"count", "per_cell", "policy", "seed", "shortfall"}, "take"
            )
        )


@dataclass(frozen=True)
class LibrarySelection:
    """Filter eligibility, then allocate from the entire eligible population."""

    filter: DesignFilter = field(default_factory=DesignFilter)
    take: Take | None = None

    def __post_init__(self) -> None:
        """Keep predicates and allocation policies independently typed."""
        if not isinstance(self.filter, DesignFilter):
            msg = "library selection requires DesignFilter"
            raise TypeError(msg)
        if self.take is not None and not isinstance(self.take, Take):
            msg = "library selection take must be Take"
            raise TypeError(msg)

    def to_dict(self) -> dict[str, object]:
        """Encode the nested predicate under the enclosing schema declaration."""
        predicate = self.filter.to_dict()
        predicate.pop("schema")
        return {
            "schema": SELECTION_SCHEMA,
            "filter": predicate,
            **({"take": self.take.to_dict()} if self.take is not None else {}),
        }

    @classmethod
    def from_dict(cls, value: object) -> LibrarySelection:
        """Reject unknown fields and nested declarations instead of overriding them."""
        data = object_fields(value, {"schema", "filter", "take"}, "library selection")
        if data.pop("schema", None) != SELECTION_SCHEMA:
            msg = "unsupported library-selection schema"
            raise ValueError(msg)
        predicate = object_fields(
            data.get("filter", {}),
            {"design_ids", "cells", "part_ids", "groups", "metrics"},
            "selection filter",
        )
        return cls(
            DesignFilter.from_dict(
                {"schema": "dense_arrays.design-filter.v1", **predicate}
            ),
            Take.from_dict(data["take"]) if "take" in data else None,
        )


def resolve_quotas(take: Take | None, cells: Iterable[str]) -> dict[str, int | None]:
    """Resolve every declared cell, including cells allocated zero designs."""
    if take is None or take.per_cell is None:
        return {"total": None if take is None else take.count}
    cells = set(cells)
    result = {}
    for label, count in take.per_cell.items():
        matches = {c for c in cells if label in {c, c.split("/")[1]}}
        if len(matches) != 1:
            kind = "unknown" if not matches else "ambiguous"
            msg = f"{kind} cell {label!r}; use a full run/cell reference"
            raise ValueError(msg)
        cell = next(iter(matches))
        if cell in result:
            msg = f"multiple quota labels resolve to the same cell {cell!r}"
            raise ValueError(msg)
        result[cell] = count
    return result
