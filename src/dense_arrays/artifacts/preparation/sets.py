"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/sets.py

Per-recipe accounting and candidate coordinates in one prepared pool.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, replace
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import object_fields, required_text
from dense_arrays.artifacts.preparation.records import COUNTS, PoolAccounting

if TYPE_CHECKING:
    from dense_arrays.parts.candidates import Candidate

SET_ACCOUNTING_SCHEMA = "dense_arrays.pool_set_accounting.v1"


@dataclass(frozen=True)
class SetAccounting:
    """Independent targets, stopping conditions and reconciled aggregate counts."""

    recipes: Mapping[str, PoolAccounting]

    def __post_init__(self) -> None:
        """Require named non-nested accounting for every executed recipe."""
        if not isinstance(self.recipes, Mapping) or not self.recipes:
            msg = "set accounting requires a nonempty recipe mapping"
            raise ValueError(msg)
        for name, item in self.recipes.items():
            required_text(name, field_name="recipe ID")
            if not isinstance(item, PoolAccounting):
                msg = "each recipe requires independent PoolAccounting"
                raise TypeError(msg)
        object.__setattr__(self, "recipes", MappingProxyType(dict(self.recipes)))

    @property
    def counts(self) -> Mapping[str, int]:
        """Sum stages after each recipe has independently reconciled its counts."""
        return MappingProxyType(
            {
                key: sum(item.counts[key] for item in self.recipes.values())
                for key in COUNTS
            }
        )

    @property
    def requested_retention(self) -> int:
        """Sum targets without allowing one recipe to satisfy another."""
        return sum(item.requested_retention for item in self.recipes.values())

    @property
    def candidate_budget(self) -> int:
        """Sum independently enforced candidate limits."""
        return sum(item.candidate_budget for item in self.recipes.values())

    @property
    def rejections(self) -> Mapping[str, int]:
        """Aggregate reason occurrences while retaining their per-recipe context."""
        result = Counter()
        for item in self.recipes.values():
            result.update(item.rejections)
        return MappingProxyType(dict(result))

    @property
    def state(self) -> str:
        """Complete only when every independent target was met without errors."""
        return (
            "completed"
            if all(item.state == "completed" for item in self.recipes.values())
            else "incomplete"
        )

    @property
    def stop_reason(self) -> str:
        """Direct the reader to individual recipe stopping conditions."""
        return "recipe_limits"

    @property
    def retention(self) -> None:
        """There is no global retention ranking across scoring models."""
        return None

    def to_dict(self) -> dict[str, object]:
        """Publish per-recipe evidence with independently reproducible totals."""
        return {
            "schema": SET_ACCOUNTING_SCHEMA,
            "counts": dict(self.counts),
            "requested_retention": self.requested_retention,
            "candidate_budget": self.candidate_budget,
            "stop_reason": self.stop_reason,
            "rejections": dict(self.rejections),
            "recipes": [
                {"id": name, "accounting": item.to_dict()}
                for name, item in self.recipes.items()
            ],
        }

    @classmethod
    def from_dict(cls, value: object) -> SetAccounting:
        """Reject inconsistent totals, duplicated names or nested set accounting."""
        fields = {
            "schema",
            "counts",
            "requested_retention",
            "candidate_budget",
            "stop_reason",
            "rejections",
            "recipes",
        }
        data = object_fields(value, fields, "set accounting")
        if (
            set(data) != fields
            or data["schema"] != SET_ACCOUNTING_SCHEMA
            or not isinstance(data["recipes"], list)
        ):
            msg = "unsupported or incomplete set accounting"
            raise ValueError(msg)
        recipes = {}
        for raw in data["recipes"]:
            item = object_fields(raw, {"id", "accounting"}, "recipe accounting")
            if set(item) != {"id", "accounting"} or item["id"] in recipes:
                msg = "recipe accounting requires unique IDs and complete records"
                raise ValueError(msg)
            recipes[item["id"]] = PoolAccounting.from_dict(item["accounting"])
        result = cls(recipes)
        if result.to_dict() != data:
            msg = "aggregate accounting disagrees with recipe counts"
            raise ValueError(msg)
        return result


def read_accounting(value: object) -> PoolAccounting | SetAccounting:
    """Decode exactly the accounting family declared by saved evidence."""
    owner = (
        SetAccounting
        if isinstance(value, dict) and value.get("schema") == SET_ACCOUNTING_SCHEMA
        else PoolAccounting
    )
    return owner.from_dict(value)


def qualify_candidate(candidate: Candidate, recipe_id: str, offset: int) -> Candidate:
    """Give local evidence a global row and stable recipe-qualified part ID."""
    return replace(
        candidate,
        index=candidate.index + offset,
        part=replace(candidate.part, part_id=f"{recipe_id}/{candidate.part.part_id}"),
        representative=None
        if candidate.representative is None
        else candidate.representative + offset,
        recipe_id=recipe_id,
        recipe_index=candidate.index,
    )


def local_candidate(candidate: Candidate, recipe_id: str, offset: int) -> Candidate:
    """Recover local coordinates only after checking the exact recipe binding."""
    index = candidate.index - offset
    if (
        candidate.recipe_id != recipe_id
        or candidate.recipe_index != index
        or candidate.part.part_id != f"{recipe_id}/candidate_{index}"
    ):
        msg = "candidate recipe origin disagrees with its recorded position"
        raise ValueError(msg)
    return replace(
        candidate,
        index=index,
        part=replace(candidate.part, part_id=f"candidate_{index}"),
        representative=None
        if candidate.representative is None
        else candidate.representative - offset,
        recipe_id=None,
        recipe_index=None,
    )
