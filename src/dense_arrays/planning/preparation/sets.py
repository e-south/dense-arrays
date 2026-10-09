"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/preparation/sets.py

Resolve independent sampled recipes without drawing or sharing quotas.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import mutable_json, object_fields, semantic_digest
from dense_arrays.parts import PreparationSet

from .requests import preparation_from_dict, preparation_to_dict
from .sampled import SampledPreparation, resolve_sampled

if TYPE_CHECKING:
    from pathlib import Path

SET_PLAN_SCHEMA = "dense_arrays.preparation_plan.v3"
SET_POLICY = "independent_recipes.v1"


@dataclass(frozen=True, repr=False)
class SetPreparation:
    """Ordered recipes with independently resolved models, policies and budgets."""

    recipes: Mapping[str, SampledPreparation]
    sequence_collisions: str = "error"
    core_collisions: str = field(default="preserve", kw_only=True)

    def __post_init__(self) -> None:
        """Require sampled evidence and validate names through the request owner."""
        if not isinstance(self.recipes, Mapping) or any(
            not isinstance(item, SampledPreparation) for item in self.recipes.values()
        ):
            msg = "resolved preparation set requires sampled recipes"
            raise TypeError(msg)
        PreparationSet(
            {name: item.request for name, item in self.recipes.items()},
            self.sequence_collisions,
            core_collisions=self.core_collisions,
        )
        object.__setattr__(self, "recipes", MappingProxyType(dict(self.recipes)))

    @property
    def request(self) -> PreparationSet:
        """Expose the effective requests without reopening their inputs."""
        return PreparationSet(
            {name: item.request for name, item in self.recipes.items()},
            self.sequence_collisions,
            core_collisions=self.core_collisions,
        )

    @property
    def identity_count(self) -> int:
        """Charge every embedded model and recipe identity against read limits."""
        return sum(1 + item.identity_count for item in self.recipes.values())

    @property
    def plan_id(self) -> str:
        """Bind recipe names, order and complete independent plan identities."""
        return semantic_digest(
            {
                "schema": SET_PLAN_SCHEMA,
                "policy": SET_POLICY,
                "sequence_collisions": self.sequence_collisions,
                **(
                    {"core_collisions": self.core_collisions}
                    if self.core_collisions != "preserve"
                    else {}
                ),
                "recipes": [
                    {"id": name, "plan_id": item.plan_id}
                    for name, item in self.recipes.items()
                ],
            }
        )

    @property
    def preview(self) -> MappingProxyType:
        """Report per-recipe limits and totals without predicting retained yield."""
        previews = [(name, item.preview) for name, item in self.recipes.items()]
        return MappingProxyType(
            {
                **{
                    key: sum(item[key] for _, item in previews)
                    for key in (
                        "source_parts",
                        "source_motifs",
                        "candidate_budget",
                        "candidate_bases_bound",
                        "requested_retention",
                    )
                },
                "retained_parts": None,
                "sequence_collisions": self.sequence_collisions,
                **(
                    {"core_collisions": self.core_collisions}
                    if self.core_collisions != "preserve"
                    else {}
                ),
                "retained_count_status": "unknown",
                "required_tools": sorted(
                    {tool for _, item in previews for tool in item["required_tools"]}
                ),
                "recipes": tuple({"id": name, **dict(item)} for name, item in previews),
            }
        )

    def verify_inputs(self) -> None:
        """Verify all bound inputs before any recipe starts."""
        for item in self.recipes.values():
            item.verify_inputs()

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Serialize complete child plans and reproducible publication order."""
        return {
            "schema": SET_PLAN_SCHEMA,
            "plan_id": self.plan_id,
            "policy": SET_POLICY,
            "request": preparation_to_dict(self.request, base=base),
            "recipes": [
                {"id": name, "plan": item.to_dict(base=base)}
                for name, item in self.recipes.items()
            ],
            "preview": mutable_json(dict(self.preview)),
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> SetPreparation:
        """Read saved models without discovering inputs or invoking scorers."""
        fields = {"schema", "plan_id", "policy", "request", "recipes", "preview"}
        data = object_fields(value, fields, "preparation set plan")
        if (
            set(data) != fields
            or data["schema"] != SET_PLAN_SCHEMA
            or not isinstance(data["recipes"], list)
        ):
            msg = "unsupported or incomplete preparation set plan"
            raise ValueError(msg)
        recipes = {}
        for raw in data["recipes"]:
            item = object_fields(raw, {"id", "plan"}, "resolved recipe")
            if set(item) != {"id", "plan"} or item["id"] in recipes:
                msg = "resolved recipes require unique IDs and plans"
                raise ValueError(msg)
            recipes[item["id"]] = SampledPreparation.from_dict(item["plan"], base=base)
        request = preparation_from_dict(data["request"], base=base)
        if not isinstance(request, PreparationSet):
            msg = "resolved set requires a preparation set request"
            raise TypeError(msg)
        result = cls(
            recipes,
            request.sequence_collisions,
            core_collisions=request.core_collisions,
        )
        if result.to_dict(base=base) != data:
            msg = "preparation set identity or resolved fields disagree"
            raise ValueError(msg)
        return result


def resolve_set(request: PreparationSet) -> SetPreparation:
    """Preflight every recipe without executing candidate work."""
    return SetPreparation(
        {name: resolve_sampled(item) for name, item in request.recipes.items()},
        request.sequence_collisions,
        core_collisions=request.core_collisions,
    )
