"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/pools/diversity.py

Recorded sequential MMR distances, scoped to each preparation recipe.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays._record_validation import digest, integer, object_fields
from dense_arrays.artifacts.preparation.sets import SetAccounting
from dense_arrays.parts.scoring.configuration import finite
from dense_arrays.planning.preparation.sets import SetPreparation

if TYPE_CHECKING:
    from collections.abc import Mapping

    from dense_arrays.artifacts.preparation.reading import VerifiedPreparation
    from dense_arrays.artifacts.preparation.records import PoolAccounting
    from dense_arrays.artifacts.reading import ReadBudget

_CONTRACT = {
    "schema": "dense_arrays.mmr_diversity.v1",
    "policy": "greedy_mmr.v1",
    "distance": "pwm_tolerant_hamming",
}


def summarize_diversity(evidence: VerifiedPreparation) -> list[dict[str, object]]:
    """Copy verified choice distances without recomputing a diversity metric."""
    source = evidence.plan.resolved
    recipes = (
        source.recipes.items()
        if isinstance(source, SetPreparation)
        else ((None, source),)
    )
    reports = {}
    for name, recipe in recipes:
        if recipe.request.retain.mmr is None:
            continue
        evidence.budget.retain(9)
        reports[name] = {
            **_CONTRACT,
            "recipe_id": name,
            "model_id": recipe.source.motif.model_id,
            "scoring_id": recipe.source.scoring.binding_id,
            "choices": [],
        }
    for candidate in evidence.candidates:
        report = reports.get(candidate.recipe_id)
        if report is not None and candidate.retained:
            evidence.budget.retain(3)
            report["choices"].append(
                {
                    "rank": candidate.rank,
                    "nearest_distance": candidate.selection.nearest_distance,
                }
            )
    for report in reports.values():
        report["choices"].sort(key=lambda choice: choice["rank"])
        report["status"] = _status(len(report["choices"]))
    return list(reports.values())


def validate_diversity(
    report: Mapping[str, object],
    accounting: PoolAccounting | SetAccounting,
    budget: ReadBudget,
) -> None:
    """Validate complete recipe-local choice populations in a portable report."""
    if "diversity" not in report:
        return
    value = report["diversity"]
    if not isinstance(value, list) or not value:
        msg = "diversity requires a nonempty ordered list of MMR recipe records"
        raise TypeError(msg)
    recipes = (
        accounting.recipes.items()
        if isinstance(accounting, SetAccounting)
        else ((None, accounting),)
    )
    expected = [(name, item) for name, item in recipes if item.retention is not None]
    if len(value) != len(expected):
        msg = "diversity recipes disagree with recorded MMR policies"
        raise ValueError(msg)
    for raw, (name, item) in zip(value, expected, strict=True):
        _validate_recipe(raw, name, item, budget)


def _validate_recipe(
    raw: object, name: str | None, item: PoolAccounting, budget: ReadBudget
) -> None:
    fields = set(_CONTRACT) | {
        "recipe_id",
        "model_id",
        "scoring_id",
        "status",
        "choices",
    }
    budget.retain(len(fields))
    data = object_fields(raw, fields, "MMR diversity")
    if set(data) != fields or any(data[k] != v for k, v in _CONTRACT.items()):
        msg = "unsupported or incomplete MMR diversity record"
        raise ValueError(msg)
    if data["recipe_id"] != name:
        msg = "diversity recipe identity or order disagrees with its population"
        raise ValueError(msg)
    for key in ("model_id", "scoring_id"):
        digest(data[key], field_name=f"diversity.{key}")
    if (
        item.score_bands is not None
        and data["scoring_id"] != item.score_bands["scoring_id"]
    ):
        msg = "diversity scorer disagrees with the recipe score bands"
        raise ValueError(msg)
    choices = data["choices"]
    if not isinstance(choices, list) or len(choices) != item.counts["retained"]:
        msg = "diversity choice count disagrees with retained parts"
        raise ValueError(msg)
    budget.retain(3 * len(choices))
    if data["status"] != _status(len(choices)):
        msg = "diversity status disagrees with its recorded comparisons"
        raise ValueError(msg)
    for rank, choice in enumerate(choices, 1):
        _validate_choice(choice, rank)


def _validate_choice(value: object, rank: int) -> None:
    choice = object_fields(value, {"rank", "nearest_distance"}, "MMR diversity choice")
    if set(choice) != {"rank", "nearest_distance"}:
        msg = "incomplete diversity choice"
        raise ValueError(msg)
    integer(choice["rank"], field_name="diversity.rank", minimum=1)
    if choice["rank"] != rank:
        msg = "diversity ranks must be consecutive in selection order"
        raise ValueError(msg)
    distance = choice["nearest_distance"]
    if rank == 1:
        if distance is not None:
            msg = "first MMR choice has no earlier core comparison"
            raise ValueError(msg)
    elif distance is None or finite(distance, "diversity.nearest_distance") < 0:
        msg = "later MMR choices require a finite nonnegative distance"
        raise ValueError(msg)


def _status(count: int) -> str:
    return "empty" if count == 0 else "singleton" if count == 1 else "recorded"
