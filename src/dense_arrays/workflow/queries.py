"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/queries.py

Translate CLI query options into the same typed Python predicates.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from pathlib import Path

from dense_arrays.parts import PartFilter
from dense_arrays.reporting import (
    AttemptFilter,
    CandidateFilter,
    DesignFilter,
    PlanFilter,
)
from dense_arrays.workflow.inputs import read_selection


def query_filter(  # noqa: PLR0913 - one shared set of CLI selection options
    view: str,
    *,
    selection: Path | None = None,
    design_id: list[str] | None = None,
    cell: list[str] | None = None,
    part_id: list[str] | None = None,
    group: list[str] | None = None,
    attempt_id: list[int] | None = None,
    outcome: list[str] | None = None,
    plan_id: list[str] | None = None,
    candidate_index: list[int] | None = None,
    recipe_id: list[str] | None = None,
    reason: list[str] | None = None,
) -> CandidateFilter | PartFilter | AttemptFilter | DesignFilter | PlanFilter | None:
    """Reject mixed predicates and file/flag precedence instead of guessing."""
    flags = any(
        value is not None
        for value in (
            design_id,
            cell,
            part_id,
            group,
            attempt_id,
            outcome,
            plan_id,
            candidate_index,
            recipe_id,
            reason,
        )
    )
    if selection is not None:
        if flags:
            msg = "--selection is exclusive with filter flags"
            raise ValueError(msg)
        return read_selection(selection)
    if not flags:
        return None
    if (
        view == "candidates"
        or candidate_index is not None
        or reason is not None
        or recipe_id is not None
    ):
        return _candidate_filter(
            view,
            candidate_index,
            outcome,
            reason,
            recipes=recipe_id,
            incompatible=(design_id, cell, part_id, group, attempt_id, plan_id),
        )
    return _record_filter(
        view,
        design_id=design_id,
        cell=cell,
        part_id=part_id,
        group=group,
        attempt_id=attempt_id,
        outcome=outcome,
        plan_id=plan_id,
    )


def _record_filter(  # noqa: PLR0913 - resolved CLI selector groups
    view: str,
    *,
    design_id: list[str] | None,
    cell: list[str] | None,
    part_id: list[str] | None,
    group: list[str] | None,
    attempt_id: list[int] | None,
    outcome: list[str] | None,
    plan_id: list[str] | None,
) -> PartFilter | AttemptFilter | DesignFilter | PlanFilter:
    """Translate remaining selectors for their existing record population."""
    if plan_id is not None:
        if view not in {"plan", "plans", "request"} or any(
            v is not None
            for v in (design_id, cell, part_id, group, attempt_id, outcome)
        ):
            msg = (
                "--plan-id applies only to plan/plans/request "
                "and cannot mix other filters"
            )
            raise ValueError(msg)
        return PlanFilter(plan_ids=tuple(plan_id))
    if (
        view in {"design", "designs", "sequences", "placements", "quality", "selection"}
        and attempt_id is None
        and outcome is None
    ):
        return DesignFilter(
            design_ids=tuple(design_id or ()),
            cells=tuple(cell or ()),
            part_ids=tuple(part_id or ()),
            groups=tuple(group or ()),
        )
    if (
        view in {"attempts", "diagnostics"}
        and design_id is None
        and part_id is None
        and group is None
    ):
        return AttemptFilter(
            attempt_ids=tuple(attempt_id or ()),
            cells=tuple(cell or ()),
            outcomes=tuple(outcome or ()),
        )
    if view == "parts" and all(
        value is None for value in (design_id, cell, attempt_id, outcome)
    ):
        return PartFilter(part_ids=tuple(part_id or ()), groups=tuple(group or ()))
    msg = f"filter flags do not apply to view {view!r}"
    raise ValueError(msg)


def _candidate_filter(  # noqa: PLR0913 - resolved candidate selector groups
    view: str,
    indices: list[int] | None,
    outcomes: list[str] | None,
    reasons: list[str] | None,
    *,
    recipes: list[str] | None,
    incompatible: tuple[object, ...],
) -> CandidateFilter:
    """Keep preparation flags scoped to their declared candidate population."""
    if view != "candidates" or any(v is not None for v in incompatible):
        msg = (
            "candidate queries require view='candidates' and only "
            "indices, outcomes, reasons or recipe IDs"
        )
        raise ValueError(msg)
    return CandidateFilter(
        indices=tuple(indices or ()),
        outcomes=tuple(outcomes or ()),
        reasons=tuple(reasons or ()),
        recipes=tuple(recipes or ()),
    )
