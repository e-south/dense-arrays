"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/objectives.py

Recount packing preferences from committed proposal geometry, without solving.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.objectives import PackingObjective
from dense_arrays.planning.batches.bindings import offered_batch

if TYPE_CHECKING:
    from dense_arrays.artifacts.batches import BatchDecision
    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Attempt
    from dense_arrays.planning import GenerationPlan


class ObjectiveHistory:
    """Retain at most one offered-batch counter per cell, under the reader budget."""

    def __init__(self, plans: dict[str, GenerationPlan], budget: ReadBudget) -> None:
        """Reserve each cell's maximum persistent usage state before allocation."""
        budget.retain(
            sum(
                len(p.request.parts)
                for p in plans.values()
                if p.request.packing_preference
            )
        )
        self._states: dict[str, tuple[tuple[object, object], PackingObjective]] = {}

    def observe(
        self, attempt: Attempt, plan: GenerationPlan, decision: BatchDecision | None
    ) -> None:
        """Check pre-solve weights, then count this attempt's recorded packed parts."""
        recorded = attempt.evidence.get("packing_objective")
        expected = (
            plan.request.packing_preference is not None
            and "solver_status" in attempt.evidence
        )
        if not expected:
            if recorded is not None:
                msg = "packing objective is not applicable to this attempt"
                raise ValueError(msg)
            return
        batch = (
            decision.batch
            if decision is not None
            else offered_batch(
                plan.request,
                attempt.evidence.get("batch_id"),
                index=attempt.evidence.get("batch_index"),
            )
        )
        context = (
            attempt.evidence.get("batch_index"),
            attempt.evidence.get("batch_id"),
        )
        prior = self._states.get(attempt.cell_id)
        state = (
            prior[1]
            if prior is not None and prior[0] == context
            else PackingObjective(
                dict.fromkeys(
                    batch.part_ids
                    if batch is not None
                    else (p.part_id for p in plan.request.parts),
                    0,
                ),
                0,
            )
        )
        if recorded != state.to_dict():
            msg = "packing objective disagrees with prior proposed packings"
            raise ValueError(msg)
        if attempt.candidate is not None:
            usage = dict(state.usage)
            for placement in attempt.candidate.packed.placements:
                usage[placement.feature_id] += 1
            state = PackingObjective(usage, state.proposed_packings + 1)
        self._states[attempt.cell_id] = (context, state)
