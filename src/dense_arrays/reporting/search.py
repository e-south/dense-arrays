"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/search.py

Bind recorded search methods and finite proposal histories to saved plans.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.search import HeuristicEvidence

if TYPE_CHECKING:
    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Attempt
    from dense_arrays.planning import GenerationPlan


class SearchHistory:
    """Retain one consumed-proposal marker per greedy cell under the reader cap."""

    def __init__(self, plans: dict[str, GenerationPlan], budget: ReadBudget) -> None:
        """Reserve bounded state before storing cell/batch identities."""
        budget.retain(sum(p.request.search == "greedy" for p in plans.values()))
        self._proposals: dict[str, tuple[object, object]] = {}

    def observe(self, attempt: Attempt, plan: GenerationPlan) -> None:
        """Reject method substitution and repeated proposals within one batch."""
        heuristic = attempt.evidence.get("heuristic")
        if plan.request.search == "exact":
            if heuristic is not None:
                msg = "recorded search method differs from exact plan"
                raise ValueError(msg)
            return
        if "solver_status" in attempt.evidence or (
            heuristic is None
            and attempt.outcome
            not in {"in_progress", "error", "interrupted_unresolved"}
        ):
            msg = "recorded search method differs from greedy plan"
            raise ValueError(msg)
        if heuristic is None:
            return
        observed = HeuristicEvidence.from_dict(heuristic)
        if observed.status == "candidate":
            context = (
                attempt.evidence.get("batch_index"),
                attempt.evidence.get("batch_id"),
            )
            if self._proposals.get(attempt.cell_id) == context:
                msg = "only one greedy proposal is supported within each offered batch"
                raise ValueError(msg)
            self._proposals[attempt.cell_id] = context
