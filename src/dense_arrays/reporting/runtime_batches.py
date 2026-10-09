"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/runtime_batches.py

Verify recorded runtime selection context against the committed attempt prefix.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import closing
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays.planning import FeedbackSnapshot
from dense_arrays.reporting.readers import RecordView, read_records

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.batches import BatchDecision
    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Attempt
    from dense_arrays.planning import GenerationPlan
    from dense_arrays.reporting.summary import RunSummary


@dataclass
class RuntimeBatches:
    """Bound registry plus independently recounted per-cell selection observations."""

    decisions: dict[tuple[str, int], BatchDecision]
    by_id: dict[tuple[str, str], BatchDecision] = field(init=False)
    seen: set = field(default_factory=set)
    used: dict = field(default_factory=lambda: defaultdict(Counter))
    failed: dict = field(default_factory=lambda: defaultdict(Counter))
    enumerated: set = field(default_factory=set)

    def __post_init__(self) -> None:
        """Index membership for accepted designs without deriving new selections."""
        self.by_id = {(d.cell_id, d.batch.batch_id): d for d in self.decisions.values()}
        if len(self.by_id) != len(self.decisions):
            msg = "runtime batch identities are not unique within cells"
            raise ValueError(msg)

    def observe(self, attempt: Attempt, plan: GenerationPlan) -> BatchDecision | None:
        """Check each new snapshot before including the current attempt's outcome."""
        policy = plan.request.resampling
        if policy is None:
            return None
        key = (attempt.cell_id, attempt.evidence.get("batch_index"))
        decision = self.decisions.get(key)
        if decision is None or decision.batch.batch_id != attempt.evidence.get(
            "batch_id"
        ):
            msg = "attempt has no matching committed batch decision"
            raise ValueError(msg)
        if key not in self.seen:
            if decision.after_attempt != attempt.attempt_id - 1:
                msg = "batch decision does not match its observed attempt prefix"
                raise ValueError(msg)
            expected = (
                None
                if policy.feedback is None
                else FeedbackSnapshot(
                    policy.feedback,
                    self.used[attempt.cell_id],
                    self.failed[attempt.cell_id],
                )
            )
            if decision.batch.feedback != expected:
                msg = "batch feedback does not match committed observations"
                raise ValueError(msg)
            self.seen.add(key)
        self._count(attempt, decision)
        return decision

    def _count(self, attempt: Attempt, decision: BatchDecision) -> None:
        key = (decision.cell_id, decision.index)
        if decision.batch.feedback is not None:
            if attempt.outcome == "accepted":
                if attempt.candidate is None or attempt.candidate.final is None:
                    msg = "runtime accepted feedback requires final placement evidence"
                    raise ValueError(msg)
                self.used[attempt.cell_id].update(
                    p.feature_id for p in attempt.candidate.final.placements
                )
            elif attempt.outcome == "rejected" or (
                attempt.evidence.get("solver_status") == "infeasible"
                and key not in self.enumerated
            ):
                self.failed[attempt.cell_id].update(decision.batch.part_ids)
        if attempt.candidate is not None:
            self.enumerated.add(key)

    def design_plan(
        self, cell: str, batch_id: str | None, plan: GenerationPlan
    ) -> GenerationPlan:
        """Bind accepted geometry to its recorded runtime selection."""
        if plan.request.resampling is None:
            return plan
        decision = self.by_id.get((cell, batch_id))
        if decision is None:
            msg = "accepted design has no runtime batch decision"
            raise ValueError(msg)
        return decision.bind(plan)

    def finish(self) -> None:
        """Reject membership records with no corresponding first reservation."""
        if self.seen != self.decisions.keys():
            msg = "committed batch decisions have no corresponding attempts"
            raise ValueError(msg)


def read_runtime_batches(
    path: Path,
    summary: RunSummary,
    plans: dict[str, GenerationPlan],
    budget: ReadBudget,
) -> RuntimeBatches:
    """Load bounded decisions once, checking their origin and declared count."""
    decisions = {}
    expected = summary.batch_count or 0
    if summary.batch_count is None:
        if any(p.request.resampling is not None for p in plans.values()):
            msg = "runtime plan requires a batch registry manifest"
            raise ValueError(msg)
        return RuntimeBatches(decisions)
    query = RecordView(
        path, summary.revision, "batches", expected + 1, run_id=summary.run_id
    )
    with closing(read_records(query, budget)) as records:
        for decision in records:
            if decision.run_id != summary.run_id or decision.cell_id not in plans:
                msg = "batch decision origin does not match its run"
                raise ValueError(msg)
            decision.bind(plans[decision.cell_id])
            key = (decision.cell_id, decision.index)
            if key in decisions:
                msg = "duplicate batch decision index"
                raise ValueError(msg)
            feedback = decision.batch.feedback
            budget.retain(
                2
                + len(decision.batch.part_ids)
                + (0 if feedback is None else len(feedback.used) + len(feedback.failed))
            )
            decisions[key] = decision
    if len(decisions) != expected:
        msg = "batch records do not reconcile with the committed batch count"
        raise ValueError(msg)
    # Feedback recounts retain at most two maps of the already bounded eligible IDs.
    budget.retain(
        sum(
            2 * len(p.request.parts)
            for p in plans.values()
            if p.request.resampling is not None
        )
    )
    return RuntimeBatches(decisions)
