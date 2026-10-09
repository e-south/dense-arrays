"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/batch_cursor.py

Per-cell batch frontiers and committed-feedback accumulation.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

from dense_arrays.artifacts.batches import BatchDecision
from dense_arrays.artifacts.search import batch_search_finished
from dense_arrays.generation.batches.progress import BatchWork
from dense_arrays.generation.batches.sampling import sample_batch
from dense_arrays.planning import CandidateBatch, FeedbackSnapshot

if TYPE_CHECKING:
    from dense_arrays.artifacts.records import Attempt
    from dense_arrays.generation.packing import PackingEngine
    from dense_arrays.planning import GenerationPlan
    from dense_arrays.planning.batches import BatchSchedule, Resampling


@dataclass
class BatchCursor:
    """Advance finite selections while keeping independent committed cell usage."""

    original: GenerationPlan
    index: int = 1
    attempts: int = 0
    accepted: int = 0
    enumerated: bool = False
    ended: bool = False
    optimizer: PackingEngine | None = None
    concrete: GenerationPlan | None = None

    batch: CandidateBatch | None = None
    used: Counter = field(default_factory=Counter)
    failed: Counter = field(default_factory=Counter)

    @property
    def policy(self) -> BatchSchedule | Resampling | None:
        """Return the declared prepared or runtime batch effort policy."""
        return self.original.request.schedule or self.original.request.resampling

    def select(self, cell: str) -> None:
        """Draw only a new runtime frontier, using previous committed observations."""
        policy = self.original.request.resampling
        if policy is not None and self.batch is None:
            self.batch = sample_batch(
                self.original,
                policy.sampling,
                stream=f"{cell}/batch/{self.index}",
                feedback=None
                if policy.feedback is None
                else FeedbackSnapshot(policy.feedback, self.used, self.failed),
            )
            self.concrete = replace(
                self.original,
                request=self.original.request.with_changes(
                    resampling=None, batch=self.batch
                ),
            )

    def decision(self, run_id: str, cell: str, after: int) -> BatchDecision | None:
        """Attach membership only to its first durable reservation."""
        if self.original.request.resampling is None or self.attempts:
            return None
        return BatchDecision(
            run_id, cell, self.original.plan_id, self.index, after, self.batch
        )

    def restore(self, decision: BatchDecision) -> None:
        """Load recorded membership without invoking the sampler."""
        self.index = decision.index
        self.attempts, self.accepted, self.enumerated, self.ended = 0, 0, False, False
        self.batch = decision.batch
        self.optimizer = None
        self.concrete = decision.bind(self.original)

    @property
    def plan(self) -> GenerationPlan:
        """Bind the current prepared selection; runtime selections are explicit."""
        if self.concrete is None:
            schedule = self.original.request.schedule
            self.concrete = (
                self.original
                if schedule is None
                else replace(
                    self.original,
                    request=self.original.request.with_changes(
                        schedule=None, batch=schedule.batches[self.index - 1]
                    ),
                )
            )
        return self.concrete

    @property
    def work(self) -> BatchWork | None:
        """Reserve logical batch ordinals independently of backend calls."""
        if self.policy is None:
            return None
        return BatchWork(
            self.index, self.attempts + 1, self.enumerated, self.policy.on_unproven
        )

    def observe(self, attempt: Attempt) -> None:
        """Include one committed outcome in the cell and batch frontier."""
        index = attempt.evidence.get("batch_index", 1)
        if self.index != index:
            self.index, self.attempts, self.accepted, self.enumerated, self.ended = (
                index,
                0,
                0,
                False,
                False,
            )
            self.optimizer = self.concrete = None
        if self.original.request.resampling is not None:
            if attempt.outcome == "accepted":
                self.used.update(
                    p.feature_id for p in attempt.candidate.final.placements
                )
            elif attempt.outcome == "rejected" or (
                attempt.evidence.get("solver_status") == "infeasible"
                and not self.enumerated
            ):
                self.failed.update(self.batch.part_ids)
        self.attempts += 1
        self.accepted += attempt.outcome == "accepted"
        self.enumerated |= attempt.candidate is not None
        schedule = self.policy
        self.ended = batch_search_finished(
            attempt.evidence, on_unproven=schedule.on_unproven if schedule else "stop"
        ) or (
            schedule is not None
            and (
                self.attempts >= schedule.attempts_per_batch
                or (
                    schedule.accepted_per_batch is not None
                    and self.accepted >= schedule.accepted_per_batch
                )
            )
        )

    def advance(self) -> bool:
        """Move only after exhaustion or a local cap; never replenish limits."""
        if not self.ended:
            return True
        schedule = self.policy
        maximum = (
            self.original.request.resampling.max_batches
            if self.original.request.resampling is not None
            else (len(schedule.batches) if schedule is not None else 1)
        )
        if schedule is None or self.index == maximum:
            return False
        self.index += 1
        self.attempts, self.accepted, self.enumerated, self.ended = 0, 0, False, False
        self.optimizer = self.concrete = self.batch = None
        return True
