"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/batch_accounting.py

Independently reconcile scheduled batch order and per-batch attempt limits.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass

from dense_arrays.artifacts.records import Attempt
from dense_arrays.artifacts.search import HeuristicEvidence, batch_search_finished
from dense_arrays.planning import BatchSchedule, Resampling


@dataclass
class BatchAccounting:
    """Check committed attempts against a schedule, retaining one frontier."""

    schedule: BatchSchedule | Resampling | None
    index: int = 1
    attempts: int = 0
    accepted: int = 0
    closed: bool = False
    failed: bool = False

    def observe(self, record: Attempt) -> None:
        """Reject skipped/reopened batches, invalid ordinals and unearned retries."""
        if self.schedule is None:
            return
        index = record.evidence.get("batch_index")
        ordinal = record.evidence.get("batch_attempt")
        if self.closed and not self.failed:
            self.index += 1
            self.attempts = 0
            self.accepted = 0
            self.closed = False
        self.attempts += 1
        self.accepted += record.outcome == "accepted"
        if (
            self.failed
            or index != self.index
            or ordinal != self.attempts
            or self.index
            > (
                self.schedule.max_batches
                if isinstance(self.schedule, Resampling)
                else len(self.schedule.batches)
            )
            or self.attempts > self.schedule.attempts_per_batch
        ):
            msg = "batch order or attempt accounting violates the declared schedule"
            raise ValueError(msg)
        status = record.evidence.get("solver_status")
        finished = batch_search_finished(
            record.evidence, on_unproven=self.schedule.on_unproven
        )
        self.failed = record.outcome == "error" or (
            not finished
            and status
            not in {
                None,
                "optimal",
                "infeasible",
            }
        )
        if "heuristic" in record.evidence:
            self.failed |= (
                HeuristicEvidence.from_dict(record.evidence["heuristic"]).status
                == "time_limit"
            )
        self.closed = (
            finished
            or self.attempts == self.schedule.attempts_per_batch
            or self.accepted == self.schedule.accepted_per_batch
        )
