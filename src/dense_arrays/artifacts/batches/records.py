"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/batches/records.py

Bind one sampled membership to its cell and committed attempt prefix.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    integer,
    object_fields,
    required_text,
)
from dense_arrays.planning import CandidateBatch

if TYPE_CHECKING:
    from dense_arrays.planning import PlanEvidence

DECISION_SCHEMA = "dense_arrays.batch_decision.v1"


@dataclass(frozen=True)
class BatchDecision:
    """A once-recorded selection; after_attempt names its observed global prefix."""

    run_id: str
    cell_id: str
    plan_id: str
    index: int
    after_attempt: int
    batch: CandidateBatch

    def __post_init__(self) -> None:
        """Require complete provenance and positive logical batch coordinates."""
        for name in ("run_id", "cell_id"):
            required_text(getattr(self, name), field_name=f"batch decision.{name}")
        digest(self.plan_id, field_name="batch decision.plan_id")
        integer(self.index, field_name="batch decision.index", minimum=1)
        integer(
            self.after_attempt, field_name="batch decision.after_attempt", minimum=0
        )
        if not isinstance(self.batch, CandidateBatch):
            msg = "batch decision requires CandidateBatch"
            raise TypeError(msg)

    def bind(self, plan: "PlanEvidence") -> "PlanEvidence":
        """Validate saved context and bind geometry without sampling again."""
        policy = plan.request.resampling
        if (
            plan.plan_id != self.plan_id
            or policy is None
            or self.index > policy.max_batches
            or self.batch.sampling != policy.sampling
            or self.batch.stream != f"{self.cell_id}/batch/{self.index}"
            or (None if self.batch.feedback is None else self.batch.feedback.policy)
            != policy.feedback
        ):
            msg = "batch decision does not match its resampling plan"
            raise ValueError(msg)
        return replace(
            plan, request=plan.request.with_changes(resampling=None, batch=self.batch)
        )

    def to_dict(self) -> dict[str, object]:
        """Encode immutable membership once, separate from attempts and designs."""
        return {
            "schema": DECISION_SCHEMA,
            "run_id": self.run_id,
            "cell_id": self.cell_id,
            "plan_id": self.plan_id,
            "index": self.index,
            "after_attempt": self.after_attempt,
            "batch": self.batch.to_dict(),
        }

    @classmethod
    def from_dict(cls, value: object) -> "BatchDecision":
        """Reject incomplete and unsupported decision records."""
        keys = {
            "schema",
            "run_id",
            "cell_id",
            "plan_id",
            "index",
            "after_attempt",
            "batch",
        }
        data = object_fields(value, keys, "batch decision")
        if set(data) != keys or data.pop("schema") != DECISION_SCHEMA:
            msg = "unsupported or incomplete batch decision"
            raise ValueError(msg)
        data["batch"] = CandidateBatch.from_dict(data["batch"])
        return cls(**data)
