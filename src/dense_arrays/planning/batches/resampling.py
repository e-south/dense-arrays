"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/batches/resampling.py

Finite runtime batch selection with explicit local search and feedback policy.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass

from dense_arrays._record_validation import integer, object_fields

from .feedback import FeedbackPolicy
from .models import BatchSampling, validate_unproven_policy

RESAMPLING_SCHEMA = "dense_arrays.resampling.v1"


@dataclass(frozen=True)
class Resampling:
    """Select a new offered batch at each exhaustion or local search boundary."""

    sampling: BatchSampling
    max_batches: int
    attempts_per_batch: int
    accepted_per_batch: int | None = None
    feedback: FeedbackPolicy | None = None
    on_unproven: str = "stop"

    def __post_init__(self) -> None:
        """Require finite effort and typed policies; global limits still apply."""
        validate_unproven_policy(self.on_unproven)
        if not isinstance(self.sampling, BatchSampling):
            msg = "resampling.sampling must be BatchSampling"
            raise TypeError(msg)
        for name in ("max_batches", "attempts_per_batch"):
            integer(getattr(self, name), field_name=f"resampling.{name}", minimum=1)
        if self.accepted_per_batch is not None:
            integer(
                self.accepted_per_batch,
                field_name="resampling.accepted_per_batch",
                minimum=1,
            )
        if self.feedback is not None and not isinstance(self.feedback, FeedbackPolicy):
            msg = "resampling.feedback must be FeedbackPolicy"
            raise TypeError(msg)

    def to_dict(self) -> dict[str, object]:
        """Persist finite effort and versioned selection meaning without drawing."""
        return {
            "schema": RESAMPLING_SCHEMA,
            "sampling": self.sampling.to_dict(),
            "max_batches": self.max_batches,
            "attempts_per_batch": self.attempts_per_batch,
            "accepted_per_batch": self.accepted_per_batch,
            "feedback": None if self.feedback is None else self.feedback.to_dict(),
            **({"on_unproven": self.on_unproven} if self.on_unproven != "stop" else {}),
        }

    @classmethod
    def from_dict(cls, value: object) -> "Resampling":
        """Reject unknown versions and incomplete persisted policies."""
        keys = {
            "schema",
            "sampling",
            "max_batches",
            "attempts_per_batch",
            "accepted_per_batch",
            "feedback",
        }
        data = object_fields(value, keys | {"on_unproven"}, "resampling")
        if not keys <= data.keys() or data.pop("schema") != RESAMPLING_SCHEMA:
            msg = "unsupported or incomplete resampling policy"
            raise ValueError(msg)
        data["sampling"] = BatchSampling.from_dict(data["sampling"])
        if data["feedback"] is not None:
            data["feedback"] = FeedbackPolicy.from_dict(data["feedback"])
        return cls(**data)
