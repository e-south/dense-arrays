"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/batches/schedules.py

Ordered replay selections with explicit per-batch search allowances.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from collections.abc import Mapping
from dataclasses import dataclass

from dense_arrays._record_validation import integer, object_fields, records

from .models import CandidateBatch, validate_unproven_policy
from .resampling import RESAMPLING_SCHEMA, Resampling

SCHEDULE_SCHEMA = "dense_arrays.batch_schedule.v1"


@dataclass(frozen=True)
class BatchSchedule:
    """Advance after exhaustion or a local cap; targets stay cell-owned."""

    batches: tuple[CandidateBatch, ...]
    attempts_per_batch: int
    accepted_per_batch: int | None = None
    on_unproven: str = "stop"

    def __post_init__(self) -> None:
        """Require a finite nonempty sequence from one eligible collection."""
        validate_unproven_policy(self.on_unproven)
        object.__setattr__(
            self,
            "batches",
            records(self.batches, CandidateBatch, field_name="schedule.batches"),
        )
        integer(self.attempts_per_batch, field_name="attempts_per_batch", minimum=1)
        if self.accepted_per_batch is not None:
            integer(self.accepted_per_batch, field_name="accepted_per_batch", minimum=1)
        if not self.batches or len({b.collection_id for b in self.batches}) != 1:
            msg = "a batch schedule requires nonempty batches from one collection"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Persist the full replay order and its finite attempt policy."""
        return {
            "schema": SCHEDULE_SCHEMA,
            "batches": [b.to_dict() for b in self.batches],
            "attempts_per_batch": self.attempts_per_batch,
            **({"on_unproven": self.on_unproven} if self.on_unproven != "stop" else {}),
            **(
                {"accepted_per_batch": self.accepted_per_batch}
                if self.accepted_per_batch is not None
                else {}
            ),
        }

    @classmethod
    def from_dict(cls, value: object) -> "BatchSchedule":
        """Load exact selections without sampling or changing their order."""
        keys = {"schema", "batches", "attempts_per_batch"}
        data = object_fields(
            value, keys | {"accepted_per_batch", "on_unproven"}, "batch schedule"
        )
        if not keys <= data.keys() or data.pop("schema") != SCHEDULE_SCHEMA:
            msg = "unsupported or incomplete batch schedule"
            raise ValueError(msg)
        if not isinstance(data["batches"], list):
            msg = "schedule.batches must be an ordered array"
            raise TypeError(msg)
        return cls(
            tuple(CandidateBatch.from_dict(b) for b in data["batches"]),
            data["attempts_per_batch"],
            accepted_per_batch=data.get("accepted_per_batch"),
            on_unproven=data.get("on_unproven", "stop"),
        )


def selection_from_dict(value: object) -> CandidateBatch | BatchSchedule | Resampling:
    """Dispatch explicit single-batch or schedule schemas at the matrix boundary."""
    if isinstance(value, Mapping) and value.get("schema") == RESAMPLING_SCHEMA:
        return Resampling.from_dict(value)
    if isinstance(value, Mapping) and value.get("schema") == SCHEDULE_SCHEMA:
        return BatchSchedule.from_dict(value)
    return CandidateBatch.from_dict(value)
