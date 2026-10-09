"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/batches/bindings.py

Resolve persisted batch references without drawing selections again.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer

from .models import CandidateBatch

if TYPE_CHECKING:
    from dense_arrays.planning.models import DesignSpec


def maximum_offered_parts(request: "DesignSpec") -> int:
    """Count the largest declared search input without drawing or building a model."""
    if request.schedule is not None:
        return max(len(batch.part_ids) for batch in request.schedule.batches)
    if request.resampling is not None:
        return request.resampling.sampling.size
    return len(request.parts) if request.batch is None else len(request.batch.part_ids)


def offered_batch(
    request: "DesignSpec", batch_id: str | None = None, *, index: int | None = None
) -> CandidateBatch | None:
    """Resolve a design or attempt to its exact declared offered membership."""
    if request.schedule is None:
        if index is not None:
            msg = "batch index requires a declared schedule"
            raise ValueError(msg)
        batch = request.batch
        if batch_id is not None and (batch is None or batch.batch_id != batch_id):
            msg = "batch identity does not match its offered plan"
            raise ValueError(msg)
        return batch
    if index is not None:
        integer(index, field_name="batch_index", minimum=1)
        if index > len(request.schedule.batches):
            msg = "batch index exceeds the declared schedule"
            raise ValueError(msg)
        candidates = (request.schedule.batches[index - 1],)
    else:
        candidates = request.schedule.batches
    for batch in candidates:
        if batch.batch_id == batch_id:
            return batch
    msg = "batch identity does not match its declared schedule"
    raise ValueError(msg)


def membership_size(request: "DesignSpec") -> int:
    """Count retained batch membership entries for bounded readers."""
    if request.schedule is not None:
        return sum(batch_state_size(b) for b in request.schedule.batches)
    return 0 if request.batch is None else batch_state_size(request.batch)


def batch_state_size(batch: CandidateBatch) -> int:
    """Count membership and observation entries retained by one selection."""
    feedback = batch.feedback
    return len(batch.part_ids) + (
        0 if feedback is None else len(feedback.used) + len(feedback.failed)
    )


def encoded_batch_size(value: object) -> int:
    """Bound membership and feedback maps before constructing typed state."""
    if value is None:
        return 0
    if not isinstance(value, dict) or not isinstance(value.get("part_ids"), list):
        msg = "prepared batch requires an ordered part identity array"
        raise TypeError(msg)
    size = len(value["part_ids"])
    if "feedback" in value:
        feedback = value["feedback"]
        if not isinstance(feedback, dict) or any(
            not isinstance(feedback.get(name), dict) for name in ("used", "failed")
        ):
            msg = "batch feedback requires observation maps"
            raise TypeError(msg)
        size += len(feedback["used"]) + len(feedback["failed"])
    return size
