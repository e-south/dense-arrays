"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/batches.py

Prepare portable executable plans by freezing offered batches per active cell.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import replace
from pathlib import Path

from dense_arrays._record_validation import integer
from dense_arrays.generation.batches.sampling import sample_batch
from dense_arrays.planning import (
    BatchSampling,
    BatchSchedule,
    CandidateBatch,
    GenerationPlan,
    MatrixPlan,
)


def prepare_batches(  # noqa: PLR0913 - explicit shared preparation settings
    plan: GenerationPlan | MatrixPlan,
    policy: BatchSampling,
    out: Path,
    *,
    batch_count: int = 1,
    attempts_per_batch: int | None = None,
    accepted_per_batch: int | None = None,
) -> GenerationPlan | MatrixPlan:
    """Verify sources once and publish the complete selection without generation."""
    if not isinstance(policy, BatchSampling):
        msg = (
            "prepare requires an explicit BatchSampling policy for a generation plan; "
            "use run to generate designs"
        )
        raise TypeError(msg)
    integer(batch_count, field_name="batch_count", minimum=1)
    if attempts_per_batch is not None:
        integer(attempts_per_batch, field_name="attempts_per_batch", minimum=1)
    if accepted_per_batch is not None:
        integer(accepted_per_batch, field_name="accepted_per_batch", minimum=1)
    if (
        batch_count > 1 or accepted_per_batch is not None
    ) and attempts_per_batch is None:
        msg = "multiple batches or an accepted cap require attempts_per_batch"
        raise ValueError(msg)
    plan.verify_inputs()
    if isinstance(plan, MatrixPlan):
        if plan.request.batches:
            msg = (
                "plan already contains prepared batches; revise the request explicitly"
            )
            raise ValueError(msg)
        batches = {
            c.cell_id: _selection(
                c.plan,
                policy,
                stream=c.cell_id,
                count=batch_count,
                attempts=attempts_per_batch,
                accepted=accepted_per_batch,
            )
            for c in plan.cells
            if c.active
        }
        base = replace(
            plan.base,
            inputs=(),
            embedded_input_digests=plan.base.evidence.input_digests,
        )
        prepared = MatrixPlan(
            replace(
                plan.request,
                batches=batches,
                sources={
                    cell: replace(source, locations=None)
                    for cell, source in plan.request.sources.items()
                },
            ),
            base,
        )
    else:
        if plan.request.batch is not None or plan.request.schedule is not None:
            msg = (
                "plan already contains a prepared batch; revise the request explicitly"
            )
            raise ValueError(msg)
        selection = _selection(
            plan,
            policy,
            stream="default",
            count=batch_count,
            attempts=attempts_per_batch,
            accepted=accepted_per_batch,
        )
        prepared = replace(
            plan,
            request=plan.request.with_changes(
                batch=selection if isinstance(selection, CandidateBatch) else None,
                schedule=selection if isinstance(selection, BatchSchedule) else None,
            ),
            inputs=(),
            embedded_input_digests=plan.evidence.input_digests,
        )
    prepared.write(out)
    return prepared


def _selection(  # noqa: PLR0913 - cell context and explicit batch bounds
    plan: GenerationPlan,
    policy: BatchSampling,
    *,
    stream: str,
    count: int,
    attempts: int | None,
    accepted: int | None,
) -> CandidateBatch | BatchSchedule:
    if count == 1 and attempts is None:
        return sample_batch(plan, policy, stream=stream)
    batches = tuple(
        sample_batch(plan, policy, stream=f"{stream}/batch/{i}")
        for i in range(1, count + 1)
    )
    return BatchSchedule(
        batches, attempts_per_batch=attempts, accepted_per_batch=accepted
    )
