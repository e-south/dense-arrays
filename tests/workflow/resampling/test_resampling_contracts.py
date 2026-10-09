"""Bounded runtime selection is explicit and portable in generation plans.

Author: Eric J. South.
"""

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.planning.serialization import request_from_dict, request_to_dict


def source():
    return planning.DesignSpec(
        parts=(parts.Part("a", "AAA"), parts.Part("b", "CCC")),
        length=planning.Length(maximum=3),
    )


def test_resampling_policy_survives_saved_plan_and_request():
    policy = planning.Resampling(
        planning.BatchSampling(1, seed=7),
        max_batches=4,
        attempts_per_batch=2,
        accepted_per_batch=1,
        feedback=planning.FeedbackPolicy(failure_alpha=2),
    )
    request = source().with_changes(resampling=policy)
    plan = da.plan(request)
    assert request_from_dict(request_to_dict(request)) == request
    assert (
        planning.GenerationPlan.from_dict(plan.to_dict()).request.resampling == policy
    )
    assert plan.plan_id != da.plan(source()).plan_id
    assert "resampling" not in request_to_dict(source())


@pytest.mark.parametrize(
    "field,value",
    [
        ("max_batches", 0),
        ("max_batches", True),
        ("attempts_per_batch", -1),
        ("accepted_per_batch", 0),
        ("feedback", {}),
        ("sampling", {}),
    ],
)
def test_resampling_rejects_implicit_or_unbounded_policy(field: str, value: object):
    values = {
        "sampling": planning.BatchSampling(1),
        "max_batches": 3,
        "attempts_per_batch": 2,
    }
    values[field] = value
    with pytest.raises((ValueError, TypeError)):
        planning.Resampling(**values)


def test_resampling_is_exclusive_with_prepared_membership():
    plan = da.plan(source())
    batch = planning.CandidateBatch(("a",), plan.collection_id)
    policy = planning.Resampling(
        planning.BatchSampling(1), max_batches=2, attempts_per_batch=1
    )
    with pytest.raises(ValueError, match="mutually exclusive"):
        source().with_changes(batch=batch, resampling=policy)
    with pytest.raises(ValueError, match="mutually exclusive"):
        source().with_changes(
            schedule=planning.BatchSchedule((batch,), 1), resampling=policy
        )


def test_invalid_runtime_eligibility_fails_in_planning():
    policy = planning.Resampling(
        planning.BatchSampling(3), max_batches=2, attempts_per_batch=1
    )
    with pytest.raises(ValueError, match="batch size"):
        da.plan(source().with_changes(resampling=policy))
    policy = planning.Resampling(
        planning.BatchSampling(1, unique_cores=True),
        max_batches=2,
        attempts_per_batch=1,
    )
    with pytest.raises(ValueError, match="group"):
        da.plan(source().with_changes(resampling=policy))


def test_matrix_cells_can_override_the_base_sampling_policy():
    base = source().with_changes(
        resampling=planning.Resampling(planning.BatchSampling(1), 2, 1)
    )
    matrix = planning.MatrixSpec(
        base,
        {"x": {"a": planning.Variant(), "b": planning.Variant()}},
        planning.Allocation(per_cell=1),
        2,
    )
    initial = da.plan(matrix)
    cell = initial.cells[0]
    specific = planning.Resampling(planning.BatchSampling(2), 4, 2)
    changed = da.plan(matrix.with_changes(batches={cell.cell_id: specific}))
    assert changed.cells[0].plan.request.resampling == specific
    assert changed.cells[1].plan.request.resampling == base.resampling
    batch = planning.CandidateBatch(("a",), cell.plan.collection_id)
    frozen = da.plan(matrix.with_changes(batches={cell.cell_id: batch}))
    assert frozen.cells[0].plan.request.resampling is None
    assert frozen.cells[0].plan.request.batch == batch
    assert (
        planning.MatrixSpec.from_dict(
            matrix.with_changes(batches={cell.cell_id: specific}).to_dict()
        ).batches[cell.cell_id]
        == specific
    )


def test_preview_describes_bounded_sampling_without_selecting():
    policy = planning.Resampling(planning.BatchSampling(1), 4, 2, accepted_per_batch=1)
    plan = da.plan(source().with_changes(resampling=policy))
    assert plan.preview["offered_parts"] == 1
    assert plan.preview["max_batches"] == 4
    assert plan.preview["path_variables"] == 6
