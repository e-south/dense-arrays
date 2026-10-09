"""Runtime selections and their first reservation share one committed boundary.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.batches import BatchDecision
from dense_arrays.artifacts.store import RunWriter, checked_payload, create_run
from dense_arrays.generation.batches.sampling import sample_batch


def recipe():
    return da.plan(
        planning.DesignSpec(
            (parts.Part("a", "AAA"), parts.Part("b", "CCC")),
            planning.Length(maximum=3),
            resampling=planning.Resampling(planning.BatchSampling(1), 3, 2),
        )
    )


def decision(writer: RunWriter, plan: planning.GenerationPlan):

    return BatchDecision(
        writer.handle.run_id,
        "default",
        plan.plan_id,
        1,
        0,
        sample_batch(plan, plan.request.resampling.sampling, stream="default/batch/1"),
    )


def test_batch_commit_rolls_back_with_failed_attempt_reservation(tmp_path: Path):
    plan = recipe()
    with create_run(plan, tmp_path / "run") as writer:
        record = decision(writer, plan)
        writer.connection.execute(
            "CREATE TRIGGER fail_attempt BEFORE INSERT ON attempts "
            "BEGIN SELECT RAISE(ABORT, 'injected failure'); END"
        )
        before = dict(writer.state)
        with pytest.raises(Exception, match="injected failure"):
            writer.reserve(
                active_seconds=0,
                batch_id=record.batch.batch_id,
                batch_index=1,
                batch_attempt=1,
                batch_decision=record,
            )
        assert writer.state == before
        assert (
            writer.connection.execute("SELECT count(*) FROM batches").fetchone()[0] == 0
        )
        assert (
            writer.connection.execute("SELECT count(*) FROM commits").fetchone()[0] == 1
        )
        writer.connection.execute("DROP TRIGGER fail_attempt")
        writer.reserve(
            active_seconds=0,
            batch_id=record.batch.batch_id,
            batch_index=1,
            batch_attempt=1,
            batch_decision=record,
        )
        assert writer.state["batch_count"] == 1
        saved = checked_payload(
            writer.connection.execute("SELECT payload,digest FROM batches").fetchone()
        )
        assert saved == record.to_dict()
        assert type(record).from_dict(saved) == record
        writer.publish(
            1, "interrupted_unresolved", {"code": "interrupted"}, active_seconds=0
        )
        writer.reserve(
            active_seconds=0,
            batch_id=record.batch.batch_id,
            batch_index=1,
            batch_attempt=2,
        )
        assert writer.state["batch_count"] == 1


def test_wrong_batch_binding_is_rejected_without_mutation(tmp_path: Path):
    plan = recipe()
    with create_run(plan, tmp_path / "run") as writer:
        record = decision(writer, plan)
        with pytest.raises(ValueError, match="batch decision"):
            writer.reserve(
                active_seconds=0,
                batch_id="0" * 64,
                batch_index=1,
                batch_attempt=1,
                batch_decision=record,
            )
        assert writer.state["counts"]["started"] == 0
        assert writer.state["batch_count"] == 0
