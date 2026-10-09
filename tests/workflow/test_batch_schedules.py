"""Ordered batches advance under local limits without changing cell targets.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app
from dense_arrays.workflow.inputs import read_source


def source() -> planning.GenerationPlan:
    return da.plan(
        planning.DesignSpec(
            parts=[
                parts.Part("a", "AAA"),
                parts.Part("b", "CCC"),
                parts.Part("c", "GGG"),
            ],
            length=planning.Length(maximum=3),
            strands="single",
            target=planning.Target(3),
        )
    )


def schedule(*members: tuple[str, ...], attempts: int = 2) -> planning.GenerationPlan:
    plan = source()
    selections = tuple(
        planning.CandidateBatch(ids, plan.collection_id, stream=f"step{i}")
        for i, ids in enumerate(members)
    )
    return da.plan(
        plan.request.with_changes(
            schedule=planning.BatchSchedule(selections, attempts_per_batch=attempts)
        )
    )


def test_schedule_advances_after_exhaustion_and_keeps_cell_uniqueness(tmp_path: Path):
    plan = schedule(("a",), ("a", "b"), ("c",), attempts=3)
    result = da.run(plan, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == 3
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    assert {a.evidence["batch_index"] for a in attempts} == {1, 2, 3}
    assert summary.counts["duplicate"] == 1
    assert attempts[0].evidence["batch_attempt"] == 1
    for record in da.inspect(result, view="designs", all=True).records():
        assert record.plan_id == plan.plan_id
        assert record.batch_id in {b.batch_id for b in plan.request.schedule.batches}
    da.export(result, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == 3


def test_per_batch_attempt_cap_stops_without_claiming_global_exhaustion(tmp_path: Path):
    plan = schedule(("a", "b"), ("c",), attempts=1)
    result = da.run(plan, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 2
    assert summary.counts["started"] == 2
    assert summary.termination_reason == "batch_schedule_exhausted"
    assert summary.state == "stopped"


def test_matrix_schedule_keeps_independent_progress_and_shared_budget(tmp_path: Path):
    plan = schedule(("a",), ("b",), ("c",), attempts=1)
    request = planning.MatrixSpec(
        base=plan.request.with_changes(
            target=planning.Target(), schedule=None, limits=planning.Limits(attempts=5)
        ),
        axes={"x": {"a": planning.Variant(), "b": planning.Variant()}},
        allocation=planning.Allocation(per_cell=3),
        max_cells=2,
        batches={"x=a": plan.request.schedule, "x=b": plan.request.schedule},
    )
    result = da.run(request, out=tmp_path / "matrix")
    summary = da.inspect(result, verify=True)
    assert summary.counts["started"] == summary.accepted == 5
    assert summary.termination_reason == "attempt_limit"
    records = list(da.inspect(result, view="attempts", all=True).records())
    assert [a.cell_id for a in records] == ["x=a", "x=b", "x=a", "x=b", "x=a"]
    assert [a.evidence["batch_index"] for a in records] == [1, 1, 2, 2, 3]


@pytest.mark.parametrize("boundary", ["reserve", "publish"])
def test_resume_restores_batch_frontier_without_restarting_old_models(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, boundary: str
):
    plan = schedule(("a",), ("b",), ("c",), attempts=1)
    original = getattr(RunWriter, boundary)

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> object:
        result = original(self, *args, **kwargs)
        if self.state["counts"]["started"] == 2:
            raise KeyboardInterrupt
        return result

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, boundary, interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(plan, out=tmp_path / "run")
    before = da.inspect(tmp_path / "run", verify=True)
    assert before.resumable
    result = da.run(resume=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.counts["started"] == 3
    assert summary.accepted == (2 if boundary == "reserve" else 3)
    records = list(da.inspect(result, view="attempts", all=True).records())
    assert [a.evidence["batch_index"] for a in records] == [1, 2, 3]


def test_prepare_schedule_python_cli_and_plan_round_trip(tmp_path: Path):
    base = source()
    base.write(tmp_path / "source.json")
    result = da.prepare(
        base,
        sampling=planning.BatchSampling(2, seed=7),
        batch_count=3,
        attempts_per_batch=2,
        out=tmp_path / "schedule.json",
    )
    assert len(result.request.schedule.batches) == 3
    assert result.request.batch is None
    assert read_source(tmp_path / "schedule.json").plan_id == result.plan_id
    cli = CliRunner().invoke(
        app,
        [
            "prepare",
            str(tmp_path / "source.json"),
            "--batch-size",
            "2",
            "--batch-seed",
            "7",
            "--batch-count",
            "3",
            "--attempts-per-batch",
            "2",
            "--out",
            str(tmp_path / "cli.json"),
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["plan_id"] == result.plan_id


@pytest.mark.parametrize("cells", [1, 2])
@pytest.mark.parametrize("boundary", ["reserve", "publish"])
def test_resume_after_exhausted_batch_preserves_committed_prefix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cells: int, boundary: str
):
    plan = schedule(("a",), ("b",), ("c",), attempts=2)
    request = (
        plan
        if cells == 1
        else planning.MatrixSpec(
            base=plan.request.with_changes(target=planning.Target(), schedule=None),
            axes={"x": {"a": planning.Variant(), "b": planning.Variant()}},
            allocation=planning.Allocation(per_cell=3),
            max_cells=2,
            batches={"x=a": plan.request.schedule, "x=b": plan.request.schedule},
        )
    )
    original = getattr(RunWriter, boundary)

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> object:
        result = original(self, *args, **kwargs)
        if self.state["counts"]["started"] == 2 * cells + 1:
            raise KeyboardInterrupt
        return result

    path = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, boundary, interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request, out=path)
    assert da.inspect(path, verify=True).resumable
    before = [
        a.to_dict() for a in da.inspect(path, view="attempts", all=True).records()
    ]
    assert sum(a["outcome"] == "no_candidate" for a in before) == cells
    result = da.run(resume=path)
    summary = da.inspect(result, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == 3 * cells
    assert summary.counts["started"] == 5 * cells
    after = [
        a.to_dict() for a in da.inspect(result, view="attempts", all=True).records()
    ]
    assert after[: len(before)] == before


def test_resume_after_final_batch_exhaustion_closes_without_another_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> object:
        result = publish(self, *args, **kwargs)
        if self.state["counts"]["started"] == 4:
            raise KeyboardInterrupt
        return result

    path = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(schedule(("a",), ("b",), attempts=2), out=path)
    result = da.run(resume=path)
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 2
    assert summary.counts["started"] == 4
    assert summary.state == "stopped"
    assert summary.termination_reason == "batch_schedule_exhausted"


@pytest.mark.parametrize(
    "status", ["backend_error", "invalid_result", "unknown", "unproven"]
)
def test_interruption_does_not_make_terminal_solver_outcome_resumable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, status: str
):
    from dense_arrays.artifacts.recovery import RecoveryError  # noqa: PLC0415
    from dense_arrays.solver import SolveReport, SolveStatus  # noqa: PLC0415

    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> object:
        publish(self, *args, **kwargs)
        raise KeyboardInterrupt

    path = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        patch.setattr(
            da.Optimizer,
            "solve_report",
            lambda *_a, **_kw: SolveReport(SolveStatus(status), None),
        )
        with pytest.raises(KeyboardInterrupt):
            da.run(schedule(("a",), ("b",), attempts=2), out=path)
    before = (path / "run.sqlite3").read_bytes()
    with pytest.raises(RecoveryError, match=r"cannot be resumed|not recoverable"):
        da.run(resume=path)
    assert (path / "run.sqlite3").read_bytes() == before


@pytest.mark.parametrize("operation", ["plan", "inspect"])
def test_cli_preview_explains_schedule_effort(tmp_path: Path, operation: str):
    plan = schedule(("a",), ("b",), ("c",), attempts=2)
    plan.write(tmp_path / "plan.json")
    args = [operation, str(tmp_path / "plan.json")]
    if operation == "inspect":
        args.extend(["--view", "plan"])
    response = CliRunner().invoke(app, args)
    assert response.exit_code == 0, response.output
    assert (
        "Batches: 3; maximum offered parts: 1; attempts per batch: 2" in response.stdout
    )
    if operation == "plan":
        assert "at most 2 path variables" in response.stdout


def test_schedule_accounting_rejects_rewritten_batch_ordinals(tmp_path: Path):
    import sqlite3  # noqa: PLC0415

    from dense_arrays._record_validation import (  # noqa: PLC0415
        canonical_json,
        semantic_digest,
    )

    result = da.run(schedule(("a",), ("b",), ("c",), attempts=1), out=tmp_path / "run")
    with sqlite3.connect(result.path / "run.sqlite3") as connection:
        rows = connection.execute(
            "SELECT revision,payload FROM attempts WHERE attempt=2"
        ).fetchall()
        for revision, raw in rows:
            value = json.loads(raw)
            value["evidence"]["batch_attempt"] = 2
            connection.execute(
                "UPDATE attempts SET payload=?,digest=? WHERE attempt=2 AND revision=?",
                (canonical_json(value), semantic_digest(value), revision),
            )
    with pytest.raises(ValueError, match="batch"):
        da.inspect(result, verify=True)


def test_schedule_membership_is_charged_to_reader_limits():
    from dense_arrays import reporting  # noqa: PLC0415

    plan = schedule(("a",), ("b",), ("c",), attempts=1)
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(plan, view="plan", read_limits=reporting.ReadLimits(identities=3))


def test_matrix_resumes_with_different_batch_frontiers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    plan = schedule(("a",), ("b",), ("c",), attempts=1)
    request = planning.MatrixSpec(
        base=plan.request.with_changes(target=planning.Target(), schedule=None),
        axes={"x": {"a": planning.Variant(), "b": planning.Variant()}},
        allocation=planning.Allocation(per_cell=3),
        max_cells=2,
        batches={"x=a": plan.request.schedule, "x=b": plan.request.schedule},
    )
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> str:
        result = publish(self, *args, **kwargs)
        if self.state["counts"]["started"] == 3:
            raise KeyboardInterrupt
        return result

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request, out=tmp_path / "run")
    before = da.inspect(tmp_path / "run", verify=True)
    assert before.resumable
    result = da.run(resume=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 6
    records = list(da.inspect(result, view="attempts", all=True).records())
    assert [(a.cell_id, a.evidence["batch_index"]) for a in records] == [
        ("x=a", 1),
        ("x=b", 1),
        ("x=a", 2),
        ("x=b", 2),
        ("x=a", 3),
        ("x=b", 3),
    ]
    da.export(result, view="request", out=tmp_path / "request.json")
    assert (
        da.plan(read_source(tmp_path / "request.json")).plan_id
        == da.plan(request).plan_id
    )


def test_backend_failure_never_advances_to_another_batch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.optimizer import Optimizer  # noqa: PLC0415
    from dense_arrays.solver import SolveReport, SolveStatus  # noqa: PLC0415

    calls = []

    def fail(_self: Optimizer, **_kwargs: object) -> SolveReport:
        calls.append(1)
        return SolveReport(SolveStatus.BACKEND_ERROR, None, detail="fixture failure")

    monkeypatch.setattr(Optimizer, "solve_report", fail)
    result = da.run(schedule(("a",), ("b",), attempts=2), out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.state == "failed"
    assert summary.counts["error"] == 1
    assert calls == [1]


def test_rejected_candidates_consume_batch_attempt_allowance(tmp_path: Path):
    plan = schedule(("a",), ("b", "c"), attempts=1)
    request = plan.request.with_changes(
        requirements=(planning.Avoid("no_a", patterns=("AAA",)),)
    )
    result = da.run(request, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.counts["rejected"] == 1
    assert summary.accepted == 1
    assert summary.counts["started"] == 2
    assert summary.termination_reason == "batch_schedule_exhausted"


@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_batch_schedule_limits_are_explicit_positive_integers(value: object):
    base = source()
    batch = planning.CandidateBatch(("a",), base.collection_id)
    with pytest.raises((TypeError, ValueError), match="attempts_per_batch"):
        planning.BatchSchedule((batch,), attempts_per_batch=value)


def test_preparation_rejects_missing_schedule_limit_before_output(tmp_path: Path):
    with pytest.raises(ValueError, match="attempts_per_batch"):
        da.prepare(
            source(),
            sampling=planning.BatchSampling(1),
            batch_count=2,
            out=tmp_path / "invalid",
        )
    assert not (tmp_path / "invalid").exists()
    planned = schedule(("a",), ("b",))
    with pytest.raises(ValueError, match="already"):
        da.prepare(planned, sampling=planning.BatchSampling(1), out=tmp_path / "again")
    assert not (tmp_path / "again").exists()


def test_schedule_requires_known_collections_and_exclusive_selection(tmp_path: Path):
    base = source()
    unknown = planning.CandidateBatch(("a",), "0" * 64)
    with pytest.raises(ValueError, match="collection"):
        da.run(
            base.request.with_changes(schedule=planning.BatchSchedule((unknown,), 1)),
            out=tmp_path / "unknown",
        )
    assert not (tmp_path / "unknown").exists()
    selection = planning.CandidateBatch(("a",), base.collection_id)
    with pytest.raises(ValueError, match="exclusive"):
        base.request.with_changes(
            batch=selection, schedule=planning.BatchSchedule((selection,), 1)
        )


def test_prepared_replay_order_matches_pinned_artifact_iterator(tmp_path: Path):
    fixture = json.loads(
        (
            Path(__file__).parents[1]
            / "fixtures/workflow/densegen-batch-replay-v1.json"
        ).read_text()
    )
    base = da.plan(
        planning.DesignSpec(
            parts=[
                parts.Part(f"p{entry['index']}", entry["sequences"][0])
                for entry in fixture["selected"]
            ],
            length=planning.Length(maximum=3),
            strands="single",
            target=planning.Target(4),
        )
    )
    batches = tuple(
        planning.CandidateBatch((f"p{entry['index']}",), base.collection_id)
        for entry in fixture["selected"]
    )
    result = da.run(
        base.request.with_changes(schedule=planning.BatchSchedule(batches, 1)),
        out=tmp_path / "run",
    )
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    assert [a.evidence["batch_index"] for a in attempts] == fixture[
        "constraint_callback_indices"
    ]
    designs = list(da.inspect(result, view="designs", all=True).records())
    assert [d.realized.sequence for d in designs] == [
        entry["sequences"][0] for entry in fixture["selected"]
    ]
    assert "exhausted" in fixture["exhausted"]
    assert (
        da.inspect(result, verify=True).termination_reason == "batch_schedule_exhausted"
    )


def test_schedule_rechecks_deadline_after_model_construction(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    import time  # noqa: PLC0415
    from types import SimpleNamespace  # noqa: PLC0415

    from dense_arrays.workflow import schedules  # noqa: PLC0415

    offset = [0]

    def construct(*_args: object, **_kwargs: object) -> object:
        offset[0] += 1000
        return object()

    monkeypatch.setattr(
        schedules,
        "time",
        SimpleNamespace(monotonic=lambda: time.monotonic() + offset[0]),
    )
    monkeypatch.setattr(schedules, "_restore_optimizer", construct)
    result = da.run(schedule(("a",), ("b",), attempts=1), out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.counts["started"] == 0
    assert summary.termination_reason == "active_time_limit"


def test_invalid_batch_reservation_does_not_advance_commit_state(tmp_path: Path):
    from dense_arrays.artifacts.store import create_run  # noqa: PLC0415

    planned = schedule(("a",), attempts=1)
    with create_run(planned, tmp_path / "run") as writer:
        revision = writer.state["revision"]
        with pytest.raises((TypeError, ValueError), match="batch_index"):
            writer.reserve(
                active_seconds=0,
                batch_id=planned.request.schedule.batches[0].batch_id,
                batch_index=0,
                batch_attempt=1,
            )
        assert writer.state["revision"] == revision
        assert writer.state["counts"]["started"] == 0
