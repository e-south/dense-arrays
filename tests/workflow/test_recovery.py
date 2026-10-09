"""Continue measured interruptions without resetting accepted records or budgets.

Author: Eric J. South.
"""

import json
import sqlite3
import time
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app


def request(*, attempts: int = 20):
    return planning.DesignSpec(
        parts=[parts.Part("a", "AAA"), parts.Part("b", "CCC"), parts.Part("c", "GGG")],
        length=planning.Length(maximum=6),
        strands="single",
        target=planning.Target(count=3),
        limits=planning.Limits(attempts=attempts),
    )


def interrupt_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, attempts: int = 20
):
    solve = da.Optimizer.solve_report
    calls = 0

    def interrupt(
        self: da.Optimizer, *args: object, **kwargs: object
    ) -> da.solver.SolveReport:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise KeyboardInterrupt
        return solve(self, *args, **kwargs)

    out = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(da.Optimizer, "solve_report", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request(attempts=attempts), out=out)
    return out


def test_resume_preserves_prefix_and_charges_interrupted_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = interrupt_run(tmp_path, monkeypatch)
    before = da.inspect(path, verify=True)
    assert before.accepted == 1
    assert before.counts["interrupted_unresolved"] == 1
    assert before.resumable
    with da.inspect(path, view="designs").records() as records:
        first = next(records).to_dict()
    lock_inode = (path / ".writer.lock").stat().st_ino
    result = da.run(resume=path)
    after = da.inspect(result, verify=True)
    assert after.run_id == before.run_id
    assert after.state == "completed"
    assert after.accepted == 3
    assert after.counts["started"] == 4
    assert after.counts["interrupted_unresolved"] == 1
    assert after.active_seconds > before.active_seconds
    with da.inspect(result, view="designs").records() as records:
        designs = list(records)
    assert designs[0].to_dict() == first
    assert len({d.sequence_id for d in designs}) == 3
    assert (path / ".writer.lock").stat().st_ino == lock_inode
    bound = da.artifacts.RunHandle(path, before.run_id, before.revision)
    assert da.inspect(bound, verify=True).to_dict() == before.to_dict()


def test_completed_resume_is_verified_noop_through_python_and_cli(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    result = da.run(request(), out=tmp_path / "run")
    before = (result.path / "run.sqlite3").read_bytes()

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("completed resume built a model")

    monkeypatch.setattr(da.Optimizer, "build_model", forbidden)
    assert da.run(resume=result.path) == result
    response = CliRunner().invoke(app, ["run", "--resume", str(result.path), "--json"])
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout)["state"] == "completed"
    assert (result.path / "run.sqlite3").read_bytes() == before


def test_resume_cannot_reset_exhausted_attempts_or_accept_overrides(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = interrupt_run(tmp_path, monkeypatch, attempts=2)
    before = (path / "run.sqlite3").read_bytes()
    assert da.inspect(path).resumable is False
    with pytest.raises(ValueError, match=r"attempt.*budget"):
        da.run(resume=path)
    with pytest.raises(ValueError, match="exclusive"):
        da.run(request(), resume=path)
    with pytest.raises(ValueError, match="exclusive"):
        da.run(resume=path, out=tmp_path / "other")
    response = CliRunner().invoke(
        app, ["run", "--resume", str(path), "--count", "5", "--json"]
    )
    assert response.exit_code == 2
    assert "exclusive" in response.stdout
    assert (path / "run.sqlite3").read_bytes() == before


def test_plan_preview_preserves_recorded_capability_and_new_plans_advertise_resume():
    plan = da.plan(request())
    assert plan.preview["resume_supported"] is True
    previous = plan.to_dict()
    previous["preview"]["resume_supported"] = False
    restored = planning.GenerationPlan.from_dict(previous)
    assert restored.plan_id == plan.plan_id
    assert restored.to_dict() == previous


@pytest.mark.parametrize("preview", [False, []])
def test_saved_preview_requires_an_object(preview: object):
    record = da.plan(request()).to_dict()
    record["preview"] = preview
    with pytest.raises(TypeError, match=r"preview.*object"):
        planning.GenerationPlan.from_dict(record)


def test_stale_inputs_and_corrupt_candidates_fail_before_resume_writes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = interrupt_run(tmp_path, monkeypatch)
    with sqlite3.connect(path / "run.sqlite3") as db:
        revision, payload = db.execute(
            "SELECT revision,payload FROM attempts WHERE attempt=1 "
            "ORDER BY revision DESC LIMIT 1"
        ).fetchone()
        value = json.loads(payload)
        value["evidence"]["candidate"]["packed"]["placements"][0]["feature_id"] = (
            "unknown"
        )
        db.execute(
            "UPDATE attempts SET payload=?,digest=? WHERE revision=?",
            (canonical_json(value), semantic_digest(value), revision),
        )
    before = (path / "run.sqlite3").read_bytes()
    with pytest.raises(ValueError, match="unknown part"):
        da.run(resume=path)
    assert (path / "run.sqlite3").read_bytes() == before

    source = tmp_path / "parts.csv"
    source.write_text("part_id,sequence\na,AAA\n")
    completed = da.run(
        planning.DesignSpec(
            parts=parts.PartTable(source, "csv"), length=planning.Length(maximum=3)
        ),
        out=tmp_path / "completed",
    )
    before = (completed.path / "run.sqlite3").read_bytes()
    source.write_text("part_id,sequence\na,CCC\n")
    with pytest.raises(ValueError, match="changed since planning"):
        da.run(resume=completed.path)
    assert (completed.path / "run.sqlite3").read_bytes() == before


def test_repeated_interruptions_consume_one_shared_time_allowance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = interrupt_run(tmp_path, monkeypatch)
    before = da.inspect(path)

    def interrupt(*_args: object, **_kwargs: object) -> None:
        raise KeyboardInterrupt

    with monkeypatch.context() as patch:
        patch.setattr(da.Optimizer, "solve_report", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(resume=path)
    interrupted = da.inspect(path, verify=True)
    assert interrupted.counts["started"] == 3
    assert interrupted.counts["interrupted_unresolved"] == 2
    assert interrupted.active_seconds > before.active_seconds
    with sqlite3.connect(path / "run.sqlite3") as db:
        revision, payload = db.execute(
            "SELECT revision,payload FROM commits ORDER BY revision DESC LIMIT 1"
        ).fetchone()
        state = json.loads(payload)
        state["active_seconds"] = 300.0
        state["resumable"] = False
        db.execute(
            "UPDATE commits SET payload=?,digest=? WHERE revision=?",
            (canonical_json(state), semantic_digest(state), revision),
        )
    saved = (path / "run.sqlite3").read_bytes()
    with pytest.raises(ValueError, match="active-time budget"):
        da.run(resume=path)
    assert (path / "run.sqlite3").read_bytes() == saved


def test_budget_expiry_during_replay_releases_reader_before_terminal_commit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = interrupt_run(tmp_path, monkeypatch)
    before = da.inspect(path)
    baseline = time.monotonic()
    ticks = iter((baseline, baseline + 301, baseline + 302, baseline + 303))
    monkeypatch.setattr(
        "dense_arrays.workflow.execution.time",
        SimpleNamespace(monotonic=lambda: next(ticks)),
    )
    result = da.run(resume=path)
    after = da.inspect(result, verify=True)
    assert after.accepted == before.accepted
    assert after.counts == before.counts
    assert after.termination_reason == "active_time_limit"
    assert after.resumable is False


@pytest.mark.parametrize("outcome", ["duplicate", "rejected"])
def test_resume_restores_nonaccepted_packing_exclusions(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, outcome: str
):
    spec = planning.DesignSpec(
        parts=[
            parts.Part("a", "AAA"),
            parts.Part("b", "AAA" if outcome == "duplicate" else "CCC"),
        ],
        length=planning.Length(maximum=3 if outcome == "duplicate" else 6),
        strands="single",
        target=planning.Target(count=3 if outcome == "duplicate" else 1),
        requirements=[]
        if outcome == "duplicate"
        else [
            planning.Fixed("a-fixed", "a", "forward"),
            planning.Fixed("b-fixed", "b", "forward"),
            planning.Avoid("no-A", patterns=("A",), strands="forward"),
        ],
    )
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> str:
        observed = publish(self, *args, **kwargs)
        if observed == outcome:
            raise KeyboardInterrupt
        return observed

    path = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(spec, out=path)
    before = da.inspect(path, verify=True)
    assert before.resumable
    after = da.inspect(da.run(resume=path), verify=True)
    assert after.counts["started"] == 3
    assert after.counts["duplicate" if outcome == "duplicate" else "rejected"] == (
        1 if outcome == "duplicate" else 2
    )
    assert after.termination_reason == "batch_exhausted"
