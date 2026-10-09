"""Matrix continuation preserves cell frontiers and one shared effort allowance.

Author: Eric J. South.
"""

import time
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.records import RunHandle
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app
from dense_arrays.workflow import matrices


def matrix(*, attempts: int = 20) -> planning.MatrixSpec:
    return planning.MatrixSpec(
        base=planning.DesignSpec(
            parts=[
                parts.Part("a", "AAA"),
                parts.Part("b", "CCC"),
                parts.Part("c", "GGG"),
            ],
            length=planning.Length(maximum=6),
            strands="single",
            limits=planning.Limits(attempts=attempts),
        ),
        axes={"x": {name: planning.Variant() for name in ("a", "b", "inactive")}},
        allocation=planning.Allocation(counts={"x=a": 2, "x=b": 2, "x=inactive": 0}),
        max_cells=3,
    )


def interrupt_after_design(
    path: Path, spec: planning.MatrixSpec, monkeypatch: pytest.MonkeyPatch
) -> None:
    publish = RunWriter.publish

    def interrupted(self: RunWriter, *args: object, **kwargs: object) -> str:
        observed = publish(self, *args, **kwargs)
        if observed == "accepted":
            raise KeyboardInterrupt
        return observed

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupted)
        with pytest.raises(KeyboardInterrupt):
            da.run(spec, out=path)


def test_resume_keeps_prefix_cell_uniqueness_and_round_robin_position(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = tmp_path / "matrix"
    interrupt_after_design(path, matrix(), monkeypatch)
    before = da.inspect(path, verify=True)
    assert before.resumable
    prefix = [d.to_dict() for d in da.inspect(path, view="designs", all=True).records()]
    inode = (path / ".writer.lock").stat().st_ino
    result = da.run(resume=path)
    after = da.inspect(result, verify=True)
    assert after.run_id == before.run_id
    assert after.state == "completed"
    assert after.accepted == after.target == after.counts["started"] == 4
    assert after.active_seconds > before.active_seconds
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    assert [a.cell_id for a in attempts] == ["x=a", "x=b", "x=a", "x=b"]
    assert [a.evidence["cell_attempt"] for a in attempts] == [1, 1, 2, 2]
    designs = list(da.inspect(result, view="designs", all=True).records())
    assert [d.to_dict() for d in designs[: len(prefix)]] == prefix
    for name in ("x=a", "x=b"):
        assert len({d.sequence_id for d in designs if d.cell_id == name}) == 2
    assert after.cells["x=inactive"].counts["started"] == 0
    assert (path / ".writer.lock").stat().st_ino == inode
    bound = RunHandle(path, before.run_id, before.revision)
    assert da.inspect(bound, verify=True).to_dict() == before.to_dict()


def test_cell_exhaustion_is_committed_with_its_attempt_before_interruption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    spec = matrix()
    spec = spec.with_changes(
        base=spec.base.with_changes(parts=[parts.Part("a", "AAA")]),
    )
    publish = RunWriter.publish

    def interrupted(self: RunWriter, *args: object, **kwargs: object) -> str:
        observed = publish(self, *args, **kwargs)
        if observed == "no_candidate":
            raise KeyboardInterrupt
        return observed

    path = tmp_path / "matrix"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupted)
        with pytest.raises(KeyboardInterrupt):
            da.run(spec, out=path)
    before = da.inspect(path, verify=True)
    assert before.resumable
    assert before.cells["x=a"].termination_reason == "batch_exhausted"
    closed = before.cells["x=a"].to_dict()
    after = da.inspect(da.run(resume=path), verify=True)
    assert after.state == "stopped"
    assert after.termination_reason == "cells_exhausted"
    assert after.counts["started"] == 4
    assert after.cells["x=a"].to_dict() == closed
    assert after.cells["x=b"].termination_reason == "batch_exhausted"


@pytest.mark.parametrize("outcome", ["duplicate", "rejected"])
def test_replay_keeps_nonaccepted_exclusions_separate_between_cells(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, outcome: str
):
    spec = matrix()
    duplicate = outcome == "duplicate"
    spec = spec.with_changes(
        base=spec.base.with_changes(
            parts=[
                parts.Part("a", "AAA"),
                parts.Part("b", "AAA" if duplicate else "CCC"),
            ],
            length=planning.Length(maximum=3 if duplicate else 6),
            requirements=[]
            if duplicate
            else [
                planning.Fixed("a-fixed", "a", "forward"),
                planning.Fixed("b-fixed", "b", "forward"),
                planning.Avoid("no-A", patterns=("A",), strands="forward"),
            ],
        ),
        allocation=planning.Allocation(per_cell=3 if duplicate else 1),
    )
    publish = RunWriter.publish

    def interrupted(self: RunWriter, *args: object, **kwargs: object) -> str:
        result = publish(self, *args, **kwargs)
        if result == outcome:
            raise KeyboardInterrupt
        return result

    path = tmp_path / "matrix"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupted)
        with pytest.raises(KeyboardInterrupt):
            da.run(spec, out=path)
    assert da.inspect(path).resumable
    after = da.inspect(da.run(resume=path), verify=True)
    assert after.counts["started"] == 9
    assert after.termination_reason == "cells_exhausted"
    assert all(
        c.counts[outcome] == (1 if duplicate else 2) for c in after.cells.values()
    )
    assert all(c.counts["no_candidate"] == 1 for c in after.cells.values())


def test_unresolved_attempt_keeps_cell_ordinal_and_shared_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    solve = da.Optimizer.solve_report
    calls = 0

    def interrupted(
        self: da.Optimizer, *args: object, **kwargs: object
    ) -> da.solver.SolveReport:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise KeyboardInterrupt
        return solve(self, *args, **kwargs)

    path = tmp_path / "matrix"
    with monkeypatch.context() as patch:
        patch.setattr(da.Optimizer, "solve_report", interrupted)
        with pytest.raises(KeyboardInterrupt):
            da.run(matrix(attempts=4), out=path)
    before = da.inspect(path, verify=True)
    assert before.counts["interrupted_unresolved"] == 1
    after = da.inspect(da.run(resume=path), verify=True)
    assert after.termination_reason == "attempt_limit"
    assert after.counts["started"] == 4
    assert after.accepted == 3
    records = list(da.inspect(path, view="attempts", all=True).records())
    assert [a.cell_id for a in records] == ["x=a", "x=b", "x=a", "x=b"]
    assert [a.evidence["cell_attempt"] for a in records] == [1, 1, 2, 2]
    assert after.counts["interrupted_unresolved"] == 1
    saved = (path / "run.sqlite3").read_bytes()
    with pytest.raises(ValueError, match="attempt budget"):
        da.run(resume=path)
    assert (path / "run.sqlite3").read_bytes() == saved


def test_interruption_at_resume_commit_is_measured_and_recoverable(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = tmp_path / "matrix"
    interrupt_after_design(path, matrix(), monkeypatch)
    before = da.inspect(path)
    commit = RunWriter._commit  # noqa: SLF001 - exact native publication boundary

    def interrupted(self: RunWriter, state: dict, **kwargs: object) -> None:
        commit(self, state, **kwargs)
        if state["state"] == "running" and not kwargs:
            raise KeyboardInterrupt

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "_commit", interrupted)
        with pytest.raises(KeyboardInterrupt):
            da.run(resume=path)
    after = da.inspect(path, verify=True)
    assert after.state == "stopped"
    assert after.termination_reason == "interrupted"
    assert after.resumable
    assert after.counts == before.counts
    assert after.active_seconds > before.active_seconds
    assert da.inspect(da.run(resume=path), verify=True).accepted == 4


def test_completed_matrix_resume_verifies_without_model_or_byte_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    result = da.run(matrix(), out=tmp_path / "matrix")
    saved = (result.path / "run.sqlite3").read_bytes()

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("completed matrix resume built a solver")

    monkeypatch.setattr(da.Optimizer, "build_model", forbidden)
    assert da.run(resume=result.path) == result
    response = CliRunner().invoke(app, ["run", "--resume", str(result.path), "--json"])
    assert response.exit_code == 0, response.output
    assert (result.path / "run.sqlite3").read_bytes() == saved


def test_replay_time_exhaustion_closes_reader_before_terminal_commit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    path = tmp_path / "matrix"
    interrupt_after_design(path, matrix(), monkeypatch)
    before = da.inspect(path)
    now = time.monotonic()
    clock = [now]
    restore = matrices._restore_optimizer  # noqa: SLF001 - charge model restoration

    def expired(*args: object, **kwargs: object) -> da.Optimizer:
        optimizer = restore(*args, **kwargs)
        clock[0] += 301
        return optimizer

    monkeypatch.setattr(matrices, "time", SimpleNamespace(monotonic=lambda: clock[0]))
    monkeypatch.setattr(matrices, "_restore_optimizer", expired)
    after = da.inspect(da.run(resume=path), verify=True)
    assert after.termination_reason == "active_time_limit"
    assert after.counts == before.counts
    assert after.active_seconds > 300
    assert not after.resumable


def test_replay_deadline_crossing_does_not_pass_negative_solver_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    path = tmp_path / "matrix"
    interrupt_after_design(path, matrix(), monkeypatch)
    before = da.inspect(path)
    baseline = time.monotonic()
    ticks = iter((baseline + 299, baseline + 301, baseline + 302))
    monkeypatch.setattr(
        matrices, "time", SimpleNamespace(monotonic=lambda: next(ticks, baseline + 303))
    )
    after = da.inspect(da.run(resume=path), verify=True)
    assert after.termination_reason == "active_time_limit"
    assert after.counts == before.counts
