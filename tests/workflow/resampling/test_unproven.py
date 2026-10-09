"""Explicit unproven-search continuation preserves evidence and finite batch bounds.

Author: Eric J. South.
"""

from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app
from dense_arrays.reporting.batch_accounting import BatchAccounting
from dense_arrays.solver import SolveReport, SolveStatus


def request(
    kind: str, *, batches: int = 2
) -> planning.DesignSpec | planning.MatrixSpec:
    base = planning.DesignSpec(
        [parts.Part("a", "AAA")],
        planning.Length(maximum=3),
        strands="single",
    )
    if kind == "schedule":
        bound = da.plan(base)
        batch = planning.CandidateBatch(("a",), bound.collection_id)
        base = base.with_changes(
            schedule=planning.BatchSchedule(
                (batch,) * batches,
                attempts_per_batch=10,
                on_unproven="next_batch",
            )
        )
    else:
        base = base.with_changes(
            resampling=planning.Resampling(
                planning.BatchSampling(1),
                max_batches=batches,
                attempts_per_batch=10,
                on_unproven="next_batch",
            )
        )
    if kind == "matrix":
        return planning.MatrixSpec(
            base,
            {"x": {"a": planning.Variant(), "b": planning.Variant()}},
            planning.Allocation(per_cell=1),
            max_cells=2,
        )
    return base


@pytest.mark.parametrize("kind", ["schedule", "resampling", "matrix"])
@pytest.mark.parametrize("status", ["unknown", "unproven"])
def test_next_batch_after_unproven_search_preserves_exact_acceptance(
    kind: str,
    status: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    original = da.Optimizer.solve_report
    calls = 0

    def solve(self: da.Optimizer, **kwargs: object) -> SolveReport:
        nonlocal calls
        calls += 1
        return (
            SolveReport(SolveStatus(status), None)
            if calls == 1
            else original(self, **kwargs)
        )

    monkeypatch.setattr(da.Optimizer, "solve_report", solve)
    run = da.run(request(kind), out=tmp_path / "run")
    summary = da.inspect(run, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == (2 if kind == "matrix" else 1)
    attempts = list(da.inspect(run, view="attempts", all=True).records())
    assert attempts[0].outcome == "no_candidate"
    assert attempts[0].candidate is None
    assert attempts[0].evidence["solver_status"] == status
    assert attempts[0].evidence["proof_scope"] is None
    cell = [a for a in attempts if a.cell_id == attempts[0].cell_id]
    assert [a.evidence["batch_index"] for a in cell] == [1, 2]
    assert cell[-1].evidence["solver_status"] == "optimal"
    da.export(run, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == summary.accepted
    policy = request("resampling").resampling
    strict = BatchAccounting(replace(policy, on_unproven="stop"))
    strict.observe(cell[0])
    with pytest.raises(ValueError, match="batch order"):
        strict.observe(cell[1])


@pytest.mark.parametrize(
    "status,expected",
    [("unknown", 2), ("unproven", 2), ("backend_error", 1), ("invalid_result", 1)],
)
def test_unproven_continuation_never_hides_errors_or_replenishes_limits(
    status: str,
    expected: int,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        da.Optimizer,
        "solve_report",
        lambda *_a, **_kw: SolveReport(SolveStatus(status), None),
    )
    run = da.run(request("resampling"), out=tmp_path / "run")
    summary = da.inspect(run, verify=True)
    assert summary.accepted == 0
    assert summary.counts["started"] == expected
    assert summary.termination_reason == (
        "batch_limit" if expected == 2 else f"solver_{status}"
    )


@pytest.mark.parametrize("kind", ["schedule", "resampling", "matrix"])
def test_resume_after_committed_unproven_outcome_advances_once(
    kind: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> None:
        publish(self, *args, **kwargs)
        raise KeyboardInterrupt

    with monkeypatch.context() as patch:
        patch.setattr(
            da.Optimizer,
            "solve_report",
            lambda *_a, **_kw: SolveReport(SolveStatus.UNKNOWN, None),
        )
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request(kind), out=tmp_path / "run")
    before = list(da.inspect(tmp_path / "run", view="attempts", all=True).records())
    assert da.inspect(tmp_path / "run", verify=True).resumable
    run = da.run(resume=tmp_path / "run")
    assert da.inspect(run, verify=True).state == "completed"
    after = list(da.inspect(run, view="attempts", all=True).records())
    assert after[:1] == before
    cell = [a for a in after if a.cell_id == before[0].cell_id]
    assert [a.evidence["batch_index"] for a in cell] == [1, 2]


def test_unproven_policy_is_explicit_and_default_wire_remains_unchanged():
    default = planning.Resampling(planning.BatchSampling(1), 2, 10)
    assert "on_unproven" not in default.to_dict()
    assert planning.Resampling.from_dict(default.to_dict()).on_unproven == "stop"
    with pytest.raises(ValueError, match="on_unproven"):
        replace(default, on_unproven="accept")
    opted = replace(default, on_unproven="next_batch")
    assert planning.Resampling.from_dict(opted.to_dict()) == opted


def test_unproven_preview_and_global_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    spec = request("resampling").with_changes(limits=planning.Limits(attempts=1))
    da.export(spec, out=tmp_path / "request.json")
    response = CliRunner().invoke(app, ["plan", str(tmp_path / "request.json")])
    assert response.exit_code == 0, response.output
    assert "Unproven search: next batch" in response.stdout
    monkeypatch.setattr(
        da.Optimizer,
        "solve_report",
        lambda *_a, **_kw: SolveReport(SolveStatus.UNKNOWN, None),
    )
    run = da.run(spec, out=tmp_path / "run")
    summary = da.inspect(run, verify=True)
    assert summary.termination_reason == "attempt_limit"
    assert summary.counts["started"] == 1
    assert summary.batch_count == 1
