"""Greedy workflow admission, honest evidence and finite offered-batch search.

Author: Eric J. South.
"""

import json
import sqlite3
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.records import Attempt, RunHandle
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app


def request():
    return planning.DesignSpec(
        [parts.Part("a", "AAA"), parts.Part("b", "CCC")],
        planning.Length(maximum=6),
        strands="single",
        search="greedy",
    )


def attempts(path: Path | RunHandle):
    with da.inspect(path, view="attempts", all=True).records() as records:
        return list(records)


def test_greedy_python_and_cli_have_no_solver_or_optimality_claim(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def forbidden(*_a: object, **_kw: object) -> None:
        msg = "greedy must not construct a solver"
        raise AssertionError(msg)

    monkeypatch.setattr(da.Optimizer, "build_model", forbidden)
    plan = da.plan(request())
    assert plan.preview["solver"] is None
    assert plan.preview["path_variables"] == 0
    saved = tmp_path / "plan.json"
    plan.write(saved)
    result = da.run(plan, out=tmp_path / "python")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 1
    assert summary.producer.solver is None
    record = attempts(result)[0]
    assert record.evidence["heuristic"]["status"] == "candidate"
    assert "solver_status" not in record.evidence
    assert record.evidence["proof_scope"] is None
    assert len(record.candidate.packed.placements) == 2
    cli = CliRunner().invoke(
        app, ["run", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert cli.exit_code == 0, cli.output
    assert (
        attempts(tmp_path / "cli")[0].evidence["heuristic"]
        == record.evidence["heuristic"]
    )


def test_greedy_exhaustion_is_neither_infeasibility_nor_optimal_enumeration(
    tmp_path: Path,
):
    result = da.run(
        request().with_changes(target=planning.Target(2)), out=tmp_path / "run"
    )
    report = da.inspect(result, verify=True)
    assert report.accepted == 1
    assert report.state == "stopped"
    assert report.termination_reason == "heuristic_exhausted"
    assert attempts(result)[-1].evidence["heuristic"]["status"] == "exhausted"


@pytest.mark.parametrize(
    "constraint", ["counts", "fixed", "spacing", "coverage", "preference", "exact"]
)
def test_greedy_rejects_unsupported_semantics_before_execution(constraint: str):
    changes = {
        "counts": {
            "requirements": [
                planning.Occurrences(
                    "count", parts.PartSelector(part_ids=("a",)), min=1
                )
            ]
        },
        "fixed": {"requirements": [planning.Fixed("fixed", "a", "forward")]},
        "spacing": {
            "requirements": [planning.Spacing("space", "a", "b", min=0, max=0)]
        },
        "coverage": {
            "requirements": [planning.GroupCoverage("coverage", ("R",), min=1)]
        },
        "preference": {"packing_preference": "underused_parts"},
        "exact": {"length": planning.Length(exact=6), "assembly": planning.Assembly()},
    }[constraint]
    with pytest.raises(ValueError, match="greedy"):
        da.plan(request().with_changes(**changes))


def test_greedy_keeps_padding_and_final_screening(tmp_path: Path):
    spec = request().with_changes(
        length=planning.Length(exact=8),
        assembly=planning.Assembly(planning.Padding("right", 2)),
        requirements=[planning.Avoid("forbidden", ("AAA",))],
    )
    result = da.run(spec, out=tmp_path / "run")
    report = da.inspect(result, verify=True)
    assert report.accepted == 0
    assert report.counts["rejected"] == 1
    assert report.termination_reason == "heuristic_exhausted"
    assert len(attempts(result)[0].candidate.final.sequence) == 8


@pytest.mark.parametrize("mode", ["single", "schedule", "resampling", "matrix"])
def test_greedy_recovery_preserves_search_consumption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
):
    spec = request().with_changes(target=planning.Target(2))
    if mode == "schedule":
        base = da.plan(spec)
        batches = tuple(
            planning.CandidateBatch((p,), base.collection_id, stream=p)
            for p in ("a", "b")
        )
        spec = spec.with_changes(
            schedule=planning.BatchSchedule(batches, attempts_per_batch=3)
        )
    elif mode == "resampling":
        spec = spec.with_changes(
            resampling=planning.Resampling(
                planning.BatchSampling(1, seed=7), max_batches=10, attempts_per_batch=3
            )
        )
    elif mode == "matrix":
        spec = planning.MatrixSpec(
            spec.with_changes(target=planning.Target()),
            axes={
                "setting": {"first": planning.Variant(), "second": planning.Variant()}
            },
            allocation=planning.Allocation(per_cell=1),
            max_cells=2,
        )
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> None:
        publish(self, *args, **kwargs)
        raise KeyboardInterrupt

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(spec, out=tmp_path / "resumed")
    prefix = attempts(tmp_path / "resumed")[0].to_dict()
    resumed = da.run(resume=tmp_path / "resumed")
    direct = da.run(spec, out=tmp_path / "direct")
    a, b = attempts(resumed), attempts(direct)
    assert a[0].to_dict() == prefix
    assert [r.evidence["heuristic"] for r in a] == [r.evidence["heuristic"] for r in b]
    assert da.inspect(resumed, verify=True).accepted == (1 if mode == "single" else 2)


@pytest.mark.parametrize("method", ["exact", "greedy"])
def test_verification_binds_search_evidence_to_declared_method(
    tmp_path: Path, method: str
):
    result = da.run(request().with_changes(search=method), out=tmp_path / "run")
    with sqlite3.connect(result.path / "run.sqlite3") as connection:
        revision, payload = connection.execute(
            "SELECT revision,payload FROM attempts WHERE attempt=1 "
            "ORDER BY revision DESC LIMIT 1"
        ).fetchone()
        data = json.loads(payload)
        evidence = data["evidence"]
        if method == "exact":
            evidence.pop("solver_status")
            evidence.pop("backend_status")
            evidence["proof_scope"] = None
            evidence["termination_reason"] = "heuristic_candidate"
            evidence["heuristic"] = {
                "schema": "dense_arrays.heuristic.v1",
                "method": "greedy_multistart.v1",
                "status": "candidate",
            }
        else:
            evidence.pop("heuristic")
            evidence.update(
                solver_status="optimal",
                backend_status=0,
                proof_scope="offered_packing_model",
                termination_reason="unknown",
            )
        connection.execute(
            "UPDATE attempts SET payload=?,digest=? WHERE attempt=1 AND revision=?",
            (canonical_json(data), semantic_digest(data), revision),
        )
    with pytest.raises(ValueError, match="search method"):
        da.inspect(result, verify=True)


def test_verification_rejects_repeated_greedy_proposals_within_batch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.generation.heuristic import GreedySearch  # noqa: PLC0415

    monkeypatch.setattr(GreedySearch, "forbid", lambda *_a: None)
    spec = request().with_changes(
        target=planning.Target(2), limits=planning.Limits(attempts=3)
    )
    result = da.run(spec, out=tmp_path / "run")
    with pytest.raises(ValueError, match="one greedy proposal"):
        da.inspect(result, verify=True)


def test_greedy_time_stop_is_explicit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    from dense_arrays.generation import heuristic  # noqa: PLC0415

    def stopped(*_a: object, **_kw: object) -> None:
        raise TimeoutError

    monkeypatch.setattr(heuristic, "realize_greedy", stopped)
    result = da.run(request(), out=tmp_path / "run")
    report = da.inspect(result, verify=True)
    assert report.accepted == 0
    assert report.termination_reason == "heuristic_time_limit"
    assert report.producer.solver is None


def test_cli_plan_identifies_greedy_search(tmp_path: Path):
    saved = tmp_path / "plan.json"
    da.plan(request()).write(saved)
    for args in (["plan", str(saved)], ["inspect", str(saved), "--view", "plan"]):
        result = CliRunner().invoke(app, args)
        assert result.exit_code == 0, result.output
        assert "greedy" in result.output
        assert "unproven" in result.output
        assert "CBC" not in result.output


@pytest.mark.parametrize(
    "change", ["missing_candidate", "solver", "proof", "termination"]
)
def test_greedy_attempt_requires_coherent_method_evidence(change: str):
    evidence = {
        "heuristic": {
            "schema": "dense_arrays.heuristic.v1",
            "method": "greedy_multistart.v1",
            "status": "exhausted",
        },
        "proof_scope": None,
        "termination_reason": "heuristic_exhausted",
    }
    if change == "missing_candidate":
        evidence["heuristic"]["status"] = "candidate"
        evidence["termination_reason"] = "heuristic_candidate"
    elif change == "solver":
        evidence["solver_status"] = "infeasible"
    elif change == "proof":
        evidence["proof_scope"] = "offered_packing_model"
    else:
        evidence["termination_reason"] = "solver_infeasible"
    with pytest.raises(ValueError, match="heuristic"):
        Attempt(1, "default", "no_candidate", evidence)


def test_resume_after_committed_greedy_exhaustion_advances_schedule(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    base = da.plan(request().with_changes(target=planning.Target(2)))
    spec = base.request.with_changes(
        schedule=planning.BatchSchedule(
            tuple(
                planning.CandidateBatch((p,), base.collection_id, stream=p)
                for p in ("a", "b")
            ),
            attempts_per_batch=3,
        )
    )
    publish = RunWriter.publish

    def interrupt(
        self: RunWriter, attempt: int, *args: object, **kwargs: object
    ) -> None:
        publish(self, attempt, *args, **kwargs)
        if attempt == 2:
            raise KeyboardInterrupt

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(spec, out=tmp_path / "run")
    result = da.run(resume=tmp_path / "run")
    assert da.inspect(result, verify=True).accepted == 2
    assert [r.evidence["batch_index"] for r in attempts(result)] == [1, 1, 2]
