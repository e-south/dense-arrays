"""Execution provenance survives inspection, snapshots and portable handoffs.

Author: Eric J. South.
"""

import json
import shutil
from importlib.metadata import version
from pathlib import Path

import pytest
from ortools.linear_solver import pywraplp
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import semantic_digest
from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.pool_records import PoolSummary
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.store import create_run, latest, reader
from dense_arrays.cli import app
from dense_arrays.solver import SolverIdentity
from dense_arrays.workflow import execution


def request():
    return planning.DesignSpec(
        (parts.Part("a", "AAA"),), planning.Length(maximum=3), strands="single"
    )


def test_run_records_actual_runtime_and_backend_without_rechecking_at_inspection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run = da.run(request(), out=tmp_path / "run")
    summary = da.inspect(run, verify=True)
    producer = summary.producer
    assert producer.package_version == version("dense-arrays")
    assert producer.ortools_version == version("ortools")
    assert producer.python_version
    assert producer.python_implementation
    assert producer.solver.name == "CBC"
    assert (
        producer.solver.version == pywraplp.Solver.CreateSolver("CBC").SolverVersion()
    )
    saved = producer.to_dict()
    assert saved["schema"] == "dense_arrays.producer.v1"

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("inspection attempted to discover current producer or solver")

    monkeypatch.setattr(Producer, "capture", forbidden)
    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", forbidden)
    assert da.inspect(run).producer == producer
    cli = CliRunner().invoke(app, ["inspect", str(run.path), "--json"])
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["producer"] == saved
    human = CliRunner().invoke(app, ["inspect", str(run.path)])
    assert human.exit_code == 0, human.output
    assert "Producer: Dense Arrays" in human.stdout
    assert producer.solver.version in human.stdout
    assert (
        da.inspect(run, view="quality").to_dict()["source_runs"][0]["producer"] == saved
    )
    da.export(summary, out=tmp_path / "summary.json")
    assert json.loads((tmp_path / "summary.json").read_text())["producer"] == saved


def test_failed_model_build_does_not_claim_a_backend(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def unavailable(*_args: object, **_kwargs: object) -> None:
        msg = "test backend unavailable"
        raise da.SolverBackendError(msg)

    monkeypatch.setattr(execution, "build_optimizer", unavailable)
    with pytest.raises(da.RunExecutionError) as failed:
        da.run(request(), out=tmp_path / "failed")
    assert isinstance(failed.value.__cause__, da.SolverBackendError)
    assert failed.value.run.path == tmp_path / "failed"
    summary = da.inspect(tmp_path / "failed", verify=True)
    assert summary.state == "failed"
    assert summary.producer.package_version == version("dense-arrays")
    assert summary.producer.solver is None
    assert summary.counts["started"] == 0


@pytest.mark.parametrize("error", [ValueError, OSError])
def test_cli_execution_failure_reports_committed_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error: type[Exception]
):
    def failed_build(*_args: object, **_kwargs: object) -> None:
        msg = "model construction failed"
        raise error(msg)

    monkeypatch.setattr(execution, "build_optimizer", failed_build)
    out = tmp_path / "failed"
    result = CliRunner().invoke(
        app, ["run", "--motif", "AAA", "--length", "3", "--out", str(out), "--json"]
    )
    assert result.exit_code == 4, result.output
    data = json.loads(result.stdout)
    assert data["code"] == "execution_error"
    assert data["artifact"] == str(out)
    assert "inspect" in result.stderr
    assert str(out) in result.stderr
    assert da.inspect(out, verify=True).state == "failed"


def test_preflight_failure_remains_invalid_input(tmp_path: Path):
    result = CliRunner().invoke(
        app,
        [
            "run",
            "--motif",
            "AAA",
            "--length",
            "0",
            "--out",
            str(tmp_path / "bad"),
            "--json",
        ],
    )
    assert result.exit_code == 2, result.output
    assert json.loads(result.stdout)["artifact"] is None
    assert not (tmp_path / "bad").exists()


def test_interrupted_run_names_its_actual_committed_destination(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def interrupted(*_args: object, **_kwargs: object) -> None:
        raise KeyboardInterrupt

    monkeypatch.setattr(execution, "build_optimizer", interrupted)
    out = tmp_path / "interrupted"
    result = CliRunner().invoke(
        app, ["run", "--motif", "AAA", "--length", "3", "--out", str(out)]
    )
    assert result.exit_code == 130, result.output
    assert str(out) in result.stderr
    assert "inspectable" in result.stderr
    assert da.inspect(out, verify=True).state == "stopped"


def test_pool_and_bundle_preserve_producer_with_no_invented_solver(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\n")
    pool = da.prepare(
        parts.PreparationSpec(parts.PartTable(table, "csv")), out=tmp_path / "pool"
    )
    summary = da.inspect(pool, verify=True)
    assert summary.producer.package_version == version("dense-arrays")
    assert summary.producer.solver is None
    # The producer is execution context, not part of prepared collection identity.
    legacy_data = summary.to_dict()
    legacy_data.pop("producer")
    legacy_data["schema"] = "dense_arrays.pool.v1"
    legacy = PoolSummary.from_dict(legacy_data)
    assert legacy.producer is None
    assert legacy.to_dict() == legacy_data
    plan = da.inspect(pool, view="plan")
    assert pool.pool_id == semantic_digest(
        {"schema": "dense_arrays.pool.v1", "preparation_plan_id": plan.plan_id}
    )
    run = da.run(request(), out=tmp_path / "run")
    producer = da.inspect(run).producer.to_dict()
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    shutil.rmtree(run.path)
    assert (
        da.inspect(bundle, verify=True).to_dict()["source_runs"][0]["producer"]
        == producer
    )
    assert (
        da.inspect(bundle, view="quality").to_dict()["source_runs"][0]["producer"]
        == producer
    )


def test_legacy_and_invalid_producer_records_have_explicit_compatibility(
    tmp_path: Path,
):
    run = da.run(request(), out=tmp_path / "run")
    with reader(run.path) as connection:
        manifest = latest(connection)
    assert manifest["schema"] == "dense_arrays.run.v2"
    value = manifest["producer"]
    assert Producer.from_dict(value).to_dict() == value
    with pytest.raises(ValueError, match="schema"):
        Producer.from_dict({**value, "schema": "dense_arrays.producer.v999"})
    invalid = {**value, "solver": {"name": "CBC", "version": ""}}
    with pytest.raises((TypeError, ValueError), match="version"):
        Producer.from_dict(invalid)
    incomplete = {k: v for k, v in manifest.items() if k != "producer"}
    with pytest.raises(ValueError, match="incomplete"):
        reporting.RunSummary.from_manifest(incomplete)
    legacy = reporting.RunSummary.from_manifest(
        {**incomplete, "schema": "dense_arrays.run.v1"}
    )
    assert legacy.producer is None
    assert "producer" not in legacy.to_dict()
    quality = da.inspect(run, view="quality").to_dict()
    quality["source_runs"][0]["producer"]["schema"] = "dense_arrays.producer.v999"
    with pytest.raises(ValueError, match="producer schema"):
        reporting.QualitySnapshot.from_dict(quality)


def test_solver_binding_is_explicit_and_revision_bound(tmp_path: Path):
    optimizer = da.Optimizer(["AAA"], 3)
    assert optimizer.solver_identity is None
    optimizer.build_model()
    identity = optimizer.solver_identity
    with create_run(da.plan(request()), tmp_path / "run") as writer:
        original = da.inspect(writer.handle)
        old = RunHandle(writer.handle.path, writer.handle.run_id, original.revision)
        assert original.producer.solver is None
        writer.bind_solver(identity)
        bound = da.inspect(writer.handle)
        assert bound.producer.solver == identity
        assert da.inspect(old).producer.solver is None
        writer.bind_solver(identity)
        assert da.inspect(writer.handle).revision == bound.revision
        with pytest.raises(ValueError, match="cannot change"):
            writer.bind_solver(SolverIdentity("SCIP", "different"))
        assert da.inspect(writer.handle) == bound
