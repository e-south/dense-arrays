"""Curated preparation shares strict import and preview contracts with generation.

Author: Eric J. South.
"""

import json
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest
from ortools.linear_solver import pywraplp
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.cli import app


def test_part_requests_import_without_storage_or_execution():
    program = """
import importlib.abc
import sys
class RejectRuntime(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.startswith(("dense_arrays.artifacts", "dense_arrays.workflow")):
            raise RuntimeError(f"domain imported runtime owner: {fullname}")
sys.meta_path.insert(0, RejectRuntime())
from dense_arrays import parts
assert parts.Part("a", "AAA").sequence == "AAA"
assert str(parts.PoolSource("pool").path) == "pool"
"""
    completed = subprocess.run(  # noqa: S603 - current interpreter and fixed test program
        [sys.executable, "-c", program], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize("source", ["", "source: null\n", "source: []\n"])
def test_malformed_preparation_source_has_a_field_diagnostic(
    tmp_path: Path, source: str
):
    request = tmp_path / "prepare.yaml"
    request.write_text("schema: dense_arrays.prepare.v1\n" + source)
    result = CliRunner().invoke(app, ["plan", str(request)])
    assert result.exit_code == 2, result.output
    assert "source" in result.stderr
    assert "Traceback" not in result.stderr


def test_pool_verification_rejects_changed_retained_ordinal(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\n")
    pool = da.prepare(
        parts.PreparationSpec(source=parts.PartTable(table, "csv")),
        out=tmp_path / "pool",
    )
    connection = sqlite3.connect(pool.path / "pool.sqlite3")
    try:
        payload = json.loads(
            connection.execute("SELECT payload FROM parts").fetchone()[0]
        )
        payload["ordinal"] = 2
        with connection:
            connection.execute(
                "UPDATE parts SET ordinal=2,payload=?,digest=?",
                (canonical_json(payload), semantic_digest(payload)),
            )
    finally:
        connection.close()
    with pytest.raises(ValueError, match=r"ordinal|retained"):
        da.inspect(pool, verify=True)


def test_preparation_preview_filters_without_creating_output_or_solver(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCCC,A\nc,GGG,B\n")
    request = parts.PreparationSpec(
        source=parts.PartTable(table, "csv"),
        retain=parts.Retention(
            select=parts.PartFilter(
                groups=("A",), metrics={"length": reporting.Range(max=3)}
            )
        ),
    )

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("preparation preview allocated a solver")

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", forbidden)
    plan = da.plan(request)
    assert plan.preview["source_parts"] == 3
    assert plan.preview["retained_parts"] == 1
    assert plan.preview["required_tools"] == ()
    assert len(repr(plan)) < 300
    assert list(tmp_path.iterdir()) == [table]
    saved = tmp_path / "prepare.plan.json"
    plan.write(saved)

    assert (
        planning.PreparationPlan.from_dict(json.loads(saved.read_text())).plan_id
        == plan.plan_id
    )


def test_part_filter_rejects_unknown_identity_and_unavailable_metrics(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\n")
    for selection in (
        parts.PartFilter(groups=("unknown",)),
        parts.PartFilter(metrics={"fimo_score": reporting.Range(min=0)}),
    ):
        with pytest.raises(ValueError, match=r"unknown|unavailable"):
            da.plan(
                parts.PreparationSpec(
                    source=parts.PartTable(table, "csv"),
                    retain=parts.Retention(select=selection),
                )
            )


def test_prepared_pool_survives_source_removal_move_and_filtered_reuse(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,B\nc,GGG,A\n")
    pool = da.prepare(
        parts.PreparationSpec(source=parts.PartTable(table, "csv")),
        out=tmp_path / "pool",
    )
    assert da.inspect(pool, verify=True).retained_parts == 3
    view = da.inspect(
        pool, view="parts", select=parts.PartFilter(groups=("A",)), all=True
    )
    with view.records() as rows:
        assert [row.part_id for row in rows] == ["a", "c"]
    with view.records() as rows:
        assert next(rows).part_id == "a"
    table.unlink()
    moved = tmp_path / "moved-pool"
    shutil.move(pool.path, moved)
    assert da.inspect(moved, verify=True).pool_id == pool.pool_id
    design = planning.DesignSpec(
        parts=parts.PoolSource(pool=moved, select=parts.PartFilter(groups=("B",))),
        length=planning.Length(maximum=3),
        strands="single",
    )
    run = da.run(design, out=tmp_path / "run")
    assert da.inspect(run, verify=True).accepted == 1
    with da.inspect(run, view="designs").records() as rows:
        record = next(rows)
    assert record.realized.sequence == "CCC"
    assert record.realized.placements[0].feature_id == "b"


def test_curated_prepare_cli_and_python_share_the_resolved_plan(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,B\n")
    source = tmp_path / "prepare.yaml"
    source.write_text(
        "schema: dense_arrays.prepare.v1\n"
        "source: {kind: table, table: parts.csv, format: csv}\n"
        "retain: {select: {groups: [A]}}\n"
    )
    runner = CliRunner()
    plan_file = tmp_path / "prepare.plan.json"
    preview = runner.invoke(
        app, ["plan", str(source), "--out", str(plan_file), "--json"]
    )
    assert preview.exit_code == 0, preview.output
    assert json.loads(preview.stdout)["preview"]["retained_parts"] == 1
    pool = tmp_path / "cli-pool"
    made = runner.invoke(app, ["prepare", str(plan_file), "--out", str(pool), "--json"])
    assert made.exit_code == 0, made.output
    report = da.inspect(pool, verify=True)
    assert report.plan_id == json.loads(preview.stdout)["plan_id"]
    rows = runner.invoke(
        app,
        ["inspect", str(pool), "--view", "parts", "--all", "--group", "A", "--json"],
    )
    assert rows.exit_code == 0, rows.output
    assert json.loads(rows.stdout)["records"][0]["part"]["part_id"] == "a"


def test_declared_part_filter_matches_python_and_rejects_mixed_flags(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCCC,A\nc,GGG,B\n")
    pool = da.prepare(
        parts.PreparationSpec(source=parts.PartTable(table, "csv")),
        out=tmp_path / "pool",
    )
    selected = parts.PartFilter(
        groups=("A",), metrics={"length": reporting.Range(max=3)}
    )
    predicate = tmp_path / "filter.json"
    predicate.write_text(json.dumps(selected.to_dict()))
    runner = CliRunner()
    args = ["inspect", str(pool.path), "--view", "parts", "--selection", str(predicate)]
    result = runner.invoke(app, [*args, "--all", "--json"])
    assert result.exit_code == 0, result.output
    assert [r["part"]["part_id"] for r in json.loads(result.stdout)["records"]] == ["a"]
    conflict = runner.invoke(app, [*args, "--group", "A", "--json"])
    assert conflict.exit_code == 2
    assert "exclusive" in conflict.stderr


def test_pool_integrity_failure_has_exit_four_and_a_json_error(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\n")
    pool = da.prepare(
        parts.PreparationSpec(source=parts.PartTable(table, "csv")),
        out=tmp_path / "pool",
    )
    connection = sqlite3.connect(pool.path / "pool.sqlite3")
    try:
        with connection:
            connection.execute("UPDATE manifest SET digest='invalid'")
    finally:
        connection.close()
    result = CliRunner().invoke(app, ["inspect", str(pool.path), "--verify", "--json"])
    assert result.exit_code == 4, result.output
    assert json.loads(result.stdout)["code"] == "artifact_integrity"
    assert "checksum" in result.stderr


def test_preparation_rejects_stale_input_and_wrong_executor_before_output(
    tmp_path: Path,
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\n")
    resolved = da.plan(parts.PreparationSpec(source=parts.PartTable(table, "csv")))
    with pytest.raises(TypeError, match="use prepare"):
        da.run(resolved, out=tmp_path / "wrong-run")
    generation = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        )
    )
    with pytest.raises(TypeError, match="use run"):
        da.prepare(generation, out=tmp_path / "wrong-pool")
    table.write_text("part_id,sequence\na,CCC\n")
    with pytest.raises(ValueError, match="changed"):
        da.prepare(resolved, out=tmp_path / "stale-pool")
    assert list(tmp_path.iterdir()) == [table]


def test_pool_summary_uses_only_manifest_and_never_loads_parts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    pool = da.prepare(
        parts.PreparationSpec(source=parts.PartTable(table, "csv")),
        out=tmp_path / "pool",
    )
    connect = sqlite3.connect
    statements = []

    def instrumented(*args: object, **kwargs: object) -> sqlite3.Connection:
        connection = connect(*args, **kwargs)
        connection.set_trace_callback(statements.append)
        return connection

    monkeypatch.setattr(sqlite3, "connect", instrumented)
    assert da.inspect(pool).retained_parts == 2
    assert any("FROM manifest" in query for query in statements)
    assert not any(
        "FROM parts" in query or "FROM preparation" in query for query in statements
    )


@pytest.mark.parametrize("field", ["schema", "policies", "unknown_field"])
def test_preparation_plan_rejects_unsupported_wire_contracts(
    tmp_path: Path, field: str
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\n")
    plan = da.plan(parts.PreparationSpec(source=parts.PartTable(table, "csv")))
    wire = plan.to_dict()
    wire[field] = "future"
    with pytest.raises(ValueError, match=r"unsupported|unknown"):
        planning.PreparationPlan.from_dict(wire)


def test_saved_preparation_plan_has_human_preview_and_collision_preflight(
    tmp_path: Path,
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\n")
    resolved = da.plan(parts.PreparationSpec(source=parts.PartTable(table, "csv")))
    saved = tmp_path / "plan.json"
    resolved.write(saved)
    runner = CliRunner()
    preview = runner.invoke(app, ["plan", str(saved)])
    assert preview.exit_code == 0, preview.output
    assert "1 / 1 supplied parts retained" in preview.stdout
    destination = tmp_path / "existing"
    destination.mkdir()
    made = runner.invoke(
        app, ["prepare", str(tmp_path / "absent.yaml"), "--out", str(destination)]
    )
    assert made.exit_code == 2
    assert "already exists" in made.stderr
