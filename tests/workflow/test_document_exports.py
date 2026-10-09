"""Portable JSON handoffs keep editable inputs separate from report envelopes.

Author: Eric J. South.
"""

import io
import json
import sqlite3
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.cli import app
from dense_arrays.workflow.inputs import read_source

from .test_collections import library
from .test_quality import shortfall_run


@pytest.mark.parametrize("kind", ["design", "preparation", "set", "matrix"])
def test_inspect_editable_request_before_sources_exist(tmp_path: Path, kind: str):
    source = parts.PartTable(tmp_path / "unavailable.csv", "csv")
    design = planning.DesignSpec(source, planning.Length(maximum=6))
    preparation = parts.PreparationSpec(source)
    request = {
        "design": design,
        "preparation": preparation,
        "set": parts.PreparationSet(
            {
                "one": parts.PreparationSpec(
                    parts.PWMArtifact(tmp_path / "unavailable.json"),
                    sampling=parts.Sampling(planning.Length(exact=6)),
                    budget=parts.CandidateBudget(10),
                    scoring=parts.FimoScoring(),
                    retain=parts.Retention(
                        count=2, policy="top_score", rank_by="best_hit_score"
                    ),
                )
            }
        ),
        "matrix": planning.MatrixSpec(
            design,
            {"mode": {"one": planning.Variant()}},
            planning.Allocation(per_cell=1),
            max_cells=1,
        ),
    }[kind]
    da.export(request, out=tmp_path / "request.json")
    assert da.inspect(request, view="request").request == request
    report = da.inspect(tmp_path / "request.json", view="request")
    assert report.request == request
    result = CliRunner().invoke(
        app, ["inspect", str(tmp_path / "request.json"), "--view", "request", "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == report.to_dict()
    with pytest.raises(TypeError, match="resolved"):
        da.inspect(tmp_path / "request.json", view="plan")
    assert not source.table.exists()


def test_request_inspection_counts_added_requirement_identities(tmp_path: Path):
    request = planning.MatrixSpec(
        planning.DesignSpec([parts.Part("a", "AAA")], planning.Length(maximum=3)),
        {
            "mode": {
                "one": planning.Variant(
                    add_requirements=(planning.Fixed("fixed", "a", "forward"),)
                )
            }
        },
        planning.Allocation(per_cell=1),
        max_cells=1,
    )
    da.export(request, out=tmp_path / "request.json")
    for value in (request, tmp_path / "request.json"):
        with pytest.raises(reporting.ReadLimitError, match="identities"):
            da.inspect(
                value, view="request", read_limits=reporting.ReadLimits(identities=2)
            )
        assert (
            da.inspect(
                value, view="request", read_limits=reporting.ReadLimits(identities=3)
            ).request
            == request
        )


def test_request_export_is_editable_and_plan_preserves_relative_bindings(
    tmp_path: Path,
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    run = da.run(
        planning.DesignSpec(
            parts=parts.PartTable(table, "csv"),
            length=planning.Length(maximum=6),
            strands="single",
        ),
        out=tmp_path / "run",
    )
    handoff = tmp_path / "handoff"
    receipt = da.export(run, view="request", out=handoff / "request.json")
    assert receipt.records == 1
    value = json.loads((handoff / "request.json").read_text())
    assert value["schema"] == "dense_arrays.design.v1"
    assert "records" not in value
    request = da.inspect(run, view="request").request
    assert isinstance(request, planning.DesignSpec)
    recovered = read_source(handoff / "request.json")
    assert reporting.RequestReport(recovered).to_dict(base=handoff) == value
    assert reporting.RequestReport(request).to_dict(base=handoff) == value
    original = da.inspect(run, view="plan")
    da.export(original, out=handoff / "plan.json")
    wire = json.loads((handoff / "plan.json").read_text())
    assert not Path(wire["inputs"][0]["path"]).is_absolute()
    restored = read_source(handoff / "plan.json")
    assert restored.plan_id == original.plan_id
    assert restored.inputs[0].path.resolve() == original.inputs[0].path.resolve()
    assert restored.inputs[0].sha256 == original.inputs[0].sha256
    assert da.inspect(da.run(restored, out=tmp_path / "replayed")).accepted == 1
    cli = CliRunner().invoke(
        app, ["export", str(run.path), "--view", "request", "--out", "-"]
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == reporting.RequestReport(request).to_dict(
        base=Path.cwd()
    )
    inspected = CliRunner().invoke(
        app, ["inspect", str(run.path), "--view", "request", "--json"]
    )
    assert inspected.exit_code == 0, inspected.output
    assert json.loads(inspected.stdout) == reporting.RequestReport(request).to_dict()
    assert (
        da.inspect(
            da.run(read_source(handoff / "request.json"), out=tmp_path / "edited")
        ).accepted
        == 1
    )
    table.unlink()
    with pytest.raises(FileNotFoundError):
        da.run(recovered, out=tmp_path / "missing-input")
    assert not (tmp_path / "missing-input").exists()


def test_report_exports_preserve_scope_and_keep_receipts_out_of_data(tmp_path: Path):
    run = shortfall_run(tmp_path)
    selected = reporting.DesignFilter(groups=("A",))
    report = da.inspect(run, view="quality", select=selected, limit=2)
    stream = io.StringIO()
    receipt = da.export(report, out=stream)
    assert not stream.closed
    assert receipt.view == "quality"
    assert json.loads(stream.getvalue()) == report.to_dict()
    assert receipt.sources[0]["run_id"] == run.run_id
    cli = CliRunner().invoke(
        app,
        [
            "export",
            str(run.path),
            "--view",
            "quality",
            "--group",
            "A",
            "--limit",
            "2",
            "--out",
            "-",
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == report.to_dict()
    assert "Read cost:" in cli.stderr
    assert "dense_arrays.export_receipt.v1" in cli.stderr
    for view in ("summary", "diagnostics"):
        output = tmp_path / f"{view}.json"
        da.export(run, view=view, out=output)
        assert json.loads(output.read_text()) == da.inspect(run, view=view).to_dict()
        with pytest.raises(FileExistsError):
            da.export(run, view=view, out=output)
    with pytest.raises(ValueError, match="all"):
        da.export(report, all=True, out=tmp_path / "ambiguous.json")


def test_document_exports_fail_before_publication_or_unsupported_reads(tmp_path: Path):
    with pytest.raises(ValueError, match="CSV/TSV"):
        da.export(
            tmp_path / "missing",
            view="quality",
            format="csv",
            out=tmp_path / "invalid.csv",
        )
    run = shortfall_run(tmp_path)
    output = tmp_path / "uncreated" / "report.json"
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            run, view="quality", out=output, read_limits=reporting.ReadLimits(records=1)
        )
    assert not output.exists()
    assert not output.parent.exists()


def test_extension_request_export_does_not_drop_exclusions(tmp_path: Path):
    parent = library(tmp_path / "parent")
    plan = da.plan(
        planning.ExtensionSpec(
            planning.ParentRun(parent.path), 1, planning.Limits(attempts=10), 19
        )
    )
    editable = tmp_path / "editable.json"
    da.export(plan, view="request", out=editable)
    assert da.plan(read_source(editable)).exclusions == plan.exclusions
    path = tmp_path / "extension.plan.json"
    da.export(plan, out=path)
    assert read_source(path).parent == plan.parent
    comparison = da.inspect(parent, view="plan", compare=plan)
    stream = io.StringIO()
    da.export(comparison, out=stream)
    assert json.loads(stream.getvalue()) == comparison.to_dict()
    cli = CliRunner().invoke(
        app,
        [
            "export",
            str(parent.path),
            "--view",
            "plan",
            "--compare",
            str(path),
            "--out",
            "-",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == comparison.to_dict()


def test_preparation_request_and_plan_exports_keep_source_paths_relative(
    tmp_path: Path,
):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,B\n")
    pool = da.prepare(
        parts.PreparationSpec(
            parts.PartTable(table, "csv"),
            retain=parts.Retention(select=parts.PartFilter(groups=("A",))),
        ),
        out=tmp_path / "pool",
    )
    plan = da.inspect(pool, view="plan")
    assert isinstance(plan, planning.PreparationPlan)
    for view in ("request", "plan"):
        output = tmp_path / "handoff" / f"{view}.json"
        da.export(pool, view=view, out=output)
        value = json.loads(output.read_text())
        request = value["request"] if view == "plan" else value
        assert not Path(request["source"]["table"]).is_absolute()
        restored = read_source(output)
        prepared = da.prepare(restored, out=tmp_path / f"prepared-{view}")
        assert da.inspect(prepared).retained_parts == 1
        if view == "plan":
            assert restored.plan_id == plan.plan_id


def test_editable_source_requests_export_without_opening_inputs(tmp_path: Path):
    request = planning.DesignSpec(
        parts=parts.PartTable(tmp_path / "not-yet-supplied.csv", "csv"),
        length=planning.Length(maximum=6),
    )
    output = tmp_path / "handoff" / "request.json"
    da.export(request, out=output)
    value = json.loads(output.read_text())
    assert value["parts"]["table"] == "../not-yet-supplied.csv"
    restored = read_source(output)
    assert restored.parts.table.resolve() == request.parts.table
    cli = CliRunner().invoke(
        app, ["export", str(output), "--view", "request", "--out", "-"]
    )
    assert cli.exit_code == 0, cli.output
    assert (
        Path(json.loads(cli.stdout)["parts"]["table"]).resolve() == request.parts.table
    )
    extension = planning.ExtensionSpec(
        planning.ParentRun(tmp_path / "future-parent"),
        2,
        planning.Limits(attempts=10),
        73,
    )
    extension_path = tmp_path / "handoff" / "extension.json"
    da.export(extension, out=extension_path)
    restored = read_source(extension_path)
    assert restored.parent.run.resolve() == extension.parent.run
    assert restored.additional == 2
    assert restored.seed == 73


def test_plan_exports_enforce_identity_limits_and_artifact_binding(tmp_path: Path):
    run = library(tmp_path / "run")
    resolved = da.inspect(run, view="plan")
    for source in (run, resolved):
        with pytest.raises(reporting.ReadLimitError, match="identities"):
            da.export(
                source,
                view="plan",
                read_limits=reporting.ReadLimits(identities=1),
                out=tmp_path / "uncreated" / "plan.json",
            )
    assert not (tmp_path / "uncreated").exists()
    connection = sqlite3.connect(run.path / "run.sqlite3")
    try:
        payload = json.loads(
            connection.execute(
                "SELECT payload FROM commits ORDER BY revision DESC LIMIT 1"
            ).fetchone()[0]
        )
        payload["plan_id"] = "0" * 64
        connection.execute(
            "UPDATE commits SET payload=?, digest=? WHERE revision=?",
            (canonical_json(payload), semantic_digest(payload), payload["revision"]),
        )
        connection.commit()
    finally:
        connection.close()
    with pytest.raises(ArtifactIntegrityError, match="plan identity"):
        da.export(run, view="plan", out=tmp_path / "invalid.json")
    (run.path / "pool.sqlite3").touch()
    with pytest.raises(ValueError, match="both"):
        da.export(run.path, view="plan", out=tmp_path / "ambiguous.json")
