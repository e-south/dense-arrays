"""CLI and Python converge before validation and preserve stream/exit semantics.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app


def test_yaml_plan_matches_python_and_saved_plan_runs(tmp_path: Path):
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,B\n")
    source = tmp_path / "design.yaml"
    source.write_text(
        "schema: dense_arrays.design.v1\n"
        "parts: {table: parts.csv, format: csv}\n"
        "length: {maximum: 6}\n"
        "requirements:\n"
        "  - {id: groups, kind: group_coverage, groups: [A, B], min: 2}\n"
        "target: {count: 1}\n"
    )
    expected = da.plan(
        planning.DesignSpec(
            parts=parts.PartTable(table, "csv"),
            length=planning.Length(maximum=6),
            requirements=[
                planning.GroupCoverage(id="groups", groups=("A", "B"), min=2)
            ],
        )
    )
    runner = CliRunner()
    preview = runner.invoke(app, ["plan", str(source), "--json"])
    assert preview.exit_code == 0, preview.output
    assert json.loads(preview.stdout)["plan_id"] == expected.plan_id
    saved = tmp_path / "plan.json"
    assert runner.invoke(app, ["plan", str(source), "--out", str(saved)]).exit_code == 0
    out = tmp_path / "run"
    generated = runner.invoke(app, ["run", str(saved), "--out", str(out), "--json"])
    assert generated.exit_code == 0, generated.output
    assert json.loads(generated.stdout)["accepted"] == 1
    inspected = runner.invoke(app, ["inspect", str(out), "--verify", "--json"])
    assert inspected.exit_code == 0, inspected.output
    assert json.loads(inspected.stdout)["verified"] is True


def test_inline_shortfall_uses_exit_three_and_clean_json(tmp_path: Path):
    result = CliRunner().invoke(
        app,
        [
            "run",
            "--motif",
            "AAA",
            "--length",
            "3",
            "--count",
            "3",
            "--strands",
            "single",
            "--out",
            str(tmp_path / "run"),
            "--json",
        ],
    )
    assert result.exit_code == 3, result.output
    report = json.loads(result.stdout)
    assert report["accepted"] == 1
    assert report["state"] == "stopped"
    assert "Traceback" not in result.stderr


@pytest.mark.parametrize(
    "extra", ["target: {count: 1, count: 2}\n", "unknown_field: 3\n"]
)
def test_malformed_request_fails_without_output(tmp_path: Path, extra: str):
    source = tmp_path / "request.yaml"
    source.write_text(
        "schema: dense_arrays.design.v1\n"
        "parts: [{part_id: a, sequence: AAA}]\nlength: {maximum: 3}\n" + extra
    )
    out = tmp_path / "run"
    result = CliRunner().invoke(app, ["run", str(source), "--out", str(out), "--json"])
    assert result.exit_code == 2
    report = json.loads(result.stdout)
    assert report["schema"] == "dense_arrays.error.v1"
    assert report["code"] == "invalid_input"
    assert report["exit_code"] == 2
    assert "Traceback" not in result.stderr
    assert not out.exists()


def test_exact_assembly_preview_and_cli_match_python(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "ACG")],
        length=planning.Length(exact=5),
        strands="single",
        assembly=planning.Assembly(
            padding=planning.Padding(side="right", max_trials=5)
        ),
        requirements=[planning.GC("gc", scope="sequence", min=0.2, max=1)],
    )
    expected = da.plan(request)
    source = tmp_path / "request.json"
    source.write_text(json.dumps(expected.to_dict()["request"]))
    runner = CliRunner()
    preview = runner.invoke(app, ["plan", str(source)])
    assert preview.exit_code == 0, preview.output
    assert "length exact 5" in preview.stdout
    assert "None" not in preview.stdout
    out = tmp_path / "run"
    result = runner.invoke(app, ["run", str(source), "--out", str(out), "--json"])
    assert result.exit_code == 0, result.output
    assert da.inspect(out, verify=True).plan_id == expected.plan_id
