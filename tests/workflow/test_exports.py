"""Portable exports contain only data and preserve native design joins.

Author: Eric J. South.
"""

import csv
import io
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.cli import app


def library(tmp_path: Path):
    return da.run(
        planning.DesignSpec(
            parts=[
                parts.Part("a", "AAA", group="A"),
                parts.Part("b", "CCC", group="B"),
            ],
            length=planning.Length(maximum=6),
            strands="single",
            target=planning.Target(count=2),
        ),
        out=tmp_path / "run",
    )


def test_sequence_placement_and_fasta_exports_share_exact_design_references(
    tmp_path: Path,
):
    run = library(tmp_path)
    csv_path = tmp_path / "sequences.csv"
    receipt = da.export(run, view="sequences", all=True, format="csv", out=csv_path)
    assert receipt.records == 2
    sequences = list(csv.DictReader(io.StringIO(csv_path.read_text())))
    assert set(receipt.design_refs) == {r["design_ref"] for r in sequences}
    assert {r["sequence"] for r in sequences} == {"AAACCC", "CCCAAA"}
    assert all(r["schema"] == "dense_arrays.sequence_record.v1" for r in sequences)
    stream = io.StringIO()
    placements_receipt = da.export(
        run, view="placements", all=True, format="tsv", out=stream
    )
    assert not stream.closed
    placements = list(csv.DictReader(io.StringIO(stream.getvalue()), delimiter="\t"))
    assert placements_receipt.records == 4
    assert {p["design_ref"] for p in placements} == set(receipt.design_refs)
    assert all(p["core_start"] == "" and p["core_end"] == "" for p in placements)
    result = CliRunner().invoke(
        app,
        [
            "export",
            str(run.path),
            "--view",
            "sequences",
            "--all",
            "--format",
            "fasta",
            "--out",
            "-",
        ],
    )
    assert result.exit_code == 0, result.output
    assert result.stdout.count(">") == 2
    assert "Read cost:" in result.stderr
    assert "Read cost:" not in result.stdout
    for row in sequences:
        assert (
            f">{row['design_ref']} sequence_id={row['sequence_id']}\n"
            f"{row['sequence']}\n" in result.stdout
        )


def test_export_requires_explicit_scope_and_publishes_files_atomically(tmp_path: Path):
    run = library(tmp_path)
    out = tmp_path / "library.json"
    with pytest.raises(ValueError, match="all"):
        da.export(run, view="designs", format="json", out=out)
    assert not out.exists()
    with pytest.raises(ValueError, match=r"FASTA|fasta"):
        da.export(run, view="placements", all=True, format="fasta", out=out)
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            run,
            view="sequences",
            all=True,
            format="json",
            out=out,
            read_limits=reporting.ReadLimits(records=1),
        )
    assert not out.exists()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["run"]
    selected = reporting.DesignFilter(metrics={"gc_fraction": reporting.Range(min=0.9)})
    receipt = da.export(
        run, view="sequences", select=selected, all=True, format="csv", out=out
    )
    assert receipt.records == 0
    assert receipt.sources[0]["record_schema"] == "dense_arrays.sequence_record.v1"
    assert out.read_text().startswith("schema,design_ref,")
    before = out.read_bytes()
    with pytest.raises(FileExistsError):
        da.export(run, view="sequences", all=True, format="csv", out=out)
    assert out.read_bytes() == before


def test_json_receipt_is_separate_and_stream_errors_use_stderr(
    tmp_path: Path,
):
    run = library(tmp_path)
    output = tmp_path / "sequence.json"
    result = CliRunner().invoke(
        app,
        [
            "export",
            str(run.path),
            "--view",
            "sequences",
            "--all",
            "--format",
            "json",
            "--out",
            str(output),
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    receipt = json.loads(result.stdout)
    assert receipt["schema"] == "dense_arrays.export_receipt.v1"
    assert receipt["records"] == 2
    data = json.loads(output.read_text())
    assert data["schema"] == "dense_arrays.record_export.v1"
    assert data["sources"][0]["run_id"] == run.run_id
    assert len(data["records"]) == 2
    partial = CliRunner().invoke(
        app,
        [
            "export",
            str(run.path),
            "--view",
            "sequences",
            "--all",
            "--format",
            "fasta",
            "--out",
            "-",
            "--max-read-records",
            "1",
            "--json",
        ],
    )
    assert partial.exit_code == 4
    assert partial.stdout.startswith(">")
    assert "dense_arrays.error" not in partial.stdout
    assert "read_limit" in partial.stderr
