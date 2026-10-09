"""Selected collections carry the evidence needed for inspection after relocation.

Author: Eric J. South.
"""

import json
import os
import shutil
import sqlite3
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.cli import app

from .test_collections import library


def test_plan_evidence_preserves_identity_without_executable_locations(tmp_path: Path):
    # Identities produced by the previously qualified package, before PlanEvidence.
    original = da.plan(
        planning.DesignSpec(
            [parts.Part("a", "AAA"), parts.Part("b", "CCC")],
            planning.Length(maximum=6),
            strands="single",
            target=planning.Target(count=2),
        )
    )
    assert (
        original.plan_id
        == "ff96988f19b52f1533df319db6e0d627f53864178a6dba7963ee689f46962ac2"  # pragma: allowlist secret  # noqa: E501
    )
    assert (
        original.collection_id
        == "ae68021f7e619db5b3596da4cf6924504fe177b42fd29115c8b4dde9b78ffacb"  # pragma: allowlist secret  # noqa: E501
    )
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    resolved = da.plan(
        planning.DesignSpec(parts.PartTable(table, "csv"), planning.Length(maximum=6))
    )
    evidence = resolved.evidence
    wire = evidence.to_dict()
    assert wire["schema"] == "dense_arrays.plan_evidence.v1"
    assert str(tmp_path) not in json.dumps(wire)
    assert evidence.plan_id == resolved.plan_id
    assert evidence.collection_id == resolved.collection_id
    assert evidence.request == resolved.request
    assert planning.PlanEvidence.from_dict(wire) == evidence
    wire["content"]["request"]["target"]["count"] += 1
    with pytest.raises(ValueError, match="digest"):
        planning.PlanEvidence.from_dict(wire)


def test_selected_bundle_survives_missing_sources_and_preserves_scope(tmp_path: Path):
    parent = da.run(
        planning.DesignSpec(
            [parts.Part(b, b * 3, group=b) for b in "ACGT"],
            planning.Length(maximum=12),
            strands="single",
            target=planning.Target(count=2),
        ),
        out=tmp_path / "parent",
    )
    child = da.run(
        planning.ExtensionSpec(
            planning.ParentRun(parent.path), 2, planning.Limits(attempts=10), 19
        ),
        out=tmp_path / "child",
    )
    with da.inspect(child, view="designs", all=True).records() as designs:
        selected = reporting.DesignFilter(design_ids=(next(designs).reference,))
    source = [parent, child, parent]
    expected = [
        d.to_dict()
        for d in da.inspect(source, view="designs", all=True, select=selected).records()
    ]
    receipt = da.export(
        source,
        view="designs",
        select=selected,
        all=True,
        format="bundle",
        out=tmp_path / "exported",
    )
    assert receipt.records == 1
    shutil.rmtree(parent.path)
    shutil.rmtree(child.path)
    moved = tmp_path / "moved"
    (tmp_path / "exported").rename(moved)
    summary = da.inspect(moved, verify=True)
    assert summary.designs == 1
    assert summary.verified
    assert summary.to_dict()["scope"] == "selected_collection"
    assert len(summary.to_dict()["source_runs"]) == 2
    assert [
        d.to_dict() for d in da.inspect(moved, view="designs", all=True).records()
    ] == expected
    placements = list(da.inspect(moved, view="placements", all=True).records())
    assert len(placements) == 4
    assert {p.design_ref for p in placements} == set(receipt.design_refs)
    fasta = tmp_path / "sequences.fasta"
    da.export(moved, view="sequences", all=True, format="fasta", out=fasta)
    assert fasta.read_text().count(">") == 1
    cli = CliRunner().invoke(app, ["inspect", str(moved), "--verify", "--json"])
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == summary.to_dict()


def test_bundle_readers_filter_page_and_reexport_the_same_identities(tmp_path: Path):
    run = library(tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    query = da.inspect(bundle, view="placements", limit=1)
    with query.records() as records:
        first = list(records)
        cursor = records.next_cursor
    rest = list(da.inspect(bundle, view="placements", all=True, after=cursor).records())
    assert [p.to_dict() for p in first + rest] == [
        p.to_dict() for p in da.inspect(run, view="placements", all=True).records()
    ]
    selected = reporting.DesignFilter(design_ids=(first[0].design_ref,))
    assert [
        d.to_dict()
        for d in da.inspect(bundle, view="designs", all=True, select=selected).records()
    ] == [
        d.to_dict()
        for d in da.inspect(
            [bundle, run], view="designs", all=True, select=selected
        ).records()
    ]
    child = tmp_path / "smaller"
    da.export(bundle, all=True, select=selected, format="bundle", out=child)
    assert da.inspect(child, verify=True).designs == 1
    with pytest.raises(ValueError, match="unknown"):
        list(
            da.inspect(
                bundle,
                view="designs",
                select=reporting.DesignFilter(groups=("absent",)),
            ).records()
        )
    with pytest.raises(reporting.ReadLimitError):
        da.inspect(bundle, verify=True, read_limits=reporting.ReadLimits(records=1))
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            run,
            all=True,
            format="bundle",
            read_limits=reporting.ReadLimits(records=1),
            out=tmp_path / "limited",
        )
    assert not (tmp_path / "limited").exists()


def test_bundle_rejects_invalid_evidence_and_preserves_collision_owner(tmp_path: Path):
    run = library(tmp_path / "run")
    destination = tmp_path / "bundle"

    def publish() -> da.artifacts.ExportReceipt | None:
        try:
            return da.export(run, all=True, format="bundle", out=destination)
        except FileExistsError:
            return None

    with ThreadPoolExecutor(max_workers=2) as workers:
        receipts = list(workers.map(lambda _: publish(), range(2)))
    assert sum(r is not None for r in receipts) == 1
    assert da.inspect(destination, verify=True).designs == 2
    with sqlite3.connect(destination / "bundle.sqlite3") as connection:
        connection.execute("CREATE TABLE tampered(value TEXT)")
    with pytest.raises(ArtifactIntegrityError, match="checksum"):
        da.inspect(destination, verify=True)
    with sqlite3.connect(run.path / "run.sqlite3") as connection:
        value = json.loads(
            connection.execute(
                "SELECT payload FROM designs WHERE ordinal=1"
            ).fetchone()[0]
        )
        value["requirements"] = [{"id": "fabricated", "observed": 1, "passed": True}]
        connection.execute(
            "UPDATE designs SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(value), semantic_digest(value)),
        )
    with pytest.raises(ValueError, match="requirement"):
        da.export(run, all=True, format="bundle", out=tmp_path / "invalid")
    assert not (tmp_path / "invalid").exists()


def test_bundle_uses_resolved_evidence_and_allows_empty_selections(tmp_path: Path):
    table = tmp_path / "original" / "parts.csv"
    table.parent.mkdir()
    table.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    run = da.run(
        planning.DesignSpec(parts.PartTable(table, "csv"), planning.Length(maximum=6)),
        out=tmp_path / "run",
    )
    shutil.rmtree(table.parent)
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    assert str(table).encode() not in (bundle / "bundle.sqlite3").read_bytes()
    assert da.inspect(bundle, verify=True).designs == 1
    empty = tmp_path / "empty"
    da.export(
        bundle,
        all=True,
        format="bundle",
        select=reporting.DesignFilter(metrics={"length": reporting.Range(min=7)}),
        out=empty,
    )
    assert da.inspect(empty, verify=True).designs == 0
    assert list(da.inspect(empty, view="designs", all=True).records()) == []


def test_interrupted_bundle_is_visibly_incomplete_and_not_overwritten(tmp_path: Path):
    run = library(tmp_path / "run")
    output = tmp_path / "interrupted"
    script = """
import os, sys
import dense_arrays as da
import dense_arrays.reporting.exporting.bundles as exporter
def interrupted(*args, **kwargs):
    os._exit(17)
exporter.write_new = interrupted
da.export(sys.argv[1], all=True, format='bundle', out=sys.argv[2])
"""
    environment = dict(os.environ)
    environment.pop("__PYVENV_LAUNCHER__", None)
    result = subprocess.run(  # noqa: S603 - fixed crash-injection fixture and argv paths
        [sys.executable, "-c", script, str(run.path), str(output)],
        env=environment,
        capture_output=True,
        timeout=30,
        check=False,
    )
    assert result.returncode == 17, result.stderr
    with pytest.raises(ArtifactIntegrityError, match="incomplete"):
        da.inspect(output)
    with pytest.raises(FileExistsError):
        da.export(run, all=True, format="bundle", out=output)
    assert (output / ".bundle-pending").exists()
