"""Additional library generation preserves parent bytes and ancestor uniqueness.

Author: Eric J. South.
"""

import hashlib
import json
import sqlite3
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app


def parent_library(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part(name, name * 3, group=name) for name in "ACGT"],
        length=planning.Length(maximum=12),
        strands="single",
        target=planning.Target(count=12),
        limits=planning.Limits(attempts=8),
    )
    return da.run(request, out=tmp_path / "parent")


def sequences(run: RunHandle):
    return {
        row.sequence_id for row in da.inspect(run, view="sequences", all=True).records()
    }


def test_additional_generation_preserves_stopped_parent_and_ancestor_exclusions(
    tmp_path: Path,
):
    parent = parent_library(tmp_path)
    before = hashlib.sha256((parent.path / "run.sqlite3").read_bytes()).hexdigest()
    summary = da.inspect(parent, verify=True)
    assert (summary.accepted, summary.target, summary.state) == (8, 12, "stopped")
    request = planning.ExtensionSpec(
        parent=planning.ParentRun(parent.path),
        additional=4,
        limits=planning.Limits(attempts=100),
        seed=19,
    )
    resolved = da.plan(request)
    assert resolved.request.target.count == 4
    assert resolved.parent.run_id == parent.run_id
    assert resolved.parent.revision == summary.revision
    assert resolved.preview["excluded_sequences"] == 8
    assert planning.GenerationPlan.from_dict(resolved.to_dict()) == resolved
    comparison = da.inspect(parent, view="plan", compare=resolved)
    assert {"target", "limits", "seed", "exclusions", "lineage"} <= set(
        comparison.changed_fields
    )
    assert {"parts", "requirements", "length", "policies"} <= set(
        comparison.unchanged_fields
    )
    child = da.run(resolved, out=tmp_path / "child")
    child_summary = da.inspect(child, verify=True)
    assert (child_summary.accepted, child_summary.target, child_summary.state) == (
        4,
        4,
        "completed",
    )
    assert not sequences(parent).intersection(sequences(child))
    combined = list(
        da.inspect([parent, child, parent], view="sequences", all=True).records()
    )
    assert len(combined) == 12
    assert len({row.sequence_id for row in combined}) == 12
    quality = da.inspect([parent, child, parent], view="quality").to_dict()
    assert quality["attainment"] is None
    assert quality["selection"]["designs"] == 12
    assert quality["selection"]["distinct_sequences"] == 12
    assert [r["state"] for r in quality["source_runs"]] == ["stopped", "completed"]
    assert [r["attainment"]["shortfall"] for r in quality["source_runs"]] == [4, 0]
    assert quality["supply"]["eligible_parts"] == 4
    assert {
        row.collection_id
        for row in da.inspect([parent, child], view="placements", all=True).records()
    } == {resolved.collection_id}
    duplicates = [
        r
        for r in da.inspect(child, view="attempts", all=True).records()
        if r.outcome == "duplicate"
    ]
    assert duplicates
    assert all(
        r.evidence["code"] == "parent_duplicate"
        and r.evidence["matched_design_ref"].startswith(parent.run_id + "/")
        for r in duplicates
    )
    grandplan = da.plan(
        planning.ExtensionSpec(
            parent=planning.ParentRun(child.path),
            additional=4,
            limits=planning.Limits(attempts=100),
            seed=20,
        )
    )
    assert grandplan.preview["excluded_sequences"] == 12
    grandchild = da.run(grandplan, out=tmp_path / "grandchild")
    assert da.inspect(grandchild, verify=True).accepted == 4
    assert len(sequences(parent) | sequences(child) | sequences(grandchild)) == 16
    assert (
        hashlib.sha256((parent.path / "run.sqlite3").read_bytes()).hexdigest() == before
    )
    assert da.inspect(parent).state == "stopped"


def test_extension_cli_plan_round_trip_and_truthful_shortfall(tmp_path: Path):
    parent = parent_library(tmp_path)
    source = tmp_path / "extend.json"
    source.write_text(
        json.dumps(
            {
                "schema": "dense_arrays.extension.v1",
                "parent": {"run": "parent"},
                "additional": 4,
                "limits": {"attempts": 1, "active_seconds": 300, "solver_seconds": 30},
                "seed": 19,
            }
        )
    )
    plan_path = tmp_path / "extension.plan.json"
    planned = CliRunner().invoke(
        app, ["plan", str(source), "--out", str(plan_path), "--json"]
    )
    assert planned.exit_code == 0, planned.output
    assert json.loads(planned.stdout)["preview"]["excluded_sequences"] == 8
    compared = CliRunner().invoke(
        app,
        [
            "inspect",
            str(parent.path),
            "--view",
            "plan",
            "--compare",
            str(plan_path),
            "--json",
        ],
    )
    assert compared.exit_code == 0, compared.output
    assert "exclusions" in json.loads(compared.stdout)["changed_fields"]
    run = CliRunner().invoke(
        app, ["run", str(plan_path), "--out", str(tmp_path / "child"), "--json"]
    )
    assert run.exit_code == 3, run.output
    value = json.loads(run.stdout)
    assert value["target"] == 4
    assert value["accepted"] < 4
    assert value["termination_reason"] == "attempt_limit"


def test_extension_rejects_live_parent_and_requires_explicit_new_effort_and_seed(
    tmp_path: Path,
):
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        )
    )
    with (
        create_run(plan, tmp_path / "live"),
        pytest.raises(ValueError, match="terminal"),
    ):
        da.plan(
            planning.ExtensionSpec(
                parent=planning.ParentRun(tmp_path / "live"),
                additional=1,
                limits=planning.Limits(),
                seed=5,
            )
        )
    for field in ("seed", "limits"):
        args = {
            "parent": planning.ParentRun(tmp_path / "live"),
            "additional": 1,
            "limits": planning.Limits(),
            "seed": 5,
        }
        del args[field]
        with pytest.raises(TypeError):
            planning.ExtensionSpec(**args)


def test_frozen_extension_runs_without_parent_paths_and_catches_tampering(
    tmp_path: Path,
):
    parent = parent_library(tmp_path)
    resolved = da.plan(
        planning.ExtensionSpec(
            parent=planning.ParentRun(parent.path),
            additional=1,
            limits=planning.Limits(attempts=100),
            seed=4,
        )
    )
    parent.path.rename(tmp_path / "moved-parent")
    saved = planning.GenerationPlan.from_dict(resolved.to_dict())
    child = da.run(saved, out=tmp_path / "child")
    assert da.inspect(child, verify=True).accepted == 1
    value = resolved.to_dict()
    value["parent"]["exclusions"][0]["sequence_id"] = "0" * 64
    with pytest.raises(ValueError, match="digest"):
        planning.GenerationPlan.from_dict(value)


def test_extension_plan_read_caps_and_malformed_requests_fail_clearly(tmp_path: Path):
    parent = parent_library(tmp_path)
    request = planning.ExtensionSpec(
        parent=planning.ParentRun(parent.path),
        additional=1,
        limits=planning.Limits(),
        seed=1,
    )
    with pytest.raises(reporting.ReadLimitError, match="records"):
        da.plan(request, read_limits=reporting.ReadLimits(records=2))
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.plan(request, read_limits=reporting.ReadLimits(identities=1))
    for missing in ("parent", "limits", "seed"):
        value = {
            "schema": "dense_arrays.extension.v1",
            "parent": {"run": str(parent.path)},
            "additional": 1,
            "limits": {"attempts": 10},
            "seed": 1,
        }
        del value[missing]
        path = tmp_path / f"missing-{missing}.json"
        path.write_text(json.dumps(value))
        response = CliRunner().invoke(app, ["plan", str(path), "--json"])
        assert response.exit_code == 2, response.exception
        assert "invalid_input" in response.stderr


def test_extension_verification_rejects_false_duplicate_lineage(tmp_path: Path):
    parent = parent_library(tmp_path)
    child = da.run(
        planning.ExtensionSpec(
            parent=planning.ParentRun(parent.path),
            additional=1,
            limits=planning.Limits(attempts=100),
            seed=7,
        ),
        out=tmp_path / "child",
    )
    with sqlite3.connect(child.path / "run.sqlite3") as connection:
        attempt, revision, payload = connection.execute(
            "SELECT attempt,revision,payload FROM attempts "
            "WHERE json_extract(payload,'$.outcome')='duplicate' LIMIT 1"
        ).fetchone()
        value = json.loads(payload)
        value["evidence"]["matched_design_ref"] = "other/default/d0"
        connection.execute(
            "UPDATE attempts SET payload=?,digest=? WHERE attempt=? AND revision=?",
            (canonical_json(value), semantic_digest(value), attempt, revision),
        )
    with pytest.raises(ValueError, match="parent duplicate"):
        da.inspect(child, verify=True)
