"""Verification budgets cover the complete declared snapshot, not a displayed page.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app


def simple_run(tmp_path: Path):
    return da.run(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")],
            length=planning.Length(maximum=3),
            strands="single",
        ),
        out=tmp_path / "run",
    )


def test_run_verification_charges_plan_attempts_and_designs_to_one_budget(
    tmp_path: Path,
):
    run = simple_run(tmp_path)
    preview = da.inspect(run, read_limits=reporting.ReadLimits(records=3))
    assert preview.cost.mode == "manifest"
    assert preview.verification_cost.records_estimate == 3
    assert preview.verification_cost.limits.records == 3
    verified = da.inspect(run, verify=True, read_limits=reporting.ReadLimits(records=3))
    assert verified.verified
    assert verified.verification.records_checked == 3
    assert verified.verification.bytes_checked > 0
    assert verified.verification.boundary == "committed_plan_attempts_designs"
    with pytest.raises(reporting.ReadLimitError, match="records"):
        da.inspect(run, verify=True, read_limits=reporting.ReadLimits(records=2))
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(run, verify=True, read_limits=reporting.ReadLimits(identities=1))
    response = CliRunner().invoke(
        app, ["inspect", str(run.path), "--verify", "--max-read-records", "2", "--json"]
    )
    assert response.exit_code == 4
    assert json.loads(response.stdout)["code"] == "read_limit"
    assert response.stderr.index("Read cost:") < response.stderr.index("read_limit")


def test_pool_verification_counts_preparation_and_retained_records(tmp_path: Path):
    source = tmp_path / "parts.csv"
    source.write_text("part_id,sequence\na,AAA\nb,CCC\n")
    pool = da.prepare(
        parts.PreparationSpec(parts.PartTable(source, "csv")), out=tmp_path / "pool"
    )
    preview = da.inspect(pool, read_limits=reporting.ReadLimits(records=3))
    assert preview.verification_cost.records_estimate == 3
    result = da.inspect(pool, verify=True, read_limits=reporting.ReadLimits(records=3))
    assert result.verification.records_checked == 3
    assert result.verification.boundary == "preparation_and_retained_parts"
    with pytest.raises(reporting.ReadLimitError, match="records"):
        da.inspect(pool, verify=True, read_limits=reporting.ReadLimits(records=2))
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(pool, verify=True, read_limits=reporting.ReadLimits(identities=1))


def test_revision_bound_handle_keeps_preview_and_verification_on_same_prefix(
    tmp_path: Path,
):
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        )
    )
    with create_run(plan, tmp_path / "run") as writer:
        first = da.inspect(writer.handle)
        bound = RunHandle(
            writer.handle.path, writer.handle.run_id, revision=first.revision
        )
        attempt = writer.reserve(active_seconds=0)
        writer.publish(attempt, "rejected", {}, active_seconds=0)
        assert da.inspect(writer.handle).revision > first.revision
        checked = da.inspect(
            bound, verify=True, read_limits=reporting.ReadLimits(records=1)
        )
        assert checked.revision == first.revision
        assert checked.counts["started"] == 0
        assert checked.verification.records_checked == 1
