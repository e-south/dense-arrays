"""Portable plan evidence is discoverable and individually addressable.

Author: Eric J. South.
"""

import json
import shutil
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.cli import app


def bundled_plans(tmp_path: Path, *, shared: bool = False):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA"), parts.Part("b", "CCC")],
        length=planning.Length(maximum=6),
        strands="single",
        target=planning.Target(count=1),
    )
    left = da.run(request, out=tmp_path / "left")
    right = da.run(replace(request, seed=0 if shared else 19), out=tmp_path / "right")
    expected = [da.inspect(run, view="plan").evidence for run in (left, right)]
    bundle = tmp_path / "bundle"
    da.export([left, right], all=True, format="bundle", out=bundle)
    shutil.rmtree(left.path)
    shutil.rmtree(right.path)
    return bundle, expected


def test_single_semantic_plan_is_unambiguous_across_two_origins(tmp_path: Path):
    bundle, expected = bundled_plans(tmp_path, shared=True)
    assert expected[0] == expected[1]
    assert da.inspect(bundle, view="plan") == expected[0]
    assert list(da.inspect(bundle, view="plans", all=True).records()) == [expected[0]]


def test_tampered_plan_is_rejected_before_scalar_export(tmp_path: Path):
    bundle, expected = bundled_plans(tmp_path)
    with sqlite3.connect(bundle / "bundle.sqlite3") as connection:
        value = json.loads(
            connection.execute(
                "SELECT payload FROM plans WHERE plan_id=?", (expected[0].plan_id,)
            ).fetchone()[0]
        )
        value["content"]["request"]["seed"] = 734
        connection.execute(
            "UPDATE plans SET payload=?,digest=? WHERE plan_id=?",
            (canonical_json(value), semantic_digest(value), expected[0].plan_id),
        )
    output = tmp_path / "corrupt-plan.json"
    with pytest.raises(ArtifactIntegrityError, match="digest"):
        da.export(
            bundle,
            view="plan",
            select=reporting.PlanFilter(plan_ids=(expected[0].plan_id,)),
            out=output,
        )
    assert not output.exists()


def test_bundle_plans_page_and_export_contained_evidence(tmp_path: Path):
    bundle, expected = bundled_plans(tmp_path)
    query = da.inspect(bundle, view="plans", limit=1)
    assert len(repr(query)) < 180
    assert query.cost.records_estimate == 2
    with query.records() as records:
        prefix = list(records)
        cursor = records.next_cursor
    suffix = list(da.inspect(bundle, view="plans", all=True, after=cursor).records())
    assert prefix + suffix == expected
    output = tmp_path / "plans.json"
    receipt = da.export(bundle, view="plans", all=True, out=output)
    assert receipt.records == 2
    assert json.loads(output.read_text())["records"] == [p.to_dict() for p in expected]
    response = CliRunner().invoke(
        app, ["inspect", str(bundle), "--view", "plans", "--all", "--json"]
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout)["records"] == [p.to_dict() for p in expected]
    human = CliRunner().invoke(app, ["inspect", str(bundle), "--view", "plans"])
    assert human.exit_code == 0, human.output
    assert human.stdout.count("Plan evidence ") == 2
    assert '"request":' not in human.stdout


def test_exact_plan_selection_and_comparison_match_python_and_cli(tmp_path: Path):
    bundle, expected = bundled_plans(tmp_path)
    with pytest.raises(ValueError, match="multiple plans"):
        da.inspect(bundle, view="plan")
    selected = reporting.PlanFilter(plan_ids=(expected[0].plan_id,))
    evidence = da.inspect(bundle, view="plan", select=selected)
    assert evidence == expected[0]
    assert len(repr(evidence)) < 180
    output = tmp_path / "one-plan.json"
    da.export(bundle, view="plan", select=selected, out=output)
    assert da.inspect(output, view="plan") == evidence
    difference = da.inspect(bundle, view="plan", select=selected, compare=expected[1])
    assert difference.changed_fields == ("seed",)
    other = tmp_path / "other-plan.json"
    da.export(expected[1], out=other)
    cli = CliRunner().invoke(
        app,
        [
            "inspect",
            str(bundle),
            "--view",
            "plan",
            "--plan-id",
            evidence.plan_id,
            "--compare",
            str(other),
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == difference.to_dict()
    cli_file = tmp_path / "cli-plan.json"
    result = CliRunner().invoke(
        app,
        [
            "export",
            str(bundle),
            "--view",
            "plan",
            "--plan-id",
            evidence.plan_id,
            "--out",
            str(cli_file),
        ],
    )
    assert result.exit_code == 0, result.output
    assert cli_file.read_bytes() == output.read_bytes()
    predicate = tmp_path / "filter.json"
    predicate.write_text(json.dumps(selected.to_dict()))
    assert reporting.PlanFilter.from_dict(json.loads(predicate.read_text())) == selected
    filtered = CliRunner().invoke(
        app,
        [
            "inspect",
            str(bundle),
            "--view",
            "plans",
            "--selection",
            str(predicate),
            "--json",
        ],
    )
    assert filtered.exit_code == 0, filtered.output
    assert json.loads(filtered.stdout)["records"] == [evidence.to_dict()]
    comparison_file = tmp_path / "comparison.json"
    da.export(
        bundle, view="plan", select=selected, compare=expected[1], out=comparison_file
    )
    assert json.loads(comparison_file.read_text()) == difference.to_dict()


def test_plan_views_reject_wrong_filters_unknown_ids_and_budget_exhaustion(
    tmp_path: Path,
):
    bundle, expected = bundled_plans(tmp_path)
    with pytest.raises(TypeError, match="PlanFilter"):
        da.inspect(bundle, view="plans", select=reporting.DesignFilter(groups=("A",)))
    with pytest.raises(ValueError, match="unknown plan"):
        list(
            da.inspect(
                bundle, view="plans", select=reporting.PlanFilter(plan_ids=("0" * 64,))
            ).records()
        )
    with pytest.raises(ValueError, match="multiple plans"):
        da.inspect(
            bundle,
            view="plan",
            select=reporting.PlanFilter(plan_ids=tuple(p.plan_id for p in expected)),
        )
    with pytest.raises(reporting.ReadLimitError):
        list(
            da.inspect(
                bundle,
                view="plans",
                all=True,
                read_limits=reporting.ReadLimits(records=1),
            ).records()
        )
    output = tmp_path / "limited.json"
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            bundle,
            view="plans",
            all=True,
            read_limits=reporting.ReadLimits(records=1),
            out=output,
        )
    assert not output.exists()
    with pytest.raises(ValueError, match="all=True"):
        da.export(bundle, view="plans", out=tmp_path / "implicit.json")
