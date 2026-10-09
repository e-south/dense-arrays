"""Descriptive quality differences retain populations and metric compatibility.

Author: Eric J. South.
"""

import json
import shutil
from pathlib import Path

import pytest
from ortools.linear_solver import pywraplp
from PIL import Image
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import reporting
from dense_arrays.cli import app

from .test_quality import shortfall_run


def test_quality_comparison_keeps_selected_and_search_denominators_distinct(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    run = shortfall_run(tmp_path)
    before = da.inspect(run, view="quality", limit=1)
    after = da.inspect(
        run, view="quality", select=reporting.DesignFilter(groups=("A",)), limit=1
    )

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("comparison invoked a solver")

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", forbidden)
    report = da.inspect(before, view="quality", compare=after)
    assert isinstance(report, reporting.QualityComparison)
    assert len(repr(report)) < 180
    metrics = {m.path: m for m in report.metrics}
    gc = metrics[("composition", "gc_fraction", "mean")]
    assert (gc.before, gc.after, gc.delta) == (0.5, 0.5, 0)
    assert (gc.before_denominator, gc.after_denominator) == (8, 4)
    assert gc.status == "comparable"
    selected = metrics[("selection", "designs")]
    assert (selected.before, selected.after, selected.delta) == (8, 4, -4)
    search = metrics[("search", "attempt_counts", "accepted")]
    assert (search.before, search.after, search.delta) == (8, 8, 0)
    assert (search.before_denominator, search.after_denominator) == (10, 10)
    data = report.to_dict()
    assert data["schema"] == "dense_arrays.quality_comparison.v1"
    assert data["mode"] == "descriptive"
    assert data["before"]["selected_designs"] == 8
    assert data["after"]["selected_designs"] == 4
    assert data["before"]["source_runs"][0]["attainment"]["shortfall"] == 4
    assert data["after"]["source_runs"][0]["attainment"]["shortfall"] == 4
    assert data["population_changed"]


def test_cli_quality_comparison_and_export_share_metric_records(tmp_path: Path):
    before = shortfall_run(tmp_path / "before")
    after = shortfall_run(tmp_path / "after")
    report = da.inspect(before, view="quality", compare=after)
    cli = CliRunner().invoke(
        app,
        [
            "inspect",
            str(before.path),
            "--view",
            "quality",
            "--compare",
            str(after.path),
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == report.to_dict()
    output = tmp_path / "comparison.json"
    exported = CliRunner().invoke(
        app,
        [
            "export",
            str(before.path),
            "--view",
            "quality",
            "--compare",
            str(after.path),
            "--out",
            str(output),
        ],
    )
    assert exported.exit_code == 0, exported.output
    assert json.loads(output.read_text()) == report.to_dict()


def test_saved_quality_reports_compare_after_sources_are_removed_and_mark_versions(
    tmp_path: Path,
):
    run = shortfall_run(tmp_path)
    before = da.inspect(run, view="quality", limit=1)
    snapshot = reporting.QualitySnapshot.from_report(before)
    native_metrics = da.inspect(before, view="quality", compare=before).to_dict()[
        "metrics"
    ]
    assert isinstance(snapshot, reporting.QualitySnapshot)
    assert snapshot.to_dict() == before.to_dict()
    old = tmp_path / "before.json"
    da.export(snapshot, out=old)
    value = snapshot.to_dict()
    value["policy"] = "library_composition.v999"
    newer = reporting.QualitySnapshot.from_dict(value)
    newer_path = tmp_path / "newer.json"
    da.export(newer, out=newer_path)
    shutil.rmtree(run.path)
    assert (
        da.inspect(old, view="quality", compare=old).to_dict()["metrics"]
        == native_metrics
    )
    compared = da.inspect(old, view="quality", compare=newer_path)
    assert all(m.status == "incomparable" for m in compared.metrics)
    assert all(
        m.delta is None and m.reason == "metric_policy_mismatch"
        for m in compared.metrics
    )
    assert compared.to_dict()["before"]["evidence"] == "saved_report"
    both_unknown = da.inspect(newer, view="quality", compare=newer)
    assert all(m.reason == "unsupported_metric_policy" for m in both_unknown.metrics)
    cli = CliRunner().invoke(
        app,
        [
            "inspect",
            str(old),
            "--view",
            "quality",
            "--compare",
            str(newer_path),
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == compared.to_dict()
    malformed = snapshot.to_dict()
    malformed["composition"]["gc_fraction"]["count"] = 100
    with pytest.raises(ValueError, match="population"):
        reporting.QualitySnapshot.from_dict(malformed)
    assert snapshot.to_dict()["policy"] == "library_composition.v2"


def test_quality_comparison_preserves_empty_and_missing_search_values(tmp_path: Path):
    run = shortfall_run(tmp_path / "native")
    full = da.inspect(run, view="quality")
    empty = da.inspect(
        run, view="quality", select=reporting.DesignFilter(groups=("C",))
    )
    comparison = da.inspect(full, view="quality", compare=empty)
    gc = next(
        m
        for m in comparison.metrics
        if m.path == ("composition", "gc_fraction", "mean")
    )
    assert (gc.before, gc.after, gc.delta) == (0.5, None, None)
    assert gc.status == "unavailable"
    assert gc.reason == "empty_population"
    assert gc.after_denominator == 0
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    portable = da.inspect(run, view="quality", compare=bundle)
    metric = next(
        m
        for m in portable.metrics
        if m.path == ("search", "attempt_counts", "accepted")
    )
    assert (metric.before, metric.after, metric.delta) == (8, None, None)
    assert metric.reason == "search_history_not_included"
    assert portable.to_dict()["after"]["source_runs"][0]["attainment"]["accepted"] == 8
    assert all(m.delta == 0 for m in portable.metrics if m.path[0] == "composition")
    other = shortfall_run(tmp_path / "other")
    partial = da.inspect([other, bundle], view="quality", compare=[other, run])
    metric = next(
        m for m in partial.metrics if m.path == ("search", "attempt_counts", "accepted")
    )
    assert (metric.before, metric.after, metric.delta) == (8, 16, None)
    assert metric.reason == "partial_search_history"


def test_comparison_honors_shared_and_prebound_read_caps_before_publication(
    tmp_path: Path,
):
    run = shortfall_run(tmp_path)
    output = tmp_path / "uncreated" / "comparison.json"
    limited = da.inspect(
        run, view="quality", read_limits=reporting.ReadLimits(records=1)
    )
    with pytest.raises(reporting.ReadLimitError):
        da.export(da.inspect(limited, view="quality", compare=run), out=output)
    assert not output.parent.exists()
    comparison = da.inspect(
        run, view="quality", compare=run, read_limits=reporting.ReadLimits(records=37)
    )
    assert comparison.cost.records_estimate == 38
    with pytest.raises(reporting.ReadLimitError):
        da.export(comparison, out=output)
    assert not output.parent.exists()
    # Pagination belongs to each report's usage tables, not to aggregate comparisons.
    with pytest.raises(ValueError, match="limit"):
        da.inspect(run, view="quality", compare=run, limit=1)


def test_saved_quality_inspection_and_export_preserve_native_document(tmp_path: Path):
    run = shortfall_run(tmp_path)
    report = da.inspect(run, view="quality", limit=1)
    path = tmp_path / "quality.json"
    da.export(report, out=path)
    shutil.rmtree(run.path)
    restored = da.inspect(path, view="quality")
    assert isinstance(restored, reporting.QualitySnapshot)
    assert restored.to_dict() == report.to_dict()
    cli = CliRunner().invoke(app, ["inspect", str(path), "--view", "quality", "--json"])
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == report.to_dict()
    output = tmp_path / "copy.json"
    da.export(path, view="quality", out=output)
    assert output.read_bytes() == path.read_bytes()
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(path, view="quality", read_limits=reporting.ReadLimits(identities=1))


def test_saved_metric_validation_and_scope_controls_fail_clearly(tmp_path: Path):
    report = da.inspect(shortfall_run(tmp_path), view="quality")
    malformed = report.to_dict()
    malformed["composition"]["gc_fraction"] = {}
    with pytest.raises(ValueError, match="missing"):
        reporting.QualitySnapshot.from_dict(malformed)
    malformed = report.to_dict()
    malformed["search"] = {"availability": "complete"}
    with pytest.raises(ValueError, match="missing"):
        reporting.QualitySnapshot.from_dict(malformed)
    malformed = report.to_dict()
    malformed["composition"]["made_up"] = dict(malformed["composition"]["gc_fraction"])
    with pytest.raises(ValueError, match="composition"):
        reporting.QualitySnapshot.from_dict(malformed)
    malformed = report.to_dict()
    malformed["part_usage"][0]["occurrences"] = -1
    with pytest.raises(ValueError, match="occurrences"):
        reporting.QualitySnapshot.from_dict(malformed)
    with pytest.raises(TypeError, match="saved selection"):
        da.inspect(
            report,
            view="quality",
            compare=report,
            select=reporting.LibrarySelection(take=reporting.Take(count=2)),
        )
    filtered = da.inspect(
        report.query.path, view="quality", select=reporting.DesignFilter(groups=("A",))
    )
    other = da.inspect(
        report.query.path, view="quality", select=reporting.DesignFilter(groups=("B",))
    )
    data = da.inspect(filtered, view="quality", compare=other).to_dict()
    assert data["population_changed"] is None
    assert data["scope_changed"] is True


def test_saved_quality_rendering_uses_recorded_metrics_and_rejects_unknown_policy(
    tmp_path: Path,
):
    run = shortfall_run(tmp_path)
    snapshot = reporting.QualitySnapshot.from_report(da.inspect(run, view="quality"))
    shutil.rmtree(run.path)
    output = tmp_path / "quality.png"
    receipt = da.render(snapshot, view="library-quality", out=output)
    assert receipt.records == 8
    with Image.open(output) as image:
        assert json.loads(image.info["DenseArraysReport"]) == snapshot.to_dict()
    unknown = snapshot.to_dict()
    unknown["policy"] = "library_composition.v999"
    invalid = tmp_path / "invalid.png"
    with pytest.raises(ValueError, match="metric policy"):
        da.render(
            reporting.QualitySnapshot.from_dict(unknown),
            view="library-quality",
            out=invalid,
        )
    assert not invalid.exists()


@pytest.mark.parametrize(
    "metric,value",
    [
        ("gc_fraction", -0.1),
        ("gc_fraction", 1.1),
        ("density", -0.1),
        ("density", 1.1),
        ("length", 0),
        ("length", 1.5),
        ("packed_span", -1),
        ("packed_span", 1.5),
        ("placement_count", 0),
        ("placement_count", 1.5),
        ("padding_length", -1),
        ("padding_length", 0.5),
        ("compression", 0),
    ],
)
def test_saved_quality_rejects_impossible_metric_values(
    tmp_path: Path, metric: str, value: float
):
    data = da.inspect(shortfall_run(tmp_path), view="quality").to_dict()
    count = data["selection"]["designs"]
    data["composition"][metric] = {
        "count": count,
        "denominator": "selected_designs",
        "histogram": [{"value": value, "count": count}],
        "min": value,
        "max": value,
        "mean": value,
    }
    with pytest.raises((TypeError, ValueError), match=metric):
        reporting.QualitySnapshot.from_dict(data)


def test_saved_quality_reconciles_placements_and_distinct_sequences(tmp_path: Path):
    report = da.inspect(shortfall_run(tmp_path), view="quality")
    no_sequences = report.to_dict()
    no_sequences["selection"]["distinct_sequences"] = 0
    with pytest.raises(ValueError, match="distinct"):
        reporting.QualitySnapshot.from_dict(no_sequences)
    placements = report.to_dict()
    count = placements["selection"]["designs"]
    placements["composition"]["placement_count"] = {
        "count": count,
        "denominator": "selected_designs",
        "histogram": [{"value": 2, "count": count}],
        "min": 2,
        "max": 2,
        "mean": 2,
    }
    with pytest.raises(ValueError, match=r"placement_count.*denominator"):
        reporting.QualitySnapshot.from_dict(placements)


def test_unknown_quality_policy_is_not_given_current_metric_semantics(tmp_path: Path):
    data = da.inspect(shortfall_run(tmp_path), view="quality").to_dict()
    data["policy"] = "different_metric_units.v1"
    count = data["selection"]["designs"]
    data["composition"]["gc_fraction"].update(
        min=200,
        max=200,
        mean=200,
        histogram=[{"value": 200, "count": count}],
    )
    snapshot = reporting.QualitySnapshot.from_dict(data)
    comparison = da.inspect(snapshot, view="quality", compare=snapshot)
    assert all(metric.status == "incomparable" for metric in comparison.metrics)
