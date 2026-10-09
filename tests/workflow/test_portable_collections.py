"""Native and portable evidence share one ordered collection query.

Author: Eric J. South.
"""

import io
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
from dense_arrays.playback.quality import quality_figure


def library(path: Path):
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
        out=path,
    )


def test_mixed_union_keeps_first_identity_and_pages_scalar_projections(tmp_path: Path):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    first = next(da.inspect(left, view="designs").records())
    bundle = tmp_path / "bundle"
    da.export(
        [left, right],
        all=True,
        select=reporting.DesignFilter(design_ids=(first.reference,)),
        format="bundle",
        out=bundle,
    )
    sources = [bundle, right, left, bundle]
    expected = [
        first,
        *da.inspect(right, view="designs", all=True).records(),
        *[
            d
            for d in da.inspect(left, view="designs", all=True).records()
            if d.reference != first.reference
        ],
    ]
    query = da.inspect(sources, view="designs", all=True)
    assert [d.to_dict() for d in query.records()] == [d.to_dict() for d in expected]
    assert query.sources[0]["bundle_id"] == da.inspect(bundle).bundle_id
    assert query.sources[1]["run_id"] == right.run_id
    with da.inspect(sources, view="placements", limit=3).records() as records:
        prefix = list(records)
        cursor = records.next_cursor
    suffix = list(
        da.inspect(sources, view="placements", all=True, after=cursor).records()
    )
    assert len(prefix + suffix) == 8
    with pytest.raises(ValueError, match="cursor"):
        da.inspect(list(reversed(sources)), view="placements", after=cursor)
    out = io.StringIO()
    receipt = da.export(sources, view="sequences", all=True, format="fasta", out=out)
    assert receipt.design_refs == tuple(d.reference for d in expected)
    cli = CliRunner().invoke(
        app,
        [
            "export",
            str(bundle),
            str(right.path),
            str(left.path),
            str(bundle),
            "--view",
            "sequences",
            "--all",
            "--format",
            "fasta",
            "--out",
            "-",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert cli.stdout == out.getvalue()


def test_mixed_filters_and_selection_share_namespaces_and_contained_evidence(
    tmp_path: Path,
):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    bundle = tmp_path / "bundle"
    da.export(left, all=True, format="bundle", out=bundle)
    sources = [bundle, right]
    with pytest.raises(ValueError, match="ambiguous"):
        list(
            da.inspect(
                sources,
                view="designs",
                select=reporting.DesignFilter(cells=("default",)),
            ).records()
        )
    chosen = da.inspect(
        sources,
        view="designs",
        all=True,
        select=reporting.DesignFilter(cells=(f"{left.run_id}/default",), groups=("A",)),
    )
    assert len(list(chosen.records())) == 2
    panel = da.inspect(
        sources,
        view="selection",
        select=reporting.LibrarySelection(
            take=reporting.Take(
                per_cell={f"{r.run_id}/default": 1 for r in (left, right)},
                policy="random",
                seed=23,
            )
        ),
    )
    assert panel.selected == 2
    selected = tmp_path / "selected"
    da.export(sources, select=panel, format="bundle", out=selected)
    assert da.inspect(selected, verify=True).designs == 2
    assert [
        r.reference for r in da.inspect(selected, view="designs", all=True).records()
    ] == list(panel.references())
    filter_file = tmp_path / "filter.json"
    filter_file.write_text(json.dumps(reporting.DesignFilter(groups=("A",)).to_dict()))
    cli = CliRunner().invoke(
        app,
        [
            "inspect",
            str(bundle),
            str(right.path),
            "--view",
            "designs",
            "--selection",
            str(filter_file),
            "--all",
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert len(json.loads(cli.stdout)["records"]) == 4


def test_bundle_quality_uses_included_designs_and_labels_missing_attempts(
    tmp_path: Path,
):
    run = library(tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(
        run,
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
        format="bundle",
        out=bundle,
    )
    shutil.rmtree(run.path)
    report = da.inspect(bundle, view="quality")
    value = report.to_dict()
    assert value["schema"] == "dense_arrays.quality.v3"
    assert value["selection"]["designs"] == 1
    assert value["source_runs"][0]["attainment"]["accepted"] == 2
    assert value["source_runs"][0]["included_designs"] == 1
    assert value["source_runs"][0]["search"] is None
    assert value["search"]["availability"] == "not_included"
    assert value["search"]["attempt_counts"] is None
    assert value["composition"]["length"]["count"] == 1
    figure = quality_figure(value)
    outcomes = next(axis for axis in figure.axes if axis.get_label() == "outcomes")
    assert not outcomes.axison
    assert not outcomes.patches
    assert "not included" in outcomes.texts[0].get_text()
    figure.clear()
    image = tmp_path / "quality.png"
    da.render(bundle, view="library-quality", out=image)
    assert image.is_file()
    response = CliRunner().invoke(
        app, ["inspect", str(bundle), "--view", "quality", "--json"]
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout) == value


def test_mixed_quality_counts_each_design_once_and_reports_available_histories(
    tmp_path: Path,
):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    bundle = tmp_path / "bundle"
    da.export(
        left,
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
        format="bundle",
        out=bundle,
    )
    partial = da.inspect([bundle, right], view="quality").to_dict()
    assert partial["selection"]["designs"] == 3
    assert partial["search"]["availability"] == "partial"
    assert partial["search"]["attempt_counts"]["accepted"] == 2
    assert partial["search"]["unavailable_source_refs"] == [
        f"{left.run_id}/{da.inspect(left).revision}"
    ]
    complete = da.inspect([bundle, right, left, bundle], view="quality").to_dict()
    assert complete["selection"]["designs"] == 4
    assert complete["search"]["availability"] == "complete"
    assert complete["search"]["attempt_counts"]["accepted"] == 4
    panel = da.inspect(
        [bundle, right],
        view="selection",
        select=reporting.LibrarySelection(take=reporting.Take(count=1)),
    )
    selected = da.inspect([bundle, right], view="quality", select=panel).to_dict()
    assert selected["selection"]["designs"] == 1
    assert selected["search"]["availability"] == "partial"


def test_missing_bundle_rows_cannot_be_hidden_by_other_sources(tmp_path: Path):
    run = library(tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    with sqlite3.connect(bundle / "bundle.sqlite3") as connection:
        connection.execute("DELETE FROM designs WHERE ordinal=1")
    for source in (bundle, [run, bundle]):
        with pytest.raises(ArtifactIntegrityError, match="contained designs"):
            da.inspect(source, view="quality").to_dict()
        with pytest.raises(ArtifactIntegrityError, match="contained designs"):
            list(da.inspect(source, view="designs", all=True).records())
    with pytest.raises(ArtifactIntegrityError, match="contained designs"):
        list(da.inspect([run, bundle], view="designs", all=True).records())
    destination = tmp_path / "corrupt.json"
    with pytest.raises(ArtifactIntegrityError, match="contained designs"):
        da.export([run, bundle], all=True, format="json", out=destination)
    assert not destination.exists()


def test_quality_unknown_read_estimate_and_missing_history_text(tmp_path: Path):
    run = library(tmp_path / "run")
    report = da.inspect(run, view="quality")
    unestimated = replace(report, query=replace(report.query, source_records=None))
    assert unestimated.cost.records_estimate is None
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    result = CliRunner().invoke(app, ["inspect", str(bundle), "--view", "quality"])
    assert result.exit_code == 0, result.output
    assert "Search history: not included" in result.stdout
    assert "2 included; 2 selected" in result.stdout
    assert "outcomes: null" not in result.stdout


def test_empty_bundle_quality_and_explicit_read_caps(tmp_path: Path):
    run = library(tmp_path / "run")
    bundle = tmp_path / "empty"
    da.export(
        run,
        select=reporting.LibrarySelection(take=reporting.Take(count=0)),
        format="bundle",
        out=bundle,
    )
    value = da.inspect(bundle, view="quality").to_dict()
    assert value["selection"]["designs"] == 0
    assert value["source_runs"][0]["attainment"]["accepted"] == 2
    assert value["source_runs"][0]["included_designs"] == 0
    assert value["search"]["attempt_counts"] is None
    mixed = da.inspect([bundle, run], view="quality")
    result = mixed.to_dict()
    assert result["selection"]["designs"] == 2
    assert result["examined"] <= mixed.cost.records_estimate
    for limits in (reporting.ReadLimits(records=1), reporting.ReadLimits(identities=2)):
        with pytest.raises(reporting.ReadLimitError):
            da.inspect([bundle, run], view="quality", read_limits=limits).to_dict()


def test_conflicting_bundle_copy_is_rejected_before_filtering(tmp_path: Path):
    run = library(tmp_path / "run")
    bundle = tmp_path / "bundle"
    da.export(run, all=True, format="bundle", out=bundle)
    with sqlite3.connect(bundle / "bundle.sqlite3") as connection:
        value = json.loads(
            connection.execute(
                "SELECT payload FROM designs WHERE ordinal=1"
            ).fetchone()[0]
        )
        value["realized"]["provenance"]["note"] = "conflicting evidence"
        connection.execute(
            "UPDATE designs SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(value), semantic_digest(value)),
        )
    last = list(da.inspect(run, view="designs", all=True).records())[-1]
    with pytest.raises(ArtifactIntegrityError, match="conflicting content"):
        list(
            da.inspect(
                [run, bundle],
                view="sequences",
                all=True,
                select=reporting.DesignFilter(design_ids=(last.reference,)),
            ).records()
        )
