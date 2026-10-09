"""Selected composition stays distinct from source attainment and search effort.

Author: Eric J. South.
"""

import json
import shutil
import sqlite3
from pathlib import Path

import pytest
from PIL import Image
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts import RunHandle
from dense_arrays.cli import app
from dense_arrays.playback.quality import quality_figure

from .test_collections import library
from .test_quality import shortfall_run


def test_filtered_quality_matches_record_selection_and_preserves_attainment(
    tmp_path: Path,
):
    run = shortfall_run(tmp_path)
    selected = reporting.DesignFilter(groups=("A",))
    records = list(da.inspect(run, view="designs", select=selected, all=True).records())
    report = da.inspect(run, view="quality", select=selected, limit=2)
    value = report.to_dict()
    assert value["schema"] == "dense_arrays.quality.v3"
    with pytest.raises(ValueError, match="unsupported quality report schema"):
        quality_figure({**value, "schema": "dense_arrays.quality.v1"})
    assert value["selection"]["designs"] == len(records) == 4
    assert value["attainment"]["accepted"] == 8
    assert value["attainment"]["target"] == 12
    assert value["attainment"]["shortfall"] == 4
    assert value["composition"]["length"]["count"] == 4
    assert value["composition"]["length"]["denominator"] == "selected_designs"
    assert value["supply"]["unused_parts"] == 6
    assert value["search"]["attempt_counts"]["started"] == 10
    assert value["search"]["population"] == "all_attempts_in_source_snapshots"
    assert value["examined"] == report.cost.records_estimate
    continued = da.inspect(
        run, view="quality", select=selected, after=value["next_cursor"]
    ).to_dict()
    assert continued["selection"] == value["selection"]
    assert continued["composition"] == value["composition"]
    with pytest.raises(ValueError, match="cursor"):
        da.inspect(run, view="quality", after=value["next_cursor"])
    response = CliRunner().invoke(
        app,
        [
            "inspect",
            str(run.path),
            "--view",
            "quality",
            "--group",
            "A",
            "--limit",
            "2",
            "--json",
        ],
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout) == value


def test_combined_quality_deduplicates_records_and_preserves_source_context(
    tmp_path: Path,
):
    left, right = shortfall_run(tmp_path / "left"), shortfall_run(tmp_path / "right")
    copied = tmp_path / "copy"
    shutil.copytree(left.path, copied)
    report = da.inspect([left, right, copied], view="quality", limit=2)
    value = report.to_dict()
    assert value["attainment"] is None
    assert value["selection"]["designs"] == 16
    assert value["selection"]["distinct_sequences"] == 8
    assert len(value["source_runs"]) == 2
    assert [s["state"] for s in value["source_runs"]] == ["stopped", "stopped"]
    assert [s["attainment"]["shortfall"] for s in value["source_runs"]] == [4, 4]
    assert value["search"]["attempt_counts"]["started"] == 20
    assert value["search"]["attempt_counts"]["duplicate"] == 2
    assert value["supply"]["eligible_parts"] == 10
    assert value["part_usage"][0]["occurrences"] == 2
    assert value["part_usage"][0]["design_denominator"] == 16
    assert value["part_usage"][0]["collection_id"]
    assert value["requirements"][0]["cell_ref"].startswith(left.run_id)
    assert len(value["cells"]) == 2
    assert [cell["selected_designs"] for cell in value["cells"]] == [8, 8]
    assert value["examined"] == report.cost.records_estimate
    filtered = da.inspect(
        [left, right, copied],
        view="quality",
        select=reporting.DesignFilter(groups=("A",)),
    ).to_dict()
    assert filtered["selection"]["designs"] == 8
    assert [cell["selected_designs"] for cell in filtered["cells"]] == [4, 4]
    assert filtered["search"] == value["search"]
    continued = da.inspect(
        [left, right, copied], view="quality", after=value["next_cursor"]
    ).to_dict()
    assert continued["selection"] == value["selection"]
    assert continued["composition"] == value["composition"]
    assert len(continued["part_usage"]) == 8
    with pytest.raises(ValueError, match=r"cursor|identity"):
        da.inspect([right, left, copied], view="quality", after=value["next_cursor"])
    cli = CliRunner().invoke(
        app,
        [
            "inspect",
            str(left.path),
            str(right.path),
            str(copied),
            "--view",
            "quality",
            "--limit",
            "2",
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout) == value


def test_empty_filtered_quality_and_bounded_render_use_the_same_population(
    tmp_path: Path,
):
    run = shortfall_run(tmp_path)
    selected = reporting.DesignFilter(metrics={"length": reporting.Range(min=6)})
    value = da.inspect(run, view="quality", select=selected).to_dict()
    assert value["selection"]["designs"] == 0
    assert value["composition"]["length"]["mean"] is None
    assert value["part_usage"][0]["design_fraction"] is None
    assert value["attainment"]["accepted"] == 8
    with pytest.raises(reporting.ReadLimitError):
        da.inspect(
            [run, run], view="quality", read_limits=reporting.ReadLimits(records=1)
        ).to_dict()
    report = da.inspect(
        [run, run], view="quality", select=reporting.DesignFilter(groups=("A",))
    )
    output = tmp_path / "selected.png"
    receipt = da.render(report, view="library-quality", out=output)
    assert receipt.records == 4
    with Image.open(output) as image:
        embedded = json.loads(image.info["DenseArraysReport"])
    assert embedded == report.to_dict()
    figure = quality_figure(embedded)
    axes = {a.get_label(): a for a in figure.axes}
    assert sum(b.get_height() for b in axes["gc_fraction"].patches) == 4
    figure.clear()
    cli = CliRunner().invoke(
        app,
        [
            "render",
            str(run.path),
            str(run.path),
            "--view",
            "library-quality",
            "--group",
            "A",
            "--out",
            str(tmp_path / "cli.png"),
            "--json",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["records"] == 4


def test_quality_keeps_colliding_part_ids_in_their_collections(tmp_path: Path):
    left, right = library(tmp_path / "left"), library(tmp_path / "right", "GGG")
    value = da.inspect([left, right], view="quality").to_dict()
    assert value["supply"]["eligible_parts"] == 4
    assert value["concentration"]["highest_part_occurrence_share"] == 0.25
    assert len({row["part_ref"] for row in value["part_usage"]}) == 4
    assert [row["part_id"] for row in value["part_usage"]].count("a") == 2
    assert len(value["group_usage"]) == 2
    initial = RunHandle(left.path, left.run_id, revision=0)
    with pytest.raises(ValueError, match="one consistent revision"):
        da.inspect([initial, left], view="quality")
    with pytest.raises(ValueError, match="ambiguous"):
        da.inspect(
            [left, right],
            view="quality",
            select=reporting.DesignFilter(part_ids=("a",)),
        ).to_dict()
    with pytest.raises(ValueError, match="unknown"):
        da.inspect(
            [left, right],
            view="quality",
            select=reporting.DesignFilter(groups=("missing",)),
        ).to_dict()


def test_quality_rejects_conflicting_search_evidence_in_repeated_sources(
    tmp_path: Path,
):
    run = shortfall_run(tmp_path)
    copied = tmp_path / "copy"
    shutil.copytree(run.path, copied)
    with sqlite3.connect(copied / "run.sqlite3") as connection:
        ordinal, revision, payload = connection.execute(
            "SELECT attempt,revision,payload FROM attempts "
            "WHERE json_extract(payload,'$.outcome')='rejected'"
        ).fetchone()
        value = json.loads(payload)
        value["evidence"]["code"] = "different_reason"
        connection.execute(
            "UPDATE attempts SET payload=?,digest=? WHERE attempt=? AND revision=?",
            (canonical_json(value), semantic_digest(value), ordinal, revision),
        )
    with pytest.raises(ValueError, match="conflicting"):
        da.inspect([run, copied], view="quality").to_dict()
