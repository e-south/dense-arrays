"""Saved quality reports reconcile source attainment and available search evidence.

Module Author(s): Eric J. South
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts import RunHandle
from dense_arrays.cli import app


def _shortfall(tmp_path: Path) -> RunHandle:
    return da.run(
        planning.DesignSpec(
            [parts.Part("site", "ACGTTGCAAGTCCTGA")],
            planning.Length(maximum=16),
            strands="single",
            target=planning.Target(2),
        ),
        out=tmp_path / "run",
    )


@pytest.mark.parametrize(
    "mutation",
    [
        "completed_shortfall",
        "unknown_state",
        "accepted_history",
        "negative_time",
        "source_accepted",
        "source_time",
        "availability",
        "missing_history",
    ],
)
def test_snapshot_rejects_contradictions_before_rendering(
    tmp_path: Path, mutation: str
):
    run = _shortfall(tmp_path)
    value = da.inspect(run, view="quality").to_dict()
    source = value["source_runs"][0]
    assert source["state"] == "stopped"
    assert source["attainment"]["shortfall"] == 1
    assert value["search"]["availability"] == "complete"
    if mutation == "completed_shortfall":
        source["state"] = "completed"
    elif mutation == "unknown_state":
        source["state"] = "invented"
    elif mutation == "accepted_history":
        value["search"]["attempt_counts"]["accepted"] = 0
        value["search"]["attempt_counts"]["rejected"] += 1
    elif mutation == "negative_time":
        value["search"]["active_seconds"] = -1
    elif mutation == "source_accepted":
        source["search"]["attempt_counts"]["accepted"] = 0
        source["search"]["attempt_counts"]["rejected"] += 1
    elif mutation == "source_time":
        source["search"]["active_seconds"] += 1
    elif mutation == "availability":
        value["search"]["availability"] = "partial"
    else:
        source["search"] = None
    with pytest.raises(ValueError, match=r"source|search|completion|active_seconds"):
        reporting.QualitySnapshot.from_dict(value)
    saved = tmp_path / "quality.json"
    saved.write_text(json.dumps(value))
    output = tmp_path / "invalid.png"
    result = CliRunner().invoke(
        app, ["render", str(saved), "--view", "library-quality", "--out", str(output)]
    )
    assert result.exit_code == 2
    assert not output.exists()


def test_snapshot_keeps_partial_and_absent_histories_separate(tmp_path: Path):
    first = _shortfall(tmp_path)
    second_root = tmp_path / "second"
    second_root.mkdir()
    second = _shortfall(second_root)
    bundle = tmp_path / "bundle"
    da.export(first, all=True, format="bundle", out=bundle)
    for sources, availability in (
        (bundle, "not_included"),
        ([bundle, second], "partial"),
    ):
        value = da.inspect(sources, view="quality").to_dict()
        assert value["search"]["availability"] == availability
        assert reporting.QualitySnapshot.from_dict(value).to_dict() == value
