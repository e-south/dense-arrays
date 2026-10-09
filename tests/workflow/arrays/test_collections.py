"""Supplied arrays carry literal geometry and identities without execution claims."""

import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays.arrays import ArrayCollection, ArrayFilter
from dense_arrays.cli import app
from dense_arrays.realized import Orientation
from dense_arrays.reporting import ReadLimitError, ReadLimits


def test_collection_round_trip_preserves_equal_sequences_and_repeated_parts(
    tmp_path: Path, source: ArrayCollection
):
    published = tmp_path / "collection"
    receipt = da.export(source, view="arrays", all=True, format="bundle", out=published)
    assert receipt.records == 2
    relocated = tmp_path / "moved"
    published.rename(relocated)
    summary = da.inspect(relocated, verify=True)
    assert summary.verified
    assert summary.arrays == 2
    assert summary.placements == 4
    assert summary.parts == 3
    value = summary.to_dict()
    assert value["evidence_boundary"] == "supplied_sequences_parts_and_placements"
    assert (
        not {"run_id", "plan_id", "attempt_id", "solver", "attainment"} & value.keys()
    )
    assert value["exporter"]["package"] == "dense-arrays"
    with da.inspect(relocated, view="arrays", all=True).records() as rows:
        actual = list(rows)
    assert [row.realized for row in actual] == list(source.arrays)
    assert actual[0].sequence_id == actual[1].sequence_id
    with da.inspect(relocated, view="placements", all=True).records() as rows:
        placements = [row.to_dict() for row in rows]
    assert [p["part_id"] for p in placements] == [
        "part-a",
        "part-b",
        "part-a",
        "part-r",
    ]
    assert (
        placements[-1]["core_start"],
        placements[-1]["core_end"],
        placements[-1]["core_orientation"],
    ) == (5, 16, "forward")
    assert placements[1]["core_start"] is None
    assert not list(tmp_path.glob("**/run.sqlite3"))


def test_publication_refuses_unknown_part_and_leaves_no_partial_bundle(
    tmp_path: Path, source: ArrayCollection
):

    incomplete = ArrayCollection(parts=source.parts[:1], arrays=source.arrays)
    with pytest.raises(ValueError, match="unknown part"):
        da.export(
            incomplete, view="arrays", all=True, format="bundle", out=tmp_path / "bad"
        )
    assert not (tmp_path / "bad").exists()


def test_existing_destination_is_never_replaced(
    tmp_path: Path, source: ArrayCollection
):
    destination = tmp_path / "owned"
    destination.mkdir()
    sentinel = destination / "notes.txt"
    sentinel.write_text("keep")
    with pytest.raises(FileExistsError):
        da.export(source, view="arrays", all=True, format="bundle", out=destination)
    assert sentinel.read_text() == "keep"


def test_cursor_continues_inside_an_array_and_rejects_changed_filter(
    tmp_path: Path, source: ArrayCollection
):

    da.export(source, view="arrays", all=True, format="bundle", out=tmp_path / "data")
    with da.inspect(tmp_path / "data", view="placements", limit=2).records() as rows:
        first = list(rows)
        cursor = rows.next_cursor
    with da.inspect(
        tmp_path / "data", view="placements", all=True, after=cursor
    ).records() as rows:
        rest = list(rows)
    assert [r.placement_id for r in first + rest] == ["one", "two", "repeat", "reverse"]
    changed = da.inspect(
        tmp_path / "data",
        view="placements",
        after=cursor,
        select=ArrayFilter(array_ids=("array-2",)),
    )
    with pytest.raises(ValueError, match="cursor"), changed.records() as rows:
        list(rows)


@pytest.mark.parametrize(
    "defect", ["duplicate-array", "part-orientation", "missing-orientation"]
)
def test_invalid_supplied_evidence_never_commits(
    tmp_path: Path, source: ArrayCollection, defect: str
):

    arrays = source.arrays
    if defect == "duplicate-array":
        arrays = (arrays[0], arrays[0])
    else:
        placement = arrays[1].placements[0]
        placement = replace(
            placement,
            orientation=Orientation.FORWARD
            if defect == "part-orientation"
            else Orientation.UNSPECIFIED,
        )
        arrays = (replace(arrays[1], placements=(placement,)),)
    with pytest.raises(ValueError, match=r"duplicate|orientation|oriented part"):
        da.export(
            ArrayCollection(source.parts, arrays),
            view="arrays",
            all=True,
            format="bundle",
            out=tmp_path / "invalid",
        )
    assert not (tmp_path / "invalid").exists()


def test_native_execution_views_remain_unavailable_and_limits_are_enforced(
    tmp_path: Path, source: ArrayCollection
):

    da.export(source, view="arrays", all=True, format="bundle", out=tmp_path / "data")
    for view in ("attempts", "plan", "quality", "designs"):
        with pytest.raises(ValueError, match="unavailable"):
            da.inspect(tmp_path / "data", view=view)
    with (
        pytest.raises(ReadLimitError),
        da.inspect(
            tmp_path / "data",
            view="arrays",
            all=True,
            read_limits=ReadLimits(records=4),
        ).records() as rows,
    ):
        list(rows)


def test_cli_and_python_publish_and_export_same_supplied_records(
    tmp_path: Path, source: ArrayCollection
):

    runner = CliRunner()
    transport = tmp_path / "arrays.jsonl"
    da.export(source, view="arrays", all=True, format="jsonl", out=transport)
    result = runner.invoke(
        app,
        [
            "export",
            str(transport),
            "--view",
            "arrays",
            "--all",
            "--format",
            "bundle",
            "--out",
            str(tmp_path / "data"),
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    summary = runner.invoke(
        app, ["inspect", str(tmp_path / "data"), "--verify", "--json"]
    )
    assert summary.exit_code == 0, summary.output
    assert json.loads(summary.stdout)["arrays"] == 2
    human = runner.invoke(app, ["inspect", str(tmp_path / "data")])
    assert human.exit_code == 0
    assert "2 arrays" in human.output
    page = runner.invoke(
        app,
        [
            "inspect",
            str(tmp_path / "data"),
            "--view",
            "placements",
            "--array-id",
            "array-2",
            "--json",
        ],
    )
    assert page.exit_code == 0, page.output
    with da.inspect(
        tmp_path / "data", view="placements", select=ArrayFilter(array_ids=("array-2",))
    ).records() as rows:
        expected = [row.to_dict() for row in rows]
    assert json.loads(page.stdout)["records"] == expected
    output = tmp_path / "placements.csv"
    result = runner.invoke(
        app,
        [
            "export",
            str(tmp_path / "data"),
            "--view",
            "placements",
            "--array-id",
            "array-2",
            "--all",
            "--format",
            "csv",
            "--out",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    with output.open() as stream:
        row = next(iter(csv.DictReader(stream)))
    assert row["core_start"] == "5"
    assert row["part_id"] == "part-r"
    fasta = tmp_path / "sequences.fa"
    da.export(tmp_path / "data", view="sequences", all=True, format="fasta", out=fasta)
    assert fasta.read_text().count(">") == 2
