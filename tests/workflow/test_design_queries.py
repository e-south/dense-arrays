"""Accepted-design queries share identity and geometry across projections.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.store import reader
from dense_arrays.cli import app


def query_run(tmp_path: Path):
    return da.run(
        planning.DesignSpec(
            parts=[
                parts.Part("a", "AAA", group="A"),
                parts.Part("b", "CCC", group="B"),
            ],
            length=planning.Length(maximum=6),
            strands="single",
            target=planning.Target(count=2),
            requirements=[
                planning.Occurrences(
                    "both", parts.PartSelector(part_ids=("a", "b")), min=2, max=2
                )
            ],
        ),
        out=tmp_path / "run",
    )


def test_design_filter_and_projections_preserve_joinable_identities(tmp_path: Path):
    run = query_run(tmp_path)
    designs = list(da.inspect(run, view="designs", all=True).records())
    assert len(designs) == 2
    selected = reporting.DesignFilter(
        design_ids=(designs[0].reference, designs[1].design_id),
        cells=("default",),
        part_ids=("a",),
        groups=("B",),
        metrics={"gc_fraction": reporting.Range(min=0.5, max=0.5)},
    )
    assert reporting.DesignFilter.from_dict(selected.to_dict()) == selected
    sequences = list(
        da.inspect(run, view="sequences", select=selected, all=True).records()
    )
    assert [s.design_ref for s in sequences] == [d.reference for d in designs]
    assert {s.sequence for s in sequences} == {"AAACCC", "CCCAAA"}
    assert all(s.length == 6 and s.gc_fraction == 0.5 for s in sequences)
    placements = list(
        da.inspect(run, view="placements", select=selected, all=True).records()
    )
    assert len(placements) == 4
    by_ref = {s.design_ref: s for s in sequences}
    for p in placements:
        assert by_ref[p.design_ref].sequence[p.start : p.end] == p.sequence
        assert p.group == {"a": "A", "b": "B"}[p.part_id]
        assert p.collection_id
        assert p.orientation == "forward"
    # Placement pages can end inside a design without dropping or repeating rows.
    first = da.inspect(run, view="placements", select=selected, limit=1)
    with first.records() as rows:
        head = list(rows)
        cursor = rows.next_cursor
    tail = list(
        da.inspect(
            run, view="placements", select=selected, after=cursor, all=True
        ).records()
    )
    assert [p.to_dict() for p in head + tail] == [p.to_dict() for p in placements]


def test_design_queries_reject_unknown_selectors_and_match_cli(tmp_path: Path):
    run = query_run(tmp_path)
    selected = reporting.DesignFilter(
        groups=("A",), metrics={"packing_density": reporting.Range(min=1)}
    )
    expected = [
        r.to_dict()
        for r in da.inspect(run, view="sequences", select=selected, all=True).records()
    ]
    selection = tmp_path / "selected.json"
    selection.write_text(json.dumps(selected.to_dict()))
    response = CliRunner().invoke(
        app,
        [
            "inspect",
            str(run.path),
            "--view",
            "sequences",
            "--selection",
            str(selection),
            "--all",
            "--json",
        ],
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout)["records"] == expected
    flags = CliRunner().invoke(
        app,
        [
            "inspect",
            str(run.path),
            "--view",
            "sequences",
            "--group",
            "A",
            "--part-id",
            "b",
            "--cell",
            "default",
            "--all",
            "--json",
        ],
    )
    assert flags.exit_code == 0, flags.output
    assert json.loads(flags.stdout)["records"] == expected
    for kwargs in (
        {"design_ids": ("missing",)},
        {"part_ids": ("missing",)},
        {"groups": ("missing",)},
        {"cells": ("missing",)},
    ):
        with pytest.raises(ValueError, match="unknown"):
            list(
                da.inspect(
                    run, view="designs", select=reporting.DesignFilter(**kwargs)
                ).records()
            )
    with pytest.raises(ValueError, match="metric"):
        reporting.DesignFilter(metrics={"similarity": reporting.Range(min=0)})
    empty = reporting.DesignFilter(metrics={"gc_fraction": reporting.Range(min=0.9)})
    assert list(da.inspect(run, view="sequences", select=empty).records()) == []


def test_reverse_core_coordinates_and_identity_filtering(tmp_path: Path):
    run = da.run(
        planning.DesignSpec(
            parts=[
                parts.Part(
                    "used",
                    "ACGTTA",
                    group="A",
                    core_start=1,
                    core_end=3,
                    core_orientation="forward",
                ),
                parts.Part("unused", "ACGTTA", group="B"),
            ],
            length=planning.Length(maximum=6),
            requirements=[
                planning.Fixed("reverse", "used", "reverse"),
                planning.Occurrences(
                    "exclude", parts.PartSelector(part_ids=("unused",)), max=0
                ),
            ],
        ),
        out=tmp_path / "run",
    )
    rows = list(da.inspect(run, view="placements", all=True).records())
    assert len(rows) == 1
    row = rows[0]
    assert (row.sequence, row.orientation) == ("TAACGT", "reverse")
    assert (row.core_start, row.core_end, row.core_orientation) == (3, 5, "reverse")
    assert (
        list(
            da.inspect(
                run, view="designs", select=reporting.DesignFilter(groups=("B",))
            ).records()
        )
        == []
    )


def test_design_query_bounds_cover_plan_and_rows_and_cli_rejects_mixed_flags(
    tmp_path: Path,
):
    run = query_run(tmp_path)
    with pytest.raises(reporting.ReadLimitError, match="records"):
        list(
            da.inspect(
                run,
                view="placements",
                all=True,
                read_limits=reporting.ReadLimits(records=2),
            ).records()
        )
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        list(
            da.inspect(
                run, view="placements", read_limits=reporting.ReadLimits(identities=1)
            ).records()
        )
    response = CliRunner().invoke(
        app,
        [
            "inspect",
            str(run.path),
            "--view",
            "sequences",
            "--outcome",
            "rejected",
            "--json",
        ],
    )
    assert response.exit_code == 2
    assert "do not apply" in response.stderr


def test_unknown_design_is_a_query_error_and_bound_render_does_not_advance(
    tmp_path: Path,
):
    run = query_run(tmp_path)
    response = CliRunner().invoke(
        app,
        [
            "inspect",
            str(run.path),
            "--view",
            "designs",
            "--design-id",
            "missing",
            "--json",
        ],
    )
    assert response.exit_code == 2
    assert "invalid_input" in response.stderr

    bound = RunHandle(run.path, run.run_id, revision=0)
    assert list(da.inspect(bound, view="designs", all=True).records()) == []
    with pytest.raises(ValueError, match="exactly one"):
        da.render(bound, out=tmp_path / "initial.png")
    with reader(run.path) as connection:
        first_revision = connection.execute(
            "SELECT MIN(revision) FROM designs"
        ).fetchone()[0]
    first_accepted = RunHandle(run.path, run.run_id, revision=first_revision)
    assert da.inspect(first_accepted).accepted == 1
    receipt = da.render(first_accepted, out=tmp_path / "first.png")
    assert receipt.records == 1
    assert receipt.sources[0]["revision"] == first_revision
