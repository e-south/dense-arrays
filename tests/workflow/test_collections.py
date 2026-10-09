"""Combined libraries preserve full design identity, snapshots and scalar joins.

Author: Eric J. South.
"""

import csv
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
from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app


def library(path: Path, second: str = "CCC"):
    return da.run(
        planning.DesignSpec(
            parts=[
                parts.Part("a", "AAA", group="A"),
                parts.Part("b", second, group="B"),
            ],
            length=planning.Length(maximum=6),
            strands="single",
            target=planning.Target(count=2),
        ),
        out=path,
    )


def test_combined_views_deduplicate_full_identity_and_page_inside_placements(
    tmp_path: Path,
):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    sources = [left, right, left]
    designs = list(da.inspect(sources, view="designs", all=True).records())
    assert len(designs) == 4
    assert [d.run_id for d in designs] == [
        left.run_id,
        left.run_id,
        right.run_id,
        right.run_id,
    ]
    assert len({d.sequence_id for d in designs}) == 2
    sequences = list(da.inspect(sources, view="sequences", all=True).records())
    placements = list(da.inspect(sources, view="placements", all=True).records())
    assert len(sequences) == 4
    assert len(placements) == 8
    with da.inspect(sources, view="placements", limit=3).records() as rows:
        first = list(rows)
        after = rows.next_cursor
    rest = list(da.inspect(sources, view="placements", after=after, all=True).records())
    assert [p.to_dict() for p in first + rest] == [p.to_dict() for p in placements]
    by_ref = {s.design_ref: s.sequence for s in sequences}
    for p in placements:
        assert by_ref[p.design_ref][p.start : p.end] == p.sequence


def test_combined_export_matches_cli_and_keeps_colliding_local_ids(tmp_path: Path):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    stream = io.StringIO()
    receipt = da.export(
        [left, right, left], view="sequences", all=True, format="csv", out=stream
    )
    assert receipt.records == 4
    assert len(receipt.design_refs) == 4
    rows = list(csv.DictReader(io.StringIO(stream.getvalue())))
    assert len({r["design_ref"] for r in rows}) == 4
    response = CliRunner().invoke(
        app,
        [
            "export",
            str(left.path),
            str(right.path),
            str(left.path),
            "--view",
            "sequences",
            "--all",
            "--format",
            "csv",
            "--out",
            "-",
        ],
    )
    assert response.exit_code == 0, response.output
    assert response.stdout == stream.getvalue()
    inspected = CliRunner().invoke(
        app,
        [
            "inspect",
            str(left.path),
            str(right.path),
            "--view",
            "sequences",
            "--all",
            "--json",
        ],
    )
    assert inspected.exit_code == 0, inspected.output
    assert len(json.loads(inspected.stdout)["records"]) == 4


def test_combined_filters_reject_ambiguous_names_and_accept_full_references(
    tmp_path: Path,
):
    left, right = library(tmp_path / "left"), library(tmp_path / "right", "GGG")
    first = next(da.inspect(left, view="designs").records())
    for selected in (
        reporting.DesignFilter(design_ids=(first.design_id,)),
        reporting.DesignFilter(cells=("default",)),
        reporting.DesignFilter(part_ids=("a",)),
    ):
        with pytest.raises(ValueError, match="ambiguous"):
            list(da.inspect([left, right], view="designs", select=selected).records())
    selected = reporting.DesignFilter(
        design_ids=(first.reference,), cells=(f"{left.run_id}/default",)
    )
    assert [
        r.reference
        for r in da.inspect([left, right], view="designs", select=selected).records()
    ] == [first.reference]
    placement = next(da.inspect(left, view="placements").records())
    scoped = reporting.DesignFilter(part_ids=(f"{placement.collection_id}/a",))
    assert {
        r.run_id
        for r in da.inspect(
            [left, right], view="designs", select=scoped, all=True
        ).records()
    } == {left.run_id}
    with pytest.raises(ValueError, match="unknown"):
        list(
            da.inspect(
                [left, right],
                view="designs",
                select=reporting.DesignFilter(design_ids=("missing",)),
            ).records()
        )


def test_conflicting_full_records_fail_even_when_sequence_projection_is_equal(
    tmp_path: Path,
):
    left = library(tmp_path / "left")
    copy = tmp_path / "copy"
    shutil.copytree(left.path, copy)
    with sqlite3.connect(copy / "run.sqlite3") as connection:
        payload = connection.execute(
            "SELECT payload FROM designs WHERE ordinal=1"
        ).fetchone()[0]
        value = json.loads(payload)
        value["realized"]["provenance"]["note"] = "conflicting evidence"
        connection.execute(
            "UPDATE designs SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(value), semantic_digest(value)),
        )
    output = tmp_path / "library.fasta"
    with pytest.raises(ValueError, match="conflicting"):
        da.export([left, copy], view="sequences", all=True, format="fasta", out=output)
    assert not output.exists()


def test_combined_cursor_keeps_each_snapshot_and_rejects_scope_changes(tmp_path: Path):
    reference = library(tmp_path / "reference")
    plan = da.inspect(reference, view="plan")
    candidates = list(da.inspect(reference, view="designs", all=True).records())
    evidence = {"solver_status": "optimal", "proof_scope": "offered_packing_model"}
    with create_run(plan, tmp_path / "active") as writer:
        ordinal = writer.reserve(active_seconds=0)
        writer.publish(
            ordinal,
            "accepted",
            evidence,
            active_seconds=0,
            design=replace(candidates[0], run_id=writer.handle.run_id),
        )
        sources = [writer.handle, reference]
        first = da.inspect(sources, view="placements", limit=1)
        with first.records() as records:
            prefix = list(records)
            cursor = records.next_cursor
        ordinal = writer.reserve(active_seconds=0)
        writer.publish(
            ordinal,
            "accepted",
            evidence,
            active_seconds=0,
            design=replace(candidates[1], run_id=writer.handle.run_id),
        )
        rest = da.inspect(sources, view="placements", after=cursor, all=True)
        assert rest.sources[0]["revision"] == first.sources[0]["revision"]
        assert rest.sources[0]["revision"] < da.inspect(writer.handle).revision
        assert len(prefix + list(rest.records())) == 6
        assert (
            len(list(da.inspect(sources, view="placements", all=True).records())) == 8
        )
        for changed in ([reference, writer.handle], [writer.handle]):
            with pytest.raises(ValueError, match=r"cursor|identity"):
                da.inspect(changed, view="placements", after=cursor)
        with pytest.raises(ValueError, match="cursor"):
            da.inspect(sources, view="sequences", after=cursor)
        with pytest.raises(ValueError, match="cursor"):
            da.inspect(
                sources,
                view="placements",
                after=cursor,
                select=reporting.DesignFilter(groups=("A",)),
            )
        bound = RunHandle(
            writer.handle.path,
            writer.handle.run_id,
            revision=da.inspect(writer.handle).revision,
        )
        with pytest.raises(ValueError, match="cursor"):
            da.inspect([bound, reference], view="placements", after=cursor)


def test_combined_caps_and_missing_sources_do_not_publish(tmp_path: Path):
    left, right = library(tmp_path / "left"), library(tmp_path / "right")
    out = tmp_path / "handoff.csv"
    with pytest.raises(sqlite3.OperationalError):
        da.export(
            [left, tmp_path / "missing"],
            view="sequences",
            all=True,
            format="csv",
            out=out,
        )
    assert not out.exists()
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            [left, right],
            view="sequences",
            all=True,
            format="csv",
            out=out,
            read_limits=reporting.ReadLimits(identities=5),
        )
    assert not out.exists()
    with pytest.raises(reporting.ReadLimitError):
        da.export(
            [left, left],
            view="sequences",
            all=True,
            format="csv",
            out=out,
            read_limits=reporting.ReadLimits(records=1),
        )
    assert not out.exists()
    query = da.inspect([left, left], view="designs", all=True)
    with query.records() as records:
        assert len(list(records)) == 2
        assert records.examined == query.cost.records_estimate == 4


def test_part_identity_survives_changed_target_and_supports_qualified_filters(
    tmp_path: Path,
):
    run = library(tmp_path / "library")
    plan = da.inspect(run, view="plan")
    changed = planning.GenerationPlan(
        replace(plan.request, target=planning.Target(count=1), seed=19),
        plan.inputs,
        plan.import_report,
    )
    assert changed.plan_id != plan.plan_id
    assert changed.collection_id == plan.collection_id
    selected = reporting.DesignFilter(part_ids=(f"{plan.collection_id}/a",))
    assert (
        len(list(da.inspect(run, view="designs", all=True, select=selected).records()))
        == 2
    )
    assert (
        len(
            list(
                da.inspect(
                    [run, run], view="designs", all=True, select=selected
                ).records()
            )
        )
        == 2
    )
    assert (
        len(
            list(
                da.inspect(
                    [run, run],
                    view="designs",
                    all=True,
                    select=reporting.DesignFilter(part_ids=("a",), groups=("A",)),
                ).records()
            )
        )
        == 2
    )
