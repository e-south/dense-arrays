"""Continuation cursors bind immutable snapshots and declared query semantics.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app


def make_pool(tmp_path: Path) -> parts.PoolHandle:
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,B\nc,GGG,A\n")
    return da.prepare(
        parts.PreparationSpec(parts.PartTable(table, "csv")), out=tmp_path / "pool"
    )


def test_filtered_pool_pages_resume_after_the_last_emitted_identity(tmp_path: Path):
    pool = make_pool(tmp_path)
    selected = parts.PartFilter(groups=("A",))
    view = da.inspect(pool, view="parts", select=selected, limit=1)
    with view.records() as rows:
        assert [r.part_id for r in rows] == ["a"]
        cursor = rows.next_cursor
    assert isinstance(cursor, str)
    second = da.inspect(pool, view="parts", select=selected, limit=2, after=cursor)
    assert second.cost.records_estimate == 2
    with second.records() as rows:
        assert [r.part_id for r in rows] == ["c"]
        assert rows.examined == 2
        assert rows.next_cursor is None
    with pytest.raises(ValueError, match=r"cursor.*query"):
        da.inspect(
            pool, view="parts", select=parts.PartFilter(groups=("B",)), after=cursor
        )
    with view.records() as again:
        assert [r.part_id for r in again] == ["a"]
        assert again.next_cursor == cursor


def test_cursor_keeps_the_original_revision_when_a_run_advances(tmp_path: Path):
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        )
    )
    with create_run(plan, tmp_path / "run") as writer:
        for _ in range(2):
            ordinal = writer.reserve(active_seconds=0)
            writer.publish(
                ordinal, "rejected", {"code": "screening_rejection"}, active_seconds=0
            )
        first = da.inspect(writer.handle, view="attempts", limit=1)
        with first.records() as rows:
            assert next(rows).attempt_id == 1
            cursor = rows.next_cursor
        ordinal = writer.reserve(active_seconds=0)
        writer.publish(
            ordinal, "rejected", {"code": "screening_rejection"}, active_seconds=0
        )
        continued = da.inspect(writer.handle, view="attempts", after=cursor, all=True)
        assert continued.revision == first.revision < da.inspect(writer.handle).revision
        with continued.records() as rows:
            assert [r.attempt_id for r in rows] == [2]
            assert rows.next_cursor is None
        with pytest.raises(ValueError, match=r"cursor.*query"):
            da.inspect(writer.handle, view="designs", after=cursor)


def test_cli_cursor_matches_python_and_rejects_malformed_tokens(tmp_path: Path):
    pool = make_pool(tmp_path)
    runner = CliRunner()
    first = runner.invoke(
        app, ["inspect", str(pool.path), "--view", "parts", "--limit", "1", "--json"]
    )
    assert first.exit_code == 0, first.output
    cursor = json.loads(first.stdout)["next_cursor"]
    second = runner.invoke(
        app,
        [
            "inspect",
            str(pool.path),
            "--view",
            "parts",
            "--after",
            cursor,
            "--all",
            "--json",
        ],
    )
    assert second.exit_code == 0, second.output
    assert [r["part"]["part_id"] for r in json.loads(second.stdout)["records"]] == [
        "b",
        "c",
    ]
    invalid = runner.invoke(
        app,
        [
            "inspect",
            str(pool.path),
            "--view",
            "parts",
            "--after",
            "not-a-cursor",
            "--json",
        ],
    )
    assert invalid.exit_code == 2
    assert "cursor" in json.loads(invalid.stdout)["message"]
