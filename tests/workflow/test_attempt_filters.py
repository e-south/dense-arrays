"""Typed attempt queries preserve pagination, source scope and complete totals.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays.artifacts.store import create_run
from dense_arrays.cli import app


def attempt_run(tmp_path: Path):
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        )
    )
    with create_run(plan, tmp_path / "run") as writer:
        for outcome in ("rejected", "duplicate", "rejected"):
            ordinal = writer.reserve(active_seconds=0)
            writer.publish(ordinal, outcome, {}, active_seconds=0)
        writer.finish("stopped", "attempt_limit", active_seconds=0)
        return writer.handle


def test_attempt_filter_uses_or_within_fields_and_and_across_fields(tmp_path: Path):
    run = attempt_run(tmp_path)
    selected = reporting.AttemptFilter(
        outcomes=("rejected",), cells=(f"{run.run_id}/default",)
    )
    first = da.inspect(run, view="attempts", select=selected, limit=1)
    with first.records() as rows:
        assert [row.attempt_id for row in rows] == [1]
        cursor = rows.next_cursor
    second = da.inspect(run, view="attempts", select=selected, after=cursor, all=True)
    with second.records() as rows:
        assert [row.attempt_id for row in rows] == [3]
        assert rows.examined == 2
    narrowed = reporting.AttemptFilter(attempt_ids=(1, 2), outcomes=("duplicate",))
    with da.inspect(run, view="attempts", select=narrowed).records() as rows:
        assert [row.attempt_id for row in rows] == [2]
    for bad in (
        reporting.AttemptFilter(attempt_ids=(4,)),
        reporting.AttemptFilter(cells=("other",)),
    ):
        with pytest.raises(ValueError, match="unknown"):
            da.inspect(run, view="attempts", select=bad)
    with pytest.raises(TypeError, match="AttemptFilter"):
        da.inspect(run, view="designs", select=selected)
    with pytest.raises(ValueError, match="outcome"):
        reporting.AttemptFilter(outcomes=("maybe",))


def test_filtered_diagnostics_label_the_population_and_match_cli(tmp_path: Path):
    run = attempt_run(tmp_path)
    selected = reporting.AttemptFilter(outcomes=("rejected",))
    report = da.inspect(run, view="diagnostics", select=selected)
    value = report.to_dict()
    assert value["population"] == "matching_attempts_at_revision"
    assert value["filter"] == selected.to_dict()
    assert value["attempt_counts"]["started"] == 2
    assert value["attempt_counts"]["duplicate"] == 0
    assert value["examined"] == 4
    args = [
        "inspect",
        str(run.path),
        "--view",
        "diagnostics",
        "--outcome",
        "rejected",
        "--json",
    ]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == value
    selection = tmp_path / "filter.json"
    selection.write_text(json.dumps(selected.to_dict()))
    result = CliRunner().invoke(
        app, [*args[:-3], "--selection", str(selection), "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == value
    conflict = CliRunner().invoke(app, [*args, "--selection", str(selection)])
    assert conflict.exit_code == 2
    assert "exclusive" in conflict.stderr


def test_filter_work_bound_and_unknown_schema_fail_before_iteration(tmp_path: Path):
    run = attempt_run(tmp_path)
    selected = reporting.AttemptFilter(attempt_ids=(1, 2))
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(
            run,
            view="attempts",
            select=selected,
            read_limits=reporting.ReadLimits(identities=1),
        )
    value = selected.to_dict()
    assert reporting.AttemptFilter.from_dict(value) == selected
    value["schema"] = "dense_arrays.attempt-filter.v99"
    with pytest.raises(ValueError, match="schema"):
        reporting.AttemptFilter.from_dict(value)
