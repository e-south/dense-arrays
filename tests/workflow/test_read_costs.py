"""Artifact reading declares and enforces work independently of displayed rows.

Author: Eric J. South.
"""

import json
import shutil
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning, reporting
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.cli import app
from dense_arrays.reporting.bundles.reading import resolve_filter
from dense_arrays.workflow import operations


def test_bundle_alias_resolution_does_not_fetch_unrequested_identities():
    with sqlite3.connect(":memory:") as connection:
        connection.execute(
            "CREATE TABLE designs (local_id TEXT, design_ref TEXT UNIQUE)"
        )
        connection.execute("CREATE INDEX local_ids ON designs(local_id)")
        connection.executemany(
            "INSERT INTO designs VALUES (?,?)",
            ((f"design_{i}", f"run/cell/design_{i}") for i in range(10_000)),
        )
        fetched = []

        def observe(_cursor: object, row: tuple) -> tuple:
            fetched.append(row)
            return row

        connection.row_factory = observe
        requested = reporting.DesignFilter(design_ids=("design_0",))
        query = SimpleNamespace(
            select=requested, summary=SimpleNamespace(manifest={"source_runs": []})
        )
        budget = ReadBudget(reporting.ReadLimits(records=1, identities=3))
        resolved = resolve_filter(connection, query, {}, budget)
    assert resolved.design_ids == ("run/cell/design_0",)
    assert fetched == [("design_0", "run/cell/design_0")]
    assert budget.identities == 1


def pool(tmp_path: Path) -> parts.PoolHandle:
    """Retain a late matching row to distinguish scans from page size."""
    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence,group\na,AAA,A\nb,CCC,A\nc,GGG,B\n")
    return da.prepare(
        parts.PreparationSpec(source=parts.PartTable(table, "csv")),
        out=tmp_path / "pool",
    )


def test_filtered_page_declares_scan_and_enforces_record_work_bound(tmp_path: Path):
    source = pool(tmp_path)
    view = da.inspect(
        source,
        view="parts",
        select=parts.PartFilter(groups=("B",)),
        limit=1,
        read_limits=reporting.ReadLimits(records=2),
    )
    assert view.cost.mode == "scan"
    assert view.cost.records_estimate == 3
    assert view.cost.source_id == source.pool_id
    assert view.cost.revision == 0
    assert view.cost.limits.records == 2
    with view.records() as rows:
        assert rows.examined == 0
        with pytest.raises(reporting.ReadLimitError, match="records"):
            next(rows)
        assert rows.examined == 2
        assert rows.returned == 0


def test_exact_cap_completes_and_iteration_state_is_independent(tmp_path: Path):
    source = pool(tmp_path)
    view = da.inspect(
        source, view="parts", all=True, read_limits=reporting.ReadLimits(records=3)
    )
    assert view.cost.mode == "indexed"
    assert view.cost.projection == "parts"
    first, second = view.records(), view.records()
    with first, second:
        assert next(first).part_id == "a"
        assert first.examined == 1
        assert second.examined == 0
        first.close()
        assert list(first) == []
        assert [r.part_id for r in second] == ["a", "b", "c"]
        assert second.examined == second.returned == 3


def test_cli_announces_read_cost_and_fails_when_scan_cap_is_reached(tmp_path: Path):
    source = pool(tmp_path)
    result = CliRunner().invoke(
        app,
        [
            "inspect",
            str(source.path),
            "--view",
            "parts",
            "--group",
            "B",
            "--limit",
            "1",
            "--max-read-records",
            "2",
            "--json",
        ],
    )
    assert result.exit_code == 4, result.output
    assert "Read cost:" in result.stderr
    assert "read_limit" in result.stderr
    assert result.stderr.index("Read cost:") < result.stderr.index("read_limit")
    assert "read_limit" not in result.stdout
    completed = CliRunner().invoke(
        app,
        [
            "inspect",
            str(source.path),
            "--view",
            "parts",
            "--all",
            "--max-read-records",
            "3",
            "--json",
        ],
    )
    assert completed.exit_code == 0, completed.output
    assert len(json.loads(completed.stdout)["records"]) == 3


@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_read_limits_reject_nonpositive_or_noninteger_caps(value: object):
    with pytest.raises((ValueError, TypeError)):
        reporting.ReadLimits(records=value)


@pytest.mark.parametrize(
    "options",
    [
        {"limit": 1},
        {"all": True},
        *({"view": "designs", "limit": value} for value in (0, -1, True, 1.5)),
        {"verify": 0},
        {"all": []},
    ],
)
def test_inspection_rejects_invalid_common_options_before_source_discovery(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, options: dict
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("invalid query reached artifact discovery")

    monkeypatch.setattr(operations, "is_bundle", forbidden)
    with pytest.raises((ValueError, TypeError)):
        da.inspect(tmp_path / "missing", **options)


@pytest.mark.parametrize("options", [{"verify": 0}, {"all": []}])
def test_in_memory_plan_checks_common_option_types(options: dict):
    plan = da.plan(
        planning.DesignSpec((parts.Part("a", "AAA"),), planning.Length(maximum=3))
    )
    with pytest.raises(TypeError, match="booleans"):
        da.inspect(plan, **options)


def test_record_view_rejects_another_run_at_the_same_path(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
    )
    first = da.run(request, out=tmp_path / "first")
    second = da.run(request, out=tmp_path / "second")
    view = da.inspect(first, view="designs")
    assert view.cost.source_id == first.run_id
    shutil.move(first.path, tmp_path / "moved")
    shutil.move(second.path, first.path)
    with (
        view.records() as rows,
        pytest.raises(ArtifactIntegrityError, match="identity"),
    ):
        next(rows)


def test_early_close_releases_the_actual_reader_connection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    source = pool(tmp_path)
    view = da.inspect(source, view="parts", all=True)
    connect = sqlite3.connect
    connections = []

    def observed(*args: object, **kwargs: object) -> sqlite3.Connection:
        connection = connect(*args, **kwargs)
        connections.append(connection)
        return connection

    monkeypatch.setattr(sqlite3, "connect", observed)
    rows = view.records()
    assert connections == []
    assert next(rows).part_id == "a"
    assert len(connections) == 1
    assert connections[0].execute("SELECT 1").fetchone() == (1,)
    rows.close()
    with pytest.raises(sqlite3.ProgrammingError, match="closed"):
        connections[0].execute("SELECT 1")


def test_filter_identity_state_is_bounded_before_scanning(tmp_path: Path):
    source = pool(tmp_path)
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        da.inspect(
            source,
            view="parts",
            select=parts.PartFilter(part_ids=("a", "b")),
            read_limits=reporting.ReadLimits(identities=1),
        )


def test_record_descriptor_repr_does_not_expand_filter_identities(tmp_path: Path):
    source = pool(tmp_path)
    view = da.inspect(
        source, view="parts", select=parts.PartFilter(part_ids=("a", "b"))
    )
    assert len(repr(view)) < 200
    assert "PartFilter(" not in repr(view)


def test_run_reader_rechecks_snapshot_schema_before_emitting_records(tmp_path: Path):
    run = da.run(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA")], length=planning.Length(maximum=3)
        ),
        out=tmp_path / "run",
    )
    view = da.inspect(run, view="designs")
    connection = sqlite3.connect(run.path / "run.sqlite3")
    try:
        wire = json.loads(
            connection.execute(
                "SELECT payload FROM commits WHERE revision=?", (view.revision,)
            ).fetchone()[0]
        )
        wire["schema"] = "dense_arrays.run.v99"
        with connection:
            connection.execute(
                "UPDATE commits SET payload=?,digest=? WHERE revision=?",
                (canonical_json(wire), semantic_digest(wire), view.revision),
            )
    finally:
        connection.close()
    with view.records() as rows, pytest.raises(ArtifactIntegrityError, match="schema"):
        next(rows)
