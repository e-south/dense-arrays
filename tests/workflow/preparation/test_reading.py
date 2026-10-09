"""Saved preparation decisions remain queryable without rerunning preparation.

Author: Eric J. South.
"""

import json
import shutil
import sqlite3
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.cli import app
from dense_arrays.reporting import (
    CandidateFilter,
    PoolQualitySnapshot,
    ReadLimitError,
    ReadLimits,
)
from dense_arrays.reporting.pools import reading

from .test_sampled import background_request


def test_candidate_pages_preserve_pool_identity_and_saved_decisions(tmp_path: Path):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    view = da.inspect(pool, view="candidates", limit=2)
    assert view.cost.mode == "indexed"
    assert view.cost.records_estimate == 2
    with view.records() as records:
        first = list(records)
        cursor = records.next_cursor
        assert records.examined == 2
    assert [r.candidate.index for r in first] == [1, 2]
    assert [r.outcome for r in first] == ["retained", "duplicate_discarded"]
    assert all(r.pool_id == pool.pool_id for r in first)
    assert first[1].candidate.representative == 1
    assert first[1].candidate.rank is None
    with da.inspect(pool, view="candidates", after=cursor, all=True).records() as rows:
        assert [r.candidate.index for r in rows] == [3, 4, 5, 6]
    with pytest.raises(ValueError, match="cursor"):
        da.inspect(pool, view="parts", after=cursor)
    limited = da.inspect(
        pool, view="candidates", all=True, read_limits=ReadLimits(records=1)
    )
    with limited.records() as rows:
        assert next(rows).candidate.index == 1
        with pytest.raises(ReadLimitError, match="records"):
            next(rows)
    output = tmp_path / "candidates.json"
    da.export(pool, view="candidates", all=True, out=output)
    exported = json.loads(output.read_text())
    assert exported["sources"][0]["pool_id"] == pool.pool_id
    assert len(exported["records"]) == 6
    assert exported["records"][1] == first[1].to_dict()
    result = CliRunner().invoke(
        app,
        ["inspect", str(pool.path), "--view", "candidates", "--limit", "2", "--json"],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["records"] == [r.to_dict() for r in first]


def test_candidate_filters_are_shared_by_python_cli_and_export(tmp_path: Path):

    pool = da.prepare(background_request(), out=tmp_path / "pool")
    predicate = CandidateFilter(indices=(2, 3, 4), outcomes=("duplicate_discarded",))
    view = da.inspect(pool, view="candidates", select=predicate, limit=1)
    assert view.cost.mode == "scan"
    with view.records() as records:
        expected = [r.to_dict() for r in records]
        cursor = records.next_cursor
    assert expected[0]["candidate"]["index"] == 2
    result = CliRunner().invoke(
        app,
        [
            "inspect",
            str(pool.path),
            "--view",
            "candidates",
            "--candidate-index",
            "2",
            "--candidate-index",
            "3",
            "--candidate-index",
            "4",
            "--outcome",
            "duplicate_discarded",
            "--limit",
            "1",
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["records"] == expected
    with da.inspect(
        pool, view="candidates", select=predicate, after=cursor, all=True
    ).records() as rows:
        assert [r.candidate.index for r in rows] == [3, 4]
    with pytest.raises(ValueError, match="cursor"):
        da.inspect(pool, view="candidates", after=cursor)
    saved = tmp_path / "filter.json"
    saved.write_text(json.dumps(predicate.to_dict()))
    target = tmp_path / "selected.json"
    result = CliRunner().invoke(
        app,
        [
            "export",
            str(pool.path),
            "--view",
            "candidates",
            "--selection",
            str(saved),
            "--all",
            "--out",
            str(target),
        ],
    )
    assert result.exit_code == 0, result.output
    assert [
        r["candidate"]["index"] for r in json.loads(target.read_text())["records"]
    ] == [2, 3, 4]
    with pytest.raises(ValueError, match="unknown candidate"):
        da.inspect(pool, view="candidates", select=CandidateFilter(indices=(7,)))
    with pytest.raises(ReadLimitError, match="identities"):
        da.inspect(
            pool,
            view="candidates",
            select=predicate,
            read_limits=ReadLimits(identities=1),
        )


def test_candidate_reasons_retain_overlap_and_cannot_filter_parts(tmp_path: Path):

    pool = da.prepare(
        background_request().with_changes(
            screening=(
                planning.Avoid("aa", ("AA",)),
                planning.GC("gc", "sequence", 0.5, 1),
            )
        ),
        out=tmp_path / "pool",
    )
    query = da.inspect(
        pool, view="candidates", select=CandidateFilter(reasons=("gc",)), all=True
    )
    with query.records() as records:
        rows = list(records)
    assert len(rows) == 6
    assert all(r.outcome == "eligibility_rejected" for r in rows)
    assert all(set(r.candidate.reasons) == {"aa", "gc"} for r in rows)
    result = CliRunner().invoke(
        app,
        [
            "inspect",
            str(pool.path),
            "--view",
            "candidates",
            "--reason",
            "gc",
            "--candidate-index",
            "1",
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["records"] == [rows[0].to_dict()]
    with pytest.raises(ValueError, match=r"unknown.*reason"):
        da.inspect(pool, view="candidates", select=CandidateFilter(reasons=("typo",)))
    with pytest.raises(TypeError, match="PartFilter"):
        da.inspect(pool, view="parts", select=CandidateFilter())


@pytest.mark.parametrize("damage", ["checksum", "missing", "index"])
def test_candidate_reads_detect_damaged_evidence_before_export_publication(
    tmp_path: Path, damage: str
):

    pool = da.prepare(background_request(), out=tmp_path / "pool")
    with sqlite3.connect(pool.path / "pool.sqlite3") as connection:
        if damage == "missing":
            connection.execute("DELETE FROM candidates WHERE ordinal=6")
        else:
            value = json.loads(
                connection.execute(
                    "SELECT payload FROM candidates WHERE ordinal=2"
                ).fetchone()[0]
            )
            value["index"] = 9
            connection.execute(
                "UPDATE candidates SET payload=?,digest=? WHERE ordinal=2",
                (
                    canonical_json(value),
                    semantic_digest(value) if damage == "index" else "0" * 64,
                ),
            )
    out = tmp_path / "candidates.json"
    with pytest.raises(ValueError, match=r"digest|ordinal|count|checksum"):
        da.export(pool, view="candidates", all=True, out=out)
    assert not out.exists()


def test_curated_pools_do_not_invent_mined_candidate_evidence(tmp_path: Path):

    table = tmp_path / "parts.csv"
    table.write_text("part_id,sequence\na,AAA\n")
    pool = da.prepare(
        parts.PreparationSpec(parts.PartTable(table, format="csv")),
        out=tmp_path / "pool",
    )
    with pytest.raises(ValueError, match="no saved candidate evidence"):
        da.inspect(pool, view="candidates")


def test_saved_pool_report_reopens_after_source_removal(tmp_path: Path):

    pool = da.prepare(background_request(), out=tmp_path / "pool")
    live = da.inspect(pool, view="quality")
    expected = live.to_dict()
    saved = tmp_path / "report.json"
    da.export(live, out=saved)
    shutil.rmtree(pool.path)
    report = da.inspect(saved, view="quality")
    assert isinstance(report, PoolQualitySnapshot)
    assert report.to_dict() == expected
    assert report.cost.mode == "manifest"
    assert report.cost.records_estimate == 1
    assert da.inspect(report, view="quality") is report
    clone = tmp_path / "copy.json"
    receipt = da.export(report, out=clone)
    assert receipt.sources[0]["pool_id"] == pool.pool_id
    assert clone.read_bytes() == saved.read_bytes()
    result = CliRunner().invoke(
        app, ["inspect", str(saved), "--view", "quality", "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout) == expected
    human = CliRunner().invoke(app, ["inspect", str(saved), "--view", "quality"])
    assert human.exit_code == 0, human.output
    assert "recorded" in human.stdout.lower()
    assert "1 / 3" in human.stdout
    with pytest.raises(ValueError, match="verify"):
        da.inspect(saved, view="quality", verify=True)
    with pytest.raises(TypeError):
        da.inspect(saved, view="quality", compare=saved)


@pytest.mark.parametrize(
    "field,value",
    [
        ("state", "completed"),
        ("pool_id", "bad"),
        ("schema", "dense_arrays.pool_quality.v99"),
        ("extra", True),
        ("rejections", {"invented": 1}),
    ],
)
def test_saved_pool_report_rejects_inconsistent_or_unknown_content(
    tmp_path: Path, field: str, value: object
):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    report = da.inspect(pool, view="quality").to_dict()
    report[field] = value
    saved = tmp_path / "report.json"
    saved.write_text(json.dumps(report))
    with pytest.raises((ValueError, TypeError)):
        da.inspect(saved, view="quality")


def test_saved_pool_report_respects_state_bound_and_freezes_values(tmp_path: Path):

    pool = da.prepare(background_request(), out=tmp_path / "pool")
    data = da.inspect(pool, view="quality").to_dict()
    report = PoolQualitySnapshot.from_dict(data)
    data["counts"]["retained"] = 9
    assert report.to_dict()["counts"]["retained"] == 1
    with pytest.raises(ReadLimitError, match="identities"):
        PoolQualitySnapshot.from_dict(
            report.to_dict(), read_limits=ReadLimits(identities=1)
        )


def test_candidate_display_is_bounded_and_close_releases_reader(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    request = background_request().with_changes(
        sampling=parts.Sampling(planning.Length(exact=2000))
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    original = reading.reader
    opened = []

    @contextmanager
    def tracked(*args: object, **kwargs: object) -> Iterator[sqlite3.Connection]:
        with original(*args, **kwargs) as connection:
            opened.append(connection)
            try:
                yield connection
            finally:
                opened.remove(connection)

    monkeypatch.setattr(reading, "reader", tracked)
    query = da.inspect(pool, view="candidates", limit=1)
    rows = query.records()
    assert opened == []
    row = next(rows)
    assert len(repr(row)) < 200
    assert len(opened) == 1
    rows.close()
    assert opened == []
    result = CliRunner().invoke(
        app, ["inspect", str(pool.path), "--view", "candidates", "--limit", "1"]
    )
    assert result.exit_code == 0, result.output
    assert len(result.stdout) < 400
    assert "retained" in result.stdout
    assert "2000" in result.stdout
    assert "--json" in result.stdout


@pytest.mark.parametrize(
    "kwargs",
    [
        {"indices": (True,)},
        {"indices": (0,)},
        {"indices": (1, 1)},
        {"outcomes": ("accepted",)},
        {"reasons": ("",)},
    ],
)
def test_candidate_predicates_reject_ambiguous_values(kwargs: dict):
    with pytest.raises((ValueError, TypeError)):
        CandidateFilter(**kwargs)
