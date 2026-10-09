"""A real packing run publishes inspectable, coherent native evidence.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning


def request(*, count: int = 1, attempts: int = 20) -> planning.DesignSpec:
    """Use repeated DNA to distinguish path enumeration from sequence uniqueness."""
    return planning.DesignSpec(
        parts=[parts.Part("a", "AAA", "A"), parts.Part("b", "AAA", "A")],
        length=planning.Length(maximum=3),
        strands="single",
        target=planning.Target(count=count),
        limits=planning.Limits(attempts=attempts),
    )


def test_real_cbc_run_publishes_placements_and_reconciles_attempts(tmp_path: Path):
    result = da.run(request(), out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.state == "completed"
    assert summary.accepted == summary.target == 1
    assert summary.counts["started"] == summary.counts["accepted"] == 1
    with da.inspect(result, view="designs").records() as records:
        design = next(records)
        assert design.realized.sequence == "AAA"
        assert design.realized.placements[0].feature_id in {"a", "b"}
        assert design.realized.placements[0].start == 0
        assert design.realized.placements[0].end == 3
        assert design.plan_id == da.plan(request()).plan_id


def test_duplicate_and_exhaustion_are_not_completed_library(tmp_path: Path):
    result = da.run(request(count=3), out=tmp_path / "run")
    report = da.inspect(result, verify=True)
    assert report.state == "stopped"
    assert report.accepted == 1
    assert report.target == 3
    assert report.counts["duplicate"] == 1
    assert report.counts["no_candidate"] == 1
    assert report.counts["started"] == 3
    assert report.termination_reason == "batch_exhausted"
    assert report.resumable is False


def test_attempt_bound_stops_search_truthfully(tmp_path: Path):
    result = da.run(request(count=3, attempts=1), out=tmp_path / "run")
    report = da.inspect(result)
    assert report.termination_reason == "attempt_limit"
    assert report.accepted == 1
    assert report.counts["started"] == 1


def test_destination_collision_is_checked_before_solving(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    out = tmp_path / "run"
    out.mkdir()
    keep = out / "keep.txt"
    keep.write_text("untouched")

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("colliding destination started a solver")

    monkeypatch.setattr(da.Optimizer, "build_model", forbidden)
    with pytest.raises(FileExistsError):
        da.run(request(), out=out)
    assert keep.read_text() == "untouched"


def test_reading_does_not_reopen_source_or_run_solver(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    result = da.run(request(), out=tmp_path / "run")

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("inspection ran a solver")

    monkeypatch.setattr(da.Optimizer, "solve_report", forbidden)
    before = (result.path / "run.sqlite3").read_bytes()
    assert da.inspect(result, verify=True).accepted == 1
    assert (result.path / "run.sqlite3").read_bytes() == before


def test_missing_resume_target_is_not_created(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        da.run(resume=tmp_path / "run")
    assert not (tmp_path / "run").exists()
