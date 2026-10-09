"""Runtime batches use committed feedback and retain independent cell frontiers.

Author: Eric J. South.
"""

import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.reading import ReadLimitError, ReadLimits
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app


def request() -> planning.DesignSpec:
    return planning.DesignSpec(
        parts=tuple(
            parts.Part(k, s, group=k)
            for k, s in zip("abc", ("AAA", "CCC", "GGG"), strict=True)
        ),
        length=planning.Length(maximum=3),
        strands="single",
        target=planning.Target(3),
        resampling=planning.Resampling(
            sampling=planning.BatchSampling(1, seed=7),
            max_batches=10,
            attempts_per_batch=1,
            accepted_per_batch=1,
            feedback=planning.FeedbackPolicy(coverage_alpha=1e6, coverage_power=20),
        ),
    )


def test_runtime_feedback_changes_between_committed_batches(tmp_path: Path):
    plan = da.plan(request())
    result = da.run(plan, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 3
    assert summary.batch_count == 3
    decisions = list(da.inspect(result, view="batches", all=True).records())
    assert [d.index for d in decisions] == [1, 2, 3]
    assert [sum(d.batch.feedback.used.values()) for d in decisions] == [0, 1, 2]
    assert {d.batch.part_ids[0] for d in decisions} == {"a", "b", "c"}
    assert all(d.plan_id == plan.plan_id for d in decisions)
    da.export(result, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == 3


@pytest.mark.parametrize("boundary", ["reserve", "publish"])
def test_runtime_resume_keeps_saved_decisions_and_feedback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, boundary: str
):
    original = getattr(RunWriter, boundary)

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> object:
        result = original(self, *args, **kwargs)
        if self.state["counts"]["started"] == 2:
            raise KeyboardInterrupt
        return result

    path = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, boundary, interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request(), out=path)
    before = [d.to_dict() for d in da.inspect(path, view="batches", all=True).records()]
    assert len(before) == 2
    result = da.run(resume=path)
    summary = da.inspect(result, verify=True)
    assert summary.state == "completed"
    after = [
        d.to_dict() for d in da.inspect(result, view="batches", all=True).records()
    ]
    assert after[:2] == before
    assert summary.counts["interrupted_unresolved"] == (boundary == "reserve")


def test_batch_limit_is_not_global_infeasibility(tmp_path: Path):
    policy = planning.Resampling(
        planning.BatchSampling(1), max_batches=2, attempts_per_batch=1
    )
    result = da.run(request().with_changes(resampling=policy), out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.state == "stopped"
    assert summary.batch_count == 2
    assert summary.termination_reason == "batch_limit"


def test_matrix_resampling_is_independent_and_cli_uses_same_plan(tmp_path: Path):
    plan = da.plan(
        planning.MatrixSpec(
            base=request().with_changes(target=planning.Target()),
            axes={"x": {"a": planning.Variant(), "b": planning.Variant()}},
            allocation=planning.Allocation(per_cell=3),
            max_cells=2,
        )
    )
    plan.write(tmp_path / "plan.json")
    response = CliRunner().invoke(
        app,
        ["run", str(tmp_path / "plan.json"), "--out", str(tmp_path / "run"), "--json"],
    )
    assert response.exit_code == 0, response.output
    summary = da.inspect(tmp_path / "run", verify=True)
    assert summary.accepted == 6
    decisions = list(da.inspect(tmp_path / "run", view="batches", all=True).records())
    assert len(decisions) == 6
    assert [sum(d.batch.feedback.used.values()) for d in decisions] == [
        0,
        0,
        1,
        1,
        2,
        2,
    ]


def test_runtime_batch_reader_enforces_identity_bounds(tmp_path: Path):

    result = da.run(request(), out=tmp_path / "run")
    view = da.inspect(
        result, view="batches", all=True, read_limits=ReadLimits(identities=1)
    )
    with pytest.raises(ReadLimitError):
        list(view.records())


def test_portable_batch_pages_survive_original_removal(tmp_path: Path):

    result = da.run(request(), out=tmp_path / "run")
    da.export(result, all=True, format="bundle", out=tmp_path / "bundle")
    shutil.rmtree(tmp_path / "run")
    page = da.inspect(tmp_path / "bundle", view="batches", limit=2)
    with page.records() as rows:
        first = list(rows)
        cursor = rows.next_cursor
    last = list(
        da.inspect(tmp_path / "bundle", view="batches", after=cursor, limit=2).records()
    )
    assert [d.index for d in first + last] == [1, 2, 3]
    da.export(tmp_path / "bundle", format="bundle", all=True, out=tmp_path / "copy")
    assert da.inspect(tmp_path / "copy", verify=True).designs == 3
