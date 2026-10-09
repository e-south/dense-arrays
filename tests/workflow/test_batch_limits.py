"""Accepted-design caps limit contribution without changing targets or effort.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app
from dense_arrays.workflow.inputs import read_source


def source() -> planning.GenerationPlan:
    return da.plan(
        planning.DesignSpec(
            parts=[
                parts.Part(k, s)
                for k, s in zip("abcd", ("AAA", "CCC", "GGG", "TTT"), strict=True)
            ],
            length=planning.Length(maximum=3),
            strands="single",
            target=planning.Target(4),
        )
    )


def schedule(*members: tuple[str, ...], accepted: int = 1) -> planning.GenerationPlan:
    plan = source()
    batches = tuple(
        planning.CandidateBatch(ids, plan.collection_id, stream=f"step/{i}")
        for i, ids in enumerate(members)
    )
    return da.plan(
        plan.request.with_changes(
            schedule=planning.BatchSchedule(
                batches, attempts_per_batch=5, accepted_per_batch=accepted
            )
        )
    )


def test_accepted_cap_advances_before_enumeration_is_exhausted(tmp_path: Path):
    plan = schedule(("a", "b"), ("c", "d"))
    result = da.run(plan, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == summary.counts["started"] == 2
    assert summary.target == 4
    assert summary.termination_reason == "batch_schedule_exhausted"
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    assert [a.evidence["batch_index"] for a in attempts] == [1, 2]
    da.export(result, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == 2


def test_duplicates_and_rejections_do_not_consume_accepted_cap(tmp_path: Path):
    plan = schedule(("a",), ("a",), ("b",), ("c",))
    result = da.run(
        plan.request.with_changes(
            requirements=(planning.Avoid("no_c", patterns=("GGG",), strands="forward"),)
        ),
        out=tmp_path / "run",
    )
    summary = da.inspect(result, verify=True)
    assert summary.accepted == 2
    assert summary.counts["duplicate"] == summary.counts["rejected"] == 1
    assert summary.counts["started"] == 6
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    assert [a.evidence["batch_index"] for a in attempts] == [1, 2, 2, 3, 4, 4]


@pytest.mark.parametrize("accepted", [1, 2])
def test_resume_recounts_accepted_cap_without_resetting_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, accepted: int
):
    plan = schedule(("a", "b"), ("c", "d"), accepted=accepted)
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> str:
        result = publish(self, *args, **kwargs)
        if self.state["counts"]["started"] == 1:
            raise KeyboardInterrupt
        return result

    path = tmp_path / "run"
    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(plan, out=path)
    result = da.run(resume=path)
    summary = da.inspect(result, verify=True)
    assert summary.accepted == summary.counts["started"] == 2 * accepted
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    assert [a.evidence["batch_index"] for a in attempts] == [1] * accepted + [
        2
    ] * accepted


def test_prepare_accepted_cap_matches_cli_and_round_trips(tmp_path: Path):
    plan = source()
    plan.write(tmp_path / "source.json")
    result = da.prepare(
        plan,
        sampling=planning.BatchSampling(2),
        batch_count=3,
        attempts_per_batch=5,
        accepted_per_batch=1,
        out=tmp_path / "api.json",
    )
    response = CliRunner().invoke(
        app,
        [
            "prepare",
            str(tmp_path / "source.json"),
            "--batch-size",
            "2",
            "--batch-count",
            "3",
            "--attempts-per-batch",
            "5",
            "--accepted-per-batch",
            "1",
            "--out",
            str(tmp_path / "cli.json"),
            "--json",
        ],
    )
    assert response.exit_code == 0, response.output
    assert read_source(tmp_path / "cli.json").plan_id == result.plan_id
    assert result.preview["accepted_per_batch"] == 1
    preview = CliRunner().invoke(app, ["plan", str(tmp_path / "cli.json")])
    assert "accepted per batch: 1" in preview.stdout


@pytest.mark.parametrize("limit", [0, -1, True, 1.5])
def test_invalid_accepted_caps_fail_before_publication(tmp_path: Path, limit: object):
    with pytest.raises((TypeError, ValueError), match="accepted_per_batch"):
        da.prepare(
            source(),
            sampling=planning.BatchSampling(2),
            batch_count=2,
            attempts_per_batch=5,
            accepted_per_batch=limit,
            out=tmp_path / "plan.json",
        )
    assert not (tmp_path / "plan.json").exists()


def test_accepted_cap_requires_a_declared_attempt_allowance(tmp_path: Path):
    with pytest.raises(ValueError, match="attempts_per_batch"):
        da.prepare(
            source(),
            sampling=planning.BatchSampling(2),
            accepted_per_batch=1,
            out=tmp_path / "plan.json",
        )
    assert not (tmp_path / "plan.json").exists()


def test_uncapped_schedule_keeps_its_existing_wire_identity():
    plan = source()
    batch = planning.CandidateBatch(("a",), plan.collection_id)
    value = {
        "schema": "dense_arrays.batch_schedule.v1",
        "batches": [batch.to_dict()],
        "attempts_per_batch": 5,
    }
    restored = planning.BatchSchedule.from_dict(value)
    assert restored.accepted_per_batch is None
    assert restored.to_dict() == value


def test_matrix_cells_apply_their_own_accepted_caps(tmp_path: Path):
    capped = schedule(("a", "b"), ("c", "d"), accepted=1)
    wider = schedule(("a", "b"), ("c", "d"), accepted=2)
    request = planning.MatrixSpec(
        base=capped.request.with_changes(schedule=None, target=planning.Target()),
        axes={"x": {"a": planning.Variant(), "b": planning.Variant()}},
        max_cells=2,
        allocation=planning.Allocation(per_cell=4),
        batches={"x=a": capped.request.schedule, "x=b": wider.request.schedule},
    )
    result = da.run(request, out=tmp_path / "run")
    summary = da.inspect(result, verify=True)
    assert summary.accepted == summary.counts["started"] == 6
    assert summary.termination_reason == "cells_exhausted"
    assert summary.cells["x=a"].termination_reason == "batch_schedule_exhausted"
    assert summary.cells["x=b"].state == "completed"


def test_verification_rejects_work_past_an_accepted_cap(tmp_path: Path):
    from dense_arrays.reporting.batch_accounting import BatchAccounting  # noqa: PLC0415

    plan = schedule(("a", "b"), ("c", "d"), accepted=2)
    result = da.run(plan, out=tmp_path / "run")
    attempts = list(da.inspect(result, view="attempts", all=True).records())
    accounting = BatchAccounting(
        planning.BatchSchedule(plan.request.schedule.batches, 5, accepted_per_batch=1)
    )
    accounting.observe(attempts[0])
    with pytest.raises(ValueError, match="batch order"):
        accounting.observe(attempts[1])


def test_accepted_cap_matches_pinned_stopping_predicate(tmp_path: Path):
    fixture = json.loads(
        (
            Path(__file__).parents[1]
            / "fixtures/workflow/densegen-accepted-cap-v1.json"
        ).read_text()
    )
    for index, case in enumerate(fixture["generation_cases"]):
        plan = schedule(("a", "b"), ("c", "d"), accepted=case["cap"])
        result = da.run(
            plan.request.with_changes(target=planning.Target(case["target"])),
            out=tmp_path / str(index),
        )
        summary = da.inspect(result, verify=True)
        assert summary.accepted == case["expected_accepted"]
        assert summary.counts["started"] == case["expected_accepted"]
