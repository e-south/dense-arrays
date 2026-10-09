"""Usage preference preserves packing count and replays committed proposal history.

Author: Eric J. South.
"""

import json
import sqlite3
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.reading import ReadLimitError, ReadLimits
from dense_arrays.artifacts.store import RunWriter
from dense_arrays.cli import app
from dense_arrays.workflow.inputs import read_source


def request():
    return planning.DesignSpec(
        [parts.Part("a", "AAA"), parts.Part("b", "CCC"), parts.Part("c", "GGG")],
        planning.Length(maximum=3),
        strands="single",
        target=planning.Target(3),
        packing_preference="underused_parts",
    )


def attempts(path: Path) -> list:
    with da.inspect(path, view="attempts", all=True).records() as records:
        return list(records)


def test_preference_records_effective_weights_and_equivalent_python_cli(tmp_path: Path):
    plan = da.plan(request())
    assert plan.preview["packing_preference"] == "underused_parts"
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).plan_id == plan.plan_id
    run = da.run(plan, out=tmp_path / "python")
    records = attempts(run.path)
    assert da.inspect(run, verify=True).accepted == 3
    assert records[0].evidence["packing_objective"]["weights"] == {
        "a": 1,
        "b": 1,
        "c": 1,
    }
    chosen = records[0].candidate.packed.placements[0].feature_id
    weights = records[1].evidence["packing_objective"]["weights"]
    assert weights[chosen] == 1
    assert {weights[p] for p in weights if p != chosen} == {1 + 0.5 / 3}
    result = CliRunner().invoke(
        app, ["run", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 0, result.output
    assert [r.evidence["packing_objective"] for r in attempts(tmp_path / "cli")] == [
        r.evidence["packing_objective"] for r in records
    ]


def test_resume_restores_weights_without_counting_interrupted_attempt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    original = da.Optimizer.solve_report
    calls = 0

    def interrupt(
        self: da.Optimizer, *args: object, **kwargs: object
    ) -> da.solver.SolveReport:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise KeyboardInterrupt
        return original(self, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(da.Optimizer, "solve_report", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(request(), out=tmp_path / "resumed")
    resumed = da.run(resume=tmp_path / "resumed")
    direct = da.run(request(), out=tmp_path / "direct")
    resumed_records = [r for r in attempts(resumed.path) if r.candidate is not None]
    direct_records = attempts(direct.path)
    assert da.inspect(resumed, verify=True).accepted == 3
    assert [r.evidence["packing_objective"] for r in resumed_records] == [
        r.evidence["packing_objective"] for r in direct_records
    ]
    assert [r.candidate.packed.sequence for r in resumed_records] == [
        r.candidate.packed.sequence for r in direct_records
    ]


def test_invalid_preference_and_plain_plan_encoding():
    plain = request().with_changes(packing_preference=None)
    assert "packing_preference" not in da.plan(plain).to_dict()["request"]
    with pytest.raises(ValueError, match="packing_preference"):
        plain.with_changes(packing_preference="sequence_diversity")


def test_rejected_packings_contribute_to_proposed_usage(tmp_path: Path):
    screened = request().with_changes(
        requirements=(planning.Avoid("no_aaa", ("AAA",)),)
    )
    run = da.run(screened, out=tmp_path / "screened")
    records = attempts(run.path)
    rejected = next(i for i, r in enumerate(records) if r.outcome == "rejected")
    following = records[rejected + 1].evidence["packing_objective"]
    assert following["usage"]["a"] == 1
    assert following["proposed_packings"] == rejected + 1
    assert da.inspect(run, verify=True).accepted == 2


@pytest.mark.parametrize("change", ["rewind", "remove"])
def test_verification_rejects_objective_history_tampering(tmp_path: Path, change: str):
    run = da.run(request(), out=tmp_path / "run")
    first = attempts(run.path)[0].to_dict()["evidence"]["packing_objective"]
    with sqlite3.connect(run.path / "run.sqlite3") as connection:
        revision, payload = connection.execute(
            "SELECT revision,payload FROM attempts WHERE attempt=2 "
            "ORDER BY revision DESC LIMIT 1"
        ).fetchone()
        value = json.loads(payload)
        if change == "rewind":
            value["evidence"]["packing_objective"] = first
        else:
            del value["evidence"]["packing_objective"]
        connection.execute(
            "UPDATE attempts SET payload=?,digest=? WHERE attempt=2 AND revision=?",
            (canonical_json(value), semantic_digest(value), revision),
        )
    with pytest.raises(ValueError, match="objective"):
        da.inspect(run, verify=True)


def test_scheduled_batches_reset_only_their_own_usage(tmp_path: Path):
    plan = da.plan(request())
    batches = tuple(
        planning.CandidateBatch(ids, plan.collection_id, stream=f"step{i}")
        for i, ids in enumerate((("a", "b"), ("b", "c")))
    )
    plan = da.plan(
        plan.request.with_changes(
            schedule=planning.BatchSchedule(batches, attempts_per_batch=2)
        )
    )
    run = da.run(plan, out=tmp_path / "run")
    records = attempts(run.path)
    firsts = [r for r in records if r.evidence["batch_attempt"] == 1]
    assert len(firsts) == 2
    assert [set(r.evidence["packing_objective"]["weights"]) for r in firsts] == [
        {"a", "b"},
        {"b", "c"},
    ]
    assert all(
        r.evidence["packing_objective"]["proposed_packings"] == 0 for r in firsts
    )
    assert set(records[1].evidence["packing_objective"]["weights"].values()) == {
        1,
        1.25,
    }
    da.inspect(run, verify=True)


def test_matrix_cells_keep_independent_proposal_usage(tmp_path: Path):
    matrix = planning.MatrixSpec(
        request().with_changes(target=planning.Target()),
        axes={"setting": {"first": planning.Variant(), "second": planning.Variant()}},
        allocation=planning.Allocation(per_cell=2),
        max_cells=2,
    )
    run = da.run(matrix, out=tmp_path / "matrix")
    assert da.inspect(run, verify=True).accepted == 4
    records = attempts(run.path)
    assert [r.evidence["packing_objective"]["proposed_packings"] for r in records] == [
        0,
        0,
        1,
        1,
    ]


def test_usage_bonus_cannot_overcome_one_additional_occurrence():
    from dense_arrays.model import part_usage_weights  # noqa: PLC0415

    weights = part_usage_weights((100, 100, 0))
    assert weights == (1, 1, 1 + 0.5 / 3)
    assert weights[0] + weights[1] > weights[2]
    assert part_usage_weights((3, 3, 3)) == (1, 1, 1)


def test_usage_weights_match_released_iterator_fixture():
    from dense_arrays.model import part_usage_weights  # noqa: PLC0415

    fixture = json.loads(
        (
            Path(__file__).parents[1]
            / "fixtures/workflow/dense-arrays-usage-weights-v1.json"
        ).read_text()
    )
    for case in fixture["cases"]:
        assert part_usage_weights(case["counts"]) == tuple(case["weights"])


@pytest.mark.parametrize("mode", ["schedule", "resampling", "matrix"])
def test_committed_preference_replays_across_execution_modes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
):
    spec = request()
    if mode == "schedule":
        source = da.plan(spec)
        batches = tuple(
            planning.CandidateBatch(
                ("a", "b", "c"), source.collection_id, stream=f"batch{i}"
            )
            for i in range(2)
        )
        spec = spec.with_changes(
            schedule=planning.BatchSchedule(batches, attempts_per_batch=3)
        )
    elif mode == "resampling":
        spec = spec.with_changes(
            resampling=planning.Resampling(
                planning.BatchSampling(2, seed=7), max_batches=10, attempts_per_batch=2
            )
        )
    else:
        spec = planning.MatrixSpec(
            spec.with_changes(target=planning.Target()),
            axes={
                "setting": {"first": planning.Variant(), "second": planning.Variant()}
            },
            allocation=planning.Allocation(per_cell=2),
            max_cells=2,
        )
    publish = RunWriter.publish

    def interrupt(self: RunWriter, *args: object, **kwargs: object) -> None:
        publish(self, *args, **kwargs)
        raise KeyboardInterrupt

    with monkeypatch.context() as patch:
        patch.setattr(RunWriter, "publish", interrupt)
        with pytest.raises(KeyboardInterrupt):
            da.run(spec, out=tmp_path / "resumed")
    first = attempts(tmp_path / "resumed")[0].to_dict()
    resumed = da.run(resume=tmp_path / "resumed")
    direct = da.run(spec, out=tmp_path / "direct")
    a, b = attempts(resumed.path), attempts(direct.path)
    assert a[0].to_dict() == first
    assert [r.evidence["packing_objective"] for r in a] == [
        r.evidence["packing_objective"] for r in b
    ]
    assert (
        da.inspect(resumed, verify=True).accepted
        == da.inspect(direct, verify=True).accepted
    )


def test_plan_comparison_and_portable_bundle_preserve_preference(tmp_path: Path):
    plan = da.plan(request())
    plain = da.plan(request().with_changes(packing_preference=None))
    comparison = da.inspect(plan, compare=plain).to_dict()
    assert "packing_preference" in comparison["changed_fields"]
    run = da.run(plan, out=tmp_path / "run")
    da.export(run, all=True, format="bundle", out=tmp_path / "bundle")
    assert da.inspect(tmp_path / "bundle", verify=True).designs == 3


def test_attempt_objective_maps_obey_identity_read_limits(tmp_path: Path):
    run = da.run(request(), out=tmp_path / "run")
    view = da.inspect(
        run, view="attempts", limit=1, read_limits=ReadLimits(identities=5)
    )
    with pytest.raises(ReadLimitError), view.records() as records:
        next(records)


def test_missing_request_in_preference_evidence_fails_as_input_error():
    value = da.plan(request()).evidence.to_dict()
    del value["content"]["request"]
    with pytest.raises((ValueError, TypeError)):
        planning.PlanEvidence.from_dict(value)


@pytest.mark.parametrize("operation", ["plan", "inspect"])
@pytest.mark.parametrize("matrix", [False, True])
def test_human_preview_explains_preference_scope(
    tmp_path: Path, operation: str, matrix: bool
):
    spec = request()
    if matrix:
        spec = planning.MatrixSpec(
            spec.with_changes(target=planning.Target()),
            axes={"setting": {"first": planning.Variant()}},
            allocation=planning.Allocation(per_cell=2),
            max_cells=1,
        )
    saved = tmp_path / "plan.json"
    da.plan(spec).write(saved)
    args = [operation, str(saved)]
    if operation == "inspect":
        args += ["--view", "plan"]
    result = CliRunner().invoke(app, args)
    assert result.exit_code == 0, result.output
    assert "Packing preference: underused parts" in result.output
    assert "proposed usage within each offered batch" in result.output


@pytest.mark.parametrize(
    "status", ["unproven", "unknown", "backend_error", "invalid_result"]
)
def test_preference_preserves_nonoptimal_solver_outcomes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, status: str
):
    from dense_arrays.solver import SolveReport, SolveStatus  # noqa: PLC0415

    monkeypatch.setattr(
        da.Optimizer,
        "solve_report",
        lambda *_a, **_kw: SolveReport(SolveStatus(status), None),
    )
    run = da.run(request(), out=tmp_path / "run")
    summary = da.inspect(run, verify=True)
    assert summary.accepted == 0
    assert summary.state != "completed"
    record = attempts(run.path)[0]
    assert record.evidence["solver_status"] == status
    assert record.candidate is None
    assert record.evidence["packing_objective"]["proposed_packings"] == 0
