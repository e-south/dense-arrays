"""Sampled recipes keep planning, effort, selection and pool evidence separate.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.workflow import preparation
from dense_arrays.workflow.inputs import read_source


def background_request():
    return parts.PreparationSpec(
        source=parts.Background(base_probabilities=(1, 0, 0, 0)),
        sampling=parts.Sampling(length=planning.Length(exact=4)),
        budget=parts.CandidateBudget(candidates=6),
        uniqueness=parts.Uniqueness(key="sequence"),
        retain=parts.Retention(count=3, policy="first_eligible"),
        seed=7,
    )


def test_interrupted_preparation_explains_that_no_pool_was_committed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def interrupted(**_kwargs: object) -> None:
        raise KeyboardInterrupt

    saved = tmp_path / "prepare.json"
    da.plan(background_request()).write(saved)
    monkeypatch.setattr(preparation, "sample_sequence", interrupted)
    out = tmp_path / "pool"
    result = CliRunner().invoke(app, ["prepare", str(saved), "--out", str(out)])
    assert result.exit_code == 130, result.output
    assert "no pool was committed" in result.stderr
    assert str(out) in result.stderr
    assert "new destination" in result.stderr
    assert "committed run prefix" not in result.stderr
    with pytest.raises(ValueError, match=r"manifest|commit"):
        da.inspect(out)


def test_uniqueness_rejects_an_unsupported_cross_group_policy_in_python():
    with pytest.raises(TypeError, match="cross_group_collisions"):
        parts.Uniqueness(key="core", cross_group_collisions="error")


def test_preparation_parser_rejects_an_unsupported_cross_group_policy():
    from dense_arrays.planning.preparation.requests import (  # noqa: PLC0415
        preparation_from_dict,
        preparation_to_dict,
    )

    encoded = preparation_to_dict(background_request())
    encoded["uniqueness"]["cross_group_collisions"] = "error"
    with pytest.raises(ValueError, match="cross_group_collisions"):
        preparation_from_dict(encoded)


@pytest.mark.parametrize("key", ["sequence", "core"])
def test_preparation_uniqueness_encodes_only_its_per_recipe_key(key: str):
    from dense_arrays.planning.preparation.requests import (  # noqa: PLC0415
        preparation_from_dict,
        preparation_to_dict,
    )

    request = background_request().with_changes(
        source=parts.PWMArtifact("motif.json"),
        scoring=parts.FimoScoring(),
        uniqueness=parts.Uniqueness(key=key),
    )
    encoded = preparation_to_dict(request)
    assert encoded["uniqueness"] == {"key": key}
    assert preparation_from_dict(encoded) == request


@pytest.mark.parametrize("operation", ["plan", "prepare"])
def test_missing_scorer_has_a_structured_cli_failure(tmp_path: Path, operation: str):
    from .test_motifs import artifact  # noqa: PLC0415

    motif = tmp_path / "motif.json"
    motif.write_text(json.dumps(artifact()))
    recipe = tmp_path / "prepare.json"
    recipe.write_text(
        json.dumps(
            {
                "schema": "dense_arrays.prepare.v1",
                "source": {"kind": "pwm_artifact", "path": "motif.json"},
                "sampling": {"length": {"exact": 4}},
                "budget": {"candidates": 6},
                "scoring": {"backend": "fimo", "executable": "missing-fimo"},
                "retain": {"count": 2, "policy": "first_eligible"},
            }
        )
    )
    output = tmp_path / "output"
    result = CliRunner().invoke(
        app, [operation, str(recipe), "--out", str(output), "--json"]
    )
    assert result.exit_code == 4, result.output
    error = json.loads(result.stdout)
    assert error["schema"] == "dense_arrays.error.v1"
    assert error["code"] == "scoring_error"
    assert error["reason"] == "unavailable"
    assert "install FIMO" in error["message"]
    assert not output.exists()


def test_sampled_preview_freezes_unknown_yield_without_generating(tmp_path: Path):
    request = background_request()
    resolved = da.plan(request)
    assert isinstance(resolved, planning.PreparationPlan)
    assert resolved.preview["candidate_budget"] == 6
    assert resolved.preview["requested_retention"] == 3
    assert resolved.preview["retained_parts"] is None
    assert resolved.preview["retained_count_status"] == "unknown"
    assert resolved.preview["required_tools"] == ()
    assert list(tmp_path.iterdir()) == []
    saved = tmp_path / "recipe.plan.json"
    resolved.write(saved)
    restored = read_source(saved)
    assert restored.to_dict() == resolved.to_dict()
    assert restored.request == request
    assert len(repr(restored)) < 300
    result = CliRunner().invoke(app, ["plan", str(saved)])
    assert result.exit_code == 0, result.output
    assert "unknown" in result.stdout
    assert "6" in result.stdout
    assert "3" in result.stdout


@pytest.mark.parametrize(
    "change",
    [
        {"sampling": None},
        {"budget": None},
        {"sampling": lambda: parts.Sampling(length=planning.Length(maximum=4))},
        {
            "retain": lambda: parts.Retention(
                count=3, policy="top_score", rank_by="best_hit_score"
            )
        },
        {"uniqueness": lambda: parts.Uniqueness(key="core")},
    ],
)
def test_background_rejects_unresolvable_or_inapplicable_policies(change: dict):
    with pytest.raises((TypeError, ValueError)):
        background_request().with_changes(
            **{k: v() if callable(v) else v for k, v in change.items()}
        )


def test_background_pool_records_shortfall_and_reconciled_stages(tmp_path: Path):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    summary = da.inspect(pool, verify=True)
    assert summary.state == "incomplete"
    assert summary.source_parts == 6
    assert summary.retained_parts == 1
    report = da.inspect(pool, view="quality").to_dict()
    assert report["counts"] == {
        "processed": 6,
        "eligibility_rejected": 0,
        "eligible": 6,
        "execution_error": 0,
        "duplicate_discarded": 5,
        "eligible_unique": 1,
        "retained": 1,
        "not_selected": 0,
    }
    assert report["requested_retention"] == 3
    assert report["stop_reason"] == "candidate_budget"
    with da.inspect(pool, view="parts", all=True).records() as rows:
        records = list(rows)
    assert [r.part.sequence for r in records] == ["AAAA"]
    output = tmp_path / "quality.json"
    da.export(pool, view="quality", out=output)
    assert json.loads(output.read_text())["counts"] == report["counts"]
    paired = da.prepare(background_request(), out=tmp_path / "repeat")
    assert paired.pool_id == pool.pool_id
    design = planning.DesignSpec(
        parts.PoolSource(pool), planning.Length(maximum=4), strands="single"
    )
    with pytest.raises(ValueError, match="pool is incomplete"):
        da.run(design, out=tmp_path / "design")
    completed = da.prepare(
        background_request().with_changes(
            retain=parts.Retention(count=1, policy="first_eligible")
        ),
        out=tmp_path / "completed",
    )
    made = da.run(
        design.with_changes(parts=parts.PoolSource(completed)), out=tmp_path / "design"
    )
    assert da.inspect(made, verify=True).accepted == 1


def test_background_uses_shared_sequence_screens_and_reports_all_reasons(
    tmp_path: Path,
):
    recipe = background_request().with_changes(
        screening=(
            planning.Avoid("no_aa", ("AA",), strands="forward"),
            planning.GC("gc", "sequence", 0.5, 1),
        )
    )
    pool = da.prepare(recipe, out=tmp_path / "pool")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["counts"]["eligibility_rejected"] == 6
    assert report["rejections"] == {"no_aa": 6, "gc": 6}
    assert report["counts"]["eligible"] == 0
    assert da.inspect(pool, verify=True).retained_parts == 0


def test_background_cli_matches_python_and_returns_visible_shortfall(tmp_path: Path):
    recipe = tmp_path / "prepare.yaml"
    recipe.write_text(
        "schema: dense_arrays.prepare.v1\n"
        "source: {kind: background, base_probabilities: [1, 0, 0, 0]}\n"
        "sampling: {strategy: stochastic, length: {exact: 4}}\n"
        "budget: {candidates: 6}\n"
        "uniqueness: {key: sequence}\n"
        "retain: {count: 3, policy: first_eligible}\n"
        "seed: 7\n"
    )
    runner = CliRunner()
    plan_path = tmp_path / "plan.json"
    result = runner.invoke(
        app, ["plan", str(recipe), "--out", str(plan_path), "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["plan_id"] == da.plan(background_request()).plan_id
    target = tmp_path / "cli-pool"
    result = runner.invoke(
        app, ["prepare", str(plan_path), "--out", str(target), "--json"]
    )
    assert result.exit_code == 3, result.output
    assert json.loads(result.stdout)["state"] == "incomplete"
    inspected = runner.invoke(
        app, ["inspect", str(target), "--view", "quality", "--json"]
    )
    assert inspected.exit_code == 0, inspected.output
    assert json.loads(inspected.stdout)["counts"]["retained"] == 1


def test_pwm_preview_and_pool_share_bound_optional_scoring(tmp_path: Path):
    from .test_fimo import scorer  # noqa: PLC0415
    from .test_motifs import artifact  # noqa: PLC0415

    value = artifact()
    value["probabilities"] = [
        {"A": 1, "C": 0, "G": 0, "T": 0},
        {"A": 0, "C": 1, "G": 0, "T": 0},
    ]
    motif = tmp_path / "motif.json"
    motif.write_text(json.dumps(value))
    rows = "".join(
        f"motif\t\tcandidate_{i}\t1\t2\t+\t3\t0.0625\t\tAC\n" for i in range(6)
    )
    tool = scorer(tmp_path, rows)
    recipe = parts.PreparationSpec(
        parts.PWMArtifact(motif),
        parts.Retention(count=2, policy="top_score", rank_by="best_hit_score"),
        sampling=parts.Sampling(planning.Length(exact=2)),
        budget=parts.CandidateBudget(candidates=6),
        scoring=parts.FimoScoring(executable=tool, hit_pvalue_max=0.1),
        eligibility=parts.Eligibility(best_hit_score_min_exclusive=0),
        uniqueness=parts.Uniqueness(key="core"),
        seed=7,
    )
    plan = da.plan(recipe)
    assert plan.preview["required_tools"] == ("fimo",)
    path = tmp_path / "plan.json"
    plan.write(path)
    restored = read_source(path)
    assert restored.to_dict() == plan.to_dict()
    pool = da.prepare(restored, out=tmp_path / "pool")
    assert da.inspect(pool, verify=True).retained_parts == 1
    with da.inspect(pool, view="parts", all=True).records() as records:
        selected = next(records).part
    assert selected.core_sequence == "AC"
    assert selected.metadata["score"]["raw"] == 3
    motif.unlink()
    tool.unlink()
    assert da.inspect(pool, verify=True).retained_parts == 1
    assert da.inspect(pool, view="plan").plan_id == plan.plan_id


def test_sampling_never_selects_zero_probability_bases(monkeypatch: pytest.MonkeyPatch):
    from dense_arrays.parts import mining  # noqa: PLC0415

    class Entropy:
        def digest(self, size: int) -> bytes:
            return b"\xff" * size

    monkeypatch.setattr(mining.hashlib, "shake_256", lambda _: Entropy())
    assert (
        mining.sample_sequence(length=4, seed=0, index=1, probabilities=(1, 0, 0, 0))
        == "AAAA"
    )


def test_pool_quality_has_human_output_and_inspection_work_caps(tmp_path: Path):
    from dense_arrays.artifacts.reading import (  # noqa: PLC0415
        ReadLimitError,
        ReadLimits,
    )

    pool = da.prepare(background_request(), out=tmp_path / "pool")
    result = CliRunner().invoke(app, ["inspect", str(pool.path), "--view", "quality"])
    assert result.exit_code == 0, result.output
    assert "1 / 3" in result.stdout
    assert "5" in result.stdout
    assert "duplicate" in result.stdout
    report = da.inspect(pool, view="quality", read_limits=ReadLimits(records=2))
    assert report.cost.records_estimate == 8
    with pytest.raises(ReadLimitError, match="records"):
        report.to_dict()
    with pytest.raises(ReadLimitError, match="identities"):
        da.inspect(pool, view="quality", read_limits=ReadLimits(identities=2)).to_dict()


def test_sampling_prefix_does_not_change_with_candidate_batch_size(tmp_path: Path):
    from dataclasses import replace  # noqa: PLC0415

    recipe = background_request().with_changes(
        source=parts.Background(),
        retain=parts.Retention(count=6, policy="first_eligible"),
    )
    one = da.prepare(recipe, out=tmp_path / "one")
    two = da.prepare(
        recipe.with_changes(budget=replace(recipe.budget, batch_size=1)),
        out=tmp_path / "two",
    )
    with da.inspect(one, view="parts", all=True).records() as rows:
        a = [r.part for r in rows]
    with da.inspect(two, view="parts", all=True).records() as rows:
        b = [r.part for r in rows]
    assert a == b


@pytest.mark.parametrize("target", [None, parts.MiningTarget(eligible_unique=2)])
def test_scorer_failure_is_recorded_separately_from_eligibility(
    tmp_path: Path, target: parts.MiningTarget | None
):
    from .test_fimo import executable  # noqa: PLC0415
    from .test_motifs import source  # noqa: PLC0415

    tool = executable(
        tmp_path,
        'if [ "$1" = "--version" ]; then printf "5.5.9\\n"; else printf '
        '"database failed\\n" >&2; exit 17; fi',
    )
    recipe = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(count=2, policy="top_score", rank_by="best_hit_score"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(6, batch_size=3),
        scoring=parts.FimoScoring(executable=tool),
        mining_target=target,
    )
    plan = da.plan(recipe)
    pool = da.prepare(plan, out=tmp_path / "pool")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["stop_reason"] == "execution_error"
    assert report["counts"]["execution_error"] == 3
    assert report["counts"]["eligibility_rejected"] == 0
    assert report["counts"]["processed"] == 3
    assert da.inspect(pool, verify=True).state == "incomplete"
    saved = tmp_path / "plan.json"
    plan.write(saved)
    result = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 4, result.output
    assert json.loads(result.stdout)["preparation"]["counts"]["execution_error"] == 3


def test_scoring_work_is_rejected_during_planning_before_mining(tmp_path: Path):
    from .test_fimo import scorer  # noqa: PLC0415
    from .test_motifs import source  # noqa: PLC0415

    recipe = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(count=2, policy="top_score", rank_by="best_hit_score"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(6),
        scoring=parts.FimoScoring(
            executable=scorer(tmp_path, ""), limits=parts.ScoringLimits(windows=2)
        ),
    )
    with pytest.raises(ValueError, match="window limit"):
        da.plan(recipe)


def test_saved_candidate_tampering_is_detected_even_with_a_new_row_checksum(
    tmp_path: Path,
):
    import sqlite3  # noqa: PLC0415

    from dense_arrays._record_validation import (  # noqa: PLC0415
        canonical_json,
        semantic_digest,
    )

    pool = da.prepare(background_request(), out=tmp_path / "pool")
    with sqlite3.connect(pool.path / "pool.sqlite3") as connection:
        payload = json.loads(
            connection.execute(
                "SELECT payload FROM candidates WHERE ordinal=2"
            ).fetchone()[0]
        )
        payload["representative"] = 2
        connection.execute(
            "UPDATE candidates SET payload=?,digest=? WHERE ordinal=2",
            (canonical_json(payload), semantic_digest(payload)),
        )
    with pytest.raises(ValueError, match=r"selection|accounting"):
        da.inspect(pool, verify=True)


@pytest.mark.parametrize("target", [None, parts.MiningTarget(eligible_unique=2)])
def test_time_budget_records_only_processed_candidates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, target: parts.MiningTarget | None
):
    from dense_arrays.workflow import preparation  # noqa: PLC0415

    clock = iter((10.0, 11.0))
    monkeypatch.setattr(preparation.time, "monotonic", lambda: next(clock))
    recipe = background_request().with_changes(
        budget=parts.CandidateBudget(6, seconds=0.5), mining_target=target
    )
    pool = da.prepare(recipe, out=tmp_path / "pool")
    assert da.inspect(pool, verify=True).state == "incomplete"
    report = da.inspect(pool, view="quality").to_dict()
    assert report["counts"]["processed"] == 0
    assert report["stop_reason"] == "time_budget"


@pytest.mark.parametrize(
    "key,representatives,ranks",
    [
        ("core", [2, 2, 2, 4], [None, 1, None, 2]),
        ("sequence", [1, 2, 3, 4], [3, 1, 2, 4]),
    ],
)
def test_best_score_representative_and_rank_use_recorded_evidence(
    key: str, representatives: list[int], ranks: list[int | None]
):
    from dense_arrays.parts.candidates import Candidate  # noqa: PLC0415
    from dense_arrays.parts.retention.selection import (  # noqa: PLC0415
        select_candidates,
    )
    from dense_arrays.parts.scoring import FimoHit  # noqa: PLC0415

    def candidate(index: int, sequence: str, raw: float) -> Candidate:
        score = FimoHit(0, 2, "forward", sequence[:2], raw, 0.01, 5)
        return Candidate(
            index,
            parts.Part(
                f"c{index}",
                sequence,
                "A",
                core_start=0,
                core_end=2,
                core_orientation="forward",
                metadata={"score": score.to_dict()},
            ),
        )

    original = (
        candidate(1, "ACT", 3),
        candidate(2, "ACG", 4),
        candidate(3, "ACA", 4),
        candidate(4, "GTT", 2),
    )
    decided = select_candidates(
        original,
        parts.Uniqueness(key),
        parts.Retention(count=1, policy="top_score", rank_by="best_hit_score"),
    )
    assert [c.representative for c in decided] == representatives
    assert [c.rank for c in decided] == ranks
    assert [c.index for c in decided if c.retained] == [2]
    assert all(c.representative is None for c in original)


def test_preparation_and_verification_never_score_during_preview_or_inspection(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.workflow import preparation  # noqa: PLC0415

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("preview or inspection invoked preparation execution")

    recipe = background_request()
    pool = da.prepare(recipe, out=tmp_path / "pool")
    monkeypatch.setattr(preparation, "_batch", forbidden)
    monkeypatch.setattr(preparation, "scan_fimo", forbidden)
    assert da.plan(recipe).preview["retained_parts"] is None
    assert da.inspect(pool, verify=True).retained_parts == 1
    assert da.inspect(pool, view="quality").to_dict()["counts"]["processed"] == 6


def test_stale_pwm_input_fails_before_creating_output(tmp_path: Path):
    from .test_fimo import scorer  # noqa: PLC0415
    from .test_motifs import source  # noqa: PLC0415

    motif = source(tmp_path)
    request = parts.PreparationSpec(
        motif,
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(1),
        scoring=parts.FimoScoring(executable=scorer(tmp_path, "")),
    )
    plan = da.plan(request)
    motif.path.write_text("{}")
    target = tmp_path / "pool"
    with pytest.raises(ValueError, match="changed"):
        da.prepare(plan, out=target)
    assert not target.exists()


def test_quality_report_reuses_its_bound_source(tmp_path: Path):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    report = da.inspect(pool, view="quality")
    assert da.inspect(report, view="quality") is report


def test_builtin_score_reasons_cannot_collide_with_named_screens():
    with pytest.raises(ValueError, match="reserved"):
        background_request().with_changes(
            screening=(planning.Avoid("no_qualifying_hit", ("AA",)),)
        )


def test_execution_error_accounting_requires_execution_errors():
    from dense_arrays.artifacts.preparation.records import (  # noqa: PLC0415
        PoolAccounting,
    )

    counts = dict.fromkeys(
        (
            "processed",
            "eligibility_rejected",
            "eligible",
            "execution_error",
            "duplicate_discarded",
            "eligible_unique",
            "retained",
            "not_selected",
        ),
        0,
    )
    with pytest.raises(ValueError, match="execution"):
        PoolAccounting(counts, 1, 6, "execution_error", {})
