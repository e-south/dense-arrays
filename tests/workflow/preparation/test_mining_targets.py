"""Preparation stops on explicit eligible supply without confusing retention success.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.artifacts.preparation.records import PoolAccounting, recount
from dense_arrays.artifacts.preparation.verification import verify_decisions
from dense_arrays.cli import app
from dense_arrays.workflow.inputs import read_source


def test_unique_mining_target_stops_at_first_complete_batch(tmp_path: Path):
    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=2, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=12)),
        budget=parts.CandidateBudget(100, batch_size=3),
        mining_target=parts.MiningTarget(eligible_unique=4),
        seed=7,
    )
    plan = da.plan(request)
    assert plan.preview["mining_target"]["eligible_unique"] == 4
    pool = da.prepare(plan, out=tmp_path / "pool")
    summary = da.inspect(pool, verify=True)
    assert summary.source_parts == 6
    assert summary.retained_parts == 2
    accounting = summary.preparation
    assert accounting.stop_reason == "mining_target"
    assert accounting.mining_target["met"] is True
    assert accounting.state == "completed"


@pytest.mark.parametrize(
    "count,fraction,expected",
    [(200, 0.001, 200000), (10, 0.1, 100), (3, 0.3, 10), (2, 0.3, 7), (0, 0.1, 0)],
)
def test_fraction_targets_resolve_to_explicit_eligible_supply(
    count: int, fraction: float, expected: int
):
    target = parts.MiningTarget(max_retained_fraction=fraction)
    assert target.resolve(count) == {
        "eligible_unique": expected,
        "minimum_candidates": 0,
    }


def test_full_retention_with_unmet_mining_target_remains_incomplete(
    tmp_path: Path,
):
    request = parts.PreparationSpec(
        parts.Background(base_probabilities=(1, 0, 0, 0)),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(8, batch_size=2),
        mining_target=parts.MiningTarget(max_retained_fraction=0.25),
    )
    plan = da.plan(request)
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).plan_id == plan.plan_id
    pool = da.prepare(plan, out=tmp_path / "python")
    summary = da.inspect(pool, verify=True)
    assert summary.retained_parts == summary.preparation.requested_retention == 1
    assert summary.state == "incomplete"
    assert summary.preparation.stop_reason == "candidate_budget"
    assert summary.preparation.mining_target == {
        "eligible_unique": 4,
        "minimum_candidates": 0,
        "met": False,
    }
    cli = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert cli.exit_code == 3, cli.output
    assert json.loads(cli.stdout)["pool_id"] == pool.pool_id
    quality = da.inspect(pool, view="quality").to_dict()
    out = tmp_path / "quality.json"
    da.export(pool, view="quality", out=out)
    assert da.inspect(out, view="quality").to_dict() == quality
    human = CliRunner().invoke(app, ["inspect", str(out), "--view", "quality"])
    assert human.exit_code == 0, human.output
    assert "Mining target: 1 / 4 eligible unique; unmet" in human.stdout


def test_generation_rejects_incomplete_supply_before_creating_output(tmp_path: Path):
    pool = da.prepare(
        parts.PreparationSpec(
            parts.Background((1, 0, 0, 0)),
            parts.Retention(count=1, policy="first_eligible"),
            sampling=parts.Sampling(planning.Length(exact=4)),
            budget=parts.CandidateBudget(1),
            mining_target=parts.MiningTarget(eligible_unique=2),
        ),
        out=tmp_path / "pool",
    )
    assert da.inspect(pool, verify=True).retained_parts == 1
    request = planning.DesignSpec(parts.PoolSource(pool), planning.Length(maximum=4))
    with pytest.raises(ValueError, match="pool is incomplete"):
        da.plan(request)
    with pytest.raises(ValueError, match="pool is incomplete"):
        da.run(request, out=tmp_path / "python-run")
    assert not (tmp_path / "python-run").exists()
    recipe = tmp_path / "design.yaml"
    recipe.write_text(
        "schema: dense_arrays.design.v1\nparts:\n  pool: pool\nlength:\n  maximum: 4\n"
    )
    result = CliRunner().invoke(
        app, ["run", str(recipe), "--out", str(tmp_path / "cli-run"), "--json"]
    )
    assert result.exit_code == 2, result.output
    assert "incomplete" in result.stderr
    assert "--view quality" in result.stderr
    assert not (tmp_path / "cli-run").exists()


def test_minimum_candidate_floor_and_zero_retention_have_explicit_stops(tmp_path: Path):
    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=12)),
        budget=parts.CandidateBudget(30, batch_size=3),
        mining_target=parts.MiningTarget(eligible_unique=1, minimum_candidates=7),
    )
    pool = da.prepare(request, out=tmp_path / "floor")
    assert da.inspect(pool, verify=True).source_parts == 9
    empty = request.with_changes(
        retain=parts.Retention(count=0, policy="first_eligible"),
        mining_target=parts.MiningTarget(max_retained_fraction=0.1),
    )
    result = da.inspect(da.prepare(empty, out=tmp_path / "empty"), verify=True)
    assert result.source_parts == result.retained_parts == 0
    assert result.state == "completed"
    assert result.preparation.stop_reason == "mining_target"


def test_rejected_sequences_do_not_advance_mining_target(tmp_path: Path):
    request = parts.PreparationSpec(
        parts.Background(base_probabilities=(1, 0, 0, 0)),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(6, batch_size=2),
        screening=(planning.Avoid("no_aa", ("AA",)),),
        mining_target=parts.MiningTarget(eligible_unique=1),
    )
    summary = da.inspect(da.prepare(request, out=tmp_path / "pool"), verify=True)
    assert summary.state == "incomplete"
    assert summary.preparation.counts["eligibility_rejected"] == 6
    assert summary.preparation.counts["eligible_unique"] == 0
    assert summary.preparation.stop_reason == "candidate_budget"
    assert summary.preparation.mining_target["met"] is False


@pytest.mark.parametrize(
    "values",
    [
        {},
        {"eligible_unique": True},
        {"eligible_unique": 0},
        {"eligible_unique": 2.5},
        {"max_retained_fraction": 0},
        {"max_retained_fraction": 1.1},
        {"max_retained_fraction": float("nan")},
        {"max_retained_fraction": True},
        {"eligible_unique": 2, "max_retained_fraction": 0.1},
        {"eligible_unique": 2, "minimum_candidates": -1},
    ],
)
def test_mining_target_declarations_are_strict(values: dict):
    with pytest.raises((TypeError, ValueError)):
        parts.MiningTarget(**values)


def test_verification_rejects_candidates_after_first_attained_batch(tmp_path: Path):

    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=2, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=12)),
        budget=parts.CandidateBudget(6, batch_size=2),
        seed=7,
    )
    pool = da.prepare(request, out=tmp_path / "fixed-budget")
    with da.inspect(pool, view="candidates", all=True).records() as records:
        candidates = tuple(r.candidate for r in records)
    target = da.plan(
        request.with_changes(mining_target=parts.MiningTarget(eligible_unique=2))
    )
    accounting = recount(
        candidates,
        target=2,
        budget=6,
        stop_reason="mining_target",
        mining_target=target.resolved.mining_target,
    )
    with pytest.raises(ValueError, match=r"first.*batch"):
        verify_decisions(candidates, target.resolved, accounting)


def test_core_target_counts_unique_cores_not_unique_flanks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.workflow import preparation  # noqa: PLC0415

    from .test_fimo import scorer  # noqa: PLC0415
    from .test_motifs import artifact  # noqa: PLC0415

    motif = tmp_path / "motif.json"
    motif.write_text(json.dumps(artifact()))
    tool = scorer(
        tmp_path,
        "".join(
            f"motif\t\tcandidate_{i}\t2\t3\t+\t3\t0.0625\t\tAC\n" for i in range(2)
        ),
    )
    request = parts.PreparationSpec(
        parts.PWMArtifact(motif),
        parts.Retention(count=1, policy="top_score", rank_by="best_hit_score"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(6, batch_size=2),
        scoring=parts.FimoScoring(executable=tool, hit_pvalue_max=0.1),
        uniqueness=parts.Uniqueness(key="core"),
        mining_target=parts.MiningTarget(eligible_unique=2),
    )
    monkeypatch.setattr(
        preparation,
        "sample_sequence",
        lambda **kw: (
            "ACGT"[(kw["index"] - 1) // 4] + "AC" + "ACGT"[(kw["index"] - 1) % 4]
        ),
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    summary = da.inspect(pool, verify=True)
    assert summary.preparation.counts["eligible"] == 6
    assert summary.preparation.counts["eligible_unique"] == 1
    assert summary.preparation.counts["duplicate_discarded"] == 5
    assert summary.retained_parts == 1
    assert summary.state == "incomplete"
    motif.unlink()
    tool.unlink()
    assert da.inspect(pool, verify=True).preparation.mining_target["met"] is False


def test_sets_keep_independent_mining_targets_and_quality_snapshots(tmp_path: Path):
    recipe = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=2, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=12)),
        budget=parts.CandidateBudget(10, batch_size=2),
        mining_target=parts.MiningTarget(eligible_unique=4),
        seed=8,
    )
    request = parts.PreparationSet(
        {
            "met": recipe,
            "unmet": recipe.with_changes(
                source=parts.Background(base_probabilities=(1, 0, 0, 0)),
                retain=parts.Retention(count=1, policy="first_eligible"),
            ),
        }
    )
    pool = da.prepare(request, out=tmp_path / "set")
    summary = da.inspect(pool, verify=True)
    assert summary.state == "incomplete"
    assert summary.preparation.recipes["met"].counts["processed"] == 4
    assert summary.preparation.recipes["unmet"].counts["processed"] == 10
    assert summary.preparation.recipes["met"].mining_target["met"] is True
    saved = tmp_path / "set-quality.json"
    da.export(pool, view="quality", out=saved)
    assert (
        da.inspect(saved, view="quality").to_dict()
        == da.inspect(pool, view="quality").to_dict()
    )


def test_target_outcome_cannot_disagree_with_counts(tmp_path: Path):

    recipe = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=2, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=12)),
        budget=parts.CandidateBudget(10, batch_size=2),
        mining_target=parts.MiningTarget(eligible_unique=3),
    )
    summary = da.inspect(da.prepare(recipe, out=tmp_path / "pool"), verify=True)
    value = summary.preparation.to_dict()
    value["mining_target"]["met"] = False
    with pytest.raises(ValueError, match="disagrees"):
        PoolAccounting.from_dict(value)


def test_fraction_targets_match_pinned_population_arithmetic():
    fixture = json.loads(
        (
            Path(__file__).parents[2]
            / "fixtures/workflow/densegen-tier-targets-v1.json"
        ).read_text()
    )
    for case in fixture["cases"]:
        required = parts.MiningTarget(
            max_retained_fraction=case["max_retained_fraction"]
        ).resolve(case["retained_count"])["eligible_unique"]
        assert required == case["required_unique"]
        assert (case["eligible_unique"] >= required) is case["met"]
