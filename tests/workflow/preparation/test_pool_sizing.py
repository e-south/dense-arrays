"""Target-relative MMR pools keep admission, effort bounds and selection distinct.

Author: Eric J. South.
"""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts
from dense_arrays.artifacts.preparation.records import PoolAccounting, recount
from dense_arrays.cli import app
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.motifs import Motif
from dense_arrays.parts.retention.selection import select_candidates
from dense_arrays.workflow.inputs import read_source

from .test_mmr import candidate, flat_motif
from .test_score_bands import score_request


@pytest.mark.parametrize(
    "factor,maximum,count,expected",
    [
        (10, 100, 3, 30),
        (10, 20, 3, 20),
        (1.1, 100, 10, 11),
        (2.5, 100, 3, 8),
        (10, 100, 0, 0),
    ],
)
def test_pool_size_resolves_decimal_multiplier_under_explicit_cap(
    factor: float, maximum: int, count: int, expected: int
):
    size = parts.PoolSize(per_retained=factor, maximum=maximum)
    assert size.resolve(count) == expected
    policy = parts.MMR(size, 0.5, "score_percentile")
    assert parts.MMR.from_dict(policy.to_dict()) == policy
    assert policy.pool_limit(count) == expected


@pytest.mark.parametrize(
    "factor,maximum", [(0.5, 10), (True, 10), (float("inf"), 10), (2, True), (2, 0)]
)
def test_pool_size_requires_finite_multiplier_and_positive_integer_cap(
    factor: object, maximum: object
):
    with pytest.raises((ValueError, TypeError)):
        parts.PoolSize(per_retained=factor, maximum=maximum)


def test_relative_pool_keeps_cutoff_hard_and_matches_explicit_admission():
    original = tuple(
        candidate(i, seq, raw)
        for i, (seq, raw) in enumerate(
            (("AA", 9), ("AC", 8), ("AG", 7), ("AT", 6), ("CC", 2), ("GG", 1)), 1
        )
    )
    options = {
        "relevance_weight": 0.5,
        "score_scaling": "fraction_of_max_clipped",
        "minimum_fraction_of_max": 0.5,
    }

    def select(size: int | parts.PoolSize) -> tuple[Candidate, ...]:
        return select_candidates(
            original,
            parts.Uniqueness("core"),
            parts.Retention(
                count=2,
                policy="mmr",
                rank_by="best_hit_score",
                mmr=parts.MMR(pool_size=size, **options),
            ),
            motif=flat_motif(),
        )

    relative = select(parts.PoolSize(per_retained=10, maximum=5))
    assert relative == select(5)
    assert [c.selection.pool_status for c in relative] == ["included"] * 4 + [
        "below_score"
    ] * 2
    assert sum(c.retained for c in relative) == 2


def test_relative_pool_preview_reports_and_cli_share_resolved_limits(tmp_path: Path):

    request = score_request(tmp_path).with_changes(
        retain=parts.Retention(
            count=1,
            policy="mmr",
            rank_by="best_hit_score",
            mmr=parts.MMR(parts.PoolSize(10, 5), 0.5, "score_percentile"),
        )
    )
    plan = da.plan(request)
    assert plan.preview["retention"]["pool_limit"] == 5
    assert plan.preview["retention"]["distance_work_bound"] == 10
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).plan_id == plan.plan_id
    preview = CliRunner().invoke(app, ["plan", str(saved)])
    assert preview.exit_code == 0, preview.output
    assert "MMR choice-pool limit: 5" in preview.stdout
    pool = da.prepare(plan, out=tmp_path / "python")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["retention"]["sizing"] == {
        "policy": "retained_multiplier.v1",
        "per_retained": 10.0,
        "maximum": 5,
        "requested": 10,
        "limit": 5,
        "available": 1,
        "capped": True,
        "shortfall": 4,
        "has_choice": False,
    }
    assert report["counts"]["retained"] == 1
    assert report["state"] == "completed"
    result = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["pool_id"] == pool.pool_id
    da.export(pool, view="quality", out=tmp_path / "quality.json")
    assert da.inspect(tmp_path / "quality.json", view="quality").to_dict() == report
    human = CliRunner().invoke(app, ["inspect", str(pool.path), "--view", "quality"])
    assert human.exit_code == 0, human.output
    assert "MMR choice pool: 1 / 5" in human.stdout
    assert "has alternatives: no" in human.stdout
    request.source.path.unlink()
    request.scoring.executable.unlink()
    assert da.inspect(pool, verify=True).retained_parts == 1


def test_relative_pool_count_change_and_zero_target_resolve_without_scoring(
    tmp_path: Path,
):

    spec = score_request(tmp_path)
    retention = parts.Retention(
        count=2,
        policy="mmr",
        rank_by="best_hit_score",
        mmr=parts.MMR(parts.PoolSize(3, 10), 0.5, "score_percentile"),
    )
    assert (
        da.plan(spec.with_changes(retain=retention)).preview["retention"]["pool_limit"]
        == 6
    )
    empty = spec.with_changes(retain=replace(retention, count=0))
    assert da.plan(empty).preview["retention"]["distance_work_bound"] == 0
    pool = da.prepare(empty, out=tmp_path / "zero")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["retention"]["pool_size"] == 0
    assert report["retention"]["sizing"]["limit"] == 0
    assert report["retention"]["sizing"]["has_choice"] is False
    assert report["state"] == "completed"


FIXTURE = json.loads(
    (
        Path(__file__).parents[2] / "fixtures/workflow/densegen-mmr-pools-v1.json"
    ).read_text()
)


@pytest.mark.parametrize("case", FIXTURE["cases"])
def test_direct_admission_matches_pinned_tier_widening_pool(case: dict):
    sizing = parts.PoolSize(per_retained=10, maximum=case["maximum"])
    policy = parts.Retention(
        count=case["target"],
        policy="mmr",
        rank_by="best_hit_score",
        mmr=parts.MMR(
            sizing,
            0.5,
            "fraction_of_max_clipped",
            minimum_fraction_of_max=case["cutoff"],
        ),
    )
    model = Motif(
        "example", ((0.25,) * 4,) * 3, (0.25,) * 4, ((0,) * 4,) * 3, "fixture"
    )
    result = select_candidates(
        tuple(
            candidate(i, c["sequence"], c["raw"], c["maximum"])
            for i, c in enumerate(case["candidates"], 1)
        ),
        parts.Uniqueness("core"),
        policy,
        motif=model,
    )
    admitted = [
        c.part.sequence for c in result if c.selection.pool_status == "included"
    ]
    assert admitted == case["pool_sequences"]
    report = recount(
        result,
        target=case["target"],
        budget=len(result),
        stop_reason="candidate_budget",
        mmr=True,
        pool_sizing=sizing,
    ).to_dict()["retention"]
    assert report["sizing"]["has_choice"] is case["has_choice"]
    if case["has_choice"]:
        assert [
            c.part.sequence
            for c in sorted((c for c in result if c.retained), key=lambda c: c.rank)
        ] == case["selected"]


@pytest.mark.parametrize(
    "field,value",
    [
        ("limit", 99),
        ("requested", 99),
        ("available", 99),
        ("shortfall", 99),
        ("capped", False),
        ("has_choice", True),
        ("has_choice", 1),
        ("limit", True),
    ],
)
def test_saved_pool_sizing_rejects_inconsistent_or_untyped_evidence(
    field: str, value: object
):
    size = parts.PoolSize(10, 5)
    policy = parts.Retention(
        count=1,
        policy="mmr",
        rank_by="best_hit_score",
        mmr=parts.MMR(size, 0.5, "score_percentile"),
    )
    result = select_candidates(
        (candidate(1, "AA", 9),), parts.Uniqueness("core"), policy, motif=flat_motif()
    )
    report = recount(
        result,
        target=1,
        budget=1,
        stop_reason="candidate_budget",
        mmr=True,
        pool_sizing=size,
    ).to_dict()
    report["retention"]["sizing"][field] = value
    with pytest.raises((ValueError, TypeError), match="sizing"):
        PoolAccounting.from_dict(report)


def test_preparation_sets_keep_relative_pool_sizing_local(tmp_path: Path):
    spec = score_request(tmp_path)

    def recipe(factor: float, maximum: int) -> parts.PreparationSpec:
        return spec.with_changes(
            retain=parts.Retention(
                count=1,
                policy="mmr",
                rank_by="best_hit_score",
                mmr=parts.MMR(parts.PoolSize(factor, maximum), 0.5, "score_percentile"),
            )
        )

    plan = da.plan(
        parts.PreparationSet(
            {"small": recipe(2, 3), "large": recipe(10, 5)},
            sequence_collisions="preserve",
        )
    )
    assert [r["retention"]["pool_limit"] for r in plan.preview["recipes"]] == [2, 5]
    pool = da.prepare(plan, out=tmp_path / "set")
    report = da.inspect(pool, view="quality").to_dict()
    assert [
        r["accounting"]["retention"]["sizing"]["limit"] for r in report["recipes"]
    ] == [2, 5]
    saved = tmp_path / "quality.json"
    da.export(pool, view="quality", out=saved)
    assert da.inspect(saved, view="quality").to_dict() == report
    human = CliRunner().invoke(app, ["inspect", str(saved), "--view", "quality"])
    assert human.exit_code == 0, human.output
    assert "MMR choice pool: 1 / 2" in human.stdout
    assert "MMR choice pool: 1 / 5" in human.stdout
