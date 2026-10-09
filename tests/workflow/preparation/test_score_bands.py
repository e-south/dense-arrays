"""Empirical score bands preserve ties and name their representative population.

Author: Eric J. South.
"""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import semantic_digest
from dense_arrays.artifacts.pools import decode_preparation
from dense_arrays.artifacts.preparation.records import PoolAccounting, recount
from dense_arrays.artifacts.preparation.verification import verify_decisions
from dense_arrays.cli import app
from dense_arrays.parts.retention.selection import select_candidates
from dense_arrays.reporting import ReadLimitError, ReadLimits
from dense_arrays.workflow.inputs import read_source

from .test_mmr import candidate

SCORING_ID = semantic_digest({"scope": "controlled test score model"})


def test_score_bands_keep_boundary_ties_and_exclude_duplicates():
    original = (
        candidate(1, "AAA", 9),
        candidate(2, "CCC", 9),
        candidate(3, "GGG", 5),
        candidate(4, "TTT", 1),
        candidate(5, "AAA", 9),
        replace(candidate(6, "ACG", 8), reasons=("screen",)),
    )
    policy = parts.ScoreBands((0.25, 0.5))
    result = select_candidates(
        original,
        parts.Uniqueness("sequence"),
        parts.Retention(count=2, policy="top_score", rank_by="best_hit_score"),
        score_bands=policy,
    )
    assert [c.score_band for c in result] == [1, 1, 3, 3, None, None]
    assert [c.index for c in result if c.retained] == [1, 2]
    report = recount(
        result,
        target=2,
        budget=6,
        stop_reason="candidate_budget",
        score_bands=policy,
        scoring_id=SCORING_ID,
    ).to_dict()["score_bands"]
    assert report["population"] == "eligible_unique"
    assert report["total"] == 4
    assert report["units"] == "fimo_log2_odds"
    assert [b["count"] for b in report["bands"]] == [2, 0, 2]
    assert [b["retained"] for b in report["bands"]] == [2, 0, 0]
    assert [b["cutoff"] for b in report["bands"]] == [9, 9, None]
    assert report["bands"][1]["scores"] is None
    assert report["bands"][2]["scores"] == {"min": 1, "median": 3, "max": 5}


@pytest.mark.parametrize(
    "fractions", [(), (0,), (1,), (0.5, 0.25), (0.25, 0.25), (True,), (float("nan"),)]
)
def test_score_bands_require_explicit_ordered_interior_fractions(fractions: tuple):
    with pytest.raises((ValueError, TypeError)):
        parts.ScoreBands(fractions)


def test_score_bands_use_the_full_representative_population_before_mmr_pool_cap():
    from dense_arrays.parts.motifs import Motif  # noqa: PLC0415

    policy = parts.ScoreBands((0.5,))
    motif = Motif(
        "example", ((0.25,) * 4,) * 3, (0.25,) * 4, ((0,) * 4,) * 3, "fixture"
    )
    result = select_candidates(
        tuple(
            candidate(i, seq, score)
            for i, (seq, score) in enumerate(
                (("AAA", 9), ("CCC", 8), ("GGG", 4), ("TTT", 1)), 1
            )
        ),
        parts.Uniqueness("core"),
        parts.Retention(
            count=1,
            policy="mmr",
            rank_by="best_hit_score",
            mmr=parts.MMR(2, 0.5, "score_percentile"),
        ),
        motif=motif,
        score_bands=policy,
    )
    assert [c.score_band for c in result] == [1, 1, 2, 2]
    assert result[-1].selection.pool_status == "beyond_limit"


def test_score_bands_describe_an_empty_population():
    policy = parts.ScoreBands((0.1, 0.5))
    report = recount(
        (),
        target=0,
        budget=1,
        stop_reason="time_budget",
        score_bands=policy,
        scoring_id=SCORING_ID,
    ).to_dict()["score_bands"]
    assert report["total"] == 0
    assert all(
        b["count"] == 0 and b["scores"] is None and b["cutoff"] is None
        for b in report["bands"]
    )


def test_score_bands_require_scoring_identity():
    with pytest.raises((TypeError, ValueError), match="scoring"):
        recount(
            (),
            target=0,
            budget=1,
            stop_reason="time_budget",
            score_bands=parts.ScoreBands((0.1,)),
        )


def score_request(tmp_path: Path):

    from .test_fimo import scorer  # noqa: PLC0415
    from .test_motifs import artifact  # noqa: PLC0415

    model = artifact()
    model["probabilities"] = [
        dict(zip("ACGT", r, strict=True)) for r in ((1, 0, 0, 0), (0, 1, 0, 0))
    ]
    path = tmp_path / "motif.json"
    path.write_text(json.dumps(model))
    tool = scorer(
        tmp_path,
        "".join(f"motif\t\tcandidate_{i}\t1\t2\t+\t3\t0.01\t\tAC\n" for i in range(6)),
    )
    request = parts.PreparationSpec(
        parts.PWMArtifact(path),
        parts.Retention(count=1, policy="top_score", rank_by="best_hit_score"),
        sampling=parts.Sampling(planning.Length(exact=2)),
        budget=parts.CandidateBudget(candidates=6),
        scoring=parts.FimoScoring(executable=tool),
        score_bands=parts.ScoreBands((0.1, 0.5)),
    )
    return request


def test_score_bands_share_plan_pool_reports_and_cli(tmp_path: Path):
    request = score_request(tmp_path)
    plan = da.plan(request)
    assert plan.preview["score_bands"]["population"] == "eligible_unique"
    assert plan.preview["score_bands"]["counts"] is None
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).plan_id == plan.plan_id
    pool = da.prepare(plan, out=tmp_path / "python")
    report = da.inspect(pool, view="quality").to_dict()
    assert (
        report["score_bands"]["scoring_id"] == plan.resolved.source.scoring.binding_id
    )
    assert [b["count"] for b in report["score_bands"]["bands"]] == [1, 0, 0]
    cli = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert cli.exit_code == 0, cli.output
    assert da.inspect(tmp_path / "cli", view="quality").to_dict() == report
    human = CliRunner().invoke(app, ["inspect", str(pool.path), "--view", "quality"])
    assert human.exit_code == 0, human.output
    assert "Band 1: 1 eligible; 1 retained" in human.stdout
    assert "FIMO log2-odds" in human.stdout
    da.export(da.inspect(pool, view="quality"), out=tmp_path / "quality.json")
    assert da.inspect(tmp_path / "quality.json", view="quality").to_dict() == report


def test_score_band_queries_are_recipe_scoped_and_shared_with_cli(tmp_path: Path):
    request = score_request(tmp_path)
    plan = da.plan(
        parts.PreparationSet(
            {
                "two": request.with_changes(score_bands=parts.ScoreBands((0.5,))),
                "three": request.with_changes(
                    scoring=request.scoring.__class__(
                        executable=request.scoring.executable, pseudocount=0.2
                    )
                ),
            },
            sequence_collisions="preserve",
        )
    )
    pool = da.prepare(plan, out=tmp_path / "set")
    report = da.inspect(pool, view="quality").to_dict()
    assert "score_bands" not in report
    first, second = (r["accounting"]["score_bands"] for r in report["recipes"])
    assert len(first["bands"]) == 2
    assert len(second["bands"]) == 3
    assert first["scoring_id"] != second["scoring_id"]
    selected = da.reporting.CandidateFilter(recipes=("three",), score_bands=(1,))
    with da.inspect(pool, view="candidates", select=selected).records() as rows:
        expected = [row.to_dict() for row in rows]
    assert len(expected) == 1
    assert expected[0]["candidate"]["recipe_id"] == "three"
    predicate = tmp_path / "filter.json"
    predicate.write_text(json.dumps(selected.to_dict()))
    cli = CliRunner().invoke(
        app,
        [
            "export",
            str(pool.path),
            "--view",
            "candidates",
            "--selection",
            str(predicate),
            "--all",
            "--out",
            "-",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["records"] == expected
    for selected in (
        da.reporting.CandidateFilter(score_bands=(1,)),
        da.reporting.CandidateFilter(recipes=("two",), score_bands=(3,)),
    ):
        with pytest.raises(ValueError, match="score band"):
            da.inspect(pool, view="candidates", select=selected)


def test_score_bands_do_not_change_selected_parts_and_verify_saved_membership(
    tmp_path: Path,
):

    request = score_request(tmp_path)
    plan = da.plan(request)
    pool = da.prepare(plan, out=tmp_path / "bands")
    plain = da.prepare(request.with_changes(score_bands=None), out=tmp_path / "plain")
    with da.inspect(pool, view="parts").records() as rows:
        banded = [r.part for r in rows]
    with da.inspect(plain, view="parts").records() as rows:
        assert [r.part for r in rows] == banded
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidates = tuple(r.candidate for r in rows)
    forged = (replace(candidates[0], score_band=2), *candidates[1:])
    with pytest.raises(ValueError, match="selection or stage accounting"):
        verify_decisions(forged, plan.resolved, da.inspect(pool).preparation)


def test_score_band_reports_and_plan_decoding_respect_read_caps(tmp_path: Path):

    request = score_request(tmp_path).with_changes(
        score_bands=parts.ScoreBands(tuple(i / 100 for i in range(1, 100)))
    )
    plan = da.plan(request)
    with pytest.raises(ReadLimitError):
        decode_preparation(plan.to_dict(), max_identities=10)
    pool = da.prepare(plan, out=tmp_path / "bands")
    saved = tmp_path / "quality.json"
    da.export(pool, view="quality", out=saved)
    for source, view in ((pool, "summary"), (saved, "quality")):
        with pytest.raises(ReadLimitError):
            da.inspect(source, view=view, read_limits=ReadLimits(identities=200))


def test_score_band_report_rejects_impossible_rank_boundary():

    policy = parts.ScoreBands((0.5,))
    result = select_candidates(
        (candidate(1, "AAA", 9), candidate(2, "CCC", 1)),
        parts.Uniqueness("sequence"),
        parts.Retention(count=1, policy="top_score", rank_by="best_hit_score"),
        score_bands=policy,
    )
    report = recount(
        result,
        target=1,
        budget=2,
        stop_reason="candidate_budget",
        score_bands=policy,
        scoring_id=SCORING_ID,
    ).to_dict()
    report["score_bands"]["bands"][0]["cutoff"] = 5
    with pytest.raises(ValueError, match="score band"):
        PoolAccounting.from_dict(report)


@pytest.mark.parametrize(
    "change", ["count", "retained", "empty", "fraction", "scoring", "tie", "rank"]
)
def test_score_band_report_rejects_inconsistent_summaries(change: str):
    policy = parts.ScoreBands((0.5,))
    decided = select_candidates(
        (candidate(1, "AAA", 9), candidate(2, "CCC", 1)),
        parts.Uniqueness("sequence"),
        parts.Retention(count=1, policy="top_score", rank_by="best_hit_score"),
        score_bands=policy,
    )
    report = recount(
        decided,
        target=1,
        budget=2,
        stop_reason="candidate_budget",
        score_bands=policy,
        scoring_id=SCORING_ID,
    ).to_dict()
    bands = report["score_bands"]["bands"]
    if change == "count":
        bands[0]["count"] = 2
    elif change == "retained":
        bands[1]["retained"] = 1
    elif change == "empty":
        bands[1]["scores"] = None
    elif change == "fraction":
        bands[0]["upper_fraction"] = 1
    elif change == "scoring":
        report["score_bands"]["scoring_id"] = "unknown"
    elif change == "tie":
        bands[1]["scores"] = {"min": 9, "median": 9, "max": 9}
    elif change == "rank":
        bands[0]["upper_fraction"] = 0.9
    with pytest.raises((ValueError, TypeError), match=r"score|scoring"):
        PoolAccounting.from_dict(report)


def test_score_bands_fail_for_unscored_preparation_before_opening_files():
    with pytest.raises(ValueError, match="score bands require PWM"):
        parts.PreparationSpec(
            parts.PartTable("absent.csv", "csv"), score_bands=parts.ScoreBands((0.5,))
        )
