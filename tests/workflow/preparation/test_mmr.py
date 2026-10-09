"""MMR selection uses explicit relevance, bounded pools and saved choice evidence.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest

from dense_arrays import parts
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.motifs import Motif
from dense_arrays.parts.retention.selection import select_candidates
from dense_arrays.parts.scoring import FimoHit

FIXTURE = json.loads(
    (Path(__file__).parents[2] / "fixtures/workflow/densegen-mmr-v1.json").read_text()
)


def candidate(index: int, sequence: str, raw: float, maximum: float = 10) -> Candidate:
    hit = FimoHit(0, len(sequence), "forward", sequence, raw, 0.01, maximum)
    return Candidate(
        index,
        parts.Part(
            f"candidate_{index}",
            sequence,
            "example",
            core_start=0,
            core_end=len(sequence),
            core_orientation="forward",
            metadata={"score": hit.to_dict()},
        ),
    )


@pytest.mark.parametrize("case", FIXTURE["cases"])
def test_mmr_matches_pinned_choices_and_choice_evidence(case: dict):
    policy = parts.MMR(
        pool_size=case["pool_size"],
        relevance_weight=case["relevance_weight"],
        score_scaling=case["score_scaling"],
        minimum_fraction_of_max=case["minimum_fraction_of_max"],
    )
    motif = Motif(
        "example",
        tuple(tuple(row[b] for b in "ACGT") for row in case["probabilities"]),
        tuple(case["background"][b] for b in "ACGT"),
        ((0, 0, 0, 0),) * 2,
        "fixture",
    )
    original = tuple(
        candidate(i, r["sequence"], r["raw"], r["maximum"])
        for i, r in enumerate(case["candidates"], 1)
    )
    result = select_candidates(
        original,
        parts.Uniqueness("core"),
        parts.Retention(
            count=case["target"], policy="mmr", rank_by="best_hit_score", mmr=policy
        ),
        motif=motif,
    )
    selected = sorted((c for c in result if c.retained), key=lambda c: c.rank)
    assert [c.part.sequence for c in selected] == case["selected"]
    for chosen in selected:
        expected = case["evidence"][chosen.part.sequence]
        assert chosen.selection.utility == pytest.approx(expected["selection_utility"])
        assert (
            chosen.selection.nearest_distance == expected["nearest_selected_distance"]
        )
        assert (
            chosen.selection.nearest_similarity
            == expected["nearest_selected_similarity"]
        )
        assert Candidate.from_dict(chosen.to_dict()) == chosen
    assert (
        sum(c.selection.pool_status == "included" for c in result)
        == case["pool_size_actual"]
    )
    assert all(c.selection is None for c in original)


@pytest.mark.parametrize(
    "change",
    [
        {"pool_size": 0},
        {"pool_size": True},
        {"relevance_weight": 0},
        {"relevance_weight": float("nan")},
        {"score_scaling": "minmax"},
        {"distance": "sequence"},
        {"minimum_fraction_of_max": -1},
        {"tie_break": "random"},
        {"algorithm": "unknown"},
    ],
)
def test_mmr_policy_rejects_ambiguous_or_unbounded_settings(change: dict):
    args = {
        "pool_size": 10,
        "relevance_weight": 0.5,
        "score_scaling": "fraction_of_max_clipped",
    }
    with pytest.raises((TypeError, ValueError)):
        parts.MMR(**(args | change))


def test_mmr_plan_pool_and_cli_share_explicit_policy_and_saved_evidence(tmp_path: Path):
    from typer.testing import CliRunner  # noqa: PLC0415

    import dense_arrays as da  # noqa: PLC0415
    from dense_arrays import planning  # noqa: PLC0415
    from dense_arrays.cli import app  # noqa: PLC0415
    from dense_arrays.workflow.inputs import read_source  # noqa: PLC0415

    from .test_fimo import scorer  # noqa: PLC0415
    from .test_motifs import artifact  # noqa: PLC0415

    model = artifact()
    model["probabilities"] = [
        dict(zip("ACGT", row, strict=True)) for row in [(1, 0, 0, 0), (0, 1, 0, 0)]
    ]
    path = tmp_path / "motif.json"
    path.write_text(json.dumps(model))
    rows = "".join(
        f"motif\t\tcandidate_{i}\t1\t2\t+\t3\t0.01\t\tAC\n" for i in range(6)
    )
    tool = scorer(tmp_path, rows)
    policy = parts.MMR(
        pool_size=6, relevance_weight=0.5, score_scaling="fraction_of_max_clipped"
    )
    request = parts.PreparationSpec(
        parts.PWMArtifact(path),
        parts.Retention(count=2, policy="mmr", rank_by="best_hit_score", mmr=policy),
        sampling=parts.Sampling(planning.Length(exact=2)),
        budget=parts.CandidateBudget(candidates=6),
        scoring=parts.FimoScoring(executable=tool),
        uniqueness=parts.Uniqueness("core"),
    )
    plan = da.plan(request)
    assert plan.preview["retention"]["algorithm"] == "greedy_mmr.v1"
    assert plan.preview["retention"]["distance_work_bound"] == 24
    assert plan.to_dict()["policies"]["retention"] == "greedy_mmr.v1"
    request_path = tmp_path / "request.json"
    da.export(request, view="request", out=request_path)
    assert da.plan(read_source(request_path)).plan_id == plan.plan_id
    saved = tmp_path / "plan.json"
    plan.write(saved)
    runner = CliRunner()
    response = runner.invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "pool"), "--json"]
    )
    assert response.exit_code == 3, response.output
    pool = tmp_path / "pool"
    quality = da.inspect(pool, view="quality").to_dict()
    assert quality["retention"] == {
        "policy": "mmr",
        "pool_size": 1,
        "below_score": 0,
        "beyond_limit": 0,
    }
    assert quality["counts"]["retained"] == 1
    human = runner.invoke(app, ["inspect", str(pool), "--view", "quality"])
    assert human.exit_code == 0, human.output
    assert "MMR pool: 1" in human.stdout
    tool.unlink()
    path.unlink()
    assert da.inspect(pool, verify=True).retained_parts == 1
    import sqlite3  # noqa: PLC0415

    from dense_arrays._record_validation import (  # noqa: PLC0415
        canonical_json,
        semantic_digest,
    )
    from dense_arrays.artifacts.errors import ArtifactIntegrityError  # noqa: PLC0415

    with sqlite3.connect(pool / "pool.sqlite3") as connection:
        row = connection.execute(
            "SELECT ordinal,payload FROM candidates WHERE ordinal=1"
        ).fetchone()
        altered = json.loads(row[1])
        altered["selection"]["utility"] += 0.125
        connection.execute(
            "UPDATE candidates SET payload=?,digest=? WHERE ordinal=?",
            (canonical_json(altered), semantic_digest(altered), row[0]),
        )
    with pytest.raises(ArtifactIntegrityError, match="selection"):
        da.inspect(pool, verify=True)


def flat_motif() -> Motif:
    return Motif("example", ((0.25,) * 4,) * 2, (0.25,) * 4, ((0,) * 4,) * 2, "fixture")


def retention(
    count: int, scale: str = "fraction_of_max_clipped", **options: object
) -> parts.Retention:
    return parts.Retention(
        count=count,
        policy="mmr",
        rank_by="best_hit_score",
        mmr=parts.MMR(
            pool_size=6, relevance_weight=0.5, score_scaling=scale, **options
        ),
    )


@pytest.mark.parametrize(
    "background,second,distance",
    [
        ((0.25, 0.25, 0.25, 0.25), "AC", 1.0),
        ((0.97, 0.01, 0.01, 0.01), "CA", 0.9933858671331224),
    ],
)
def test_mmr_distance_uses_declared_background_not_variability_alone(
    background: tuple[float, ...], second: str, distance: float
):
    # The first position is invariant A; the second has maximum base entropy.
    # A common invariant base carries little information against a biased null.
    motif = Motif(
        "example",
        ((1.0, 0.0, 0.0, 0.0), (0.25, 0.25, 0.25, 0.25)),
        background,
        None,
        "test",
    )
    original = tuple(
        candidate(i, sequence, 10) for i, sequence in enumerate(("AA", "AC", "CA"), 1)
    )
    result = select_candidates(
        original, parts.Uniqueness("core"), retention(2), motif=motif
    )
    chosen = sorted((c for c in result if c.retained), key=lambda c: c.rank)
    assert [c.part.sequence for c in chosen] == ["AA", second]
    assert chosen[1].selection.nearest_distance == pytest.approx(distance)


def test_percentile_mmr_accepts_zero_maximum_without_inventing_a_score_ratio():
    original = tuple(
        candidate(i, core, 0, maximum=0) for i, core in enumerate(("AA", "AT", "TT"), 1)
    )
    result = select_candidates(
        original,
        parts.Uniqueness("core"),
        retention(2, "score_percentile"),
        motif=flat_motif(),
    )
    chosen = sorted((c for c in result if c.retained), key=lambda c: c.rank)
    assert [c.part.sequence for c in chosen] == ["AA", "TT"]
    assert all(
        c.selection.relevance == 1 and c.score.fraction_of_max is None for c in chosen
    )
    with pytest.raises(ValueError, match="positive theoretical maximum"):
        select_candidates(
            original, parts.Uniqueness("core"), retention(2), motif=flat_motif()
        )


def test_mmr_order_is_prefix_stable_when_all_pool_members_are_retained():
    original = tuple(
        candidate(i, row["sequence"], row["raw"])
        for i, row in enumerate(FIXTURE["cases"][0]["candidates"], 1)
    )
    result = select_candidates(
        original, parts.Uniqueness("core"), retention(6), motif=flat_motif()
    )
    chosen = sorted(result, key=lambda c: c.rank)
    assert [c.part.sequence for c in chosen[:2]] == ["AA", "TT"]
    assert all(c.retained and c.selection.utility is not None for c in chosen)


def test_mmr_cap_and_score_filter_do_not_change_eligibility_counts():
    from dense_arrays.artifacts.preparation.records import recount  # noqa: PLC0415

    original = tuple(
        candidate(i, row["sequence"], row["raw"])
        for i, row in enumerate(FIXTURE["cases"][0]["candidates"], 1)
    )
    policy = parts.Retention(
        count=4,
        policy="mmr",
        rank_by="best_hit_score",
        mmr=parts.MMR(
            pool_size=2,
            relevance_weight=0.5,
            score_scaling="fraction_of_max_clipped",
            minimum_fraction_of_max=0.85,
        ),
    )
    result = select_candidates(
        original, parts.Uniqueness("core"), policy, motif=flat_motif()
    )
    accounting = recount(
        result, target=4, budget=6, stop_reason="candidate_budget", mmr=True
    )
    assert accounting.state == "incomplete"
    assert accounting.counts["eligible"] == 6
    assert accounting.counts["retained"] == 2
    assert accounting.retention == {
        "policy": "mmr",
        "pool_size": 2,
        "below_score": 3,
        "beyond_limit": 1,
    }


def test_mmr_rejects_unscored_representatives_explicitly():
    with pytest.raises(ValueError, match="scored"):
        select_candidates(
            (Candidate(1, parts.Part("a", "AA")),),
            parts.Uniqueness(),
            retention(1),
            motif=flat_motif(),
        )
