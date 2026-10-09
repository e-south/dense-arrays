"""Pool reports expose saved, bounded MMR choice evidence per recipe.

Author: Eric J. South.
"""

from pathlib import Path

import pytest

import dense_arrays as da
from dense_arrays import parts, planning, reporting

from .test_fimo import scorer
from .test_motifs import source


def mmr_pool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, count: int = 3, recipes: int = 1
):
    import dense_arrays.workflow.preparation as execution  # noqa: PLC0415

    sequences = ("AC", "AA", "CC")
    monkeypatch.setattr(
        execution, "sample_sequence", lambda **kw: sequences[kw["index"] - 1]
    )
    rows = "".join(
        f"motif\t\tcandidate_{i}\t1\t2\t+\t{3 - i}\t0.00001\t\t{seq}\n"
        for i, seq in enumerate(sequences)
    )
    request = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(
            count=count,
            policy="mmr",
            rank_by="best_hit_score",
            mmr=parts.MMR(
                pool_size=3,
                relevance_weight=0.5,
                score_scaling="fraction_of_max_clipped",
            ),
        ),
        sampling=parts.Sampling(planning.Length(exact=2)),
        scoring=parts.FimoScoring(executable=scorer(tmp_path, rows)),
        budget=parts.CandidateBudget(candidates=3),
        uniqueness=parts.Uniqueness("core"),
    )
    prepared = (
        request
        if recipes == 1
        else parts.PreparationSet(
            {f"recipe_{i}": request for i in range(recipes)},
            sequence_collisions="preserve",
        )
    )
    plan = da.plan(prepared)
    pool = da.prepare(plan, out=tmp_path / "pool")
    return pool, plan


def test_live_report_reuses_verified_candidates_and_keeps_saved_choice_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.artifacts.preparation import verification  # noqa: PLC0415

    pool, plan = mmr_pool(tmp_path, monkeypatch)
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidates = tuple(row.candidate for row in rows)
    selected = sorted((c for c in candidates if c.retained), key=lambda c: c.rank)
    reads = []
    original = verification.read_candidates

    def counted(*args: object) -> object:
        reads.append(1)
        return original(*args)

    monkeypatch.setattr(verification, "read_candidates", counted)
    (tmp_path / "fimo fixture").unlink()
    (tmp_path / "motif.json").unlink()
    data = da.inspect(pool, view="quality").to_dict()
    assert data["diversity"] == [
        {
            "schema": "dense_arrays.mmr_diversity.v1",
            "recipe_id": None,
            "policy": "greedy_mmr.v1",
            "distance": "pwm_tolerant_hamming",
            "model_id": plan.resolved.source.motif.model_id,
            "scoring_id": plan.resolved.source.scoring.binding_id,
            "status": "recorded",
            "choices": [
                {"rank": c.rank, "nearest_distance": c.selection.nearest_distance}
                for c in selected
            ],
        }
    ]
    assert reads == [1]
    assert data["pool_id"] == pool.pool_id
    assert selected[0].selection.nearest_distance is None
    assert all(c.selection.nearest_distance > 0 for c in selected[1:])
    assert reporting.PoolQualitySnapshot.from_dict(data).to_dict() == data
    legacy = {k: v for k, v in data.items() if k != "diversity"}
    assert reporting.PoolQualitySnapshot.from_dict(legacy).to_dict() == legacy
    reporting.PoolQualitySnapshot.from_dict(
        legacy, read_limits=reporting.ReadLimits(identities=35)
    )
    with pytest.raises(reporting.ReadLimitError, match="identities"):
        reporting.PoolQualitySnapshot.from_dict(
            data, read_limits=reporting.ReadLimits(identities=35)
        )


@pytest.mark.parametrize("recipes", [1, 2])
def test_mmr_verification_enforces_one_shared_pair_allowance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, recipes: int
):
    pool, _ = mmr_pool(tmp_path, monkeypatch, recipes=recipes)
    with pytest.raises(reporting.ReadLimitError, match="pairs"):
        da.inspect(
            pool,
            view="quality",
            read_limits=reporting.ReadLimits(pairs=9 * recipes - 1),
        ).to_dict()
    assert (
        da.inspect(
            pool, view="quality", read_limits=reporting.ReadLimits(pairs=9 * recipes)
        ).to_dict()["counts"]["retained"]
        == 3 * recipes
    )


def test_pair_limit_is_admitted_before_selection_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    pool, _ = mmr_pool(tmp_path, monkeypatch)

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("pair allowance must be checked before selection replay")

    monkeypatch.setattr(
        "dense_arrays.artifacts.preparation.verification.select_candidates", forbidden
    )
    with pytest.raises(reporting.ReadLimitError, match="pairs"):
        da.inspect(
            pool, view="quality", read_limits=reporting.ReadLimits(pairs=8)
        ).to_dict()


@pytest.mark.parametrize("count,status", [(0, "empty"), (1, "singleton")])
def test_empty_and_first_choice_have_no_invented_distance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, count: int, status: str
):
    pool, _ = mmr_pool(tmp_path, monkeypatch, count=count)
    data = da.inspect(pool, view="quality").to_dict()
    record = data["diversity"][0]
    assert record["status"] == status
    assert record["choices"] == (
        [] if not count else [{"rank": 1, "nearest_distance": None}]
    )
    assert reporting.PoolQualitySnapshot.from_dict(data).to_dict() == data


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema", "dense_arrays.mmr_diversity.v2"),
        ("policy", "unknown"),
        ("distance", "hamming"),
        ("recipe_id", "unknown"),
        ("model_id", "bad"),
        ("scoring_id", "bad"),
        ("status", "singleton"),
        ("choices", []),
    ],
)
def test_saved_diversity_rejects_mismatched_contracts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, field: str, value: object
):
    pool, _ = mmr_pool(tmp_path, monkeypatch)
    data = da.inspect(pool, view="quality").to_dict()
    data["diversity"][0][field] = value
    with pytest.raises((ValueError, TypeError)):
        reporting.PoolQualitySnapshot.from_dict(data)


@pytest.mark.parametrize(
    "index,value",
    [(0, 0.0), (1, None), (1, True), (1, -0.1), (1, float("nan")), (1, float("inf"))],
)
def test_saved_distances_reject_invented_or_invalid_numbers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, index: int, value: object
):
    pool, _ = mmr_pool(tmp_path, monkeypatch)
    data = da.inspect(pool, view="quality").to_dict()
    data["diversity"][0]["choices"][index]["nearest_distance"] = value
    with pytest.raises((ValueError, TypeError)):
        reporting.PoolQualitySnapshot.from_dict(data)


def test_recipe_diversity_keeps_all_mmr_populations_separate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    pool, _ = mmr_pool(tmp_path, monkeypatch, recipes=2)
    data = da.inspect(pool, view="quality").to_dict()
    assert [row["recipe_id"] for row in data["diversity"]] == ["recipe_0", "recipe_1"]
    assert reporting.PoolQualitySnapshot.from_dict(data).to_dict() == data
    data["diversity"].reverse()
    with pytest.raises(ValueError, match="recipe"):
        reporting.PoolQualitySnapshot.from_dict(data)
