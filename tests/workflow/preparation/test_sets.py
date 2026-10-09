"""Independent preparation recipes publish one pool with per-recipe evidence.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.workflow.inputs import read_source


def background(base: str, *, candidates: int, count: int = 1):
    return parts.PreparationSpec(
        parts.Background(
            group=base, base_probabilities=tuple(int(b == base) for b in "ACGT")
        ),
        parts.Retention(count=count, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=3)),
        budget=parts.CandidateBudget(candidates),
        seed=7,
    )


def test_set_preserves_recipe_budgets_origins_and_python_cli_pool_identity(
    tmp_path: Path,
):
    request = parts.PreparationSet(
        {
            "adenine": background("A", candidates=3),
            "cytosine": background("C", candidates=4),
        }
    )
    plan = da.plan(request)
    assert plan.preview["candidate_budget"] == 7
    assert plan.preview["requested_retention"] == 2
    assert [r["id"] for r in plan.preview["recipes"]] == ["adenine", "cytosine"]
    saved = tmp_path / "set.plan.json"
    plan.write(saved)
    assert read_source(saved).to_dict() == plan.to_dict()
    pool = da.prepare(plan, out=tmp_path / "python")
    result = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["pool_id"] == pool.pool_id
    summary = da.inspect(pool, verify=True)
    assert (summary.source_parts, summary.retained_parts, summary.state) == (
        7,
        2,
        "completed",
    )
    quality = da.inspect(pool, view="quality").to_dict()
    assert [r["accounting"]["counts"]["processed"] for r in quality["recipes"]] == [
        3,
        4,
    ]
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidates = [row.candidate for row in rows]
    assert [c.index for c in candidates] == list(range(1, 8))
    assert [(c.recipe_id, c.recipe_index) for c in candidates] == [
        *(("adenine", i) for i in range(1, 4)),
        *(("cytosine", i) for i in range(1, 5)),
    ]
    assert [c.representative for c in candidates] == [1, 1, 1, 4, 4, 4, 4]
    with da.inspect(pool, view="parts", all=True).records() as rows:
        retained = [row.part for row in rows]
    assert [p.part_id for p in retained] == [
        "adenine/candidate_1",
        "cytosine/candidate_1",
    ]
    assert [p.sequence for p in retained] == ["AAA", "CCC"]


def test_set_request_export_and_saved_quality_preserve_each_shortfall(tmp_path: Path):
    request = parts.PreparationSet(
        {
            "complete": background("A", candidates=3),
            "short": background("C", candidates=4, count=2),
        }
    )
    saved = tmp_path / "request.json"
    da.export(request, view="request", out=saved)
    assert da.plan(read_source(saved)).plan_id == da.plan(request).plan_id
    pool = da.prepare(request, out=tmp_path / "pool")
    quality = da.inspect(pool, view="quality").to_dict()
    assert quality["state"] == "incomplete"
    assert quality["requested_retention"] == 3
    assert quality["counts"]["retained"] == 2
    quality_file = tmp_path / "quality.json"
    da.export(pool, view="quality", out=quality_file)
    assert da.inspect(quality_file, view="quality").to_dict() == quality
    cli = CliRunner().invoke(app, ["plan", str(saved)])
    assert cli.exit_code == 0, cli.output
    assert "complete: at most 3 candidates; requested retention 1" in cli.stdout
    assert "short: at most 4 candidates; requested retention 2" in cli.stdout


def test_cross_recipe_sequences_require_explicit_preservation(tmp_path: Path):
    recipes = {
        "first": background("A", candidates=3),
        "second": background("A", candidates=2),
    }
    with pytest.raises(ValueError, match=r"sequence.*recipes"):
        da.prepare(parts.PreparationSet(recipes), out=tmp_path / "collision")
    with pytest.raises(ValueError, match=r"committed|manifest"):
        da.inspect(tmp_path / "collision")
    pool = da.prepare(
        parts.PreparationSet(recipes, sequence_collisions="preserve"),
        out=tmp_path / "preserved",
    )
    summary = da.inspect(pool, verify=True)
    assert summary.retained_parts == 2
    with da.inspect(pool, view="parts", all=True).records() as rows:
        values = [row.part for row in rows]
    assert values[0].sequence == values[1].sequence
    assert values[0].part_id != values[1].part_id


def test_recipe_filter_and_human_quality_keep_local_context(tmp_path: Path):
    pool = da.prepare(
        parts.PreparationSet(
            {
                "A": background("A", candidates=3),
                "C": background("C", candidates=4, count=2),
            }
        ),
        out=tmp_path / "pool",
    )
    selected = da.reporting.CandidateFilter(recipes=("C",))
    with da.inspect(
        pool, view="candidates", all=True, select=selected
    ).records() as rows:
        expected = [row.to_dict() for row in rows]
    assert len(expected) == 4
    cli = CliRunner().invoke(
        app,
        [
            "export",
            str(pool.path),
            "--view",
            "candidates",
            "--recipe-id",
            "C",
            "--all",
            "--out",
            "-",
        ],
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["records"] == expected
    quality = CliRunner().invoke(app, ["inspect", str(pool.path), "--view", "quality"])
    assert quality.exit_code == 0, quality.output
    assert "C: 1 / 2 retained; stopped: candidate_budget" in quality.stdout
    with pytest.raises(ValueError, match="recipe"):
        da.inspect(
            pool,
            view="candidates",
            select=da.reporting.CandidateFilter(recipes=("missing",)),
        )


def test_saved_set_quality_charges_nested_accounting(tmp_path: Path):
    from dense_arrays.reporting import ReadLimitError, ReadLimits  # noqa: PLC0415

    pool = da.prepare(
        parts.PreparationSet({b: background(b, candidates=2) for b in "ACGT"}),
        out=tmp_path / "pool",
    )
    saved = tmp_path / "quality.json"
    da.export(pool, view="quality", out=saved)
    with pytest.raises(ReadLimitError):
        da.inspect(saved, view="quality", read_limits=ReadLimits(identities=30))


def test_reordering_recipes_preserves_local_draws_and_part_identities(tmp_path: Path):
    first = background("A", candidates=7).with_changes(
        source=parts.Background(group="one")
    )
    second = background("C", candidates=5).with_changes(
        source=parts.Background(group="two")
    )
    observations = []
    for order in (
        {"first": first, "second": second},
        {"second": second, "first": first},
    ):
        pool = da.prepare(
            parts.PreparationSet(order, sequence_collisions="preserve"),
            out=tmp_path / str(len(observations)),
        )
        with da.inspect(
            pool, verify=True, view="candidates", all=True
        ).records() as rows:
            observations.append(
                {
                    (r.candidate.recipe_id, r.candidate.recipe_index): (
                        r.candidate.part,
                        r.candidate.rank,
                        r.candidate.retained,
                    )
                    for r in rows
                }
            )
    assert observations[0] == observations[1]


def test_independent_recipe_failure_does_not_mask_other_results(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.parts.scoring import ScoringError  # noqa: PLC0415
    from dense_arrays.workflow import preparation  # noqa: PLC0415

    from .test_fimo import scorer  # noqa: PLC0415
    from .test_motifs import source  # noqa: PLC0415

    screen = parts.PWMExclusion(
        "exclude",
        (source(tmp_path),),
        parts.FimoScoring(executable=scorer(tmp_path, "")),
    )

    def failed_scan(*_args: object, **_kwargs: object) -> None:
        reason = "timeout"
        raise ScoringError(reason, "controlled scoring failure")

    monkeypatch.setattr(preparation, "scan_fimo", failed_scan)
    pool = da.prepare(
        parts.PreparationSet(
            {
                "failed": background("A", candidates=3).with_changes(
                    screening=(screen,)
                ),
                "good": background("C", candidates=4),
            }
        ),
        out=tmp_path / "pool",
    )
    summary = da.inspect(pool, verify=True)
    assert summary.state == "incomplete"
    assert summary.preparation.recipes["failed"].stop_reason == "execution_error"
    assert summary.preparation.recipes["good"].counts["retained"] == 1


def test_set_plan_and_origins_reject_tampering(tmp_path: Path):
    from dense_arrays.artifacts.preparation.sets import local_candidate  # noqa: PLC0415
    from dense_arrays.parts.candidates import Candidate  # noqa: PLC0415
    from dense_arrays.planning import PreparationPlan  # noqa: PLC0415

    plan = da.plan(parts.PreparationSet({"A": background("A", candidates=3)}))
    changed = plan.to_dict()
    changed["preview"]["candidate_budget"] = 4
    with pytest.raises(ValueError, match="disagree"):
        PreparationPlan.from_dict(changed)
    pool = da.prepare(plan, out=tmp_path / "pool")
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        first = next(rows).candidate.to_dict()
    first["recipe_index"] = 2
    with pytest.raises(ValueError, match="origin"):
        local_candidate(Candidate.from_dict(first), "A", 0)


def test_set_read_limit_admission_precedes_model_materialization(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.planning import PreparationPlan  # noqa: PLC0415
    from dense_arrays.reporting import ReadLimitError, ReadLimits  # noqa: PLC0415

    pool = da.prepare(
        parts.PreparationSet(
            {"A": background("A", candidates=2), "C": background("C", candidates=2)}
        ),
        out=tmp_path / "pool",
    )

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("set models were materialized before identity admission")

    monkeypatch.setattr(PreparationPlan, "from_dict", forbidden)
    with pytest.raises(ReadLimitError):
        da.inspect(pool, verify=True, read_limits=ReadLimits(identities=1))


def test_two_motifs_bind_separate_models_and_survive_source_removal(tmp_path: Path):
    from .test_fimo import scorer  # noqa: PLC0415

    path = tmp_path / "motifs.meme"
    path.write_text(
        "MEME version 5\nMOTIF first\nletter-probability matrix: w= 2\n"
        ".7 .1 .1 .1\n.1 .7 .1 .1\n\nMOTIF second\n"
        "letter-probability matrix: w= 2\n.85 .05 .05 .05\n.05 .85 .05 .05\n"
    )
    tool = scorer(tmp_path, "")
    base = parts.PreparationSpec(
        parts.PWMArtifact(path, format="meme", motif_ids=("first",)),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=6), strategy="consensus"),
        budget=parts.CandidateBudget(3),
        scoring=parts.FimoScoring(executable=tool),
    )
    request = parts.PreparationSet(
        {
            "first": base,
            "second": base.with_changes(
                source=parts.PWMArtifact(path, format="meme", motif_ids=("second",))
            ),
        }
    )
    plan = da.plan(request)
    assert (
        plan.resolved.recipes["first"].source.motif.model_id
        != plan.resolved.recipes["second"].source.motif.model_id
    )
    pool = da.prepare(plan, out=tmp_path / "pool")
    saved = tmp_path / "set.plan.json"
    plan.write(saved)
    path.unlink()
    tool.unlink()
    assert read_source(saved).plan_id == plan.plan_id
    assert da.inspect(pool, verify=True).source_parts == 6


def test_set_preflights_all_inputs_and_rolls_back_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.parts.scoring import ScoringError  # noqa: PLC0415
    from dense_arrays.workflow import preparation  # noqa: PLC0415

    from .test_motifs import source  # noqa: PLC0415

    missing = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=6)),
        budget=parts.CandidateBudget(2),
        scoring=parts.FimoScoring(executable=tmp_path / "absent"),
    )
    with pytest.raises(ScoringError):
        da.prepare(
            parts.PreparationSet(
                {"first": background("A", candidates=3), "missing": missing}
            ),
            out=tmp_path / "preflight",
        )
    assert not (tmp_path / "preflight").exists()
    original = preparation.publish

    def failed_publish(*args: object, **kwargs: object) -> None:
        original(*args, **kwargs)
        msg = "controlled publication failure"
        raise RuntimeError(msg)

    monkeypatch.setattr(preparation, "publish", failed_publish)
    with pytest.raises(RuntimeError, match="controlled"):
        da.prepare(
            parts.PreparationSet(
                {"A": background("A", candidates=2), "C": background("C", candidates=2)}
            ),
            out=tmp_path / "rollback",
        )
    with pytest.raises(ValueError, match=r"committed|manifest"):
        da.inspect(tmp_path / "rollback")


@pytest.mark.parametrize(
    "value", [{}, {"": background("A", candidates=2)}, {"A": "not a request"}]
)
def test_invalid_set_requests_fail_before_execution(value: dict):
    with pytest.raises((TypeError, ValueError)):
        parts.PreparationSet(value)
