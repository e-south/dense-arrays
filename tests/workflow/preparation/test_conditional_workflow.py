"""Conditional preparation is one native operation with portable evidence.

Author: Eric J. South.
"""

import json
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.cli import app
from dense_arrays.parts.background import compiler
from dense_arrays.workflow import preparation
from dense_arrays.workflow.inputs import read_source

from .test_fimo import scorer
from .test_motifs import source


def recipe(**changes: object):
    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=4, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=30), strategy="conditional"),
        budget=parts.CandidateBudget(8, batch_size=2),
        screening=(
            planning.GC("gc", "sequence", 29 / 30, 1),
            planning.Avoid("triples", ("CCC", "GGG")),
        ),
        seed=7,
    )
    return replace(request, **changes)


def test_conditional_preparation_python_cli_and_portable_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    request = recipe()
    resolved = da.plan(request)
    assert resolved.preview["proposal"]["strategy"] == "conditional"
    assert resolved.preview["proposal"]["constraint_ids"] == ("gc", "triples")
    assert resolved.preview["proposal"]["construction_limits"]["seconds"] == 30
    saved = tmp_path / "plan.json"
    resolved.write(saved)
    assert read_source(saved).to_dict() == resolved.to_dict()
    pool = da.prepare(resolved, out=tmp_path / "python")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["counts"]["processed"] == 8
    assert report["counts"]["eligibility_rejected"] == 0
    assert report["counts"]["retained"] == 4
    assert report["construction"]["status"] == "feasible"
    assert report["construction"]["reason"] == "completed"
    runner = CliRunner()
    made = runner.invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert made.exit_code == 0, made.output
    assert json.loads(made.stdout)["pool_id"] == pool.pool_id
    quality = runner.invoke(
        app, ["inspect", str(tmp_path / "cli"), "--view", "quality", "--json"]
    )
    assert json.loads(quality.stdout) == report
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        for row in rows:
            seq = row.candidate.part.sequence
            assert seq.count("C") + seq.count("G") >= 29
            assert "CCC" not in seq
            assert "GGG" not in seq
            assert (
                row.candidate.part.metadata["proposal"]["model_id"]
                == report["construction"]["model_id"]
            )

    def forbidden_rebuild(**_kwargs: object) -> None:
        pytest.fail("plan/inspection must not rebuild the sampling distribution")

    monkeypatch.setattr(compiler, "compile_background", forbidden_rebuild)
    assert da.plan(request).plan_id == resolved.plan_id
    assert da.inspect(pool, verify=True).state == "completed"
    snapshot = tmp_path / "quality.json"
    da.export(pool, view="quality", out=snapshot)
    assert da.inspect(snapshot, view="quality").to_dict() == report


@pytest.mark.parametrize(
    "limits,screening,status,reason",
    [
        (None, (planning.Avoid("none", tuple("ACGT")),), "infeasible", "completed"),
        ({"states": 1}, (), "limited", "states"),
        ({"mass_bits": 1}, (), "limited", "mass_bits"),
    ],
)
def test_construction_termination_does_not_fabricate_rejected_candidates(
    tmp_path: Path, limits: dict | None, screening: tuple, status: str, reason: str
):
    request = recipe(
        sampling=parts.Sampling(
            planning.Length(exact=4),
            strategy="conditional",
            limits=None if limits is None else parts.ConditionalLimits(**limits),
        ),
        screening=screening,
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["counts"]["processed"] == 0
    assert report["counts"]["eligibility_rejected"] == 0
    assert report["rejections"] == {}
    assert report["construction"]["status"] == status
    assert report["construction"]["reason"] == reason
    assert report["stop_reason"] == f"construction_{status}"
    assert da.inspect(pool, verify=True).state == "incomplete"
    saved = tmp_path / "request.json"
    da.plan(request).write(saved)
    result = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 3, result.output


def test_stochastic_background_plan_omits_conditional_evidence(
    tmp_path: Path,
):
    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=2, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(8),
        seed=7,
    )
    assert "limits" not in da.plan(request).to_dict()["request"]["sampling"]
    pool = da.prepare(request, out=tmp_path / "pool")
    assert "construction" not in da.inspect(pool, view="quality").to_dict()


def test_ranged_preview_discloses_length_conditioning_without_running_counting():
    request = recipe(
        sampling=parts.Sampling(parts.LengthRange(1, 2), strategy="conditional"),
        screening=(planning.Avoid("word", ("AA",), strands="forward"),),
    )
    assert (
        da.plan(request).preview["sampled_length"]["distribution"]
        == "uniform_prior_conditioned_on_constraints"
    )


def test_conditional_has_one_distribution_owner_and_rejects_motif_use(tmp_path: Path):

    with pytest.raises(ValueError, match="Background"):
        recipe(source=source(tmp_path), scoring=parts.FimoScoring())
    with pytest.raises(ValueError, match="conditional"):
        parts.Sampling(planning.Length(exact=4), limits=parts.ConditionalLimits())


def test_human_preview_explains_conditioning_and_resource_limits(tmp_path: Path):
    request = recipe(
        sampling=parts.Sampling(parts.LengthRange(1, 2), strategy="conditional"),
        screening=(planning.Avoid("word", ("AA",), strands="forward"),),
    )
    saved = tmp_path / "request.json"
    da.plan(request).write(saved)
    result = CliRunner().invoke(app, ["plan", str(saved)])
    assert result.exit_code == 0, result.output
    assert "conditioned on constraints" in result.stdout
    assert "Construction limits" in result.stdout
    assert "uniform integer draws" not in result.stdout


def test_human_quality_explains_unknown_construction(tmp_path: Path):
    request = recipe(
        sampling=parts.Sampling(
            planning.Length(exact=30),
            strategy="conditional",
            limits=parts.ConditionalLimits(states=1),
        )
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    result = CliRunner().invoke(app, ["inspect", str(pool.path), "--view", "quality"])
    assert result.exit_code == 0, result.output
    assert "feasibility unknown" in result.stdout
    assert "states" in result.stdout


def test_counting_is_reused_across_batches_and_draws_do_not_depend_on_batch_size(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    calls = []
    original = preparation.compile_background

    def tracked(**kwargs: object) -> object:
        calls.append(1)
        return original(**kwargs)

    monkeypatch.setattr(preparation, "compile_background", tracked)
    first = da.prepare(recipe(), out=tmp_path / "one")
    second = da.prepare(
        recipe(budget=parts.CandidateBudget(8, batch_size=4)), out=tmp_path / "two"
    )
    assert len(calls) == 2
    with da.inspect(first, view="candidates", all=True).records() as rows:
        expected = [r.candidate for r in rows]
    with da.inspect(second, view="candidates", all=True).records() as rows:
        assert [r.candidate for r in rows] == expected


def test_independent_fimo_exclusions_can_reject_conditionally_valid_parts(
    tmp_path: Path,
):
    motif = source(tmp_path)
    tool = scorer(tmp_path, "motif\t\tcandidate_0\t1\t2\t+\t3\t0.0625\t\tAA\n")
    request = recipe(
        source=parts.Background((1, 0, 0, 0)),
        sampling=parts.Sampling(planning.Length(exact=4), strategy="conditional"),
        budget=parts.CandidateBudget(2),
        screening=(
            planning.GC("gc", "sequence", 0, 0),
            parts.PWMExclusion(
                "motif_hit", (motif,), parts.FimoScoring(executable=tool)
            ),
        ),
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["construction"]["status"] == "feasible"
    assert report["rejections"] == {"motif_hit": 1}
    assert report["counts"]["retained"] == 1
    motif.path.unlink()
    tool.unlink()
    assert da.inspect(pool, verify=True).state == "incomplete"


@pytest.mark.parametrize("damage", ["model", "sequence"])
def test_verification_checks_conditional_proposals_before_pool_digest(
    tmp_path: Path, damage: str
):
    pool = da.prepare(recipe(), out=tmp_path / "pool")
    with sqlite3.connect(pool.path / "pool.sqlite3") as connection:
        value = json.loads(
            connection.execute(
                "SELECT payload FROM candidates WHERE ordinal=1"
            ).fetchone()[0]
        )
        if damage == "model":
            value["part"]["metadata"]["proposal"]["model_id"] = "0" * 64
        else:
            value["part"]["sequence"] = "A" * 30
        connection.execute(
            "UPDATE candidates SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(value), semantic_digest(value)),
        )
    with pytest.raises(ValueError, match="conditional proposal"):
        da.inspect(pool, verify=True)


def test_zero_mining_target_requires_no_counting_work(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def forbidden(**_kwargs: object) -> None:
        pytest.fail("zero target must not compile")

    monkeypatch.setattr(preparation, "compile_background", forbidden)
    request = recipe(
        retain=parts.Retention(count=0, policy="first_eligible"),
        mining_target=parts.MiningTarget(max_retained_fraction=1),
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    assert da.inspect(pool, verify=True).state == "completed"
    assert "construction" not in da.inspect(pool, view="quality").to_dict()
