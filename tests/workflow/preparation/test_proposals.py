"""PWM proposals keep intended sequence construction separate from scored hits.

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
from dense_arrays.artifacts.preparation.records import recount
from dense_arrays.artifacts.preparation.verification import verify_decisions
from dense_arrays.cli import app
from dense_arrays.parts import mining
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.motifs import Motif
from dense_arrays.workflow import preparation
from dense_arrays.workflow.inputs import read_source

from .test_fimo import scorer
from .test_motifs import source


def motif() -> Motif:
    return Motif(
        "example",
        ((0.7, 0.1, 0.1, 0.1), (0.1, 0.7, 0.1, 0.1)),
        (0.25,) * 4,
        ((2, -1, -1, -1), (-1, 3, -1, -1)),
        "synthetic-fixture",
    )


def test_consensus_and_background_have_explicit_proposal_geometry():
    model = motif()
    for index in range(1, 9):
        consensus = mining.propose_sequence(
            length=8,
            seed=7,
            index=index,
            probabilities=(0, 0, 0, 1),
            motif=model,
            strategy="consensus",
        )
        assert consensus.sequence == "T" * consensus.start + "AC" + "T" * (
            8 - consensus.end
        )
        assert consensus.end - consensus.start == 2
        background = mining.propose_sequence(
            length=8,
            seed=7,
            index=index,
            probabilities=(0, 0, 0, 1),
            motif=model,
            strategy="background",
        )
        assert background.sequence == "TTTTTTTT"
        assert background.start is None
        assert background.end is None
    ties = Motif(
        "ties",
        ((0.5, 0.5, 0, 0), (0, 0, 0.5, 0.5)),
        (0.25,) * 4,
        ((1, 1, -1, -1), (-1, -1, 1, 1)),
        "test",
    )
    assert (
        mining.propose_sequence(
            length=2,
            seed=7,
            index=1,
            probabilities=(0.25,) * 4,
            motif=ties,
            strategy="consensus",
        ).sequence
        == "AG"
    )


def test_stochastic_streams_keep_preexisting_sequence_prefixes():
    assert [
        mining.sample_sequence(length=8, seed=7, index=i, probabilities=(0.25,) * 4)
        for i in range(1, 5)
    ] == ["TTCTAGCT", "CCATCCTC", "CGGGCCCT", "CTTCGCGC"]
    assert [
        mining.sample_sequence(
            length=8, seed=7, index=i, probabilities=(0.25,) * 4, motif=motif()
        )
        for i in range(1, 5)
    ] == ["TTTACTGG", "TGTATACG", "GCTCTCTA", "ACAGAAAT"]


@pytest.mark.parametrize("error", [None, "scoring failed"])
def test_background_verification_rejects_zero_probability_bases(error: str | None):
    request = parts.PreparationSpec(
        parts.Background((1, 0, 0, 0)),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=1)),
        budget=parts.CandidateBudget(1),
    )
    candidate = Candidate(
        1,
        parts.Part("candidate_1", "C", group="background", source="background"),
        representative=1 if error is None else None,
        rank=1 if error is None else None,
        retained=error is None,
        error=error,
    )
    recorded = recount(
        (candidate,),
        target=1,
        budget=1,
        stop_reason="candidate_budget" if error is None else "execution_error",
    )
    # An all-A distribution cannot produce C, regardless of later scoring errors.
    with pytest.raises(ValueError, match="sampling support"):
        verify_decisions((candidate,), da.plan(request).resolved, recorded)


def test_saved_background_support_is_checked_before_pool_identity(tmp_path: Path):
    pool = da.prepare(
        parts.PreparationSpec(
            parts.Background((1, 0, 0, 0)),
            parts.Retention(count=1, policy="first_eligible"),
            sampling=parts.Sampling(planning.Length(exact=1)),
            budget=parts.CandidateBudget(1),
        ),
        out=tmp_path / "pool",
    )
    assert da.inspect(pool, verify=True).retained_parts == 1
    with sqlite3.connect(pool.path / "pool.sqlite3") as connection:
        value = json.loads(
            connection.execute(
                "SELECT payload FROM candidates WHERE ordinal=1"
            ).fetchone()[0]
        )
        value["part"]["sequence"] = "C"
        connection.execute(
            "UPDATE candidates SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(value), semantic_digest(value)),
        )
    with pytest.raises(ValueError, match="sampling support"):
        da.inspect(pool, verify=True)


@pytest.mark.parametrize("damage", ["core", "error"])
def test_background_verification_rejects_undeclared_scoring_evidence(damage: str):
    request = parts.PreparationSpec(
        parts.Background((1, 0, 0, 0)),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=1)),
        budget=parts.CandidateBudget(1),
    )
    candidate = Candidate(
        1,
        parts.Part(
            "candidate_1",
            "A",
            group="background",
            source="background",
            **(
                {"core_start": 0, "core_end": 1, "core_orientation": "forward"}
                if damage == "core"
                else {}
            ),
        ),
        representative=1 if damage == "core" else None,
        rank=1 if damage == "core" else None,
        retained=damage == "core",
        error="FIMO timeout" if damage == "error" else None,
    )
    recorded = recount(
        (candidate,),
        target=1,
        budget=1,
        stop_reason="candidate_budget" if damage == "core" else "execution_error",
    )
    with pytest.raises(ValueError, match=r"core|scoring"):
        verify_decisions((candidate,), da.plan(request).resolved, recorded)


@pytest.mark.parametrize("strategy", ["consensus", "background"])
def test_new_strategies_share_saved_plan_python_cli_and_recorded_decisions(
    tmp_path: Path, strategy: str
):

    tool = scorer(
        tmp_path,
        ""
        if strategy == "background"
        else "".join(
            f"motif\t\tcandidate_{i}\t1\t2\t+\t3\t0.0625\t\tAC\n" for i in range(4)
        ),
    )
    request = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(
            planning.Length(exact=2), strategy=strategy, base_probabilities=(0, 0, 0, 1)
        ),
        budget=parts.CandidateBudget(4),
        scoring=parts.FimoScoring(executable=tool),
    )
    plan = da.plan(request)
    assert plan.preview["proposal"]["strategy"] == strategy
    assert plan.preview["proposal"]["base_probabilities"] == (0, 0, 0, 1)
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).plan_id == plan.plan_id
    human = CliRunner().invoke(app, ["plan", str(saved)])
    assert human.exit_code == 0, human.output
    assert strategy in human.stdout
    assert "ACGT" in human.stdout
    pool = da.prepare(plan, out=tmp_path / "python")
    result = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == (3 if strategy == "background" else 0), result.output
    assert json.loads(result.stdout)["pool_id"] == pool.pool_id
    with da.inspect(pool, view="candidates", all=True).records() as records:
        candidates = [row.candidate for row in records]
    assert [c.part.sequence for c in candidates] == (
        ["TT"] * 4 if strategy == "background" else ["AC"] * 4
    )
    expected_reasons = ("no_qualifying_hit",) if strategy == "background" else ()
    assert all(c.reasons == expected_reasons for c in candidates)
    assert candidates[0].part.metadata["proposal"]["start"] == (
        None if strategy == "background" else 0
    )
    assert plan.to_dict()["policies"]["sampling"] != "candidate_shake256.v1"
    request.source.path.unlink()
    tool.unlink()
    assert da.inspect(pool, verify=True).retained_parts == (
        0 if strategy == "background" else 1
    )


def test_sampling_override_is_explicit_and_background_source_has_one_distribution_owner(
    tmp_path: Path,
):

    probabilities = [0, 0, 0, 1]
    sampling = parts.Sampling(
        planning.Length(exact=6), base_probabilities=probabilities
    )
    probabilities[3] = 0
    assert sampling.base_probabilities == (0, 0, 0, 1)
    with pytest.raises(ValueError, match="Background"):
        parts.PreparationSpec(
            parts.Background(),
            parts.Retention(count=1, policy="first_eligible"),
            sampling=sampling,
            budget=parts.CandidateBudget(4),
        )
    request = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=sampling,
        budget=parts.CandidateBudget(4),
        scoring=parts.FimoScoring(executable=scorer(tmp_path, "")),
    )
    plan = da.plan(request)
    assert plan.preview["proposal"]["background_source"] == "sampling"
    assert plan.resolved.source.scoring.background == (0.25,) * 4
    pool = da.prepare(plan, out=tmp_path / "pool")
    with da.inspect(pool, view="candidates", all=True).records() as records:
        for row in records:
            p = row.candidate.part
            proposal = p.metadata["proposal"]
            assert set(
                p.sequence[: proposal["start"]] + p.sequence[proposal["end"] :]
            ) <= {"T"}
    assert da.inspect(pool, verify=True).source_parts == 4


@pytest.mark.parametrize("damage", ["missing", "coordinate", "policy"])
def test_verification_checks_saved_proposal_before_pool_digest(
    tmp_path: Path, damage: str
):

    pool = da.prepare(
        parts.PreparationSpec(
            source(tmp_path),
            parts.Retention(count=0, policy="first_eligible"),
            sampling=parts.Sampling(planning.Length(exact=6), strategy="consensus"),
            budget=parts.CandidateBudget(2),
            scoring=parts.FimoScoring(executable=scorer(tmp_path, "")),
        ),
        out=tmp_path / "pool",
    )
    with sqlite3.connect(pool.path / "pool.sqlite3") as connection:
        value = json.loads(
            connection.execute(
                "SELECT payload FROM candidates WHERE ordinal=1"
            ).fetchone()[0]
        )
        metadata = value["part"]["metadata"]
        if damage == "missing":
            metadata.pop("proposal")
        elif damage == "coordinate":
            metadata["proposal"]["end"] = 99
        else:
            metadata["proposal"]["policy"] = "unknown.v1"
        connection.execute(
            "UPDATE candidates SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(value), semantic_digest(value)),
        )
    with pytest.raises(ValueError, match="proposal"):
        da.inspect(pool, verify=True)


def test_intended_placement_does_not_replace_best_hit_geometry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    monkeypatch.setattr(
        preparation, "propose_sequence", lambda **_: mining.Proposal("ACGTAC", 4, 6)
    )
    tool = scorer(tmp_path, "motif\t\tcandidate_0\t1\t2\t+\t3\t0.0625\t\tAC\n")
    pool = da.prepare(
        parts.PreparationSpec(
            source(tmp_path),
            parts.Retention(count=1, policy="first_eligible"),
            sampling=parts.Sampling(planning.Length(exact=6), strategy="consensus"),
            budget=parts.CandidateBudget(1),
            scoring=parts.FimoScoring(executable=tool),
        ),
        out=tmp_path / "pool",
    )
    with da.inspect(pool, view="parts", all=True).records() as rows:
        part = next(rows).part
    assert (part.core_start, part.core_end) == (0, 2)
    assert part.metadata["proposal"]["start"] == 4
    assert part.metadata["proposal"]["end"] == 6
    assert da.inspect(pool, verify=True).retained_parts == 1


@pytest.mark.parametrize("strategy", ["consensus", "background"])
def test_new_proposal_prefixes_are_independent_of_batch_size(
    tmp_path: Path, strategy: str
):

    request = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=12), strategy=strategy),
        budget=parts.CandidateBudget(8),
        scoring=parts.FimoScoring(executable=scorer(tmp_path, "")),
    )
    one = da.prepare(request, out=tmp_path / "one")
    two = da.prepare(
        request.with_changes(budget=parts.CandidateBudget(8, batch_size=2)),
        out=tmp_path / "two",
    )
    with da.inspect(one, view="candidates", all=True).records() as rows:
        expected = [r.candidate for r in rows]
    with da.inspect(two, view="candidates", all=True).records() as rows:
        assert [r.candidate for r in rows] == expected


def test_consensus_matches_pinned_canonical_column_fixtures():

    fixture = Path(__file__).parents[2] / "fixtures/workflow/densegen-consensus-v1.json"
    for case in json.loads(fixture.read_text())["cases"]:
        rows = tuple(tuple(r) for r in case["probabilities"])
        model = Motif(
            case["name"],
            rows,
            (0.25,) * 4,
            tuple((0, 0, 0, 0) for _ in rows),
            "fixture",
        )
        assert mining.consensus(model) == case["consensus"]
