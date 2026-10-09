"""Persist candidate bytes and placements for rejected and duplicate attempts.

Author: Eric J. South.
"""

import json
import sqlite3
from pathlib import Path
from types import SimpleNamespace

import pytest
from ortools.linear_solver import pywraplp
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts import Attempt
from dense_arrays.cli import app
from dense_arrays.generation.acceptance import evaluate, realize
from dense_arrays.generation.packing import build_optimizer, restore_packing
from dense_arrays.solution import DenseArray


def rejected_request(*, padding: bool = False):
    return planning.DesignSpec(
        parts=[parts.Part("a", "AA"), parts.Part("b", "CC")],
        length=planning.Length(exact=6) if padding else planning.Length(maximum=4),
        assembly=planning.Assembly(padding=planning.Padding(side="left", max_trials=2))
        if padding
        else None,
        strands="single",
        limits=planning.Limits(attempts=1),
        requirements=[
            planning.Fixed("first", "a", "forward"),
            planning.Fixed("second", "b", "forward"),
            planning.Avoid("no-A", patterns=("A",), strands="forward"),
        ],
    )


@pytest.mark.parametrize("padding", [False, True])
def test_rejected_candidate_is_inspectable_without_generation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, padding: bool
):
    plan = da.plan(rejected_request(padding=padding))
    run = da.run(plan, out=tmp_path / "run")
    with da.inspect(run, view="attempts").records() as records:
        attempt = next(records)
    assert attempt.outcome == "rejected"
    assert attempt.candidate is not None
    packed = attempt.candidate.packed
    final = attempt.candidate.final
    assert len(packed.sequence) == 4
    assert final is not None
    assert len(final.sequence) == (6 if padding else 4)
    assert {p.feature_id for p in packed.placements} == {"a", "b"}
    assert "assembly" not in packed.provenance
    if padding:
        assert final.provenance["assembly"]["trial"] == 2
        assert min(p.start for p in final.placements) == 2
    assert attempt.to_dict()["evidence"]["requirements"] == list(evaluate(final, plan))
    assert Attempt.from_dict(attempt.to_dict()) == attempt

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("candidate inspection attempted to regenerate evidence")

    monkeypatch.setattr(pywraplp.Solver, "CreateSolver", forbidden)
    monkeypatch.setattr("dense_arrays.generation.assembly.padding_dna", forbidden)
    assert da.inspect(run, verify=True).verified
    response = CliRunner().invoke(
        app, ["inspect", str(run.path), "--view", "attempts", "--json"]
    )
    assert response.exit_code == 0, response.output
    assert json.loads(response.stdout)["records"] == [attempt.to_dict()]


def test_duplicate_keeps_identity_specific_packing_for_exclusion_replay(tmp_path: Path):
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAA"), parts.Part("b", "AAA")],
            length=planning.Length(maximum=3),
            strands="single",
            target=planning.Target(count=3),
        )
    )
    run = da.run(plan, out=tmp_path / "run")
    with da.inspect(run, view="attempts").records() as records:
        attempts = list(records)
    assert [a.outcome for a in attempts] == ["accepted", "duplicate", "no_candidate"]
    assert attempts[-1].candidate is None
    assert {a.candidate.packed.placements[0].feature_id for a in attempts[:2]} == {
        "a",
        "b",
    }
    optimizer = build_optimizer(plan, seconds=2)
    for attempt in attempts[:2]:
        packed = restore_packing(attempt.candidate.packed, plan)
        assert packed.sequence == "AAA"
        optimizer.forbid(packed)
    assert optimizer.solve_report().status.value == "infeasible"
    assert da.inspect(run, verify=True).counts["duplicate"] == 1


def test_rejected_checks_are_verified_against_candidate_even_with_valid_checksum(
    tmp_path: Path,
):
    run = da.run(rejected_request(), out=tmp_path / "run")
    with sqlite3.connect(run.path / "run.sqlite3") as connection:
        revision, payload = connection.execute(
            "SELECT revision,payload FROM attempts ORDER BY revision DESC LIMIT 1"
        ).fetchone()
        record = json.loads(payload)
        record["evidence"]["requirements"][-1]["observed"] = []
        connection.execute(
            "UPDATE attempts SET payload=?,digest=? WHERE revision=?",
            (
                canonical_json(record),
                semantic_digest(record),
                revision,
            ),
        )
    with pytest.raises(ValueError, match="candidate requirement"):
        da.inspect(run, verify=True)


def test_legacy_attempt_has_no_invented_candidate_and_unknown_candidate_is_rejected(
    tmp_path: Path,
):
    legacy = {
        "schema": "dense_arrays.attempt.v1",
        "attempt_id": 1,
        "cell_id": "default",
        "outcome": "rejected",
        "evidence": {},
    }
    assert Attempt.from_dict(legacy).candidate is None
    assert Attempt.from_dict(legacy).to_dict() == legacy
    run = da.run(rejected_request(), out=tmp_path / "run")
    with da.inspect(run, view="attempts").records() as records:
        record = next(records).to_dict()
    record["evidence"]["candidate"]["schema"] = "dense_arrays.candidate.v99"
    with pytest.raises(ValueError, match="candidate schema"):
        Attempt.from_dict(record)


@pytest.mark.parametrize("trials", [0, 1])
def test_interrupted_assembly_retains_only_evaluated_final_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, trials: int
):
    ticks = iter([0] * trials + [float("inf")])
    monkeypatch.setattr(
        "dense_arrays.generation.assembly.time",
        SimpleNamespace(monotonic=lambda: next(ticks)),
    )
    run = da.run(rejected_request(padding=True), out=tmp_path / "run")
    with da.inspect(run, view="attempts").records() as records:
        attempt = next(records)
    assert attempt.outcome == "no_candidate"
    assert attempt.evidence["code"] == "active_time_limit"
    assert attempt.evidence["assembly_trials"] == trials
    assert attempt.candidate is not None
    assert (attempt.candidate.final is None) == (trials == 0)
    assert da.inspect(run, verify=True).termination_reason == "active_time_limit"


def test_restored_offsets_must_describe_an_exact_packing_path():
    plan = da.plan(
        planning.DesignSpec(
            parts=[parts.Part("a", "AAAA"), parts.Part("b", "AA")],
            length=planning.Length(maximum=4),
            strands="single",
        )
    )
    contained = DenseArray(["AAAA", "AA"], 4, [0, 1], [None, None])
    packed = realize(contained, plan, source_id="candidate")
    with pytest.raises(ValueError, match="packing path"):
        restore_packing(packed, plan)
