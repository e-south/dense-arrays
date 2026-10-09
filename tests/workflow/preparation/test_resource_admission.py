"""Candidate base admission precedes input reads and allocation, not inspection.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.planning.preparation import PreparationPlan
from dense_arrays.planning.preparation.requests import (
    preparation_from_dict,
    preparation_to_dict,
)
from dense_arrays.planning.preparation.sampled import SampledPreparation
from dense_arrays.planning.preparation.sets import SetPreparation
from dense_arrays.workflow.inputs import read_source

from .test_sampled import background_request


@pytest.mark.parametrize("name", ["batch_bases", "total_bases"])
@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_base_caps_require_positive_integers(name: str, value: object):
    with pytest.raises((ValueError, TypeError), match=name):
        parts.CandidateBudget(1, **{name: value})


def test_base_caps_are_keyword_only():
    with pytest.raises(TypeError):
        parts.CandidateBudget(1, None, 1, 4, 4)


@pytest.mark.parametrize(
    "batch_bases,total_bases,match",
    [(7, 40, r"8.*batch_bases.*7"), (8, 39, r"40.*total_bases.*39")],
)
def test_batch_and_total_admission_are_independent(
    batch_bases: int, total_bases: int, match: str
):
    budget = parts.CandidateBudget(
        10, batch_size=2, batch_bases=batch_bases, total_bases=total_bases
    )
    with pytest.raises(ValueError, match=match):
        budget.admit(4)


def test_exact_caps_admit_the_full_candidate_batch_and_budget(tmp_path: Path):
    budget = parts.CandidateBudget(6, batch_size=2, batch_bases=8, total_bases=24)
    budget.admit(4)
    request = background_request().with_changes(budget=budget)
    result = da.prepare(request, out=tmp_path / "pool")
    assert da.inspect(result, verify=True).source_parts == 6
    # An oversized batch_size does not invent candidates beyond the total budget.
    parts.CandidateBudget(1, batch_size=1000, batch_bases=4, total_bases=4).admit(4)


def test_length_ranges_use_the_maximum_for_admission():
    request = background_request().with_changes(
        sampling=parts.Sampling(parts.LengthRange(1, 5)),
        budget=parts.CandidateBudget(2, batch_size=2, batch_bases=9),
    )
    with pytest.raises(ValueError, match=r"10.*batch_bases.*9"):
        da.plan(request)


@pytest.mark.parametrize("source", ["background", "pwm", "screened"])
@pytest.mark.parametrize("operation", ["plan", "prepare"])
def test_enormous_declarations_fail_before_sources_or_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source: str, operation: str
):
    from dense_arrays.planning.preparation import sampled  # noqa: PLC0415

    request = background_request().with_changes(
        sampling=parts.Sampling(planning.Length(exact=10**12)),
        budget=parts.CandidateBudget(1, seconds=0.001, batch_size=1),
    )
    motif = parts.PWMArtifact(tmp_path / "absent.json")
    scoring = parts.FimoScoring(executable=tmp_path / "absent-fimo")
    if source == "pwm":
        request = request.with_changes(source=motif, scoring=scoring)
    elif source == "screened":
        request = request.with_changes(
            screening=(parts.PWMExclusion("exclude", (motif,), scoring),)
        )
    monkeypatch.setattr(
        sampled,
        "resolve_screens",
        lambda *_a: pytest.fail("input/scorer resolution preceded base admission"),
    )
    output = tmp_path / "untouched"
    call = da.plan if operation == "plan" else da.prepare
    options = {} if operation == "plan" else {"out": output}
    with pytest.raises(ValueError, match=r"1000000000000.*batch_bases.*1000000"):
        call(request, **options)
    assert not output.exists()


def test_default_budget_encoding_preserves_existing_plan_identity():
    request = background_request()
    expected = {"candidates": 6, "seconds": None, "batch_size": 1000}
    assert request.budget.to_dict() == expected
    assert preparation_to_dict(request)["budget"] == expected
    assert (
        da.plan(request).plan_id
        == (
            "43ec9244fcb90a3082dd6e7007598c8ede2510c41ee578038aea7e7c07fb1645"  # pragma: allowlist secret  # noqa: E501
        )
    )


def test_custom_caps_round_trip_with_matching_python_cli_results(tmp_path: Path):
    request = background_request().with_changes(
        budget=parts.CandidateBudget(6, batch_size=2, batch_bases=8, total_bases=24)
    )
    encoded = preparation_to_dict(request)
    assert encoded["budget"] == {
        "candidates": 6,
        "seconds": None,
        "batch_size": 2,
        "batch_bases": 8,
        "total_bases": 24,
    }
    assert preparation_from_dict(encoded) == request
    resolved = da.plan(request)
    defaults = request.with_changes(budget=parts.CandidateBudget(6, batch_size=2))
    assert resolved.plan_id != da.plan(defaults).plan_id
    saved = tmp_path / "request.json"
    da.export(request, view="request", out=saved)
    output = tmp_path / "plan.json"
    cli_plan = CliRunner().invoke(
        app, ["plan", str(saved), "--out", str(output), "--json"]
    )
    assert cli_plan.exit_code == 0, cli_plan.output
    assert read_source(output).to_dict() == resolved.to_dict()
    pool = da.prepare(resolved, out=tmp_path / "python")
    cli_prepare = CliRunner().invoke(
        app,
        ["prepare", str(output), "--out", str(tmp_path / "cli"), "--json"],
    )
    assert cli_prepare.exit_code == 3, cli_prepare.output
    assert json.loads(cli_prepare.stdout)["pool_id"] == pool.pool_id


def _saved_oversized_plan() -> PreparationPlan:
    """Represent an inspectable historical recipe without allocating candidates."""
    request = background_request().with_changes(
        sampling=parts.Sampling(planning.Length(exact=10**12)),
        budget=parts.CandidateBudget(1),
    )
    return PreparationPlan(SampledPreparation(request, request.source))


def test_saved_plan_inspection_stays_available_but_execution_rechecks_caps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    saved = _saved_oversized_plan()
    encoded = saved.to_dict()
    restored = PreparationPlan.from_dict(encoded)
    assert restored.to_dict() == encoded
    assert da.inspect(restored, view="plan").plan_id == saved.plan_id
    monkeypatch.setattr(
        SampledPreparation,
        "verify_inputs",
        lambda *_a: pytest.fail("input verification preceded execution admission"),
    )
    output = tmp_path / "untouched"
    with pytest.raises(ValueError, match="batch_bases"):
        da.prepare(restored, out=output)
    assert not output.exists()


def test_set_admits_every_recipe_before_inputs_output_or_first_mining(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dense_arrays.workflow import preparation  # noqa: PLC0415

    first = da.plan(background_request())
    restored = PreparationPlan(
        SetPreparation(
            {"first": first.resolved, "second": _saved_oversized_plan().resolved}
        )
    )
    monkeypatch.setattr(
        preparation, "_mine", lambda *_a: pytest.fail("mining preceded set admission")
    )
    monkeypatch.setattr(
        SampledPreparation,
        "verify_inputs",
        lambda *_a: pytest.fail("input verification preceded set admission"),
    )
    output = tmp_path / "untouched"
    with pytest.raises(ValueError, match="batch_bases"):
        da.prepare(restored, out=output)
    assert not output.exists()
