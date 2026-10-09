"""Motif windows retain source coordinates and bind the actual scoring model.

Author: Eric J. South.
"""

import json
from dataclasses import replace
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.parts.motifs import Motif
from dense_arrays.parts.motifs.windows import select_window
from dense_arrays.workflow import preparation
from dense_arrays.workflow.inputs import read_source

from .test_fimo import scorer
from .test_motifs import artifact


def wide_source(tmp_path: Path):
    value = artifact()
    uniform = dict.fromkeys("ACGT", 0.25)
    value["probabilities"] = [uniform, *value["probabilities"], uniform]
    zero = dict.fromkeys("ACGT", 0)
    value["log_odds"] = [zero, *value["log_odds"], zero]
    value["length"] = 4
    path = tmp_path / "wide.json"
    path.write_text(json.dumps(value))
    return parts.PWMArtifact(path)


def test_window_selects_information_not_background_score_or_label():
    motif = Motif(
        "example",
        ((0.25,) * 4, (1, 0, 0, 0), (0, 1, 0, 0), (0.25,) * 4),
        (0.7, 0.1, 0.1, 0.1),
        ((9, 8, 7, 6), (2, -1, -1, -1), (-1, 3, -1, -1), (9, 8, 7, 6)),
        "test",
    )
    selected = select_window(motif, parts.MotifWindow(length=2, background="uniform"))
    assert (selected.start, selected.end) == (1, 3)
    assert selected.information_bits == 4
    assert selected.source_information_bits == 4
    assert selected.motif.probabilities == ((1, 0, 0, 0), (0, 1, 0, 0))
    assert selected.motif.log_odds == ((2, -1, -1, -1), (-1, 3, -1, -1))
    assert selected.motif.model_id != motif.model_id
    assert selected.motif.motif_id == motif.motif_id
    assert selected.to_dict()["source_model_id"] == motif.model_id
    assert selected.to_dict()["retained_information_fraction"] == 1
    assert (
        select_window(motif, parts.MotifWindow(length=1, background="uniform")).start
        == 1
    )
    full = select_window(motif, parts.MotifWindow(length=4, background="uniform"))
    assert full.motif == motif
    assert full.information_bits == 4
    uniform = replace(motif, probabilities=((0.25,) * 4,) * 4)
    empty = select_window(uniform, parts.MotifWindow(length=2, background="uniform"))
    assert empty.start == 0
    assert empty.information_bits == 0
    assert empty.to_dict()["retained_information_fraction"] is None


@pytest.mark.parametrize("length", [0, -1, True, 2.5, "2"])
def test_window_rejects_invalid_lengths(length: object):
    with pytest.raises((ValueError, TypeError)):
        parts.MotifWindow(length=length)


def test_window_plan_pool_cli_and_source_free_evidence(tmp_path: Path):
    source = wide_source(tmp_path)
    tool = scorer(tmp_path, "")
    request = parts.PreparationSpec(
        replace(source, window=parts.MotifWindow(length=2)),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=6), strategy="consensus"),
        budget=parts.CandidateBudget(6),
        scoring=parts.FimoScoring(executable=tool),
        seed=7,
    )
    plan = da.plan(request)
    evidence = plan.preview["motif_window"]
    assert (evidence["start"], evidence["end"]) == (1, 3)
    assert evidence["policy"] == "max_relative_entropy.v1"
    assert evidence["background_source"] == "motif"
    assert evidence["background"] == [0.25] * 4
    assert plan.resolved.source.input.motif.width == 4
    assert plan.resolved.source.scoring.motif.width == 2
    assert plan.to_dict()["source"]["window"] == evidence
    path = tmp_path / "plan.json"
    plan.write(path)
    assert read_source(path).to_dict() == plan.to_dict()
    pool = da.prepare(plan, out=tmp_path / "python")
    cli = CliRunner().invoke(
        app, ["prepare", str(path), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert cli.exit_code == 0, cli.output
    assert json.loads(cli.stdout)["pool_id"] == pool.pool_id
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        for row in rows:
            candidate = row.candidate
            intended = candidate.part.metadata["proposal"]
            assert candidate.part.sequence[intended["start"] : intended["end"]] == "AC"
    source.path.unlink()
    tool.unlink()
    assert da.inspect(pool, verify=True).source_parts == 6
    assert read_source(path).preview["motif_window"] == evidence


def test_window_bounds_fail_before_scoring_preflight(tmp_path: Path):
    source = wide_source(tmp_path)
    request = parts.PreparationSpec(
        replace(source, window=parts.MotifWindow(length=5)),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=6)),
        budget=parts.CandidateBudget(1),
        scoring=parts.FimoScoring(executable=tmp_path / "absent"),
    )
    with pytest.raises(ValueError, match="exceeds"):
        da.plan(request)


@pytest.mark.parametrize(
    "field,value", [("start", 0), ("information_bits", 0), ("policy", "unknown")]
)
def test_saved_window_rejects_changed_meaning(
    tmp_path: Path, field: str, value: object
):
    request = parts.PreparationSpec(
        replace(wide_source(tmp_path), window=parts.MotifWindow(length=2)),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=2)),
        budget=parts.CandidateBudget(1),
        scoring=parts.FimoScoring(executable=scorer(tmp_path, "")),
    )
    value_dict = da.plan(request).to_dict()
    value_dict["source"]["window"][field] = value
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(value_dict))
    with pytest.raises(ValueError, match=r"window|motif"):
        read_source(path)


def test_legacy_fixed_window_fixture_preserves_coordinates_and_aligned_rows():
    fixture = Path(__file__).parents[2] / "fixtures/workflow/densegen-windows-v1.json"
    for case in json.loads(fixture.read_text())["cases"]:
        motif = Motif(
            case["name"],
            case["probabilities"],
            (0.25,) * 4,
            case["log_odds"],
            "synthetic-fixture",
        )
        result = select_window(
            motif, parts.MotifWindow(case["length"], background="uniform")
        )
        assert result.start == case["start"]
        assert result.motif.probabilities == tuple(
            map(tuple, case["selected_probabilities"])
        )
        assert result.motif.log_odds == tuple(map(tuple, case["selected_log_odds"]))
        if case["name"] == "whole":
            assert result.information_bits == 4
            assert case["information_bits"] == 0
        else:
            assert result.information_bits == pytest.approx(
                case["information_bits"], abs=1e-9
            )


def test_trimmed_exclusion_keeps_window_preview_binding_and_score_width(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    motif = replace(wide_source(tmp_path), window=parts.MotifWindow(2))
    tool = scorer(tmp_path, "motif\t\tcandidate_0\t1\t2\t+\t3\t0.0625\t\tAC\n")
    request = parts.PreparationSpec(
        parts.Background(base_probabilities=(0.5, 0.5, 0, 0)),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=2)),
        budget=parts.CandidateBudget(1),
        screening=(
            parts.PWMExclusion("exclude", (motif,), parts.FimoScoring(executable=tool)),
        ),
    )
    plan = da.plan(request)
    assert plan.to_dict()["policies"]["motif_window"] == "max_relative_entropy.v1"
    assert plan.preview["screening_window_bound"] == 4
    assert plan.preview["screening_windows"][0]["window"]["start"] == 1
    path = tmp_path / "screened.json"
    plan.write(path)
    assert read_source(path).to_dict() == plan.to_dict()
    monkeypatch.setattr(preparation, "sample_sequence", lambda **_: "AC")
    pool = da.prepare(plan, out=tmp_path / "pool")
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidate = next(rows).candidate
    assert candidate.reasons == ("exclude",)
    assert candidate.screening[0].hit.core == "AC"
    assert candidate.screening[0].hit.theoretical_max == 3
    assert da.inspect(pool, verify=True).source_parts == 1


def test_human_preview_explains_selected_coordinates_and_information(tmp_path: Path):
    request = parts.PreparationSpec(
        replace(wide_source(tmp_path), window=parts.MotifWindow(2)),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=2)),
        budget=parts.CandidateBudget(1),
        scoring=parts.FimoScoring(executable=scorer(tmp_path, "")),
    )
    path = tmp_path / "plan.json"
    da.plan(request).write(path)
    result = CliRunner().invoke(app, ["plan", str(path)])
    assert result.exit_code == 0, result.output
    assert "[1, 3)" in result.output
    assert "information" in result.output


def test_declared_background_changes_discrimination_window_without_changing_model():
    motif = Motif(
        "example",
        ((1, 0, 0, 0), (0, 1, 0, 0)),
        (0.5, 0.125, 0.125, 0.25),
        ((9, -9, -9, -9), (-9, 1, -9, -9)),
        "test",
    )
    uniform = select_window(motif, parts.MotifWindow(1, background="uniform"))
    assert uniform.start == 0
    assert uniform.information_bits == 2
    relative = select_window(motif, parts.MotifWindow(1, background="motif"))
    assert relative.start == 1
    assert relative.information_bits == 3
    assert relative.source_information_bits == 4
    assert relative.to_dict()["retained_information_fraction"] == 0.75
    assert relative.to_dict()["discarded_information_bits"] == 1
    override = select_window(
        motif, parts.MotifWindow(1, background=(0.125, 0.5, 0.25, 0.125))
    )
    assert override.start == 0
    assert override.information_bits == 3
    assert override.to_dict()["background_source"] == "explicit"
    assert override.motif.background == motif.background
    assert override.motif.log_odds == motif.log_odds[:1]
    assert parts.MotifWindow.from_dict(
        parts.MotifWindow(1).to_dict()
    ) == parts.MotifWindow(1)


@pytest.mark.parametrize(
    "background", ["scoring", (0, 0.5, 0.25, 0.25), (0.2,) * 4, (True, 0, 0, 0)]
)
def test_window_background_must_be_positive_normalized_and_explicit(background: object):
    with pytest.raises((ValueError, TypeError)):
        parts.MotifWindow(1, background=background)


def test_uniform_model_can_be_informative_against_biased_background():
    motif = Motif(
        "uniform", ((0.25,) * 4,), (0.5, 0.125, 0.125, 0.25), ((0,) * 4,), "test"
    )
    result = select_window(motif, parts.MotifWindow(1))
    assert result.information_bits == 0.25
    assert result.to_dict()["retained_information_fraction"] == 1
