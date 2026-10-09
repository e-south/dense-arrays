"""Optional PWM exclusions preserve declared thresholds and complete hit evidence.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.parts.scoring import FimoHit

from .test_fimo import scorer
from .test_motifs import artifact


def rule(tmp_path: Path, **options: object):
    return parts.PWMExclusion(
        "exclude",
        (parts.PWMArtifact(tmp_path / "motif.json"),),
        parts.FimoScoring(executable=tmp_path / "fimo"),
        **options,
    )


def test_exclusion_uses_an_explicit_scale_and_strict_threshold(tmp_path: Path):
    hit = FimoHit(0, 2, "forward", "AC", 2, 0.01, 4)
    assert rule(tmp_path).rejects(hit)
    assert not rule(tmp_path).rejects(None)
    relative = rule(
        tmp_path, reject="score_above", score_field="fraction_of_max", threshold=0.5
    )
    assert not relative.rejects(hit)
    assert rule(
        tmp_path, reject="score_above", score_field="raw", threshold=1.9
    ).rejects(hit)
    with pytest.raises(ValueError, match="positive theoretical maximum"):
        relative.rejects(FimoHit(0, 2, "forward", "AC", 0, 0.01, 0))


@pytest.mark.parametrize(
    "options",
    [
        {"reject": "unknown"},
        {"threshold": 0.5},
        {"score_field": "raw"},
        {"reject": "score_above"},
        {"reject": "score_above", "score_field": "normalized", "threshold": 0.5},
        {"reject": "score_above", "score_field": "raw", "threshold": True},
    ],
)
def test_exclusion_rejects_ambiguous_policies(tmp_path: Path, options: dict):
    with pytest.raises((TypeError, ValueError)):
        rule(tmp_path, **options)


def test_background_exclusions_round_trip_and_verify_without_scoring(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    from dense_arrays.workflow import preparation  # noqa: PLC0415
    from dense_arrays.workflow.inputs import read_source  # noqa: PLC0415

    motif = tmp_path / "motif.json"
    motif.write_text(json.dumps(artifact()))
    tool = scorer(
        tmp_path,
        "motif\t\tcandidate_0\t1\t2\t+\t3\t0.0625\t\tAC\n"
        "motif\t\tcandidate_2\t3\t4\t-\t3\t0.0625\t\tAC\n",
    )
    exclude = parts.PWMExclusion(
        "motif_hit",
        (parts.PWMArtifact(motif),),
        parts.FimoScoring(executable=tool, hit_pvalue_max=0.1),
    )
    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=2, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(candidates=4),
        screening=(
            exclude,
            planning.GC("gc", "sequence", 0.25, 1),
            planning.Avoid("literal", ("AA",), strands="forward"),
        ),
    )
    plan = da.plan(request)
    assert plan.preview["screening_motifs"] == 1
    assert plan.preview["required_tools"] == ("fimo",)
    assert plan.preview["screening_window_bound"] == 26
    saved_request = tmp_path / "request.json"
    da.export(request, view="request", out=saved_request)
    assert da.plan(read_source(saved_request)).plan_id == plan.plan_id
    saved_plan = tmp_path / "plan.json"
    plan.write(saved_plan)
    assert read_source(saved_plan).plan_id == plan.plan_id
    monkeypatch.setattr(
        preparation,
        "sample_sequence",
        lambda **kw: ("ACAC", "AAAA", "TTGT", "CCCC")[kw["index"] - 1],
    )
    pool = da.prepare(plan, out=tmp_path / "pool")
    quality = da.inspect(pool, view="quality").to_dict()
    assert quality["rejections"] == {"motif_hit": 2, "gc": 1, "literal": 1}
    assert quality["counts"]["eligibility_rejected"] == 3
    assert quality["counts"]["retained"] == 1
    human = CliRunner().invoke(app, ["inspect", str(pool.path), "--view", "quality"])
    assert "motif_hit=2" in human.stdout
    result = CliRunner().invoke(
        app, ["prepare", str(saved_plan), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 3, result.output
    assert json.loads(result.stdout)["pool_id"] == pool.pool_id
    motif.unlink()
    tool.unlink()
    assert da.inspect(pool, verify=True).retained_parts == 1
    with da.inspect(pool, view="parts", all=True).records() as records:
        assert [row.part.sequence for row in records] == ["CCCC"]


def background_request(
    tmp_path: Path, *, tool: Path | None = None
) -> parts.PreparationSpec:
    motif = tmp_path / "motif.json"
    motif.write_text(json.dumps(artifact()))
    return parts.PreparationSpec(
        parts.Background(base_probabilities=(1, 0, 0, 0)),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=4)),
        budget=parts.CandidateBudget(candidates=4, batch_size=2),
        screening=(
            parts.PWMExclusion(
                "exclude",
                (parts.PWMArtifact(motif),),
                parts.FimoScoring(executable=tool or scorer(tmp_path, "")),
            ),
        ),
    )


def test_background_preview_includes_scoring_and_each_batch_calibration(tmp_path: Path):
    plan = da.plan(background_request(tmp_path))
    assert "scoring" in plan.preview["screening_stages"]
    assert plan.preview["screening_window_bound"] == 28


def test_background_scorer_failure_is_not_an_exclusion(tmp_path: Path):
    from .test_fimo import executable  # noqa: PLC0415

    tool = executable(
        tmp_path, 'if [ "$1" = "--version" ]; then echo 5.5.9; else exit 92; fi'
    )
    pool = da.prepare(background_request(tmp_path, tool=tool), out=tmp_path / "failed")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["stop_reason"] == "execution_error"
    assert report["counts"]["execution_error"] == 2
    assert report["counts"]["eligibility_rejected"] == 0
    assert report["rejections"] == {}
    assert da.inspect(pool, verify=True).state == "incomplete"


def test_screen_sources_are_verified_before_output_ownership(tmp_path: Path):
    plan = da.plan(background_request(tmp_path))
    (tmp_path / "motif.json").write_text("{}")
    with pytest.raises(ValueError, match="changed"):
        da.prepare(plan, out=tmp_path / "untouched")
    assert not (tmp_path / "untouched").exists()


def test_multi_motif_rule_binds_every_model_and_rejects_repeated_ids(tmp_path: Path):
    from dataclasses import replace  # noqa: PLC0415

    request = background_request(tmp_path)
    first = request.screening[0]
    second = tmp_path / "second.json"
    value = artifact()
    value["motif_id"] = "second"
    second.write_text(json.dumps(value))
    request = request.with_changes(
        screening=(replace(first, motifs=(*first.motifs, parts.PWMArtifact(second))),)
    )
    plan = da.plan(request)
    assert plan.preview["screening_motifs"] == 2
    assert plan.preview["screening_window_bound"] == 56
    pool = da.prepare(plan, out=tmp_path / "multi")
    assert da.inspect(pool, verify=True).retained_parts == 1
    value["motif_id"] = "example"
    second.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="identities must be unique"):
        da.plan(request)


def test_missing_no_hit_observation_is_not_treated_as_passed(tmp_path: Path):
    import sqlite3  # noqa: PLC0415

    from dense_arrays._record_validation import (  # noqa: PLC0415
        canonical_json,
        semantic_digest,
    )
    from dense_arrays.artifacts.errors import ArtifactIntegrityError  # noqa: PLC0415

    pool = da.prepare(background_request(tmp_path), out=tmp_path / "pool")
    with sqlite3.connect(pool.path / "pool.sqlite3") as connection:
        payload = json.loads(
            connection.execute(
                "SELECT payload FROM candidates WHERE ordinal=1"
            ).fetchone()[0]
        )
        assert payload["screening"][0]["hit"] is None
        del payload["screening"]
        connection.execute(
            "UPDATE candidates SET payload=?,digest=? WHERE ordinal=1",
            (canonical_json(payload), semantic_digest(payload)),
        )
    with pytest.raises(ArtifactIntegrityError, match=r"missing.*screen observations"):
        da.inspect(pool, verify=True)


def test_undefined_exclusion_ratio_records_execution_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    from dataclasses import replace  # noqa: PLC0415

    from dense_arrays.parts.scoring import FimoResult  # noqa: PLC0415
    from dense_arrays.workflow import preparation  # noqa: PLC0415

    request = background_request(tmp_path)
    request = request.with_changes(
        screening=(
            replace(
                request.screening[0],
                reject="score_above",
                score_field="fraction_of_max",
                threshold=0.5,
            ),
        )
    )
    monkeypatch.setattr(
        preparation,
        "scan_fimo",
        lambda binding, sequences: FimoResult(
            binding.binding_id,
            (FimoHit(0, 2, "forward", "AA", 0, 0.01, 0),) * len(sequences),
            len(sequences) * 6,
            2,
            len(sequences),
            theoretical_max=0,
        ),
    )
    pool = da.prepare(request, out=tmp_path / "undefined")
    report = da.inspect(pool, view="quality").to_dict()
    assert report["counts"]["execution_error"] == 2
    assert report["counts"]["eligibility_rejected"] == 0
    assert report["rejections"] == {}


def test_screen_window_limit_fails_in_preflight(tmp_path: Path):
    from dataclasses import replace  # noqa: PLC0415

    request = background_request(tmp_path)
    screen = request.screening[0]
    request = request.with_changes(
        screening=(
            replace(
                screen,
                scoring=replace(screen.scoring, limits=parts.ScoringLimits(windows=13)),
            ),
        )
    )
    with pytest.raises(ValueError, match=r"screen batch.*window limit"):
        da.prepare(request, out=tmp_path / "untouched")
    assert not (tmp_path / "untouched").exists()


@pytest.mark.parametrize(
    "case",
    json.loads(
        (
            Path(__file__).parents[2] / "fixtures/workflow/densegen-exclusion-v1.json"
        ).read_text()
    )["cases"],
)
def test_exclusion_decisions_match_pinned_acceptance_loop(tmp_path: Path, case: dict):
    exclusion = rule(
        tmp_path,
        reject=case["reject"],
        score_field="fraction_of_max" if case["reject"] == "score_above" else None,
        threshold=case["threshold"],
    )
    accepted = []
    for record in case["records"]:
        hits = [
            None
            if h is None
            else FimoHit(0, 2, "forward", "AC", h["raw"], 0.01, h["maximum"])
            for h in record["hits"]
        ]
        if not any(exclusion.rejects(hit) for hit in hits):
            accepted.append(record["sequence"])
    assert accepted == case["accepted"]
