"""Bounded preparation length draws preserve independent candidate streams.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.parts import mining, sampling
from dense_arrays.workflow import preparation
from dense_arrays.workflow.inputs import read_source

from .test_fimo import scorer
from .test_motifs import source


@pytest.mark.parametrize(
    "minimum,maximum", [(0, 4), (4, 3), (True, 4), (2, False), (1.5, 4)]
)
def test_length_range_rejects_invalid_integer_bounds(minimum: object, maximum: object):
    with pytest.raises((ValueError, TypeError)):
        parts.LengthRange(minimum, maximum)


def test_range_preparation_shares_preview_cli_and_candidate_evidence(tmp_path: Path):
    request = parts.PreparationSpec(
        parts.Background(base_probabilities=(1, 0, 0, 0)),
        parts.Retention(count=3, policy="first_eligible"),
        sampling=parts.Sampling(parts.LengthRange(2, 4)),
        budget=parts.CandidateBudget(30),
        seed=7,
    )
    plan = da.plan(request)
    assert plan.preview["candidate_bases_bound"] == 120
    assert plan.preview["sampled_length"]["distribution"] == "uniform_integer"
    assert plan.to_dict()["policies"]["length"] == "candidate_length_shake256.v1"
    path = tmp_path / "plan.json"
    plan.write(path)
    assert read_source(path).to_dict() == plan.to_dict()
    pool = da.prepare(plan, out=tmp_path / "python")
    assert da.inspect(pool, verify=True).retained_parts == 3
    with da.inspect(pool, view="parts", all=True).records() as rows:
        assert {row.part.sequence for row in rows} == {"AA", "AAA", "AAAA"}
    result = CliRunner().invoke(
        app, ["prepare", str(path), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["pool_id"] == pool.pool_id
    changed = da.prepare(
        request.with_changes(budget=parts.CandidateBudget(30, batch_size=2)),
        out=tmp_path / "batches",
    )
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidates = [row.candidate for row in rows]
    with da.inspect(changed, view="candidates", all=True).records() as rows:
        assert [row.candidate for row in rows] == candidates
    human = CliRunner().invoke(app, ["plan", str(path)])
    assert human.exit_code == 0, human.output
    assert "2..4" in human.stdout


def test_preview_never_draws_lengths_and_degenerate_range_keeps_sequences(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=3, policy="first_eligible"),
        sampling=parts.Sampling(parts.LengthRange(5, 5)),
        budget=parts.CandidateBudget(5),
        seed=7,
    )

    def forbidden(**_: object) -> int:
        pytest.fail("preview drew candidate lengths")

    with monkeypatch.context() as scoped:
        scoped.setattr(sampling, "sample_length", forbidden)
        assert da.plan(request).preview["candidate_bases_bound"] == 25
    ranged = da.prepare(request, out=tmp_path / "range")
    exact = da.prepare(
        request.with_changes(sampling=parts.Sampling(planning.Length(exact=5))),
        out=tmp_path / "exact",
    )
    with da.inspect(ranged, view="candidates", all=True).records() as rows:
        candidates = [r.candidate for r in rows]
    with da.inspect(exact, view="candidates", all=True).records() as rows:
        assert [r.candidate for r in rows] == candidates


def test_pwm_and_exclusion_admission_use_maximum_length(tmp_path: Path):

    motif = source(tmp_path)
    tool = scorer(tmp_path, "")
    request = parts.PreparationSpec(
        motif,
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(parts.LengthRange(2, 6)),
        budget=parts.CandidateBudget(4),
        scoring=parts.FimoScoring(
            executable=tool, limits=parts.ScoringLimits(windows=41)
        ),
    )
    with pytest.raises(ValueError, match="window limit"):
        da.plan(request)
    plan = da.plan(
        request.with_changes(
            scoring=parts.FimoScoring(
                executable=tool, limits=parts.ScoringLimits(windows=42)
            )
        )
    )
    assert plan.preview["candidate_bases_bound"] == 24
    with pytest.raises(ValueError, match="shorter than selected motif") as failed:
        da.plan(request.with_changes(sampling=parts.Sampling(parts.LengthRange(1, 6))))
    assert "length 1" in str(failed.value)
    assert "width 2" in str(failed.value)
    assert "source.window" in str(failed.value)
    screened = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(parts.LengthRange(1, 6)),
        budget=parts.CandidateBudget(4, batch_size=2),
        screening=(
            parts.PWMExclusion(
                "motif",
                (motif,),
                parts.FimoScoring(
                    executable=tool, limits=parts.ScoringLimits(windows=22)
                ),
            ),
        ),
    )
    assert da.plan(screened).preview["screening_window_bound"] == 44
    pool = da.prepare(screened, out=tmp_path / "screened")
    assert da.inspect(pool, verify=True).source_parts == 4


@pytest.mark.parametrize("strategy", ["stochastic", "consensus", "background"])
def test_pwm_ranges_keep_valid_coordinates_and_recorded_lengths(
    tmp_path: Path, strategy: str
):

    request = parts.PreparationSpec(
        source(tmp_path),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(parts.LengthRange(2, 8), strategy=strategy),
        budget=parts.CandidateBudget(50),
        scoring=parts.FimoScoring(executable=scorer(tmp_path, "")),
        seed=7,
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidates = [r.candidate for r in rows]
    assert {len(c.part.sequence) for c in candidates} == set(range(2, 9))
    if strategy == "consensus":
        for candidate in candidates:
            p = candidate.part.metadata["proposal"]
            assert candidate.part.sequence[p["start"] : p["end"]] == "AC"
    request.source.path.unlink()
    request.scoring.executable.unlink()
    assert da.inspect(pool, verify=True).source_parts == 50


@pytest.mark.parametrize(
    "length",
    [
        {"exact": 4, "minimum": 2, "maximum": 6},
        {"minimum": 2},
        {"maximum": 6},
        {"minimum": 2, "maximum": 6, "distribution": "normal"},
    ],
)
def test_cli_rejects_ambiguous_or_unsupported_length_requests(
    tmp_path: Path, length: dict
):
    source = tmp_path / "request.json"
    source.write_text(
        json.dumps(
            {
                "schema": "dense_arrays.prepare.v1",
                "source": {"kind": "background"},
                "sampling": {"length": length},
                "budget": {"candidates": 3},
                "retain": {"count": 1, "policy": "first_eligible"},
            }
        )
    )
    out = tmp_path / "plan.json"
    result = CliRunner().invoke(app, ["plan", str(source), "--out", str(out), "--json"])
    assert result.exit_code != 0
    assert not out.exists()
    assert json.loads(result.stdout)["schema"] == "dense_arrays.error.v1"


def test_length_entropy_rejection_is_bounded(monkeypatch: pytest.MonkeyPatch):

    class Entropy:
        def digest(self, size: int) -> bytes:
            return b"\xff" * size

    monkeypatch.setattr(mining.hashlib, "shake_256", lambda _: Entropy())
    with pytest.raises(RuntimeError, match="bounded entropy"):
        mining.sample_length(minimum=1, maximum=3, seed=7, index=1)


def test_short_exclusion_candidates_preserve_hit_mapping(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    motif = source(tmp_path)
    tool = scorer(
        tmp_path,
        "".join(
            f"motif\t\tcandidate_{i}\t1\t2\t+\t3\t0.0625\t\tAC\n" for i in range(2)
        ),
    )
    monkeypatch.setattr(
        preparation,
        "sample_sequence",
        lambda **kwargs: {1: "A", 2: "AC", 3: "T", 4: "AC"}[kwargs["index"]],
    )
    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=2, policy="first_eligible"),
        sampling=parts.Sampling(parts.LengthRange(1, 2)),
        budget=parts.CandidateBudget(4),
        screening=(
            parts.PWMExclusion("exclude", (motif,), parts.FimoScoring(executable=tool)),
        ),
    )
    pool = da.prepare(request, out=tmp_path / "pool")
    with da.inspect(pool, view="candidates", all=True).records() as rows:
        candidates = [r.candidate for r in rows]
    assert [c.reasons for c in candidates] == [(), ("exclude",), (), ("exclude",)]
    assert [c.screening[0].hit is None for c in candidates] == [
        True,
        False,
        True,
        False,
    ]
    assert da.inspect(pool, verify=True).retained_parts == 2


def test_all_short_exclusion_candidates_do_not_invoke_a_scan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    request = parts.PreparationSpec(
        parts.Background(base_probabilities=(1, 0, 0, 0)),
        parts.Retention(count=1, policy="first_eligible"),
        sampling=parts.Sampling(parts.LengthRange(1, 1)),
        budget=parts.CandidateBudget(2),
        screening=(
            parts.PWMExclusion(
                "exclude",
                (source(tmp_path),),
                parts.FimoScoring(executable=scorer(tmp_path, "")),
            ),
        ),
    )
    plan = da.plan(request)

    def forbidden(*_: object, **__: object) -> None:
        pytest.fail("short sequences have no full-width exclusion windows")

    monkeypatch.setattr(preparation, "_scan", forbidden)
    pool = da.prepare(plan, out=tmp_path / "pool")
    assert da.inspect(pool, verify=True).retained_parts == 1
