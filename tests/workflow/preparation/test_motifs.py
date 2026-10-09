"""Neutral motif inputs preserve sampling and scoring meaning independently.

Author: Eric J. South.
"""

import hashlib
import json
import random
from dataclasses import replace
from pathlib import Path

import pytest

from dense_arrays import parts
from dense_arrays.parts.motifs import Motif, best_hit, read_artifact, score_core


def artifact() -> dict:
    return {
        "schema_version": "1.0",
        "producer": "synthetic-fixture",
        "motif_id": "example",
        "alphabet": "ACGT",
        "matrix_semantics": "probabilities",
        "background": dict.fromkeys("ACGT", 0.25),
        "probabilities": [
            dict(zip("ACGT", row, strict=True))
            for row in [(0.7, 0.1, 0.1, 0.1), (0.1, 0.7, 0.1, 0.1)]
        ],
        "log_odds": [
            dict(zip("ACGT", row, strict=True))
            for row in [(2, -1, -1, -1), (-1, 3, -1, -1)]
        ],
        "length": 2,
    }


def source(tmp_path: Path):
    path = tmp_path / "motif.json"
    path.write_text(json.dumps(artifact()))
    return parts.PWMArtifact(path, motif_ids=("example",))


def test_read_artifact_binds_content_and_detaches_mutable_rows(tmp_path: Path):

    declared = source(tmp_path)
    loaded = read_artifact(declared)
    assert loaded.sha256 == hashlib.sha256(declared.path.read_bytes()).hexdigest()
    assert loaded.motif.motif_id == "example"
    assert loaded.motif.width == 2
    assert loaded.motif.background == (0.25, 0.25, 0.25, 0.25)
    assert loaded.motif.log_odds == ((2.0, -1.0, -1.0, -1.0), (-1.0, 3.0, -1.0, -1.0))
    saved = loaded.motif.to_dict()
    restored = type(loaded.motif).from_dict(saved)
    assert restored == loaded.motif
    saved["probabilities"][0][0] = 0
    assert loaded.motif.probabilities[0][0] == 0.7
    loaded.verify()
    declared.path.write_text("{}")
    with pytest.raises(ValueError, match="changed"):
        loaded.verify()


@pytest.mark.parametrize(
    "field,value",
    [
        ("schema_version", "2.0"),
        ("alphabet", "ACGU"),
        ("matrix_semantics", "counts"),
        ("length", True),
        ("length", 2.5),
        ("length", 3),
        ("producer", ""),
        ("motif_id", 4),
    ],
)
def test_artifact_rejects_ambiguous_or_unsupported_declarations(
    tmp_path: Path, field: str, value: object
):

    declared = source(tmp_path)
    data = artifact()
    data[field] = value
    declared.path.write_text(json.dumps(data))
    with pytest.raises((ValueError, TypeError)):
        read_artifact(declared)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, "0.25", -1])
def test_non_probability_values_are_rejected(tmp_path: Path, value: object):

    declared = source(tmp_path)
    data = artifact()
    data["probabilities"][0]["A"] = value
    declared.path.write_text(json.dumps(data))
    with pytest.raises((ValueError, TypeError)):
        read_artifact(declared)


def test_missing_or_extra_bases_and_unknown_motif_selection_fail(tmp_path: Path):

    declared = source(tmp_path)
    with pytest.raises(ValueError, match="motif"):
        read_artifact(parts.PWMArtifact(declared.path, motif_ids=("missing",)))
    data = artifact()
    data["background"]["N"] = 0.1
    declared.path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="ACGT"):
        read_artifact(declared)


def test_exact_core_scoring_and_best_hit_use_declared_weights(tmp_path: Path):

    motif = read_artifact(source(tmp_path)).motif
    score = score_core(motif, "AC")
    assert score.raw == 5
    assert score.per_base == 2.5
    assert score.theoretical_max == 5
    assert score.fraction_of_max == 1
    assert score.units == "declared_log_odds"
    assert score.model_id == motif.model_id
    hit = best_hit(motif, "TTGTTT", strands="double")
    assert (hit.start, hit.end, hit.strand, hit.core) == (2, 4, "reverse", "AC")
    assert hit.score.raw == 5
    tied = best_hit(motif, "ACAC", strands="single")
    assert (tied.start, tied.end, tied.strand) == (0, 2, "forward")
    for invalid in ("A", "ACA", "AN", "ac"):
        with pytest.raises(ValueError, match=r"length|uppercase"):
            score_core(motif, invalid)
    with pytest.raises(ValueError, match="shorter"):
        best_hit(motif, "A")


def test_fraction_of_nonpositive_maximum_is_unavailable(tmp_path: Path):

    declared = source(tmp_path)
    data = artifact()
    data["log_odds"] = [dict.fromkeys("ACGT", 0) for _ in range(2)]
    declared.path.write_text(json.dumps(data))
    assert score_core(read_artifact(declared).motif, "AC").fraction_of_max is None


def test_normalization_is_stable_across_serialization(tmp_path: Path):

    declared = source(tmp_path)
    data = artifact()
    data["probabilities"][0] = {"A": 0.70007, "C": 0.10001, "G": 0.10001, "T": 0.10001}
    declared.path.write_text(json.dumps(data))
    motif = read_artifact(declared).motif
    for _ in range(5):
        assert type(motif).from_dict(motif.to_dict()) == motif
    assert sum(motif.probabilities[0]) == pytest.approx(1)


def test_model_identity_excludes_labels_and_includes_background(tmp_path: Path):

    motif = read_artifact(source(tmp_path)).motif
    assert replace(motif, motif_id="alias", producer="other").model_id == motif.model_id
    assert replace(motif, background=(0.1, 0.2, 0.3, 0.4)).model_id != motif.model_id


def test_score_and_hit_records_validate_denominators_and_geometry(tmp_path: Path):

    motif = read_artifact(source(tmp_path)).motif
    score = score_core(motif, "AC")
    hit = best_hit(motif, "TTGTTT")
    assert type(score).from_dict(score.to_dict()) == score
    assert type(hit).from_dict(hit.to_dict()) == hit
    for changes in (
        {"width": 0},
        {"raw": float("nan")},
        {"raw": 6},
        {"model_id": "bad"},
        {"units": "pvalue"},
    ):
        with pytest.raises((ValueError, TypeError)):
            replace(score, **changes)
    for changes in ({"start": -1}, {"end": 20}, {"strand": "unknown"}, {"core": "A"}):
        with pytest.raises((ValueError, TypeError)):
            replace(hit, **changes)


def test_score_numeric_overflow_is_explicit(tmp_path: Path):

    motif = read_artifact(source(tmp_path)).motif
    overflow = replace(motif, log_odds=((1e308, 0, 0, 0), (1e308, 0, 0, 0)))
    with pytest.raises(ValueError, match="numeric range"):
        score_core(overflow, "AA")


def test_serialized_model_normalization_is_idempotent():

    rng = random.Random(13)  # noqa: S311 - deterministic numeric stability fixture
    for index in range(200):
        row = [rng.random() for _ in range(4)]
        total = sum(row)
        values = tuple(v / total for v in row)
        motif = Motif(str(index), (values,), (0.25,) * 4, ((1, -1, 0, 0),), "fixture")
        assert Motif.from_dict(motif.to_dict()) == motif


def test_duplicate_json_keys_are_rejected(tmp_path: Path):

    declared = source(tmp_path)
    payload = declared.path.read_text().replace(
        '"motif_id": "example"', '"motif_id": "ignored", "motif_id": "example"'
    )
    declared.path.write_text(payload)
    with pytest.raises(ValueError, match="duplicate"):
        read_artifact(declared)
