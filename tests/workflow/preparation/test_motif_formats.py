"""Motif input adapters preserve probabilities and explicit unavailable scores.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app
from dense_arrays.parts.motifs import Motif, best_hit, read_artifact, score_core
from dense_arrays.parts.motifs.windows import select_window
from dense_arrays.workflow.inputs import read_source

from .test_fimo import scorer

MEME = """MEME version 5
ALPHABET= ACGT
strands: + -
Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF first shared_name
letter-probability matrix: alength= 4 w= 2 nsites= 3 E= 0.01
0.7 0.1 0.1 0.1
0.1 0.7 0.1 0.1

MOTIF second shared_name
letter-probability matrix: alength= 4 w= 1
0 0 0 1
"""


def test_probability_only_model_roundtrips_and_rejects_undeclared_score():
    motif = Motif("test", ((1, 0, 0, 0),), (0.25,) * 4, None, "test")
    assert motif.log_odds is None
    assert motif.to_dict()["schema"] == "dense_arrays.motif.v2"
    assert motif.to_dict()["score_units"] is None
    assert Motif.from_dict(motif.to_dict()) == motif
    assert select_window(motif, parts.MotifWindow(1)).motif == motif
    for operation in (score_core, best_hit):
        with pytest.raises(ValueError, match="score matrix"):
            operation(motif, "A")


def test_meme_selection_keeps_probabilities_without_count_rounding(tmp_path: Path):
    path = tmp_path / "motifs.meme"
    path.write_text(MEME)
    declared = parts.PWMArtifact(path, motif_ids=("first",), format="meme")
    loaded = read_artifact(declared)
    assert loaded.motif.motif_id == "first"
    assert loaded.motif.probabilities == ((0.7, 0.1, 0.1, 0.1), (0.1, 0.7, 0.1, 0.1))
    assert loaded.motif.log_odds is None
    assert loaded.motif.metadata["nsites"] == 3
    assert loaded.motif.metadata["alternate_name"] == "shared_name"
    assert loaded.format == "meme"
    assert loaded.motif.background == (0.25,) * 4
    loaded.verify()
    with pytest.raises(ValueError, match="select one"):
        read_artifact(parts.PWMArtifact(path, format="meme"))
    with pytest.raises(ValueError, match="missing"):
        read_artifact(
            parts.PWMArtifact(path, motif_ids=("shared_name",), format="meme")
        )


@pytest.mark.parametrize("input_format", ["meme", "jaspar"])
def test_imported_plan_and_pool_share_python_cli_and_saved_models(
    tmp_path: Path, input_format: str
):
    path = tmp_path / "motifs.meme"
    path.write_text(MEME if input_format == "meme" else ">first\n7 1\n1 7\n1 1\n1 1\n")
    tool = scorer(tmp_path, "")
    request = parts.PreparationSpec(
        parts.PWMArtifact(path, motif_ids=("first",), format=input_format),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=6), strategy="consensus"),
        budget=parts.CandidateBudget(8),
        scoring=parts.FimoScoring(executable=tool),
        seed=7,
    )
    plan = da.plan(request)
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).to_dict() == plan.to_dict()
    assert plan.preview["motif_import"]["format"] == input_format
    preview = CliRunner().invoke(app, ["plan", str(saved)])
    assert preview.exit_code == 0, preview.output
    assert f"Motif input: {input_format}; supplied score matrix: unavailable." in (
        preview.stdout
    )
    pool = da.prepare(plan, out=tmp_path / "python")
    result = CliRunner().invoke(
        app, ["prepare", str(saved), "--out", str(tmp_path / "cli"), "--json"]
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["pool_id"] == pool.pool_id
    path.unlink()
    tool.unlink()
    assert da.inspect(pool, verify=True).source_parts == 8


@pytest.mark.parametrize(
    "old,new",
    [
        ("MOTIF second", "MOTIF first"),
        ("ALPHABET= ACGT", "ALPHABET= ACGU"),
        ("nsites= 3", "nsites= 0"),
        ("w= 2", "w= 3"),
        ("w= 2", "w= 2 w= 2"),
        ("w= 2", "width= 2"),
        ("0.7 0.1 0.1 0.1", "0.7 0.1 0.1 0.2"),
        ("0 0 0 1", "0 0 0 nan"),
        ("A 0.25 C 0.25 G 0.25 T 0.25", "A 0.25 A 0.25 G 0.25 T 0.25"),
        ("strands: + -", "strands: -"),
    ],
)
def test_meme_rejects_ambiguous_or_malformed_records_including_unselected(
    tmp_path: Path, old: str, new: str
):
    path = tmp_path / "invalid.meme"
    path.write_text(MEME.replace(old, new))
    with pytest.raises(ValueError, match=r"MEME|motif"):
        read_artifact(parts.PWMArtifact(path, motif_ids=("first",), format="meme"))


def test_minimal_meme_defaults_are_visible_without_fabricating_counts(tmp_path: Path):
    path = tmp_path / "minimal.meme"
    path.write_text(
        "MEME version 5\nMOTIF simple\nletter-probability matrix:\n0 1 0 0\n"
    )
    motif = read_artifact(parts.PWMArtifact(path, format="meme")).motif
    assert motif.probabilities == ((0, 1, 0, 0),)
    assert motif.background == (0.25,) * 4
    assert motif.metadata["background_origin"] == "uniform_default"
    assert "nsites" not in motif.metadata


@pytest.mark.parametrize(
    "rows",
    [
        "A [ 7 1 ]\nC [ 1 7 ]\nG [ 1 1 ]\nT [ 1 1 ]\n",
        "7 1\n1 7\n1 1\n1 1\n",
    ],
)
def test_jaspar_preserves_source_counts_and_normalizes_each_position(
    tmp_path: Path, rows: str
):
    path = tmp_path / "sites.jaspar"
    path.write_text(">test.1 shared name\n" + rows)
    motif = read_artifact(parts.PWMArtifact(path, format="jaspar")).motif
    assert motif.motif_id == "test.1"
    assert motif.probabilities == ((0.7, 0.1, 0.1, 0.1), (0.1, 0.7, 0.1, 0.1))
    assert motif.log_odds is None
    assert motif.metadata["source_counts"] == ((7, 1, 1, 1), (1, 7, 1, 1))
    assert motif.metadata["column_totals"] == (10, 10)
    assert motif.metadata["name"] == "shared name"


@pytest.mark.parametrize(
    "rows",
    [
        "A [ 1 ]\nA [ 1 ]\nG [ 1 ]\nT [ 1 ]",
        "A [ 1 ]\nC [ 1 ]\nG [ 1 ]",
        "A [ 1 2 ]\nC [ 1 ]\nG [ 1 ]\nT [ 1 ]",
        "0\n0\n0\n0",
        "1\n-1\n1\n1",
        "inf\n1\n1\n1",
        "A [ 1 ]\n1\n1\n1",
        "1\n1\n1\n1\n1",
    ],
)
def test_jaspar_rejects_invalid_or_ambiguous_matrices(tmp_path: Path, rows: str):
    path = tmp_path / "invalid.jaspar"
    path.write_text(">test\n" + rows + "\n")
    with pytest.raises(ValueError, match="JASPAR"):
        read_artifact(parts.PWMArtifact(path, format="jaspar"))


def test_meme_whitespace_and_windowed_exclusion_roundtrip(tmp_path: Path):
    path = tmp_path / "motifs.meme"
    path.write_text(
        MEME.replace("MOTIF ", "MOTIF\t").replace("ALPHABET=", "ALPHABET =")
    )
    motif = parts.PWMArtifact(
        path, motif_ids=("first",), window=parts.MotifWindow(2), format="meme"
    )
    request = parts.PreparationSpec(
        parts.Background(),
        parts.Retention(count=0, policy="first_eligible"),
        sampling=parts.Sampling(planning.Length(exact=6)),
        budget=parts.CandidateBudget(3),
        screening=(
            parts.PWMExclusion(
                "exclude", (motif,), parts.FimoScoring(executable=scorer(tmp_path, ""))
            ),
        ),
    )
    plan = da.plan(request)
    saved = tmp_path / "plan.json"
    plan.write(saved)
    assert read_source(saved).to_dict() == plan.to_dict()
    pool = da.prepare(plan, out=tmp_path / "pool")
    assert da.inspect(pool, verify=True).source_parts == 3
