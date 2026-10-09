"""Preparation figures preserve recorded populations and portable report meaning.

Author: Eric J. South.
"""

import importlib
import json
import shutil
from pathlib import Path

import pytest
from PIL import Image
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, reporting
from dense_arrays._record_validation import semantic_digest
from dense_arrays.artifacts.preparation.records import recount
from dense_arrays.cli import app
from dense_arrays.parts.motifs import Motif
from dense_arrays.parts.retention.selection import select_candidates
from dense_arrays.playback.quality.preparation import preparation_figure

from .test_mmr import candidate
from .test_sampled import background_request


def test_preparation_quality_renders_live_and_detached_reports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    report = da.inspect(pool, view="quality")
    expected = report.to_dict()
    before = (pool.path / "pool.sqlite3").read_bytes()

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("rendering invoked preparation or scoring")

    monkeypatch.setattr("dense_arrays.workflow.preparation.sample_sequence", forbidden)
    monkeypatch.setattr("dense_arrays.parts.scoring.scan_fimo", forbidden)
    out = tmp_path / "live.png"
    receipt = da.render(pool, view="preparation-quality", out=out)
    assert receipt.view == "preparation-quality"
    assert receipt.records == expected["counts"]["retained"]
    assert receipt.design_refs == ()
    with Image.open(out) as image:
        assert json.loads(image.info["DenseArraysReport"]) == expected
        assert image.info["ReportDigest"] == receipt.sources[0]["report_sha256"]
    assert (pool.path / "pool.sqlite3").read_bytes() == before
    saved = tmp_path / "quality.json"
    da.export(report, out=saved)
    shutil.rmtree(pool.path)
    result = CliRunner().invoke(
        app,
        [
            "render",
            str(saved),
            "--view",
            "preparation-quality",
            "--out",
            str(tmp_path / "detached.png"),
            "--max-read-records",
            "1",
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    assert json.loads(result.stdout)["sources"] == receipt.to_dict()["sources"]
    assert "Read cost:" in result.stderr
    with Image.open(tmp_path / "detached.png") as image:
        assert json.loads(image.info["DenseArraysReport"]) == expected


def test_preparation_render_limits_and_filters_fail_before_publication(tmp_path: Path):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    out = tmp_path / "not-created" / "quality.png"
    with pytest.raises(reporting.ReadLimitError, match="records"):
        da.render(
            pool,
            view="preparation-quality",
            out=out,
            read_limits=reporting.ReadLimits(records=1),
        )
    assert not out.parent.exists()
    with pytest.raises(ValueError, match="filter"):
        da.render(
            pool,
            view="preparation-quality",
            out=out,
            select=reporting.DesignFilter(cells=("x",)),
        )
    assert not out.parent.exists()


def test_preparation_report_requires_its_named_render_view(tmp_path: Path):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    report = da.inspect(pool, view="quality")
    for view in ("design", "library-quality"):
        with pytest.raises((ValueError, TypeError), match="preparation-quality"):
            da.render(report, view=view, out=tmp_path / "wrong.png")
    assert not (tmp_path / "wrong.png").exists()


def test_preparation_set_keeps_recipe_rows_and_unavailable_metrics(tmp_path: Path):
    recipe = background_request()
    pool = da.prepare(
        parts.PreparationSet(
            {"one": recipe, "two": recipe}, sequence_collisions="preserve"
        ),
        out=tmp_path / "set",
    )
    figure = preparation_figure(da.inspect(pool, view="quality").to_dict())
    axes = {axis.get_label(): axis for axis in figure.axes}
    assert len(axes) == 6
    for index in range(2):
        assert list(axes[f"yield_{index}"].lines[0].get_xdata()) == [6, 6, 1, 1]
        assert (
            "not recorded"
            in " ".join(t.get_text() for t in axes[f"diversity_{index}"].texts).lower()
        )
        assert (
            "not declared"
            in " ".join(t.get_text() for t in axes[f"bands_{index}"].texts).lower()
        )
    figure.clear()


@pytest.mark.parametrize(
    "fractions,labels",
    [
        ((0.25, 0.5), ["0% to 25%", "25% to 50%", "50% to 100%"]),
        ((0.001, 0.002), ["0% to 0.1%", "0.1% to 0.2%", "0.2% to 100%"]),
    ],
)
def test_recorded_distances_and_tied_score_bands_use_their_exact_populations(
    fractions: tuple[float, ...], labels: list[str]
):
    motif = Motif("flat", ((0.25,) * 4,) * 3, (0.25,) * 4, None, "fixture")
    bands = parts.ScoreBands(fractions)
    original = tuple(
        candidate(i, seq, score)
        for i, (seq, score) in enumerate(
            (("AAA", 9), ("CCC", 9), ("GGG", 5), ("TTT", 1)), 1
        )
    )
    chosen = select_candidates(
        original,
        parts.Uniqueness("core"),
        parts.Retention(
            count=3,
            policy="mmr",
            rank_by="best_hit_score",
            mmr=parts.MMR(4, 0.5, "fraction_of_max_clipped"),
        ),
        motif=motif,
        score_bands=bands,
    )
    scoring_id = semantic_digest({"scorer": "fixture"})
    accounting = recount(
        chosen,
        target=3,
        budget=4,
        stop_reason="candidate_budget",
        mmr=True,
        score_bands=bands,
        scoring_id=scoring_id,
    )
    value = {
        **accounting.to_dict(),
        "schema": "dense_arrays.pool_quality.v1",
        "pool_id": semantic_digest({"pool": "fixture"}),
        "plan_id": semantic_digest({"plan": "fixture"}),
        "state": accounting.state,
        "diversity": [
            {
                "schema": "dense_arrays.mmr_diversity.v1",
                "recipe_id": None,
                "policy": "greedy_mmr.v1",
                "distance": "pwm_tolerant_hamming",
                "model_id": motif.model_id,
                "scoring_id": scoring_id,
                "status": "recorded",
                "choices": [
                    {"rank": 1, "nearest_distance": None},
                    {"rank": 2, "nearest_distance": 3.0},
                    {"rank": 3, "nearest_distance": 3.0},
                ],
            }
        ],
    }
    figure = preparation_figure(value)
    axes = {axis.get_label(): axis for axis in figure.axes}
    assert list(axes["diversity_0"].lines[0].get_xdata()) == [2, 3]
    assert list(axes["diversity_0"].lines[0].get_ydata()) == [3.0, 3.0]
    assert [bar.get_width() for bar in axes["bands_0"].patches] == [2, 0, 2, 2, 0, 1]
    assert [text.get_text() for text in axes["bands_0"].get_yticklabels()] == labels
    figure.clear()


def test_preparation_optional_renderer_fails_before_report_scan(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    original = importlib.import_module

    def missing(name: str) -> object:
        if name == "matplotlib.figure":
            raise ModuleNotFoundError(name)
        return original(name)

    monkeypatch.setattr("dense_arrays.reporting.rendering.import_module", missing)
    out = tmp_path / "absent" / "report.png"
    with pytest.raises(ValueError, match="playback dependencies"):
        da.render(
            pool,
            view="preparation-quality",
            out=out,
            read_limits=reporting.ReadLimits(records=1),
        )
    assert not out.parent.exists()


def test_preparation_render_rejects_existing_or_unsupported_destinations(
    tmp_path: Path,
):
    pool = da.prepare(background_request(), out=tmp_path / "pool")
    out = tmp_path / "report.png"
    out.write_bytes(b"existing")
    with pytest.raises(FileExistsError):
        da.render(pool, view="preparation-quality", out=out)
    assert out.read_bytes() == b"existing"
    with pytest.raises(ValueError, match="png"):
        da.render(pool, view="preparation-quality", out=tmp_path / "report.pdf")
    assert not (tmp_path / "report.pdf").exists()
