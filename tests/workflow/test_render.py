"""Native rendering consumes persisted identities without solving or screening again.

Author: Eric J. South.
"""

import builtins
import importlib
import json
from pathlib import Path

import pytest
from PIL import Image
from typer.testing import CliRunner

import dense_arrays as da
from dense_arrays import parts, planning
from dense_arrays.cli import app


def request(*, spacing: int = 0) -> planning.DesignSpec:
    """Keep geometry known so rendered constraint evidence can be checked."""
    return planning.DesignSpec(
        parts=[parts.Part("a", "AACC"), parts.Part("b", "CCGT")],
        length=planning.Length(maximum=8),
        strands="single",
        requirements=[
            planning.Fixed("a", "a", "forward"),
            planning.Fixed("b", "b", "forward"),
            planning.Spacing("gap", "a", "b", min=spacing, max=spacing),
        ]
        if spacing < 0
        else [],
    )


def test_render_python_cli_receipts_share_saved_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    run = da.run(request(), out=tmp_path / "run")
    before = (run.path / "run.sqlite3").read_bytes()

    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("render invoked generation or acceptance")

    monkeypatch.setattr(da.Optimizer, "solve_report", forbidden)
    monkeypatch.setattr("dense_arrays.generation.acceptance.evaluate", forbidden)
    python_out = tmp_path / "python.png"
    receipt = da.render(run, out=python_out)
    assert receipt.records == 1
    with Image.open(python_out) as image:
        assert image.width > 100
        assert run.run_id in image.info["Description"]
    cli_out = tmp_path / "cli.png"
    response = CliRunner().invoke(
        app, ["render", str(run.path), "--out", str(cli_out), "--json"]
    )
    assert response.exit_code == 0, response.output
    data = json.loads(response.stdout)
    assert data["sources"] == receipt.to_dict()["sources"]
    assert data["design_refs"] == receipt.to_dict()["design_refs"]
    assert (run.path / "run.sqlite3").read_bytes() == before


def test_negative_spacing_fails_before_render_publication(tmp_path: Path):
    run = da.run(request(spacing=-2), out=tmp_path / "run")
    out = tmp_path / "unwritten" / "array.png"
    with pytest.raises(ValueError, match="negative spacing"):
        da.render(run, out=out)
    assert not out.parent.exists()


def test_render_collision_does_not_touch_existing_output(tmp_path: Path):
    out = tmp_path / "keep.png"
    out.write_bytes(b"unchanged")
    with pytest.raises(FileExistsError):
        da.render(tmp_path / "missing-run", out=out)
    assert out.read_bytes() == b"unchanged"


@pytest.mark.parametrize(
    "missing", ["dense_arrays.playback.matplotlib_renderer", "matplotlib.pyplot"]
)
def test_missing_render_dependencies_fail_before_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, missing: str
):
    run = da.run(request(), out=tmp_path / "run")
    original = builtins.__import__

    def unavailable(name: str, *args: object, **kwargs: object) -> object:
        if name == missing:
            msg = "optional renderer unavailable"
            raise ModuleNotFoundError(msg)
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", unavailable)
    import_module = importlib.import_module

    def optional_module(name: str) -> object:
        if name == missing:
            msg = "optional renderer unavailable"
            raise ModuleNotFoundError(msg)
        return import_module(name)

    monkeypatch.setattr(
        "dense_arrays.reporting.rendering.import_module", optional_module
    )
    out = tmp_path / "uncreated" / "array.png"
    with pytest.raises(ValueError, match="playback dependencies"):
        da.render(run, out=out)
    assert not out.parent.exists()


def test_render_explains_saved_checks_without_raw_record_syntax(tmp_path: Path):
    request = planning.DesignSpec(
        parts=[
            parts.Part("u", "AACC"),
            parts.Part("bridge", "CCGT"),
            parts.Part("d", "GTTA"),
        ],
        length=planning.Length(exact=14),
        strands="single",
        assembly=planning.Assembly(
            padding=planning.Padding(side="right", max_trials=20)
        ),
        requirements=[
            planning.Fixed("u", "u", "forward", planning.StartWindow(max=0)),
            planning.Fixed("d", "d", "forward"),
            planning.Spacing("adjacent", "u", "d", min=0, max=0),
            planning.GC("final-gc", scope="sequence", min=0.2, max=0.8),
        ],
    )
    run = da.run(request, out=tmp_path / "run")
    out = tmp_path / "array.png"
    da.render(run, out=out)
    with Image.open(out) as image:
        evidence = image.info["PlaybackEvidence"]
    assert "final-gc" in evidence
    assert "GC" in evidence
    assert '"passed":' not in evidence
    assert "adjacent" in evidence
