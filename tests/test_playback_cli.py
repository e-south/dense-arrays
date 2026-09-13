"""Exercise playback CLI input diagnostics and output publication contracts.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import builtins
from typing import TYPE_CHECKING

import pytest
from typer.testing import CliRunner

from dense_arrays.playback.cli import app
from dense_arrays.playback.serialization import dumps_realized_array
from dense_arrays.realized import Placement, PlacementKind, RealizedArray

if TYPE_CHECKING:
    from pathlib import Path

runner = CliRunner()


def _input(tmp_path: Path) -> Path:
    realized = RealizedArray(
        "cli-example",
        "CAGCGT",
        (
            Placement("one", "one", PlacementKind.OTHER, "CAG", 0),
            Placement("two", "two", PlacementKind.OTHER, "AGC", 1),
            Placement("three", "three", PlacementKind.OTHER, "CGT", 3),
        ),
    )
    path = tmp_path / "realized.json"
    path.write_text(dumps_realized_array(realized), encoding="utf-8")
    return path


@pytest.mark.parametrize("payload", ["[]", "null", "{", '{"schema_version":"unknown"}'])
def test_invalid_input_reports_error_without_writing(tmp_path: Path, payload: str):
    source = tmp_path / "bad.json"
    source.write_text(payload, encoding="utf-8")
    output = tmp_path / "new" / "poster.png"
    result = runner.invoke(app, [str(source), "--poster", str(output)])
    assert result.exit_code != 0
    assert "Error" in result.output
    assert "Traceback" not in result.output
    assert not output.parent.exists()


def test_missing_input_has_a_readable_error(tmp_path: Path):
    result = runner.invoke(
        app, [str(tmp_path / "missing.json"), "--poster", str(tmp_path / "out.png")]
    )
    assert result.exit_code != 0
    assert "Error" in result.output
    assert "Traceback" not in result.output


def test_input_cannot_be_an_output(tmp_path: Path):
    source = _input(tmp_path)
    original = source.read_bytes()
    result = runner.invoke(app, [str(source), "--poster", str(source)])
    assert result.exit_code != 0
    assert "input" in result.output.lower()
    assert source.read_bytes() == original


def test_outputs_cannot_alias_each_other(tmp_path: Path):
    source = _input(tmp_path)
    output = tmp_path / "artifact"
    result = runner.invoke(
        app, [str(source), "--poster", str(output), "--mp4", str(output)]
    )
    assert result.exit_code != 0
    assert "output" in result.output.lower()
    assert not output.exists()


def test_existing_output_requires_explicit_replacement(tmp_path: Path):
    source = _input(tmp_path)
    output = tmp_path / "out.png"
    output.write_text("keep this", encoding="utf-8")
    result = runner.invoke(app, [str(source), "--poster", str(output)])
    assert result.exit_code != 0
    assert "replace" in result.output.lower()
    assert output.read_text(encoding="utf-8") == "keep this"


def test_success_reports_the_output_path(tmp_path: Path):
    source = _input(tmp_path)
    output = tmp_path / "out.png"
    result = runner.invoke(app, [str(source), "--poster", str(output)])
    assert result.exit_code == 0, result.output
    assert output.name in result.output
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_missing_ffmpeg_leaves_existing_poster_untouched(tmp_path: Path):
    source = _input(tmp_path)
    output = tmp_path / "out.png"
    output.write_text("keep this", encoding="utf-8")
    video = tmp_path / "video.mp4"
    result = runner.invoke(
        app,
        [str(source), "--poster", str(output), "--mp4", str(video), "--replace"],
        env={"PATH": ""},
    )
    assert result.exit_code != 0
    assert "ffmpeg" in result.output.lower()
    assert output.read_text(encoding="utf-8") == "keep this"
    assert not video.exists()


@pytest.mark.parametrize("link_kind", ["hard", "symbolic"])
def test_link_to_input_cannot_be_an_output(tmp_path: Path, link_kind: str):
    source = _input(tmp_path)
    original = source.read_bytes()
    output = tmp_path / "alias.png"
    if link_kind == "hard":
        output.hardlink_to(source)
    else:
        output.symlink_to(source)
    result = runner.invoke(app, [str(source), "--poster", str(output), "--replace"])
    assert result.exit_code == 1
    assert "Error:" in result.stderr
    assert source.read_bytes() == original


def test_explicit_replacement_publishes_poster(tmp_path: Path):
    source = _input(tmp_path)
    poster = tmp_path / "poster.png"
    poster.write_text("previous output", encoding="utf-8")
    result = runner.invoke(
        app,
        [str(source), "--poster", str(poster), "--replace"],
    )
    assert result.exit_code == 0, result.output
    assert poster.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert not list(tmp_path.glob(".dense-arrays-*"))


def test_export_format_is_explicit(tmp_path: Path):
    source = _input(tmp_path)
    result = runner.invoke(app, [str(source)])
    assert result.exit_code == 1
    assert "--poster, --mp4, or --gif" in result.stderr
    assert list(tmp_path.iterdir()) == [source]


def test_missing_playback_dependency_has_an_install_hint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    source = _input(tmp_path)
    original_import = builtins.__import__

    def without_matplotlib(name: str, *args: object, **kwargs: object) -> object:
        if name.partition(".")[0] == "matplotlib":
            message = "No module named matplotlib"
            raise ModuleNotFoundError(message)
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_matplotlib)
    output = tmp_path / "poster.png"
    result = runner.invoke(app, [str(source), "--poster", str(output)])
    assert result.exit_code == 1
    assert "playback dependencies" in result.stderr
    assert "--extra playback" in result.stderr
    assert not output.exists()
