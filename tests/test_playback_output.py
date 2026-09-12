"""Verify export staging and concurrent destination failures.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from dense_arrays.playback.output import publish_exports

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path


def test_render_failure_preserves_existing_export_set(tmp_path: Path):
    source = tmp_path / "input.json"
    source.write_text("{}", encoding="utf-8")
    outputs = {name: tmp_path / name for name in ("first.png", "second.mp4")}
    for path in outputs.values():
        path.write_text("previous output", encoding="utf-8")

    def fail_during_second_export(paths: Mapping[str, Path]) -> None:
        paths["first.png"].write_text("new poster", encoding="utf-8")
        message = "second renderer failed"
        raise RuntimeError(message)

    with pytest.raises(RuntimeError, match="second renderer failed"):
        publish_exports(source, outputs, fail_during_second_export, replace=True)

    assert all(
        path.read_text(encoding="utf-8") == "previous output"
        for path in outputs.values()
    )
    assert not list(tmp_path.glob(".dense-arrays-*"))


def test_publication_reports_partial_set_without_overwriting_new_file(tmp_path: Path):
    source = tmp_path / "input.json"
    source.write_text("{}", encoding="utf-8")
    outputs = {name: tmp_path / name for name in ("first.png", "second.mp4")}

    def render_with_concurrent_destination(paths: Mapping[str, Path]) -> None:
        for path in paths.values():
            path.write_text("rendered", encoding="utf-8")
        outputs["second.mp4"].write_text("concurrent file", encoding="utf-8")

    with pytest.raises(OSError, match="published files:") as caught:
        publish_exports(
            source, outputs, render_with_concurrent_destination, replace=False
        )

    assert str(outputs["first.png"]) in str(caught.value).split("published files:")[1]
    assert outputs["first.png"].read_text(encoding="utf-8") == "rendered"
    assert outputs["second.mp4"].read_text(encoding="utf-8") == "concurrent file"
    assert not list(tmp_path.glob(".dense-arrays-*"))
