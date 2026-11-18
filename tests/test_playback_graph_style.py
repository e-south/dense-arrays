"""Keep graph arrows thin, unstroked and independent of highlight state.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import json
from dataclasses import replace
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import pytest
from matplotlib.patches import PathPatch
from PIL import Image
from typer.testing import CliRunner

from dense_arrays.playback import (
    GraphArrowStyle,
    PlaybackDocument,
    PlaybackPresentation,
    cli,
    reconstruct_playback,
)
from dense_arrays.playback.cli import app
from dense_arrays.playback.graph_drawing import draw_graph
from dense_arrays.playback.matplotlib_renderer import render_collection_poster_png
from dense_arrays.realized import Placement, PlacementKind, RealizedArray

from .test_playback_cli import _input

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path


def _document() -> PlaybackDocument:
    array = RealizedArray(
        "test",
        "AAACCCGGG",
        tuple(
            Placement(str(i), str(i), PlacementKind.TFBS, bases, i * 3)
            for i, bases in enumerate(("AAA", "CCC", "GGG"))
        ),
    )
    return PlaybackDocument(
        reconstruct_playback(array),
        "test",
        presentation=PlaybackPresentation(graph_detail="reduced"),
    )


def test_graph_arrow_sizes_follow_type_and_never_change_with_highlight():
    seen = []
    for progress in (0, 0.5, 1):
        figure, axis = plt.subplots(figsize=(12, 5))
        try:
            draw_graph(axis, _document(), 0, progress, kmer_font_size_pt=20)
            shafts = [p for p in axis.patches if p.get_gid() == "graph-arrow-shaft"]
            heads = [p for p in axis.patches if p.get_gid() == "graph-arrow-head"]
            assert shafts
            assert heads
            assert all(p.get_linewidth() == pytest.approx(1.4) for p in shafts)
            assert all(
                isinstance(p, PathPatch)
                and p.get_linewidth() == 0
                and p.get_edgecolor()[-1] == 0
                for p in heads
            )
            seen.append([p.get_path().vertices.tolist() for p in heads])
        finally:
            plt.close(figure)
    assert seen[0] == seen[1] == seen[2]


def test_custom_arrow_sizes_and_shared_mapping_contract():
    settings = {
        "graph_detail": "reduced",
        "graph_arrows": {
            "shaft_width_em": 0.04,
            "head_length_em": 0.6,
            "head_width_em": 0.3,
        },
    }
    presentation = PlaybackPresentation.from_mapping(settings)
    assert presentation.graph_arrows == GraphArrowStyle(
        shaft_width_em=0.04, head_length_em=0.6, head_width_em=0.3
    )
    assert PlaybackPresentation.from_mapping(presentation.to_mapping()) == presentation
    figure, axis = plt.subplots(figsize=(12, 5))
    try:
        draw_graph(
            axis,
            replace(_document(), presentation=presentation),
            0,
            1,
            kmer_font_size_pt=20,
        )
        assert all(
            p.get_linewidth() == pytest.approx(0.8)
            for p in axis.patches
            if p.get_gid() == "graph-arrow-shaft"
        )
    finally:
        plt.close(figure)


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf"), True, "0.07"])
def test_invalid_arrow_dimensions_fail_before_drawing(bad: object):
    with pytest.raises((TypeError, ValueError), match="shaft_width_em"):
        GraphArrowStyle(shaft_width_em=bad)


def test_unknown_graph_style_keys_fail_before_output():
    with pytest.raises((TypeError, ValueError), match="unknown"):
        PlaybackPresentation.from_mapping({"graph_arrows": {"shaft_thickness": 2}})


def test_cli_passes_the_same_presentation_contract(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):

    config = tmp_path / "presentation.json"
    config.write_text(json.dumps({"graph_arrows": {"shaft_width_em": 0.04}}))
    seen = []

    def capture(document: PlaybackDocument, paths: Mapping[str, Path]) -> None:
        seen.append(document.presentation)
        paths["poster.png"].write_bytes(b"proof")

    monkeypatch.setattr(cli, "_render_exports", capture)
    result = CliRunner().invoke(
        app,
        [
            str(_input(tmp_path)),
            "--presentation",
            str(config),
            "--poster",
            str(tmp_path / "out.png"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert seen[0].graph_arrows.shaft_width_em == 0.04


def test_bad_cli_style_fails_without_output_directory(tmp_path: Path):
    config = tmp_path / "presentation.json"
    config.write_text('{"graph_arrows":{"shaft_width_em":0}}')
    output = tmp_path / "new" / "poster.png"
    result = CliRunner().invoke(
        app,
        [str(_input(tmp_path)), "--presentation", str(config), "--poster", str(output)],
    )
    assert result.exit_code != 0
    assert "shaft_width_em" in result.output
    assert not output.parent.exists()


def test_cli_cannot_overwrite_its_style_input(tmp_path: Path):
    config = tmp_path / "presentation.json"
    config.write_text("{}")
    result = CliRunner().invoke(
        app,
        [
            str(_input(tmp_path)),
            "--presentation",
            str(config),
            "--poster",
            str(config),
            "--replace",
        ],
    )
    assert result.exit_code != 0
    assert "presentation" in result.output.lower()
    assert config.read_text() == "{}"


def test_public_poster_reports_measured_type_for_downstream_placement(tmp_path: Path):
    target = render_collection_poster_png(
        (_document(),), tmp_path / "poster.png", dpi=180
    )
    with Image.open(target) as image:
        metrics = json.loads(image.info["PlaybackTypography"])
    assert metrics["schema"] == "dense-arrays.playback-typography/v1"
    assert 0 < metrics["nucleotide_cap_height_pt"] < metrics["graph_font_size_pt"]


def test_arrow_opacity_is_explicit_and_consistent_with_the_style():
    assert GraphArrowStyle().opacity == 0.55
    assert GraphArrowStyle.from_mapping({"opacity": 0.8}).opacity == 0.8
    for value in (0, -1, 1.1, True, "0.5"):
        with pytest.raises((TypeError, ValueError), match="opacity"):
            GraphArrowStyle(opacity=value)
