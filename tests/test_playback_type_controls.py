"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/tests/test_playback_type_controls.py

Public typography controls preserve physical size across API and CLI rendering.

Module Author(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import json
from dataclasses import replace
from typing import TYPE_CHECKING, NoReturn

import matplotlib.pyplot as plt
import pytest
from PIL import Image
from typer.testing import CliRunner

from dense_arrays.playback import PlaybackPresentation, PlaybackTypography
from dense_arrays.playback.cli import app
from dense_arrays.playback.matplotlib_renderer import render_collection_poster_png
from dense_arrays.playback.scene_drawing import draw_document, transition_frame_counts

from .test_playback_cli import _input
from .test_playback_typography import compact_document

if TYPE_CHECKING:
    from pathlib import Path


def test_explicit_type_is_shared_by_graph_strands_features_and_annotations():
    document = replace(
        compact_document(),
        label_overrides={"0": "Anchor"},
        presentation=PlaybackPresentation(
            graph_detail="reduced",
            canvas_height_inches=4,
            typography=PlaybackTypography(
                nucleotide_font_size_pt=18, annotation_scale=1
            ),
        ),
    )
    figure = plt.figure(figsize=(16, 4))
    try:
        states = []
        for index, progress in ((0, 0), (4, 1)):
            metrics = draw_document(
                document, transition_index=index, progress=progress, figure=figure
            )
            assert metrics["graph_font_size_pt"] == pytest.approx(18)
            graph = [t for t in figure.axes[0].texts if len(t.get_text()) == 16]
            bases = [t for t in figure.axes[1].texts if t.get_text() in "ACGT"]
            assert graph
            assert bases
            assert all(t.get_fontsize() == pytest.approx(18) for t in (*graph, *bases))
            labels = [
                t
                for a in figure.axes
                for t in a.texts
                if t.get_text() in {"Anchor", "5'", "3'"}
            ]
            assert labels
            assert all(t.get_fontsize() == pytest.approx(18) for t in labels)
            states.append([(t.get_position(), t.get_fontsize()) for t in bases])
        assert states[0] == states[1]
        assert (
            len(transition_frame_counts(document, figure, fps=24, seconds_per_step=0.5))
            == 5
        )
    finally:
        plt.close(figure)


@pytest.mark.parametrize("value", [0, -1, True, "22", float("nan"), float("inf")])
@pytest.mark.parametrize("field", ["nucleotide_font_size_pt", "annotation_scale"])
def test_invalid_typography_is_rejected_at_admission(field: str, value: object):
    with pytest.raises((TypeError, ValueError), match=field):
        PlaybackTypography(**{field: value})


def test_explicit_type_that_cannot_fit_preserves_existing_output(tmp_path: Path):
    document = replace(
        compact_document(),
        presentation=PlaybackPresentation(
            typography=PlaybackTypography(nucleotide_font_size_pt=80),
        ),
    )
    output = tmp_path / "poster.png"
    output.write_bytes(b"previous output")
    with pytest.raises(ValueError, match=r"nucleotide_font_size_pt.*fit"):
        render_collection_poster_png((document,), output)
    assert output.read_bytes() == b"previous output"


def test_cli_and_api_use_identical_type_and_canvas_settings(tmp_path: Path):
    raw = {
        "graph_detail": "reduced",
        "canvas_height_inches": 4,
        "typography": {"nucleotide_font_size_pt": 22, "annotation_scale": 1},
    }
    presentation = PlaybackPresentation.from_mapping(raw)
    assert PlaybackPresentation.from_mapping(presentation.to_mapping()) == presentation
    config = tmp_path / "presentation.json"
    config.write_text(json.dumps(raw))
    output = tmp_path / "poster.png"
    result = CliRunner().invoke(
        app,
        [str(_input(tmp_path)), "--presentation", str(config), "--poster", str(output)],
    )
    assert result.exit_code == 0, result.output
    with Image.open(output) as image:
        assert image.width / image.height == pytest.approx(4)
        assert json.loads(image.info["PlaybackTypography"])[
            "graph_font_size_pt"
        ] == pytest.approx(22)


def test_unknown_typography_fields_and_wrong_types_fail_before_rendering():
    for raw in ({"size": 22}, {"nucleotide_font_size_pt": "22"}, []):
        with pytest.raises((TypeError, ValueError)):
            PlaybackPresentation.from_mapping({"typography": raw})
    with pytest.raises(TypeError, match="typography"):
        PlaybackPresentation(typography={})


@pytest.mark.parametrize("value", [0, -1, True, "4", float("nan"), float("inf")])
def test_invalid_canvas_height_fails_at_admission(value: object):
    with pytest.raises((TypeError, ValueError), match="canvas_height_inches"):
        PlaybackPresentation(canvas_height_inches=value)


def test_collection_rejects_changing_canvas_before_output(tmp_path: Path):
    document = compact_document()
    taller = replace(
        document, presentation=PlaybackPresentation(canvas_height_inches=4)
    )
    output = tmp_path / "absent" / "poster.png"
    with pytest.raises(ValueError, match="same canvas_height_inches"):
        render_collection_poster_png((document, taller), output)
    assert not output.parent.exists()


def test_explicit_native_size_cannot_silently_override_a_producer(tmp_path: Path):
    document = replace(
        compact_document(),
        presentation=PlaybackPresentation(
            typography=PlaybackTypography(nucleotide_font_size_pt=18),
        ),
    )

    def never_called(_document: object, _index: object) -> NoReturn:
        msg = "callback must not run before contract admission"
        raise AssertionError(msg)

    output = tmp_path / "absent" / "poster.png"
    with pytest.raises(ValueError, match="producer callbacks own their typography"):
        render_collection_poster_png(
            (document,), output, duplex_frame_renderer=never_called
        )
    assert not output.parent.exists()


def test_larger_annotation_ratio_reserves_space_between_captioned_tracks():
    document = replace(
        compact_document(),
        label_overrides={"0": "Anchor", "1": "Reference"},
        presentation=PlaybackPresentation(
            graph_detail="none",
            graph_fraction=0,
            canvas_height_inches=5,
            typography=PlaybackTypography(
                nucleotide_font_size_pt=18, annotation_scale=2
            ),
        ),
    )
    figure = plt.figure(figsize=(16, 5))
    try:
        draw_document(document, transition_index=4, progress=1, figure=figure)
        figure.canvas.draw()
        axis = figure.axes[0]
        renderer = figure.canvas.get_renderer()
        labels = [
            t.get_window_extent(renderer)
            for t in axis.texts
            if t.get_text() in {"Anchor", "Reference"}
        ]
        bars = [p.get_window_extent(renderer) for p in axis.patches]
        assert len(labels) == 2
        assert all(not label.overlaps(bar) for label in labels for bar in bars)
    finally:
        plt.close(figure)
