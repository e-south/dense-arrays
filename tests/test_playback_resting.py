"""Complete neutral scene context and stationary progressive emphasis.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.colors import to_rgb
from matplotlib.patches import FancyBboxPatch
from PIL import Image

from dense_arrays.playback import reconstruct_playback
from dense_arrays.playback.duplex_frames import DuplexFrames, duplex_transition_frame
from dense_arrays.playback.frame_schedule import PlaybackTiming, scene_frame_schedule
from dense_arrays.playback.matplotlib_renderer import render_collection_gif
from dense_arrays.playback.presentation import PlaybackDocument
from dense_arrays.playback.scene_drawing import draw_document
from dense_arrays.playback.theme import LegendEntry, PlaybackPresentation
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray

if TYPE_CHECKING:
    from pathlib import Path


def longer_document() -> PlaybackDocument:
    sequence = "ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGT"
    realized = RealizedArray(
        source_id="test:resting",
        sequence=sequence,
        placements=tuple(
            Placement(
                str(index),
                str(index),
                PlacementKind.TFBS,
                sequence[start : start + 16],
                start,
                orientation=Orientation.REVERSE if index == 1 else Orientation.FORWARD,
            )
            for index, start in enumerate((0, 7, 14, 21))
        ),
    )
    return PlaybackDocument(
        reconstruct_playback(realized),
        title="Four overlapping motifs",
        presentation=PlaybackPresentation(
            legend_entries=(LegendEntry("selected", "Selected motif", "#365E80"),),
        ),
    )


def test_zero_lead_still_has_a_complete_resting_frame() -> None:
    timing = PlaybackTiming(fps=10, lead_seconds=0, hold_seconds=0)
    frames = tuple(scene_frame_schedule((3, 5), timing, first=True))
    assert frames[0].transition_index == 0
    assert frames[0].progress == 0


def test_scene_boundaries_allow_orientation_with_complete_resting_context() -> None:
    timing = PlaybackTiming(fps=10, lead_seconds=0, scene_transition_seconds=0.3)
    frames = tuple(scene_frame_schedule((3, 5), timing, first=False))
    assert all(frame.progress == 0 for frame in frames[:4])
    assert frames[4].progress > 0


def test_native_rest_preserves_full_scene_geometry_and_neutral_artists() -> None:
    document = longer_document()
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(document, transition_index=0, progress=0, figure=figure)
        resting_text = [
            (text.get_text(), text.get_position(), text.get_fontsize())
            for axis in figure.axes
            for text in axis.texts
        ]
        graph, duplex, _legend = figure.axes
        assert (
            len(
                [
                    text
                    for text in duplex.texts
                    if text.get_position()[1] in (-0.30, 0.30)
                    and text.get_text() in "ACGT"
                ]
            )
            == 74
        )
        boxes = [
            patch
            for axis in (graph, duplex)
            for patch in axis.patches
            if isinstance(patch, FancyBboxPatch)
        ]
        assert len(boxes) == 8
        assert all(patch.get_facecolor()[:3] == to_rgb("#D2D2D2") for patch in boxes)
        for axis in figure.axes:
            for collection in axis.collections:
                assert all(
                    np.ptp(color[:3]) == 0 for color in collection.get_facecolors()
                )
        draw_document(document, transition_index=4, progress=1, figure=figure)
        final_text = [
            (text.get_text(), text.get_position(), text.get_fontsize())
            for axis in figure.axes
            for text in axis.texts
        ]
        assert resting_text == final_text
    finally:
        plt.close(figure)


def test_producer_resting_frame_is_explicit_and_crossfade_is_stationary() -> None:
    calls = []
    rest = np.full((20, 40, 3), 255, dtype=np.uint8)
    rest[6:10, 12:24] = 210
    placed = rest.copy()
    placed[6:10, 12:24] = (40, 100, 80)

    def render(_document: PlaybackDocument, index: int | None) -> np.ndarray:
        calls.append(index)
        return rest if index is None else placed

    document = longer_document()
    frames = DuplexFrames(render, document)
    assert np.array_equal(duplex_transition_frame(frames, document, 0, 0), rest)
    middle = duplex_transition_frame(frames, document, 0, 0.5)
    assert np.array_equal(middle[0:6], rest[0:6])
    assert np.array_equal(middle[10:], rest[10:])
    assert np.all(middle[6:10, 12:24] == (125, 155, 145))
    assert None in calls
    assert 0 in calls
    assert len(calls) == 2


def test_producer_rejects_resting_shape_drift_and_negative_indices() -> None:
    def render(_document: PlaybackDocument, index: int | None) -> np.ndarray:
        return np.zeros((20 if index is None else 21, 40, 3), dtype=np.uint8)

    frames = DuplexFrames(render, longer_document())
    frames.frame(None)
    with pytest.raises(ValueError, match="shape"):
        frames.frame(0)
    with pytest.raises(IndexError, match="index"):
        frames.frame(-1)


def test_long_native_placement_tracks_do_not_overlap_or_leave_axes() -> None:
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(longer_document(), transition_index=0, progress=0, figure=figure)
        axis = figure.axes[1]
        bounds = [
            patch.get_bbox()
            for patch in axis.patches
            if isinstance(patch, FancyBboxPatch)
        ]
        assert len(bounds) == 4
        assert not any(
            left.overlaps(right)
            for index, left in enumerate(bounds)
            for right in bounds[index + 1 :]
        )
        assert all(
            axis.get_ylim()[0] <= box.y0 < box.y1 <= axis.get_ylim()[1]
            for box in bounds
        )
    finally:
        plt.close(figure)


@pytest.mark.parametrize("orientation", [Orientation.FORWARD, Orientation.REVERSE])
@pytest.mark.parametrize("with_legend", [False, True])
def test_native_captions_clear_bars_and_remain_in_fixed_scene_bounds(
    orientation: Orientation, with_legend: bool
) -> None:
    document = longer_document()
    document = replace(
        document,
        plan=replace(
            document.plan,
            steps=tuple(
                replace(step, orientation=orientation) for step in document.plan.steps
            ),
        ),
        label_overrides={str(index): f"Motif {index + 1}" for index in range(4)},
        presentation=replace(
            document.presentation,
            legend_entries=document.presentation.legend_entries if with_legend else (),
        ),
    )
    figure = plt.figure(figsize=(16, 2.4))
    states = []
    try:
        for transition, progress in ((0, 0), (4, 1)):
            draw_document(
                document, transition_index=transition, progress=progress, figure=figure
            )
            figure.canvas.draw()
            axis = figure.axes[1]
            renderer = figure.canvas.get_renderer()
            labels = [
                text.get_window_extent(renderer)
                for text in axis.texts
                if text.get_text().startswith("Motif ")
            ]
            bars = [patch.get_window_extent(renderer) for patch in axis.patches]
            bases = [
                text.get_window_extent(renderer)
                for text in axis.texts
                if text.get_text() in "ACGT"
            ]
            assert len(labels) == 4
            assert not any(label.overlaps(bar) for label in labels for bar in bars)
            assert not any(label.overlaps(base) for label in labels for base in bases)
            bounds = axis.get_window_extent(renderer)
            assert all(
                bounds.contains(label.x0, label.y0)
                and bounds.contains(label.x1, label.y1)
                for label in labels
            )
            states.append([tuple(label.bounds) for label in labels])
        assert states[0] == states[1]
    finally:
        plt.close(figure)


@pytest.mark.parametrize("orientation", [Orientation.FORWARD, Orientation.REVERSE])
@pytest.mark.parametrize("with_legend", [False, True])
def test_native_motif_glyphs_share_duplex_grid_font_and_fit_their_bars(
    orientation: Orientation, with_legend: bool
) -> None:
    document = longer_document()
    document = replace(
        document,
        plan=replace(
            document.plan,
            steps=tuple(
                replace(step, orientation=orientation) for step in document.plan.steps
            ),
        ),
        presentation=replace(
            document.presentation,
            legend_entries=document.presentation.legend_entries if with_legend else (),
        ),
    )
    figure = plt.figure(figsize=(16, 2.4))
    try:
        for transition, progress in ((0, 0), (4, 1)):
            draw_document(
                document, transition_index=transition, progress=progress, figure=figure
            )
            figure.canvas.draw()
            renderer = figure.canvas.get_renderer()
            axis = figure.axes[1]
            strand_y = -0.30 if orientation == Orientation.REVERSE else 0.30
            strand = {
                text.get_position()[0]: text
                for text in axis.texts
                if text.get_position()[1] == strand_y and text.get_text() in "ACGT"
            }
            for step, patch in zip(document.plan.steps, axis.patches, strict=True):
                glyphs = [
                    text
                    for text in axis.texts
                    if text.get_text() in "ACGT"
                    and abs(
                        text.get_position()[1]
                        - (patch.get_y() + patch.get_height() / 2)
                    )
                    < 1e-9
                    and step.start <= text.get_position()[0] < step.end
                ]
                assert len(glyphs) == 16
                for offset, text in enumerate(glyphs):
                    reference = strand[step.start + offset + 0.5]
                    assert text.get_text() == reference.get_text()
                    assert text.get_position()[0] == reference.get_position()[0]
                    assert text.get_fontsize() == reference.get_fontsize()
                    assert text.get_fontfamily() == reference.get_fontfamily()
                    glyph = text.get_window_extent(renderer)
                    expected = reference.get_window_extent(renderer)
                    assert glyph.x0 == pytest.approx(expected.x0)
                    assert glyph.x1 == pytest.approx(expected.x1)
                    # Actual cap ink, rather than phantom font descenders, is
                    # checked against cell centers and bars in the typography tests.
            assert not any(text.get_text().isdigit() for text in axis.texts)
    finally:
        plt.close(figure)


def test_progress_colors_only_completed_and_current_graph_placements() -> None:
    document = longer_document()
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(document, transition_index=1, progress=0.5, figure=figure)
        nodes = [
            patch
            for patch in figure.axes[0].patches
            if isinstance(patch, FancyBboxPatch)
        ]
        assert nodes[0].get_facecolor()[:3] == to_rgb(document.step_color(0))
        assert nodes[1].get_facecolor()[:3] != to_rgb("#D2D2D2")
        assert nodes[1].get_facecolor()[:3] != to_rgb(document.step_color(1))
        assert all(node.get_facecolor()[:3] == to_rgb("#D2D2D2") for node in nodes[2:])
    finally:
        plt.close(figure)


def test_uncovered_bases_remain_visible_gray_without_fabricated_traversal() -> None:
    document = longer_document()
    realized = RealizedArray(
        source_id="test:gap",
        sequence="ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGT",
        placements=(
            Placement("left", "left", PlacementKind.TFBS, "ACGTTGCAAGTC", 0),
            Placement("right", "right", PlacementKind.TFBS, "GCTTAGGACGT", 26),
        ),
    )
    document = replace(document, plan=reconstruct_playback(realized))
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(document, transition_index=2, progress=1, figure=figure)
        duplex = figure.axes[1]
        bases = [
            text
            for text in duplex.texts
            if text.get_text() in "ACGT" and text.get_position()[1] in (-0.30, 0.30)
        ]
        assert len(bases) == 74
        assert all(
            to_rgb(text.get_color()) == to_rgb("#D2D2D2")
            for text in bases
            if 12 <= text.get_position()[0] < 26
        )
        assert not figure.axes[0].collections
    finally:
        plt.close(figure)


def test_encoded_gif_opens_with_complete_neutral_native_scene(tmp_path: Path) -> None:
    document = longer_document()
    path = render_collection_gif(
        (document, document),
        tmp_path / "rest.gif",
        fps=5,
        seconds_per_step=0.1,
        lead_seconds=0,
        hold_seconds=0,
        scene_transition_seconds=0.2,
    )
    with Image.open(path) as image:
        first = np.asarray(image.convert("RGB"))
        assert image.size == (1600, 240)
        for index in range(image.n_frames):
            image.seek(index)
            frame = np.asarray(image.convert("RGB"))
            assert (
                np.count_nonzero(np.min(frame[10:190, 600:1500], axis=2) < 240) > 1000
            )
    right = first[10:190, 600:1500].astype(int)
    assert np.count_nonzero(np.min(right, axis=2) < 240) > 1000
    assert np.max(np.ptp(right, axis=2)) <= 3


def test_invalid_resting_frame_fails_before_output(tmp_path: Path) -> None:
    def render(_document: PlaybackDocument, index: int | None) -> np.ndarray:
        return np.zeros((0 if index is None else 20, 40, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="frame shape"):
        render_collection_gif(
            (longer_document(),),
            tmp_path / "invalid.gif",
            fps=1,
            duplex_frame_renderer=render,
        )
    assert list(tmp_path.iterdir()) == []
