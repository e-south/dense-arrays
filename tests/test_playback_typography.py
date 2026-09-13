"""Measured nucleotide spacing and canvas evidence contracts.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.textpath import TextPath
from matplotlib.transforms import Bbox

from dense_arrays.playback import reconstruct_playback
from dense_arrays.playback.duplex_frames import DuplexFrames
from dense_arrays.playback.presentation import PlaybackDocument, evidence_metadata
from dense_arrays.playback.scene_drawing import draw_document
from dense_arrays.playback.theme import LegendEntry, PlaybackPresentation
from dense_arrays.realized import (
    DeclaredConstraint,
    Orientation,
    Placement,
    PlacementKind,
    RealizedArray,
)

if TYPE_CHECKING:
    from matplotlib.text import Text


def ink_bounds(text: Text) -> Bbox:
    """Measure drawn letter outlines rather than the font's phantom descenders."""
    outline = TextPath((0, 0), text.get_text(), prop=text.get_fontproperties())
    bounds = outline.get_extents()
    x, y = text.get_transform().transform(text.get_position())
    scale = text.figure.dpi / 72
    # Text alignment offsets use layout metrics, independent of ink extents.
    window = text.get_window_extent(text.figure.canvas.get_renderer())
    if text.get_ha() == "center":
        x -= window.width / 2
    if text.get_va() == "center":
        _, _, descent = text.figure.canvas.get_renderer().get_text_width_height_descent(
            "lp", text.get_fontproperties(), ismath=False
        )
        y -= window.height / 2 - descent
    return Bbox.from_bounds(
        x + bounds.x0 * scale,
        y + bounds.y0 * scale,
        bounds.width * scale,
        bounds.height * scale,
    )


def compact_document(
    orientation: Orientation = Orientation.FORWARD,
) -> PlaybackDocument:
    sequence = "ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGTTCA"
    realized = RealizedArray(
        source_id="test:compact",
        sequence=sequence,
        placements=tuple(
            Placement(
                str(index),
                str(index),
                PlacementKind.TFBS,
                sequence[start : start + 16],
                start,
                orientation=orientation,
            )
            for index, start in enumerate((0, 8, 16, 24))
        ),
    )
    return PlaybackDocument(reconstruct_playback(realized), title="Four motifs")


@pytest.mark.parametrize("width", [16, 24])
@pytest.mark.parametrize("orientation", [Orientation.FORWARD, Orientation.REVERSE])
def test_native_cells_remain_compact_and_center_actual_ink(
    width: int, orientation: Orientation
) -> None:
    document = compact_document(orientation)
    figure = plt.figure(figsize=(width, 2.4), dpi=100)
    try:
        states = []
        for transition, progress in ((0, 0), (4, 1)):
            draw_document(
                document, transition_index=transition, progress=progress, figure=figure
            )
            figure.canvas.draw()
            axis = figure.axes[1]
            strand = [
                text
                for text in axis.texts
                if text.get_text() in "ACGT" and text.get_position()[1] == 0.3
            ]
            pitch = (
                axis.transData.transform((1, 0))[0]
                - axis.transData.transform((0, 0))[0]
            )
            assert (
                1.02 <= pitch / max(ink_bounds(text).width for text in strand) <= 1.35
            )
            assert (
                len(
                    {
                        text.get_fontsize()
                        for text in axis.texts
                        if text.get_text() in "ACGT"
                    }
                )
                == 1
            )
            graph = [
                text for text in figure.axes[0].texts if len(text.get_text()) == 16
            ]
            assert graph
            assert graph[0].get_fontsize() == pytest.approx(strand[0].get_fontsize())
            for step, patch in zip(document.plan.steps, axis.patches, strict=True):
                glyphs = [
                    text
                    for text in axis.texts
                    if text.get_text() in "ACGT"
                    and text.get_position()[1] == patch.get_y() + patch.get_height() / 2
                    and step.start <= text.get_position()[0] < step.end
                ]
                assert len(glyphs) == 16
                for text in glyphs:
                    ink = ink_bounds(text)
                    cell_center = axis.transData.transform(text.get_position())
                    assert (ink.x0 + ink.x1) / 2 == pytest.approx(
                        cell_center[0], abs=0.01
                    )
                    assert (ink.y0 + ink.y1) / 2 == pytest.approx(
                        cell_center[1], abs=0.01
                    )
                    box = patch.get_window_extent(figure.canvas.get_renderer())
                    assert box.contains(ink.x0, ink.y0)
                    assert box.contains(ink.x1, ink.y1)
            states.append(
                [
                    (text.get_position(), tuple(ink_bounds(text).bounds))
                    for text in strand
                ]
            )
        assert states[0] == states[1]
    finally:
        plt.close(figure)


@pytest.mark.parametrize("producer", [False, True])
def test_routine_authority_is_metadata_only_even_with_optional_notices(
    producer: bool,
) -> None:
    document = replace(
        compact_document(),
        presentation=PlaybackPresentation(show_authority_notice=True),
    )
    figure = plt.figure(figsize=(16, 2.4))
    try:
        frames = (
            DuplexFrames(
                lambda _doc, _index: np.zeros((100, 500, 3), dtype=np.uint8), document
            )
            if producer
            else None
        )
        draw_document(
            document,
            transition_index=4,
            progress=1,
            figure=figure,
            duplex_frames=frames,
        )
        visible = " ".join(text.get_text() for text in figure.texts)
        visible += " ".join(
            text.get_text() for axis in figure.axes for text in axis.texts
        )
        assert "reconstructed" not in visible.lower()
        assert "Unique coordinate order" not in visible
        assert (
            "Reconstructed from placements · Unique coordinate order"
            in evidence_metadata(document)
        )
    finally:
        plt.close(figure)


def test_legend_uses_same_relative_scale_as_native_annotation() -> None:
    document = replace(
        compact_document(),
        label_overrides={"0": "Anchor"},
        presentation=PlaybackPresentation(
            legend_entries=(LegendEntry("anchor", "Anchor", "#365E80"),)
        ),
    )
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(document, transition_index=4, progress=1, figure=figure)
        label = next(
            text for text in figure.axes[1].texts if text.get_text() == "Anchor"
        )
        legend = next(
            text for text in figure.axes[2].texts if text.get_text() == "Anchor"
        )
        assert label.get_fontsize() == pytest.approx(legend.get_fontsize())
    finally:
        plt.close(figure)


def test_native_distance_endpoints_share_the_compact_coordinate_grid() -> None:
    source = compact_document()
    realized = RealizedArray(
        source_id="test:compact-distance",
        sequence=source.plan.realized_sequence,
        placements=tuple(
            Placement(
                step.placement_id,
                step.placement_id,
                PlacementKind.TFBS,
                step.placement_sequence,
                step.start,
            )
            for step in source.plan.steps
        ),
        constraints=(DeclaredConstraint("spacing", "0", "3", 8, 8),),
    )
    document = replace(source, plan=reconstruct_playback(realized))
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(document, transition_index=4, progress=1, figure=figure)
        axis = figure.axes[1]
        bracket = next(
            line for line in axis.lines if line.get_gid() == "distance-bracket"
        )
        points = bracket.get_transform().transform(bracket.get_xydata())
        expected = axis.transData.transform(((16, 0), (24, 0)))
        assert points[0, 0] == pytest.approx(expected[0, 0])
        assert points[-1, 0] == pytest.approx(expected[-1, 0])
        label = next(
            text for text in axis.texts if text.get_text().startswith("spacing:")
        )
        nucleotide = next(text for text in axis.texts if text.get_text() == "A")
        assert label.get_fontsize() / nucleotide.get_fontsize() == pytest.approx(
            11.5 / 13.2
        )
    finally:
        plt.close(figure)


def test_ambiguous_nucleotide_ink_remains_inside_each_coordinate_cell() -> None:
    sequence = "ACGTWSMKRYBDHVN"
    document = PlaybackDocument(
        reconstruct_playback(
            RealizedArray(
                source_id="test:iupac",
                sequence=sequence,
                placements=(
                    Placement("iupac", "iupac", PlacementKind.OTHER, sequence, 0),
                ),
            )
        ),
        title="IUPAC bases",
    )
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(document, transition_index=1, progress=1, figure=figure)
        figure.canvas.draw()
        axis = figure.axes[1]
        for text in axis.texts:
            if len(text.get_text()) != 1:
                continue
            ink = ink_bounds(text)
            x, y = text.get_position()
            left = axis.transData.transform((x - 0.5, y))[0]
            right = axis.transData.transform((x + 0.5, y))[0]
            assert left < ink.x0 < ink.x1 < right
    finally:
        plt.close(figure)


@pytest.mark.parametrize("orientation", [Orientation.FORWARD, Orientation.REVERSE])
def test_two_native_distance_labels_clear_placement_tracks(
    orientation: Orientation,
) -> None:
    source = compact_document(orientation)
    realized = RealizedArray(
        source_id="test:bracket-clearance",
        sequence=source.plan.realized_sequence,
        placements=tuple(
            Placement(
                step.placement_id,
                step.placement_id,
                PlacementKind.TFBS,
                step.placement_sequence,
                step.start,
                orientation=orientation,
            )
            for step in source.plan.steps
        ),
        constraints=(
            DeclaredConstraint("separation", "0", "3", 8, 8),
            DeclaredConstraint("abutment", "0", "2", 0, 0),
        ),
    )
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(
            replace(source, plan=reconstruct_playback(realized)),
            transition_index=4,
            progress=1,
            figure=figure,
        )
        figure.canvas.draw()
        axis = figure.axes[1]
        renderer = figure.canvas.get_renderer()
        labels = [
            text.get_window_extent(renderer)
            for text in axis.texts
            if "required" in text.get_text()
        ]
        bars = [patch.get_window_extent(renderer) for patch in axis.patches]
        assert len(labels) == 2
        assert not labels[0].overlaps(labels[1])
        assert not any(label.overlaps(bar) for label in labels for bar in bars)
    finally:
        plt.close(figure)
