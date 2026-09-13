"""Visible evidence and producer-owned presentation contracts.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.patches import FancyArrowPatch

from dense_arrays.playback import reconstruct_playback
from dense_arrays.playback.duplex_frames import DuplexFrames
from dense_arrays.playback.graph.projection import project_explanation_graph
from dense_arrays.playback.presentation import PlaybackDocument, resolve_evidence
from dense_arrays.playback.scene_drawing import document_axes, draw_document
from dense_arrays.playback.theme import PlaybackPresentation, step_color
from dense_arrays.realized import (
    DeclaredConstraint,
    Placement,
    PlacementKind,
    RealizedArray,
)


def gapped_document() -> PlaybackDocument:
    return PlaybackDocument(
        plan=reconstruct_playback(
            RealizedArray(
                source_id="test:gap",
                sequence="AAATTTCCC",
                placements=(
                    Placement("left", "left", PlacementKind.TFBS, "AAA", 0),
                    Placement("right", "right", PlacementKind.TFBS, "CCC", 6),
                ),
                constraints=(DeclaredConstraint("spacing", "left", "right", 0, 0),),
            )
        ),
        title="Declared spacing",
        subtitle="Persisted coordinate evidence",
        label_overrides={"left": "Anchor A"},
    )


def test_layout_only_has_no_traversal_geometry() -> None:
    document = gapped_document()
    graph = project_explanation_graph(document.plan)
    assert graph.traversal_edges == ()


def test_raster_displays_order_failure_and_labels() -> None:
    figure = plt.figure(figsize=(16, 4))
    try:
        draw_document(
            gapped_document(), transition_index=1, progress=0.5, figure=figure
        )
        text = "\n".join(artist.get_text() for artist in figure.texts)
        text += "\n" + "\n".join(
            artist.get_text() for axis in figure.axes for artist in axis.texts
        )
        assert "Reconstructed from placements" in text
        assert "Layout only" in text
        assert "FAILED" in text
        assert "3 bp" in text
        assert "Anchor A" in text
        assert not any(
            isinstance(patch, FancyArrowPatch) and patch.get_linewidth() == 2.4
            for axis in figure.axes
            for patch in axis.patches
        )
    finally:
        plt.close(figure)


def test_caller_colors_and_labels_are_resolved_for_raster() -> None:

    document = replace(
        gapped_document(),
        color_overrides={"left": "#123456"},
        presentation=PlaybackPresentation(graph_detail="none", graph_fraction=0),
    )
    figure = plt.figure(figsize=(16, 4))
    try:
        draw_document(document, transition_index=1, progress=1.0, figure=figure)
        colors = [
            patch.get_facecolor() for axis in figure.axes for patch in axis.patches
        ]
        assert any(color[:3] == (18 / 255, 52 / 255, 86 / 255) for color in colors)
        assert len(figure.axes) == 1
    finally:
        plt.close(figure)


def test_generic_palette_does_not_infer_categories_from_labels() -> None:

    step = gapped_document().plan.steps[0]
    assert step_color(replace(step, label="BaeR -10 downstream"), 0) == step_color(
        step, 0
    )
    fixed = replace(step, placement_kind="fixed_element")
    assert step_color(replace(fixed, label="downstream -10"), 0) == step_color(fixed, 0)


def test_unsupported_presentation_settings_fail() -> None:

    for settings in (
        {"color_profile": "secg"},
        {"color_profile": "typo"},
        {"graph_detail": "unknown"},
        {"graph_detail": "inset"},
    ):
        with pytest.raises(ValueError, match=r"color_profile|graph_detail"):
            PlaybackPresentation(**settings)


def test_presentation_map_identity_and_color_validation() -> None:

    for fields in (
        {"label_overrides": {"unknown": "X"}},
        {"color_overrides": {"left": "url(x)"}},
    ):
        with pytest.raises(ValueError, match=r"unknown placement|colors must"):
            replace(gapped_document(), **fields)


def test_distance_bracket_control_changes_native_artists_() -> None:

    for mode, expected in (("never", 0), ("when_declared", 1), ("always", 1)):
        document = replace(
            gapped_document(),
            presentation=PlaybackPresentation(show_distance_bracket=mode),
        )
        figure = plt.figure(figsize=(16, 4))
        try:
            draw_document(document, transition_index=1, progress=1.0, figure=figure)
            brackets = [
                line
                for axis in figure.axes
                for line in axis.lines
                if line.get_gid() == "distance-bracket"
            ]
            assert len(brackets) == expected
        finally:
            plt.close(figure)


def test_optional_notices_and_required_evidence_have_distinct_controls() -> None:

    document = gapped_document()
    base = resolve_evidence(document)
    verbose = resolve_evidence(
        replace(document, presentation=PlaybackPresentation(show_authority_notice=True))
    )
    assert base.qualification == verbose.qualification
    assert base.constraints == verbose.constraints
    assert base.notices == ()
    assert verbose.notices
    figure = plt.figure(figsize=(16, 4))
    try:
        draw_document(
            replace(
                document, presentation=PlaybackPresentation(show_authority_notice=True)
            ),
            transition_index=1,
            progress=1.0,
            figure=figure,
        )
        texts = "\n".join(text.get_text() for text in figure.texts)
        assert "notices; full text in metadata" in texts
    finally:
        plt.close(figure)


def test_producer_distance_bracket_capability_prevents_duplicate_artists() -> None:

    class Producer:
        renders_distance_brackets = True

        def render(self, _document: PlaybackDocument, _index: int) -> np.ndarray:
            return np.zeros((80, 120, 3), dtype=np.uint8)

    document = gapped_document()
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(
            document,
            transition_index=1,
            progress=1.0,
            figure=figure,
            duplex_frames=DuplexFrames(Producer().render, document),
        )
        assert not any(
            line.get_gid() == "distance-bracket"
            for axis in figure.axes
            for line in axis.lines
        )
    finally:
        plt.close(figure)


def test_reduced_preserves_requested_graph_fraction_and_omits_context() -> None:

    document = replace(
        gapped_document(),
        presentation=PlaybackPresentation(graph_detail="reduced", graph_fraction=0.30),
    )
    figure = plt.figure(figsize=(16, 2.4))
    try:
        graph, duplex, _legend = document_axes(figure, document)
        assert graph is not None
        assert (
            graph.get_position().width / duplex.get_position().width
            == pytest.approx(0.3 / 0.7)
        )
        draw_document(document, transition_index=1, progress=1.0, figure=figure)
        assert not any(
            isinstance(patch, FancyArrowPatch)
            for axis in figure.axes
            for patch in axis.patches
        )
    finally:
        plt.close(figure)
