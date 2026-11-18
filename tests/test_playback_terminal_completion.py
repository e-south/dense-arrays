"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/tests/test_playback_terminal_completion.py

Finish playback with the last placed part, without a separate trip to End.

Module Author(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np
import pytest

from dense_arrays.playback.frame_schedule import PlaybackTiming, scene_frame_schedule
from dense_arrays.playback.scene_drawing import draw_document, transition_frame_counts
from dense_arrays.playback.theme import PlaybackPresentation

from .test_playback_typography import compact_document


@pytest.mark.parametrize("graph", [True, False])
def test_schedule_finishes_when_the_last_part_arrives(graph: bool):
    document = compact_document()
    if not graph:
        document = replace(
            document,
            presentation=PlaybackPresentation(graph_detail="none", graph_fraction=0),
        )
    figure = plt.figure(figsize=(16, 2.4))
    try:
        counts = transition_frame_counts(document, figure, fps=24, seconds_per_step=0.5)
        assert len(counts) == len(document.plan.steps)
        frames = tuple(scene_frame_schedule(counts, PlaybackTiming(fps=24), first=True))
        assert frames[-1].transition_index == len(document.plan.steps) - 1
        assert frames[-1].progress == 1
    finally:
        plt.close(figure)


def test_last_arrival_and_complete_poster_are_identical():
    document = compact_document()
    figure = plt.figure(figsize=(16, 2.4))
    try:
        images = []
        for index in (len(document.plan.steps) - 1, len(document.plan.steps)):
            draw_document(document, transition_index=index, progress=1, figure=figure)
            figure.canvas.draw()
            images.append(np.asarray(figure.canvas.buffer_rgba()).copy())
        assert np.array_equal(*images)
    finally:
        plt.close(figure)


def test_no_moving_marker_can_visit_end_after_sequence_completion():
    document = compact_document()
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(
            document,
            transition_index=len(document.plan.steps),
            progress=0.5,
            figure=figure,
        )
        assert not figure.axes[0].collections
    finally:
        plt.close(figure)
