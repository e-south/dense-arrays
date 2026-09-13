"""Frame validation, timing domains, and export cleanup.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import weakref
from dataclasses import replace
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
import pytest
from matplotlib.animation import PillowWriter
from PIL import Image

from dense_arrays.playback import export, matplotlib_renderer, reconstruct_playback
from dense_arrays.playback.duplex_frames import DuplexFrames
from dense_arrays.playback.frame_schedule import PlaybackTiming, scene_frame_schedule
from dense_arrays.playback.models import PlaybackNotice
from dense_arrays.playback.presentation import PlaybackDocument
from dense_arrays.playback.scene_drawing import draw_document
from dense_arrays.playback.theme import PlaybackPresentation
from dense_arrays.realized import (
    DeclaredConstraint,
    Placement,
    PlacementKind,
    RealizedArray,
)

if TYPE_CHECKING:
    from pathlib import Path


def document() -> PlaybackDocument:
    return PlaybackDocument(
        reconstruct_playback(
            RealizedArray(
                source_id="test:export",
                sequence="AAACCC",
                placements=(
                    Placement("left", "left", PlacementKind.TFBS, "AAA", 0),
                    Placement("right", "right", PlacementKind.TFBS, "CCC", 3),
                ),
            )
        ),
        title="Frame contract",
    )


@pytest.mark.parametrize("invalid", [float("nan"), float("inf"), True, 1.2])
def test_export_rejects_invalid_fps_before_output(
    tmp_path: Path, invalid: object
) -> None:
    with pytest.raises((ValueError, TypeError), match="fps"):
        matplotlib_renderer.render_collection_gif(
            (document(),), tmp_path / "out.gif", fps=invalid
        )
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "field",
    ["seconds_per_step", "lead_seconds", "hold_seconds", "scene_transition_seconds"],
)
def test_export_rejects_nonfinite_durations(tmp_path: Path, field: str) -> None:
    with pytest.raises((ValueError, TypeError), match=field):
        matplotlib_renderer.render_collection_gif(
            (document(),), tmp_path / "out.gif", **{field: float("nan")}
        )
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "invalid",
    [
        np.zeros((10, 10)),
        np.zeros((0, 10, 3), dtype=np.uint8),
        np.zeros((10, 10, 3)),
        np.zeros((10, 10, 2), dtype=np.uint8),
    ],
)
def test_invalid_raster_callback_frames_are_rejected(
    tmp_path: Path, invalid: np.ndarray
) -> None:
    with pytest.raises((ValueError, TypeError), match="frame"):
        matplotlib_renderer.render_collection_poster_png(
            (document(),),
            tmp_path / "out.png",
            duplex_frame_renderer=lambda *_: invalid,
        )
    assert list(tmp_path.iterdir()) == []


def test_callback_shape_changes_are_rejected_during_animation(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="shape"):
        matplotlib_renderer.render_collection_gif(
            (document(),),
            tmp_path / "out.gif",
            fps=2,
            duplex_frame_renderer=lambda _doc, step: np.zeros(
                (10 + step, 10, 3), dtype=np.uint8
            ),
        )
    assert list(tmp_path.iterdir()) == []


def test_failed_callback_leaves_no_figure_or_partial_output(tmp_path: Path) -> None:
    before = plt.get_fignums()

    def fail(*_args: object) -> None:
        message = "producer failed"
        raise RuntimeError(message)

    with pytest.raises(RuntimeError, match="producer failed"):
        matplotlib_renderer.render_collection_poster_png(
            (document(),), tmp_path / "out.png", duplex_frame_renderer=fail
        )
    assert plt.get_fignums() == before
    assert list(tmp_path.iterdir()) == []


def test_zero_lead_hold_schedule_has_only_transition_frames() -> None:

    timing = PlaybackTiming(
        fps=2,
        seconds_per_step=1,
        lead_seconds=0,
        hold_seconds=0,
        scene_transition_seconds=0,
    )
    frames = tuple(scene_frame_schedule((3, 5), timing, first=True, last=True))
    assert len(frames) == 8
    assert frames[0].transition_index == 0
    assert frames[-1].transition_index == 1
    assert frames[-1].progress == 1


def test_callback_count_and_final_frame_are_preserved(tmp_path: Path) -> None:

    calls = []

    def frame(_document: PlaybackDocument, index: int) -> np.ndarray:
        calls.append(index)
        return np.full(
            (80, 120, 3), (0, 40, 160) if index == 1 else 255, dtype=np.uint8
        )

    path = matplotlib_renderer.render_collection_poster_png(
        (document(),), tmp_path / "poster.png", duplex_frame_renderer=frame, dpi=30
    )
    assert calls == [1]
    with Image.open(path) as image:
        pixels = np.asarray(image)
    assert np.any(np.all(pixels[..., :3] == (0, 40, 160), axis=2))


def test_frames_load_lazily_and_release_older_images() -> None:

    calls = []

    def frame(_document: PlaybackDocument, index: int) -> np.ndarray:
        calls.append(index)
        return np.zeros((10, 10, 3), dtype=np.uint8)

    plan = reconstruct_playback(
        RealizedArray(
            source_id="test:cache",
            sequence="A" * 12,
            placements=tuple(
                Placement(str(i), str(i), PlacementKind.TFBS, "AAA", i * 3)
                for i in range(4)
            ),
        )
    )
    frames = DuplexFrames(frame, PlaybackDocument(plan, title="Cache"))
    assert calls == []
    references = [weakref.ref(frames.frame(index)) for index in range(4)]
    assert calls == [0, 1, 2, 3]
    assert sum(reference() is not None for reference in references) <= 2


def test_late_invalid_frame_preserves_prior_output_and_closes_figure(
    tmp_path: Path,
) -> None:
    path = tmp_path / "prior.gif"
    path.write_bytes(b"prior output")
    before = plt.get_fignums()

    def frame(_document: PlaybackDocument, index: int) -> np.ndarray:
        if index == 1:
            return np.zeros((10, 10), dtype=np.uint8)
        return np.zeros((10, 10, 3), dtype=np.uint8)

    with pytest.raises(ValueError, match="frame shape"):
        matplotlib_renderer.render_collection_gif(
            (document(),), path, fps=2, duplex_frame_renderer=frame
        )
    assert path.read_bytes() == b"prior output"
    assert list(tmp_path.iterdir()) == [path]
    assert plt.get_fignums() == before


def test_poster_preserves_producer_figure_size(tmp_path: Path) -> None:
    class Producer:
        preferred_figure_height_inches = 2.4

        def render(self, _document: PlaybackDocument, _index: int) -> np.ndarray:
            return np.zeros((80, 120, 3), dtype=np.uint8)

    path = matplotlib_renderer.render_collection_poster_png(
        (document(),),
        tmp_path / "size.png",
        dpi=100,
        duplex_frame_renderer=Producer().render,
    )
    with Image.open(path) as image:
        assert image.size == (1600, 240)


@pytest.mark.parametrize("invalid_from", [0, 1])
def test_invalid_producer_metric_preserves_primary_error_and_prior_output(
    tmp_path: Path, invalid_from: int
) -> None:
    class Producer:
        native_nucleotide_cap_height_px = 10.0

        def render(self, _document: PlaybackDocument, index: int) -> np.ndarray:
            if index >= invalid_from:
                self.native_nucleotide_cap_height_px = float("nan")
            return np.zeros((10, 10, 3), dtype=np.uint8)

    path = tmp_path / "prior.gif"
    path.write_bytes(b"prior output")
    before = plt.get_fignums()
    with pytest.raises(ValueError, match="native_nucleotide_cap_height_px"):
        matplotlib_renderer.render_collection_gif(
            (document(),), path, fps=1, duplex_frame_renderer=Producer().render
        )
    assert path.read_bytes() == b"prior output"
    assert list(tmp_path.iterdir()) == [path]
    assert plt.get_fignums() == before


def test_encoder_cleanup_failure_does_not_mask_late_callback_error(
    tmp_path: Path,
) -> None:
    class FailingCleanupWriter(PillowWriter):
        def finish(self) -> None:
            message = "encoder cleanup failed"
            raise RuntimeError(message)

    def frame(_document: PlaybackDocument, index: int) -> np.ndarray:
        if index == 1:
            message = "producer late failure"
            raise ValueError(message)
        return np.zeros((10, 10, 3), dtype=np.uint8)

    before = plt.get_fignums()
    with pytest.raises(ValueError, match="producer late failure"):
        export.render_animation(
            (document(),),
            tmp_path / "out.gif",
            timing=PlaybackTiming(fps=1),
            writer=FailingCleanupWriter(fps=1),
            dpi=30,
            duplex_frame_renderer=frame,
        )
    assert list(tmp_path.iterdir()) == []
    assert plt.get_fignums() == before


def crowded_document(*, failures: bool) -> PlaybackDocument:
    realized = RealizedArray(
        source_id="test:crowded",
        sequence="AAACCC",
        placements=(
            Placement("left", "left", PlacementKind.TFBS, "AAA", 0),
            Placement("right", "right", PlacementKind.TFBS, "CCC", 3),
        ),
        constraints=tuple(
            DeclaredConstraint(f"c{i}", "left", "right", 1, 1) for i in range(15)
        )
        if failures
        else (),
    )
    notices = () if failures else (PlaybackNotice("source", "explanation " * 210),)
    return replace(
        document(),
        plan=reconstruct_playback(realized, notices=notices),
        presentation=PlaybackPresentation(show_authority_notice=not failures),
    )


@pytest.mark.parametrize("failures", [True, False])
def test_crowded_evidence_preserves_size_and_complete_png_metadata(
    tmp_path: Path, failures: bool
) -> None:
    scene = crowded_document(failures=failures)
    figure = plt.figure(figsize=(16, 2.4))
    try:
        draw_document(scene, transition_index=2, progress=1.0, figure=figure)
        text = "\n".join(artist.get_text() for artist in figure.texts)
        assert "full" in text
        assert "metadata" in text
        assert "Reconstructed from placements" in text
        assert all(axis.get_position().y0 < 0.25 for axis in figure.axes)
        if failures:
            assert "FAILED 15 distance constraints" in text
            assert any(
                artist.get_text() == "15 declared distances; full results in metadata"
                for axis in figure.axes
                for artist in axis.texts
            )
            assert not any(
                line.get_gid() == "distance-bracket"
                for axis in figure.axes
                for line in axis.lines
            )
    finally:
        plt.close(figure)
    path = matplotlib_renderer.render_collection_poster_png(
        (scene,), tmp_path / "crowded.png", dpi=100
    )
    with Image.open(path) as image:
        assert image.size == (1600, 240)
        evidence = image.info["PlaybackEvidence"]
    if failures:
        assert all(
            f"FAILED c{i}: 0 bp (required 1..1 bp)" in evidence for i in range(15)
        )
    else:
        assert "explanation " * 210 in evidence


def test_gif_retains_complete_evidence_in_native_comment(tmp_path: Path) -> None:
    scene = crowded_document(failures=True)
    path = matplotlib_renderer.render_collection_gif(
        (scene,), tmp_path / "evidence.gif", fps=1
    )
    with Image.open(path) as image:
        evidence = image.info["comment"].decode("utf-8")
        assert image.size == (1600, 240)
    assert all(f"FAILED c{i}: 0 bp (required 1..1 bp)" in evidence for i in range(15))


def test_png_retains_passed_distances_when_brackets_are_summarized(
    tmp_path: Path,
) -> None:
    scene = crowded_document(failures=True)
    plan = replace(
        scene.plan,
        constraint_results=tuple(
            replace(result, min_distance_bp=0, max_distance_bp=0, passed=True)
            for result in scene.plan.constraint_results
        ),
    )
    path = matplotlib_renderer.render_collection_poster_png(
        (replace(scene, plan=plan),), tmp_path / "passed.png", dpi=30
    )
    with Image.open(path) as image:
        evidence = image.info["PlaybackEvidence"]
    assert all(f"PASSED c{i}: 0 bp (required 0..0 bp)" in evidence for i in range(15))
