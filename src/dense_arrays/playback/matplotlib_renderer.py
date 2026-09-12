"""Public raster media exports for synchronized placement playback.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .duplex_frames import DuplexFrameRenderer
from .export import render_animation, render_poster, require_documents
from .frame_schedule import PlaybackTiming, positive_integer
from .presentation import PlaybackDocument, collection_comment

if TYPE_CHECKING:
    from pathlib import Path


__all__ = (
    "DuplexFrameRenderer",
    "PlaybackDocument",
    "render_collection_gif",
    "render_collection_mp4",
    "render_collection_poster_png",
)


def render_collection_poster_png(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    dpi: int = 180,
    duplex_frame_renderer: DuplexFrameRenderer | None = None,
) -> Path:
    """Write the first scene as PNG, preserving visible authority and failures."""
    require_documents(documents)
    positive_integer(dpi, "dpi")
    return render_poster(
        documents, output_path, dpi=dpi, duplex_frame_renderer=duplex_frame_renderer
    )


def render_collection_mp4(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    fps: int = 30,
    seconds_per_step: float = 0.70,
    hold_seconds: float = 0.75,
    lead_seconds: float = 0.25,
    scene_transition_seconds: float = 0.0,
    duplex_frame_renderer: DuplexFrameRenderer | None = None,
) -> Path:
    """Encode placement scenes with FFmpeg using the shared playback clock."""
    require_documents(documents)
    timing = PlaybackTiming(
        fps, seconds_per_step, hold_seconds, lead_seconds, scene_transition_seconds
    )
    from matplotlib.animation import FFMpegWriter

    if not FFMpegWriter.isAvailable():
        msg = "FFmpeg is required for MP4 export"
        raise RuntimeError(msg)
    writer = FFMpegWriter(
        fps=fps,
        codec="h264",
        metadata={
            "title": " / ".join(document.title for document in documents),
            "comment": collection_comment(documents),
        },
        extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
    )
    return render_animation(
        documents,
        output_path,
        timing=timing,
        writer=writer,
        dpi=150,
        duplex_frame_renderer=duplex_frame_renderer,
    )


def render_collection_gif(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    fps: int = 15,
    seconds_per_step: float = 0.70,
    hold_seconds: float = 0.70,
    lead_seconds: float = 0.25,
    scene_transition_seconds: float = 0.0,
    duplex_frame_renderer: DuplexFrameRenderer | None = None,
) -> Path:
    """Encode placement scenes with Pillow using the shared playback clock."""
    require_documents(documents)
    timing = PlaybackTiming(
        fps, seconds_per_step, hold_seconds, lead_seconds, scene_transition_seconds
    )
    from .gif_writer import EvidencePillowWriter

    return render_animation(
        documents,
        output_path,
        timing=timing,
        writer=EvidencePillowWriter(
            fps=fps, metadata={"comment": collection_comment(documents)}
        ),
        dpi=100,
        duplex_frame_renderer=duplex_frame_renderer,
    )
