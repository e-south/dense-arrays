"""Own media writer lifecycle, staged output publication, and figure cleanup.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import contextlib
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from .duplex_frames import (
    DuplexFrameRenderer,
    DuplexFrames,
    preferred_figure_height_inches,
)
from .frame_schedule import PlaybackFrame, PlaybackTiming, scene_frame_schedule
from .presentation import PlaybackDocument, evidence_metadata
from .scene_drawing import draw_document, transition_frame_counts

if TYPE_CHECKING:
    from collections.abc import Iterator

    from matplotlib.animation import AbstractMovieWriter
    from matplotlib.figure import Figure

_PAPER = "#ffffff"


def require_documents(documents: tuple[PlaybackDocument, ...]) -> None:
    """Require at least one validated document before allocating output."""
    if not documents or any(
        not isinstance(item, PlaybackDocument) for item in documents
    ):
        msg = "at least one PlaybackDocument is required"
        raise ValueError(msg)


@contextlib.contextmanager
def staged_output(output_path: Path) -> Iterator[Path]:
    """Publish a complete file while preserving any prior output on failure."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=output_path.parent,
        prefix=f".{output_path.stem}-",
        suffix=output_path.suffix,
        delete=False,
    ) as temporary:
        path = Path(temporary.name)
    try:
        yield path
        path.replace(output_path)
    finally:
        path.unlink(missing_ok=True)


def render_poster(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    dpi: int,
    duplex_frame_renderer: DuplexFrameRenderer | None,
) -> Path:
    """Render the first scene with a validated producer frame set."""
    import matplotlib.pyplot as plt

    require_documents(documents)
    document = documents[0]
    frames = (
        DuplexFrames(duplex_frame_renderer, document)
        if duplex_frame_renderer is not None
        else None
    )
    height = preferred_figure_height_inches(duplex_frame_renderer)
    figure = plt.figure(figsize=(16, height), facecolor=_PAPER)
    output_path = Path(output_path)
    try:
        draw_document(
            document,
            transition_index=len(document.plan.steps),
            progress=1.0,
            figure=figure,
            duplex_frames=frames,
        )
        with staged_output(output_path) as path:
            figure.savefig(
                path,
                dpi=dpi,
                facecolor=_PAPER,
                metadata={
                    "Title": document.title,
                    "Description": document.subtitle,
                    "PlaybackEvidence": evidence_metadata(document),
                },
            )
    finally:
        plt.close(figure)
    return output_path


@contextlib.contextmanager
def _saving_writer(
    writer: AbstractMovieWriter, figure: Figure, path: Path, dpi: int
) -> Iterator[None]:
    """Finish the writer while preserving any primary frame-generation error."""
    primary_error: BaseException | None = None
    try:
        with writer.saving(figure, str(path), dpi=dpi):
            try:
                yield
            except BaseException as error:
                primary_error = error
                raise
    except BaseException:
        if primary_error is not None:
            raise primary_error from None
        raise


def _draw_frame(
    document: PlaybackDocument,
    frame: PlaybackFrame,
    figure: Figure,
    duplex_frames: DuplexFrames | None,
) -> None:
    draw_document(
        document,
        transition_index=frame.transition_index,
        progress=frame.progress,
        figure=figure,
        duplex_frames=duplex_frames,
    )


def render_animation(
    documents: tuple[PlaybackDocument, ...],
    output_path: Path,
    *,
    timing: PlaybackTiming,
    writer: AbstractMovieWriter,
    dpi: int,
    duplex_frame_renderer: DuplexFrameRenderer | None,
) -> Path:
    """Encode the shared frame schedule and close resources on every exit."""
    import matplotlib.pyplot as plt

    require_documents(documents)
    height = preferred_figure_height_inches(duplex_frame_renderer)
    figure = plt.figure(figsize=(16, height), facecolor=_PAPER)
    output_path = Path(output_path)
    try:
        counts = tuple(
            transition_frame_counts(
                document,
                figure,
                fps=timing.fps,
                seconds_per_step=timing.seconds_per_step,
            )
            for document in documents
        )
        first_frames = (
            DuplexFrames(duplex_frame_renderer, documents[0])
            if duplex_frame_renderer is not None
            else None
        )
        first_schedule = iter(scene_frame_schedule(counts[0], timing, first=True))
        _draw_frame(documents[0], next(first_schedule), figure, first_frames)
        figure.canvas.draw()
        with (
            staged_output(output_path) as path,
            _saving_writer(writer, figure, path, dpi),
        ):
            writer.grab_frame(facecolor=_PAPER)
            for index, document in enumerate(documents):
                frames = (
                    first_frames
                    if index == 0
                    else (
                        DuplexFrames(duplex_frame_renderer, document)
                        if duplex_frame_renderer is not None
                        else None
                    )
                )
                first_frames = None
                schedule = (
                    first_schedule
                    if index == 0
                    else scene_frame_schedule(
                        counts[index],
                        timing,
                        first=False,
                    )
                )
                for frame in schedule:
                    _draw_frame(document, frame, figure, frames)
                    writer.grab_frame(facecolor=_PAPER)
    finally:
        plt.close(figure)
    return output_path
