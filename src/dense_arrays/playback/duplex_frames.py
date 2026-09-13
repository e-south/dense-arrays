"""Validate producer raster frames and adapt them to playback transitions.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import math
from collections import OrderedDict
from typing import TYPE_CHECKING, Protocol

import numpy as np

from .timeline import placement_progress

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from numpy.typing import NDArray

    from .presentation import PlaybackDocument

_DUPLEX_RASTER_OVERSAMPLE = 2.0
_IMAGE_DIMENSIONS = 3
_FRAME_CACHE_SIZE = 2
_SETTLED_PROGRESS = 0.999


class DuplexFrameRenderer(Protocol):
    """Return fixed-shape uint8 RGB/RGBA scenes; None requests complete gray rest."""

    def __call__(
        self, document: PlaybackDocument, step_index: int | None
    ) -> NDArray[np.uint8]:
        """Return the complete resting scene or state after an integer placement."""
        ...


class DuplexFrames:
    """Load and validate scene frames lazily, retaining only two image snapshots."""

    def __init__(
        self, renderer: DuplexFrameRenderer, document: PlaybackDocument
    ) -> None:
        if not callable(renderer):
            msg = "duplex_frame_renderer must be callable"
            raise TypeError(msg)
        self.renderer = renderer
        self.document = document
        self._images: OrderedDict[int | None, NDArray[np.uint8]] = OrderedDict()
        self._shape: tuple[int, ...] | None = None
        owner = getattr(renderer, "__self__", renderer)
        capability = getattr(owner, "renders_distance_brackets", False)
        if not isinstance(capability, bool):
            msg = "renders_distance_brackets must be a boolean"
            raise TypeError(msg)
        self.renders_distance_brackets = capability

    @property
    def native_cap_height_px(self) -> float | None:
        """Read the producer's current nucleotide sizing metric."""
        return _duplex_native_cap_height_px(self.renderer)

    def frame(self, index: int | None) -> NDArray[np.uint8]:
        """Validate the requested image and retain only two snapshots."""
        if index is not None and (
            isinstance(index, bool)
            or not isinstance(index, int)
            or not 0 <= index < len(self.document.plan.steps)
        ):
            msg = "producer frame index is outside the playback steps"
            raise IndexError(msg)
        if index in self._images:
            self._images.move_to_end(index)
            return self._images[index]
        image = self.renderer(self.document, index)
        if not isinstance(image, np.ndarray) or image.dtype != np.uint8:
            msg = "producer frame must be a uint8 ndarray"
            raise TypeError(msg)
        if (
            image.ndim != _IMAGE_DIMENSIONS
            or image.shape[2] not in (3, 4)
            or min(image.shape[:2]) <= 0
        ):
            msg = "producer frame shape must be nonempty (height, width, 3 or 4)"
            raise ValueError(msg)
        if self._shape is not None and image.shape != self._shape:
            msg = "producer frame shape must remain constant within a scene"
            raise ValueError(msg)
        self._shape = image.shape
        snapshot = image.copy()
        snapshot.setflags(write=False)
        self._images[index] = snapshot
        while len(self._images) > _FRAME_CACHE_SIZE:
            self._images.popitem(last=False)
        return snapshot


def duplex_frame_for_axis(
    frames: DuplexFrames,
    document: PlaybackDocument,
    transition_index: int,
    progress: float,
    axis: Axes,
) -> tuple[NDArray[np.uint8], float | None]:
    """Prepare a duplex frame at its final display resolution."""
    import numpy as np
    from PIL import Image

    frame = duplex_transition_frame(
        frames,
        document,
        transition_index,
        progress,
    )
    source_height, source_width = frame.shape[:2]
    axis_bounds = axis.get_window_extent()
    scale = min(
        float(axis_bounds.width) / max(source_width, 1),
        float(axis_bounds.height) / max(source_height, 1),
        1.0,
    )
    raster_scale = min(1.0, scale * _DUPLEX_RASTER_OVERSAMPLE)
    target_width = max(1, round(source_width * raster_scale))
    target_height = max(1, round(source_height * raster_scale))
    native_cap_height = frames.native_cap_height_px
    displayed_cap_height = (
        native_cap_height * scale if native_cap_height is not None else None
    )
    if target_width == source_width and target_height == source_height:
        return frame, displayed_cap_height
    resized = Image.fromarray(frame).resize(
        (target_width, target_height),
        resample=Image.Resampling.LANCZOS,
    )
    return np.asarray(resized), displayed_cap_height


def _duplex_native_cap_height_px(renderer: DuplexFrameRenderer) -> float | None:
    owner = getattr(renderer, "__self__", None)
    value = getattr(owner, "native_nucleotide_cap_height_px", None)
    if value is None:
        return None
    numeric = float(value)
    if isinstance(value, bool) or not math.isfinite(numeric) or numeric <= 0:
        msg = "native_nucleotide_cap_height_px must be finite and positive"
        raise ValueError(msg)
    return numeric


def preferred_figure_height_inches(renderer: DuplexFrameRenderer | None) -> float:
    """Read the producer figure height metric, when declared on its bound owner."""
    owner = getattr(renderer, "__self__", None)
    value = getattr(owner, "preferred_figure_height_inches", 2.4)
    numeric = float(value)
    if isinstance(value, bool) or not 0.0 < numeric < float("inf"):
        msg = "preferred_figure_height_inches must be finite and positive"
        raise ValueError(msg)
    return numeric


def duplex_transition_frame(
    frames: DuplexFrames,
    document: PlaybackDocument,
    transition_index: int,
    progress: float,
) -> NDArray[np.uint8]:
    """Crossfade emphasis between complete scenes without translating glyphs."""
    import numpy as np

    final_index = len(document.plan.steps) - 1
    if transition_index > final_index:
        return frames.frame(final_index)
    previous_index = None if transition_index == 0 else transition_index - 1
    if progress <= 0:
        return frames.frame(previous_index)
    if progress >= _SETTLED_PROGRESS:
        return frames.frame(transition_index)
    previous = np.asarray(frames.frame(previous_index), dtype=np.float32)
    current = np.asarray(frames.frame(transition_index), dtype=np.float32)
    alpha = placement_progress(transition_index, transition_index, progress)
    return np.rint(previous * (1 - alpha) + current * alpha).astype(np.uint8)
