"""Finite timing contracts and shared media frame scheduling.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator


def positive_integer(value: int, name: str) -> None:
    """Require an actual positive integer without coercion."""
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        msg = f"{name} must be a positive integer"
        raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class PlaybackTiming:
    """Requested seconds per media phase; zero lead/hold means zero frames."""

    fps: int = 30
    seconds_per_step: float = 0.70
    hold_seconds: float = 0.75
    lead_seconds: float = 0.25
    scene_transition_seconds: float = 0.0

    def __post_init__(self) -> None:
        """Validate exact frame rate and finite duration domains."""
        positive_integer(self.fps, "fps")
        for name in (
            "seconds_per_step",
            "hold_seconds",
            "lead_seconds",
            "scene_transition_seconds",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                msg = f"{name} must be numeric"
                raise TypeError(msg)
            if (
                not math.isfinite(value)
                or value < 0
                or (name == "seconds_per_step" and value == 0)
            ):
                domain = "positive" if name == "seconds_per_step" else "non-negative"
                msg = f"{name} must be finite and {domain}"
                raise ValueError(msg)
            if not math.isfinite(self.fps * value):
                msg = f"{name} produces a nonfinite frame count"
                raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class PlaybackFrame:
    """One requested scene state and optional white transition overlay."""

    transition_index: int
    progress: float
    fade_alpha: float = 0.0


def scene_frame_schedule(
    transition_counts: tuple[int, ...],
    timing: PlaybackTiming,
    *,
    first: bool,
    last: bool,
) -> Iterator[PlaybackFrame]:
    """Yield the same lead, traversal, hold, and fade schedule for every encoder."""
    _validate_transition_counts(transition_counts)
    fade_frames = round(timing.fps * timing.scene_transition_seconds)
    if not first:
        for index in range(fade_frames):
            yield PlaybackFrame(0, 0.0, 1.0 - (index + 1) / fade_frames)
    for _ in range(round(timing.fps * timing.lead_seconds)):
        yield PlaybackFrame(0, 0.0)
    for transition, count in enumerate(transition_counts):
        for index in range(count):
            yield PlaybackFrame(transition, (index + 1) / count)
    final_transition = len(transition_counts) - 1
    for _ in range(round(timing.fps * timing.hold_seconds)):
        yield PlaybackFrame(final_transition, 1.0)
    if not last:
        for index in range(fade_frames):
            yield PlaybackFrame(final_transition, 1.0, (index + 1) / fade_frames)


def _validate_transition_counts(counts: tuple[int, ...]) -> None:
    if not counts:
        msg = "transition_counts must not be empty"
        raise ValueError(msg)
    for count in counts:
        positive_integer(count, "transition frame count")
