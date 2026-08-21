"""Dependency-free reference positions for playback renderers."""

from __future__ import annotations

import math


def radial_path_positions(count: int) -> tuple[tuple[float, float], ...]:
    """Return deterministic radial positions for ``count`` path nodes."""
    if count < 0:
        raise ValueError("count must be non-negative")
    if count == 0:
        return ()
    if count == 1:
        return ((0.5, 0.5),)
    if count == 2:
        return ((0.5, 0.33), (0.5, 0.67))
    return tuple(
        (
            0.5 + 0.30 * math.cos((-math.pi / 2.0) + 2.0 * math.pi * index / count),
            0.5 + 0.30 * math.sin((-math.pi / 2.0) + 2.0 * math.pi * index / count),
        )
        for index in range(count)
    )


def journey_path_positions(count: int) -> tuple[tuple[float, float], ...]:
    """Return compact deterministic positions for a path explanation."""
    if count < 0:
        raise ValueError("count must be non-negative")
    if count == 0:
        return ()
    if count == 1:
        return ((0.5, 0.5),)
    rows = 2 if count <= 8 else 3
    columns = (count + rows - 1) // rows
    lanes = (0.79, 0.21) if rows == 2 else (0.82, 0.50, 0.18)
    return tuple(
        (
            0.18 if columns == 1 else 0.18 + 0.64 * (index // rows) / (columns - 1),
            lanes[index % rows],
        )
        for index in range(count)
    )
