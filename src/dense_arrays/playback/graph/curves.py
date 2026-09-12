"""Evaluate quadratic centerlines, subcurves, and arc-length timing.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import itertools
import math
from functools import lru_cache
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .model import QuadraticCurve


_ARC_EPSILON = 1e-9


def quadratic_point(
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
    progress: float,
) -> tuple[float, float]:
    """Evaluate a quadratic centerline at a normalized parameter."""
    inverse = 1.0 - progress
    return (
        inverse * inverse * start[0]
        + 2.0 * inverse * progress * control[0]
        + progress * progress * end[0],
        inverse * inverse * start[1]
        + 2.0 * inverse * progress * control[1]
        + progress * progress * end[1],
    )


def quadratic_derivative(
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
    progress: float,
) -> tuple[float, float]:
    """Return the tangent vector at a normalized curve parameter."""
    inverse = 1.0 - progress
    return (
        2.0 * (inverse * (control[0] - start[0]) + progress * (end[0] - control[0])),
        2.0 * (inverse * (control[1] - start[1]) + progress * (end[1] - control[1])),
    )


def quadratic_segment_points(
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
    start_t: float,
    end_t: float,
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    """Return the three control points for an exact quadratic subcurve."""
    segment_start = quadratic_point(start, control, end, start_t)
    segment_end = quadratic_point(start, control, end, end_t)
    derivative = quadratic_derivative(start, control, end, start_t)
    duration = end_t - start_t
    segment_control = (
        segment_start[0] + derivative[0] * duration / 2.0,
        segment_start[1] + derivative[1] * duration / 2.0,
    )
    return segment_start, segment_control, segment_end


def quadratic_segment(
    curve: QuadraticCurve, start_t: float, end_t: float
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    """Return an exact subcurve of the canonical centerline."""
    start_t = max(0.0, min(1.0, start_t))
    end_t = max(start_t, min(1.0, end_t))
    return quadratic_segment_points(
        curve.motion_start,
        curve.control,
        curve.motion_end,
        start_t,
        end_t,
    )


@lru_cache(maxsize=512)
def _arc_table(
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
    samples: int = 160,
) -> tuple[float, ...]:
    points = tuple(
        quadratic_point(start, control, end, index / samples)
        for index in range(samples + 1)
    )
    cumulative = [0.0]
    for left, right in itertools.pairwise(points):
        cumulative.append(cumulative[-1] + math.dist(left, right))
    return tuple(cumulative)


def quadratic_arc_length(curve: QuadraticCurve) -> float:
    """Return sampled centerline length in scene units."""
    return _arc_table(curve.motion_start, curve.control, curve.motion_end)[-1]


def quadratic_arc_t(curve: QuadraticCurve, fraction: float) -> float:
    """Map a normalized distance to the canonical curve parameter."""
    fraction = max(0.0, min(1.0, fraction))
    cumulative = _arc_table(curve.motion_start, curve.control, curve.motion_end)
    total = cumulative[-1]
    if total <= _ARC_EPSILON:
        return fraction
    target = total * fraction
    for index, distance in enumerate(cumulative[1:], start=1):
        if distance < target:
            continue
        previous = cumulative[index - 1]
        span = max(distance - previous, 1e-9)
        local = (target - previous) / span
        samples = len(cumulative) - 1
        return ((index - 1) + local) / samples
    return 1.0
