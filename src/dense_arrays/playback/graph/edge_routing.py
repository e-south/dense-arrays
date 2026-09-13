"""Choose deterministic edge centerlines around measured node obstacles.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import math

from .curves import quadratic_segment_points
from .model import GraphEdge, GraphScene, QuadraticCurve
from .obstacles import curve_obstacle_hits, visible_interval


def _curve_for_offset(
    scene: GraphScene, edge: GraphEdge, offset: float
) -> QuadraticCurve:
    source, target = scene.position(edge.source_id), scene.position(edge.target_id)
    dx, dy = target[0] - source[0], target[1] - source[1]
    distance = max(math.hypot(dx, dy), 1e-9)
    control = (
        (source[0] + target[0]) / 2.0 - dy / distance * offset,
        (source[1] + target[1]) / 2.0 + dx / distance * offset,
    )
    visible_t_start, visible_t_end = visible_interval(
        scene,
        edge,
        source,
        control,
        target,
    )
    visible_start, visible_control, visible_end = quadratic_segment_points(
        source,
        control,
        target,
        visible_t_start,
        visible_t_end,
    )
    return QuadraticCurve(
        visible_start,
        control,
        visible_end,
        source,
        target,
        visible_control=visible_control,
        visible_t_start=visible_t_start,
        visible_t_end=visible_t_end,
    )


def _preferred_sign(edge: GraphEdge) -> float:
    key = f"{edge.source_id}>{edge.target_id}"
    return (
        1.0
        if sum((index + 1) * ord(character) for index, character in enumerate(key)) % 2
        == 0
        else -1.0
    )


def route_edge(
    scene: GraphScene,
    edge: GraphEdge,
    *,
    allow_masked: bool = False,
) -> QuadraticCurve:
    """Search deterministic curve offsets for a route around node obstacles."""
    sign = _preferred_sign(edge)
    offsets = (
        sign * 10.0,
        -sign * 10.0,
        0.0,
        sign * 20.0,
        -sign * 20.0,
        sign * 32.0,
        -sign * 32.0,
        sign * 46.0,
        -sign * 46.0,
        sign * 62.0,
        -sign * 62.0,
    )
    best_masked: tuple[int, QuadraticCurve] | None = None
    for offset in offsets:
        curve = _curve_for_offset(scene, edge, offset)
        obstacle_hits = curve_obstacle_hits(scene, edge, curve)
        if obstacle_hits == 0:
            return curve
        if best_masked is None or obstacle_hits < best_masked[0]:
            best_masked = (obstacle_hits, curve)
    if allow_masked and best_masked is not None:
        return best_masked[1]
    msg = f"no collision-free route for edge {edge.source_id!r} -> {edge.target_id!r}"
    raise ValueError(msg)
