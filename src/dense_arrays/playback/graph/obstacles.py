"""Measure node occupancy and curve intersections in scene coordinates.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from .curves import quadratic_point

if TYPE_CHECKING:
    from .model import GraphEdge, GraphScene, QuadraticCurve


def node_bounds(
    scene: GraphScene, node_id: str, inflate: float = 0.0
) -> tuple[float, float, float, float]:
    """Return the node rectangle expanded by the requested clearance."""
    x, y = scene.position(node_id)
    geometry = scene.geometry(node_id)
    return (
        x - geometry.width_pt / 2.0 - inflate,
        y - geometry.height_pt / 2.0 - inflate,
        x + geometry.width_pt / 2.0 + inflate,
        y + geometry.height_pt / 2.0 + inflate,
    )


def contains(
    bounds: tuple[float, float, float, float], point: tuple[float, float]
) -> bool:
    """Test whether a point lies inside or on a rectangle."""
    return bounds[0] <= point[0] <= bounds[2] and bounds[1] <= point[1] <= bounds[3]


def intersects(
    left: tuple[float, float, float, float], right: tuple[float, float, float, float]
) -> bool:
    """Test whether two rectangles overlap with positive area."""
    return not (
        left[2] <= right[0]
        or right[2] <= left[0]
        or left[3] <= right[1]
        or right[3] <= left[1]
    )


def inside_node(scene: GraphScene, node_id: str, point: tuple[float, float]) -> bool:
    """Test circular terminals or rectangular placement bounds."""
    center = scene.position(node_id)
    geometry = scene.geometry(node_id)
    if scene.node(node_id).terminal:
        return math.dist(center, point) <= geometry.width_pt / 2.0
    return (
        abs(point[0] - center[0]) <= geometry.width_pt / 2.0
        and abs(point[1] - center[1]) <= geometry.height_pt / 2.0
    )


def visible_interval(
    scene: GraphScene,
    edge: GraphEdge,
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
) -> tuple[float, float]:
    """Find the centerline interval outside its endpoint nodes."""
    samples = 160
    start_t = None
    for index in range(1, samples + 1):
        candidate = index / samples
        point = quadratic_point(start, control, end, candidate)
        if inside_node(scene, edge.source_id, point):
            continue
        low = (index - 1) / samples
        high = candidate
        for _ in range(14):
            middle = (low + high) / 2.0
            if inside_node(
                scene,
                edge.source_id,
                quadratic_point(start, control, end, middle),
            ):
                low = middle
            else:
                high = middle
        start_t = high
        break
    end_t = None
    for index in range(samples - 1, -1, -1):
        candidate = index / samples
        point = quadratic_point(start, control, end, candidate)
        if inside_node(scene, edge.target_id, point):
            continue
        low = candidate
        high = (index + 1) / samples
        for _ in range(14):
            middle = (low + high) / 2.0
            if inside_node(
                scene,
                edge.target_id,
                quadratic_point(start, control, end, middle),
            ):
                high = middle
            else:
                low = middle
        end_t = low
        break
    if start_t is None or end_t is None or start_t >= end_t:
        msg = (
            f"edge has no visible centerline span: {edge.source_id!r} "
            f"-> {edge.target_id!r}"
        )
        raise ValueError(msg)
    return start_t, end_t


def curve_obstacle_hits(
    scene: GraphScene, edge: GraphEdge, curve: QuadraticCurve
) -> int:
    """Count sampled centerline points that intersect other nodes."""
    excluded = {edge.source_id, edge.target_id}
    obstacles = tuple(
        node_bounds(scene, node.node_id, scene.layout_spec.route_clearance_pt)
        for node in scene.graph.nodes
        if node.node_id not in excluded
    )
    return sum(
        any(
            contains(
                obstacle,
                quadratic_point(
                    curve.motion_start,
                    curve.control,
                    curve.motion_end,
                    index / 80.0,
                ),
            )
            for obstacle in obstacles
        )
        for index in range(1, 80)
    )


def curve_intersects_bounds(
    curve: QuadraticCurve,
    bounds: tuple[float, float, float, float],
) -> bool:
    """Test sampled visible curve points against a rectangle."""
    control = curve.visible_control or curve.control
    return any(
        contains(
            bounds,
            quadratic_point(
                curve.visible_start,
                control,
                curve.visible_end,
                index / 80.0,
            ),
        )
        for index in range(81)
    )


def inside_viewport(
    scene: GraphScene,
    bounds: tuple[float, float, float, float],
    *,
    safety_pt: float = 3.0,
) -> bool:
    """Test whether a rectangle fits within the padded viewport."""
    inset = scene.layout_spec.viewport.padding_pt + safety_pt
    return (
        bounds[0] >= inset
        and bounds[1] >= inset
        and bounds[2] <= scene.layout_spec.viewport.width_pt - inset
        and bounds[3] <= scene.layout_spec.viewport.height_pt - inset
    )
