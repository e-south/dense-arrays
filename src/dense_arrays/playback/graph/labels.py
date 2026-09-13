"""Place edge-owned cost labels while avoiding nodes and other labels.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass
from typing import TYPE_CHECKING

from .curves import quadratic_derivative, quadratic_point
from .geometry import KMER_FONT_FAMILY
from .obstacles import curve_intersects_bounds, intersects, node_bounds

if TYPE_CHECKING:
    from .model import GraphEdge, GraphScene, QuadraticCurve

EDGE_LABEL_FONT_SIZE_PT = 12.0


@dataclass(frozen=True, slots=True)
class LabelObstacles:
    """Occupied label rectangles and other curves to avoid during placement."""

    labels: tuple[tuple[float, float, float, float], ...]
    curves: tuple[QuadraticCurve, ...]


def _label_extent(
    text: str, font_size: float = EDGE_LABEL_FONT_SIZE_PT
) -> tuple[float, float]:
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    bounds = TextPath(
        (0, 0), text, prop=FontProperties(family=KMER_FONT_FAMILY, size=font_size)
    ).get_extents()
    return float(bounds.width) + 4.0, float(bounds.height) + 3.0


def place_edge_aware_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied: LabelObstacles,
    font_size: float,
) -> tuple[tuple[float, float], tuple[float, float, float, float], tuple[float, float]]:
    """Find an edge-owned label near the curve, then in adjacent free space."""
    local = _near_curve_label(scene, edge, curve, occupied, font_size)
    if local is not None:
        return local
    return _offset_label(scene, edge, curve, occupied, font_size)


def _near_curve_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied: LabelObstacles,
    font_size: float,
) -> (
    tuple[tuple[float, float], tuple[float, float, float, float], tuple[float, float]]
    | None
):
    width, height = _label_extent(str(edge.added_bases), font_size)
    obstacles = tuple(
        node_bounds(scene, node.node_id, 1.5) for node in scene.graph.nodes
    )
    visible_control = curve.visible_control or curve.control
    best_candidate = None
    for progress_index, progress in enumerate(
        (0.50, 0.38, 0.62, 0.28, 0.72, 0.20, 0.80, 0.12, 0.88)
    ):
        point = quadratic_point(
            curve.visible_start, visible_control, curve.visible_end, progress
        )
        tx, ty = quadratic_derivative(
            curve.visible_start,
            visible_control,
            curve.visible_end,
            progress,
        )
        length = max(math.hypot(tx, ty), 1e-9)
        normal = (-ty / length, tx / length)
        for offset_index, offset in enumerate(
            (
                0.0,
                5.0,
                -5.0,
                9.0,
                -9.0,
                13.0,
                -13.0,
                18.0,
                -18.0,
                24.0,
                -24.0,
                30.0,
                -30.0,
            )
        ):
            center = (
                point[0] + normal[0] * offset,
                point[1] + normal[1] * offset,
            )
            bounds = (
                center[0] - width / 2.0,
                center[1] - height / 2.0,
                center[0] + width / 2.0,
                center[1] + height / 2.0,
            )
            expanded = (
                bounds[0] - 1.5,
                bounds[1] - 1.5,
                bounds[2] + 1.5,
                bounds[3] + 1.5,
            )
            if any(intersects(bounds, obstacle) for obstacle in obstacles):
                continue
            if any(intersects(bounds, bounds_used) for bounds_used in occupied.labels):
                continue
            curve_crossings = sum(
                curve_intersects_bounds(other, expanded) for other in occupied.curves
            )
            rank = (curve_crossings, abs(offset), progress_index, offset_index)
            if best_candidate is None or rank < best_candidate[0]:
                best_candidate = (rank, center, bounds, point)
            if curve_crossings == 0:
                return center, bounds, point
    if best_candidate is not None:
        return best_candidate[1], best_candidate[2], best_candidate[3]
    return None


def _offset_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied: LabelObstacles,
    font_size: float,
) -> tuple[tuple[float, float], tuple[float, float, float, float], tuple[float, float]]:
    width, height = _label_extent(str(edge.added_bases), font_size)
    obstacles = tuple(
        node_bounds(scene, node.node_id, 1.5) for node in scene.graph.nodes
    )
    visible_control = curve.visible_control or curve.control
    search_offsets = (
        (0.0, -18.0),
        (0.0, 18.0),
        (-18.0, 0.0),
        (18.0, 0.0),
        (-24.0, -18.0),
        (24.0, -18.0),
        (-24.0, 18.0),
        (24.0, 18.0),
        (0.0, -30.0),
        (0.0, 30.0),
        (-36.0, 0.0),
        (36.0, 0.0),
    )
    fallback_candidate = None
    for progress_index, progress in enumerate((0.50, 0.35, 0.65, 0.20, 0.80)):
        anchor = quadratic_point(
            curve.visible_start,
            visible_control,
            curve.visible_end,
            progress,
        )
        for offset_index, (dx, dy) in enumerate(search_offsets):
            center = (anchor[0] + dx, anchor[1] + dy)
            bounds = (
                center[0] - width / 2.0,
                center[1] - height / 2.0,
                center[0] + width / 2.0,
                center[1] + height / 2.0,
            )
            if (
                bounds[0] < 1.0
                or bounds[1] < 1.0
                or bounds[2] > scene.layout_spec.viewport.width_pt - 1.0
                or bounds[3] > scene.layout_spec.viewport.height_pt - 1.0
            ):
                continue
            if any(intersects(bounds, obstacle) for obstacle in obstacles):
                continue
            if any(intersects(bounds, bounds_used) for bounds_used in occupied.labels):
                continue
            expanded = (
                bounds[0] - 1.5,
                bounds[1] - 1.5,
                bounds[2] + 1.5,
                bounds[3] + 1.5,
            )
            curve_crossings = sum(
                curve_intersects_bounds(other, expanded) for other in occupied.curves
            )
            rank = (
                curve_crossings,
                math.hypot(dx, dy),
                progress_index,
                offset_index,
            )
            if fallback_candidate is None or rank < fallback_candidate[0]:
                fallback_candidate = (rank, center, bounds, anchor)
            if curve_crossings == 0:
                return center, bounds, anchor
    if fallback_candidate is not None:
        return (
            fallback_candidate[1],
            fallback_candidate[2],
            fallback_candidate[3],
        )
    msg = (
        f"no edge-owned label position for edge {edge.source_id!r} "
        f"-> {edge.target_id!r}"
    )
    raise ValueError(msg)


def place_leader_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied: LabelObstacles,
    font_size: float,
) -> tuple[tuple[float, float], tuple[float, float, float, float], tuple[float, float]]:
    """Place a label in free scene space and bind it to its edge with a leader."""
    text = str(edge.added_bases)
    width, height = _label_extent(text, font_size)
    padding = scene.layout_spec.viewport.padding_pt + 3.0
    viewport_width = scene.layout_spec.viewport.width_pt
    viewport_height = scene.layout_spec.viewport.height_pt
    anchor = quadratic_point(
        curve.motion_start,
        curve.control,
        curve.motion_end,
        0.5,
    )
    x_min = padding + width / 2.0
    x_max = viewport_width - padding - width / 2.0
    y_min = padding + height / 2.0
    y_max = viewport_height - padding - height / 2.0
    candidates = []
    for y_index in range(9):
        y_value = y_min + (y_max - y_min) * y_index / 8.0
        for x_index in range(19):
            x_value = x_min + (x_max - x_min) * x_index / 18.0
            bounds = (
                x_value - width / 2.0,
                y_value - height / 2.0,
                x_value + width / 2.0,
                y_value + height / 2.0,
            )
            if any(intersects(bounds, bounds_used) for bounds_used in occupied.labels):
                continue
            if any(
                intersects(bounds, node_bounds(scene, node.node_id, inflate=2.0))
                for node in scene.graph.nodes
            ):
                continue
            crossings = sum(
                curve_intersects_bounds(other, bounds) for other in occupied.curves
            )
            distance = math.hypot(x_value - anchor[0], y_value - anchor[1])
            rank = (
                distance + crossings * 6.0,
                crossings,
                abs(y_value - viewport_height / 2.0),
                x_index,
                y_index,
            )
            candidates.append((rank, (x_value, y_value), bounds))
    if not candidates:
        msg = (
            f"no leader-label position for edge {edge.source_id!r} "
            f"-> {edge.target_id!r}"
        )
        raise ValueError(msg)
    candidates.sort(key=operator.itemgetter(0))
    _, center, bounds = candidates[0]
    return center, bounds, anchor
