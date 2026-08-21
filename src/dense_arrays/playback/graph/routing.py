"""Obstacle-aware deterministic curve and cost-label routing."""

import itertools
import math
import operator
from functools import lru_cache

from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

from .geometry import KMER_FONT_FAMILY
from .model import GraphEdge, GraphRoutes, GraphScene, QuadraticCurve, RoutedEdge

EDGE_LABEL_FONT_SIZE_PT = 12.0


def quadratic_point(start, control, end, progress: float) -> tuple[float, float]:
    inverse = 1.0 - progress
    return (
        inverse * inverse * start[0]
        + 2.0 * inverse * progress * control[0]
        + progress * progress * end[0],
        inverse * inverse * start[1]
        + 2.0 * inverse * progress * control[1]
        + progress * progress * end[1],
    )


def _quadratic_derivative(start, control, end, progress: float) -> tuple[float, float]:
    inverse = 1.0 - progress
    return (
        2.0 * (inverse * (control[0] - start[0]) + progress * (end[0] - control[0])),
        2.0 * (inverse * (control[1] - start[1]) + progress * (end[1] - control[1])),
    )


def _quadratic_segment_points(start, control, end, start_t: float, end_t: float):
    segment_start = quadratic_point(start, control, end, start_t)
    segment_end = quadratic_point(start, control, end, end_t)
    derivative = _quadratic_derivative(start, control, end, start_t)
    duration = end_t - start_t
    segment_control = (
        segment_start[0] + derivative[0] * duration / 2.0,
        segment_start[1] + derivative[1] * duration / 2.0,
    )
    return segment_start, segment_control, segment_end


def quadratic_segment(curve: QuadraticCurve, start_t: float, end_t: float):
    """Return an exact subcurve of the canonical centerline."""
    start_t = max(0.0, min(1.0, start_t))
    end_t = max(start_t, min(1.0, end_t))
    return _quadratic_segment_points(
        curve.motion_start,
        curve.control,
        curve.motion_end,
        start_t,
        end_t,
    )


@lru_cache(maxsize=512)
def _arc_table(start, control, end, samples: int = 160):
    points = tuple(
        quadratic_point(start, control, end, index / samples)
        for index in range(samples + 1)
    )
    cumulative = [0.0]
    for left, right in itertools.pairwise(points):
        cumulative.append(cumulative[-1] + math.dist(left, right))
    return tuple(cumulative)


def quadratic_arc_length(curve: QuadraticCurve) -> float:
    return _arc_table(curve.motion_start, curve.control, curve.motion_end)[-1]


def quadratic_arc_t(curve: QuadraticCurve, fraction: float) -> float:
    """Map a normalized distance to the canonical curve parameter."""
    fraction = max(0.0, min(1.0, fraction))
    cumulative = _arc_table(curve.motion_start, curve.control, curve.motion_end)
    total = cumulative[-1]
    if total <= 1e-9:
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


def _node_bounds(
    scene: GraphScene, node_id: str, inflate: float = 0.0
) -> tuple[float, float, float, float]:
    x, y = scene.position(node_id)
    geometry = scene.geometry(node_id)
    return (
        x - geometry.width_pt / 2.0 - inflate,
        y - geometry.height_pt / 2.0 - inflate,
        x + geometry.width_pt / 2.0 + inflate,
        y + geometry.height_pt / 2.0 + inflate,
    )


def _contains(bounds, point) -> bool:
    return bounds[0] <= point[0] <= bounds[2] and bounds[1] <= point[1] <= bounds[3]


def _intersects(left, right) -> bool:
    return not (
        left[2] <= right[0]
        or right[2] <= left[0]
        or left[3] <= right[1]
        or right[3] <= left[1]
    )


def _boundary_point(
    scene: GraphScene, node_id: str, center, toward
) -> tuple[float, float]:
    geometry = scene.geometry(node_id)
    node = scene.node(node_id)
    dx, dy = toward[0] - center[0], toward[1] - center[1]
    distance = math.hypot(dx, dy)
    if distance <= 1e-9:
        return center
    if node.terminal:
        scale = geometry.width_pt / 2.0 / distance
    else:
        scale_x = math.inf if abs(dx) <= 1e-9 else geometry.width_pt / 2.0 / abs(dx)
        scale_y = math.inf if abs(dy) <= 1e-9 else geometry.height_pt / 2.0 / abs(dy)
        scale = min(scale_x, scale_y)
    return (center[0] + dx * scale, center[1] + dy * scale)


def _inside_node(scene: GraphScene, node_id: str, point: tuple[float, float]) -> bool:
    center = scene.position(node_id)
    geometry = scene.geometry(node_id)
    if scene.node(node_id).terminal:
        return math.dist(center, point) <= geometry.width_pt / 2.0
    return (
        abs(point[0] - center[0]) <= geometry.width_pt / 2.0
        and abs(point[1] - center[1]) <= geometry.height_pt / 2.0
    )


def _visible_interval(
    scene: GraphScene,
    edge: GraphEdge,
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
) -> tuple[float, float]:
    samples = 160
    start_t = None
    for index in range(1, samples + 1):
        candidate = index / samples
        point = quadratic_point(start, control, end, candidate)
        if _inside_node(scene, edge.source_id, point):
            continue
        low = (index - 1) / samples
        high = candidate
        for _ in range(14):
            middle = (low + high) / 2.0
            if _inside_node(
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
        if _inside_node(scene, edge.target_id, point):
            continue
        low = candidate
        high = (index + 1) / samples
        for _ in range(14):
            middle = (low + high) / 2.0
            if _inside_node(
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
        raise ValueError(
            f"edge has no visible centerline span: {edge.source_id!r} -> {edge.target_id!r}"
        )
    return start_t, end_t


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
    visible_t_start, visible_t_end = _visible_interval(
        scene,
        edge,
        source,
        control,
        target,
    )
    visible_start, visible_control, visible_end = _quadratic_segment_points(
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


def _curve_obstacle_hits(
    scene: GraphScene, edge: GraphEdge, curve: QuadraticCurve
) -> int:
    excluded = {edge.source_id, edge.target_id}
    obstacles = tuple(
        _node_bounds(scene, node.node_id, scene.layout_spec.route_clearance_pt)
        for node in scene.graph.nodes
        if node.node_id not in excluded
    )
    return sum(
        any(
            _contains(
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


def _preferred_sign(edge: GraphEdge) -> float:
    key = f"{edge.source_id}>{edge.target_id}"
    return (
        1.0
        if sum((index + 1) * ord(character) for index, character in enumerate(key)) % 2
        == 0
        else -1.0
    )


def _route_edge(
    scene: GraphScene,
    edge: GraphEdge,
    *,
    allow_masked: bool = False,
) -> QuadraticCurve:
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
        obstacle_hits = _curve_obstacle_hits(scene, edge, curve)
        if obstacle_hits == 0:
            return curve
        if best_masked is None or obstacle_hits < best_masked[0]:
            best_masked = (obstacle_hits, curve)
    if allow_masked and best_masked is not None:
        return best_masked[1]
    raise ValueError(
        f"no collision-free route for edge {edge.source_id!r} -> {edge.target_id!r}"
    )


def _label_extent(
    text: str, font_size: float = EDGE_LABEL_FONT_SIZE_PT
) -> tuple[float, float]:
    bounds = TextPath(
        (0, 0), text, prop=FontProperties(family=KMER_FONT_FAMILY, size=font_size)
    ).get_extents()
    return float(bounds.width) + 4.0, float(bounds.height) + 3.0


def _place_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied_labels: list[tuple[float, float, float, float]],
):
    width, height = _label_extent(str(edge.added_bases))
    node_bounds = tuple(
        _node_bounds(scene, node.node_id, 1.5) for node in scene.graph.nodes
    )
    for progress in (0.50, 0.40, 0.60, 0.30, 0.70, 0.22, 0.78):
        visible_control = curve.visible_control or curve.control
        point = quadratic_point(
            curve.visible_start, visible_control, curve.visible_end, progress
        )
        tx = 2.0 * (1.0 - progress) * (
            visible_control[0] - curve.visible_start[0]
        ) + 2.0 * progress * (curve.visible_end[0] - visible_control[0])
        ty = 2.0 * (1.0 - progress) * (
            visible_control[1] - curve.visible_start[1]
        ) + 2.0 * progress * (curve.visible_end[1] - visible_control[1])
        length = max(math.hypot(tx, ty), 1e-9)
        normal = (-ty / length, tx / length)
        for offset in (0.0, 7.0, -7.0, 12.0, -12.0, 18.0, -18.0):
            center = (point[0] + normal[0] * offset, point[1] + normal[1] * offset)
            bounds = (
                center[0] - width / 2.0,
                center[1] - height / 2.0,
                center[0] + width / 2.0,
                center[1] + height / 2.0,
            )
            if any(_intersects(bounds, obstacle) for obstacle in node_bounds) or any(
                _intersects(bounds, occupied) for occupied in occupied_labels
            ):
                continue
            return center, bounds
    raise ValueError(
        f"no collision-free label position for edge {edge.source_id!r} -> {edge.target_id!r}"
    )


def _curve_intersects_bounds(
    curve: QuadraticCurve,
    bounds: tuple[float, float, float, float],
) -> bool:
    control = curve.visible_control or curve.control
    return any(
        _contains(
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


def _label_inside_viewport(
    scene: GraphScene,
    bounds: tuple[float, float, float, float],
    *,
    safety_pt: float = 3.0,
) -> bool:
    inset = scene.layout_spec.viewport.padding_pt + safety_pt
    return (
        bounds[0] >= inset
        and bounds[1] >= inset
        and bounds[2] <= scene.layout_spec.viewport.width_pt - inset
        and bounds[3] <= scene.layout_spec.viewport.height_pt - inset
    )


def _place_edge_aware_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied_labels: list[tuple[float, float, float, float]],
    other_curves: tuple[QuadraticCurve, ...],
    font_size: float,
):
    width, height = _label_extent(str(edge.added_bases), font_size)
    node_bounds = tuple(
        _node_bounds(scene, node.node_id, 1.5) for node in scene.graph.nodes
    )
    visible_control = curve.visible_control or curve.control
    best_candidate = None
    for progress_index, progress in enumerate(
        (0.50, 0.38, 0.62, 0.28, 0.72, 0.20, 0.80, 0.12, 0.88)
    ):
        point = quadratic_point(
            curve.visible_start, visible_control, curve.visible_end, progress
        )
        tx, ty = _quadratic_derivative(
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
            if any(_intersects(bounds, obstacle) for obstacle in node_bounds):
                continue
            if any(_intersects(bounds, occupied) for occupied in occupied_labels):
                continue
            curve_crossings = sum(
                _curve_intersects_bounds(other, expanded) for other in other_curves
            )
            rank = (curve_crossings, abs(offset), progress_index, offset_index)
            if best_candidate is None or rank < best_candidate[0]:
                best_candidate = (rank, center, bounds, point)
            if curve_crossings == 0:
                return center, bounds, point
    if best_candidate is not None:
        return best_candidate[1], best_candidate[2], best_candidate[3]
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
            if any(_intersects(bounds, obstacle) for obstacle in node_bounds):
                continue
            if any(_intersects(bounds, occupied) for occupied in occupied_labels):
                continue
            expanded = (
                bounds[0] - 1.5,
                bounds[1] - 1.5,
                bounds[2] + 1.5,
                bounds[3] + 1.5,
            )
            curve_crossings = sum(
                _curve_intersects_bounds(other, expanded) for other in other_curves
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
    raise ValueError(
        f"no edge-owned label position for edge {edge.source_id!r} -> {edge.target_id!r}"
    )


def _place_leader_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied_labels: list[tuple[float, float, float, float]],
    other_curves: tuple[QuadraticCurve, ...],
    font_size: float,
):
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
            if any(_intersects(bounds, occupied) for occupied in occupied_labels):
                continue
            if any(
                _intersects(bounds, _node_bounds(scene, node.node_id, inflate=2.0))
                for node in scene.graph.nodes
            ):
                continue
            crossings = sum(
                _curve_intersects_bounds(other, bounds) for other in other_curves
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
        raise ValueError(
            f"no leader-label position for edge {edge.source_id!r} -> {edge.target_id!r}"
        )
    candidates.sort(key=operator.itemgetter(0))
    _, center, bounds = candidates[0]
    return center, bounds, anchor


@lru_cache(maxsize=128)
def route_graph_scene(scene: GraphScene) -> GraphRoutes:
    traversal_curves = []
    for edge in scene.graph.traversal_edges:
        curve = _route_edge(scene, edge, allow_masked=True)
        traversal_curves.append((edge, curve))
    context = []
    for edge in scene.display_context_edges:
        try:
            curve = _route_edge(
                scene,
                edge,
                allow_masked=edge.relation_kind == "declared_constraint",
            )
        except ValueError:
            # Context edges are a presentation subset of the complete semantic
            # relation set. An unroutable context curve may be omitted without
            # weakening or inventing the realized traversal claim.
            continue
        context.append(RoutedEdge(edge, curve))
    all_curves = tuple(curve for _, curve in traversal_curves) + tuple(
        routed.curve for routed in context
    )
    occupied_labels = []
    traversal = []
    for edge_index, (edge, curve) in enumerate(traversal_curves):
        label_position = label_bounds = None
        label_font_size = None
        label_anchor = None
        if edge.added_bases is not None:
            other_curves = all_curves[:edge_index] + all_curves[edge_index + 1 :]
            for candidate_font_size in (EDGE_LABEL_FONT_SIZE_PT, 11.0, 10.0, 9.0):
                try:
                    label_position, label_bounds, label_anchor = (
                        _place_edge_aware_label(
                            scene,
                            edge,
                            curve,
                            occupied_labels,
                            other_curves,
                            candidate_font_size,
                        )
                    )
                except ValueError:
                    continue
                if not _label_inside_viewport(scene, label_bounds):
                    label_position = None
                    label_bounds = None
                    label_anchor = None
                    continue
                label_font_size = candidate_font_size
                break
            if label_position is None or label_bounds is None:
                for candidate_font_size in (EDGE_LABEL_FONT_SIZE_PT, 11.0, 10.0, 9.0):
                    try:
                        label_position, label_bounds, label_anchor = (
                            _place_leader_label(
                                scene,
                                edge,
                                curve,
                                occupied_labels,
                                other_curves,
                                candidate_font_size,
                            )
                        )
                    except ValueError:
                        continue
                    label_font_size = candidate_font_size
                    break
            if label_position is None or label_bounds is None:
                raise ValueError(
                    f"no collision-free cost label for edge {edge.source_id!r} -> {edge.target_id!r}"
                )
            occupied_labels.append(label_bounds)
        traversal.append(
            RoutedEdge(
                edge,
                curve,
                label_position,
                label_bounds,
                label_font_size,
                label_anchor,
            )
        )
    return GraphRoutes(tuple(context), tuple(traversal))


def edge_curve(scene: GraphScene, edge: GraphEdge) -> QuadraticCurve:
    return _route_edge(scene, edge)
