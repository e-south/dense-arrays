"""Deterministic, symmetry-scored layout for publication playback graphs."""

from __future__ import annotations

import itertools
import math
import operator
from dataclasses import dataclass

from .layout_candidates import generate_layout_candidates, internal_edge_pairs
from .model import (
    END_NODE_ID,
    START_NODE_ID,
    ExplanationGraph,
    GraphLayoutSpec,
    GraphNode,
    GraphPosition,
    GraphScene,
)


def _rotate(
    raw: dict[str, tuple[float, float]], angle: float, reflected: bool
) -> dict[str, tuple[float, float]]:
    center_x = sum(point[0] for point in raw.values()) / len(raw)
    center_y = sum(point[1] for point in raw.values()) / len(raw)
    cosine, sine = math.cos(angle), math.sin(angle)
    rotated: dict[str, tuple[float, float]] = {}
    for node_id, (x_value, y_value) in raw.items():
        centered_x = x_value - center_x
        centered_y = y_value - center_y
        centered_x = -centered_x if reflected else centered_x
        rotated[node_id] = (
            cosine * centered_x - sine * centered_y,
            sine * centered_x + cosine * centered_y,
        )
    return rotated


def _internal_box(spec: GraphLayoutSpec) -> tuple[float, float, float, float]:
    start_width = spec.geometry(START_NODE_ID).width_pt
    end_width = spec.geometry(END_NODE_ID).width_pt
    return (
        spec.viewport.padding_pt + start_width + spec.terminal_gap_pt,
        spec.viewport.width_pt
        - spec.viewport.padding_pt
        - end_width
        - spec.terminal_gap_pt,
        spec.viewport.padding_pt,
        spec.viewport.height_pt - spec.viewport.padding_pt,
    )


def _fit_isotropic(
    raw: dict[str, tuple[float, float]],
    internal: tuple[GraphNode, ...],
    spec: GraphLayoutSpec,
) -> dict[str, list[float]]:
    left, right, bottom, top = _internal_box(spec)
    max_half_width = max(node.width for node in internal) / 2.0
    max_half_height = max(node.height for node in internal) / 2.0
    center_left, center_right = left + max_half_width, right - max_half_width
    center_bottom, center_top = bottom + max_half_height, top - max_half_height
    if center_left >= center_right or center_bottom >= center_top:
        raise ValueError("graph viewport is too small for measured node geometry")

    raw_left = min(point[0] for point in raw.values())
    raw_right = max(point[0] for point in raw.values())
    raw_bottom = min(point[1] for point in raw.values())
    raw_top = max(point[1] for point in raw.values())
    raw_width = max(raw_right - raw_left, 1e-9)
    raw_height = max(raw_top - raw_bottom, 1e-9)
    scale = 0.92 * min(
        (center_right - center_left) / raw_width,
        (center_top - center_bottom) / raw_height,
    )
    raw_center_x = (raw_left + raw_right) / 2.0
    raw_center_y = (raw_bottom + raw_top) / 2.0
    target_x = (center_left + center_right) / 2.0
    target_y = (center_bottom + center_top) / 2.0
    return {
        node_id: [
            target_x + (point[0] - raw_center_x) * scale,
            target_y + (point[1] - raw_center_y) * scale,
        ]
        for node_id, point in raw.items()
    }


def _symmetric_collision_relaxation(
    positions: dict[str, list[float]],
    internal: tuple[GraphNode, ...],
    spec: GraphLayoutSpec,
) -> bool:
    left, right, bottom, top = _internal_box(spec)
    target_x = sum(point[0] for point in positions.values()) / len(positions)
    target_y = sum(point[1] for point in positions.values()) / len(positions)
    nodes = tuple(sorted(internal, key=lambda node: node.node_id))
    for _ in range(280):
        moved = False
        for index, first in enumerate(nodes):
            for second in nodes[index + 1 :]:
                first_position = positions[first.node_id]
                second_position = positions[second.node_id]
                delta_x = second_position[0] - first_position[0]
                delta_y = second_position[1] - first_position[1]
                required_x = (first.width + second.width) / 2.0 + spec.node_clearance_pt
                required_y = (
                    first.height + second.height
                ) / 2.0 + spec.node_clearance_pt
                overlap_x = required_x - abs(delta_x)
                overlap_y = required_y - abs(delta_y)
                if overlap_x <= 0.0 or overlap_y <= 0.0:
                    continue
                moved = True
                if overlap_x / required_x < overlap_y / required_y:
                    direction = 1.0 if delta_x > 0 else -1.0
                    if abs(delta_x) < 1e-9:
                        direction = 1.0 if first.node_id < second.node_id else -1.0
                    displacement = overlap_x / 2.0 + 0.05
                    first_position[0] -= direction * displacement
                    second_position[0] += direction * displacement
                else:
                    direction = 1.0 if delta_y > 0 else -1.0
                    if abs(delta_y) < 1e-9:
                        direction = 1.0 if first.node_id < second.node_id else -1.0
                    displacement = overlap_y / 2.0 + 0.05
                    first_position[1] -= direction * displacement
                    second_position[1] += direction * displacement

        current_x = sum(point[0] for point in positions.values()) / len(positions)
        current_y = sum(point[1] for point in positions.values()) / len(positions)
        for node in nodes:
            point = positions[node.node_id]
            point[0] += target_x - current_x
            point[1] += target_y - current_y
            point[0] = min(
                right - node.width / 2.0, max(left + node.width / 2.0, point[0])
            )
            point[1] = min(
                top - node.height / 2.0, max(bottom + node.height / 2.0, point[1])
            )
        if not moved:
            return True
    return False


def _segment_crosses(
    first: tuple[tuple[float, float], tuple[float, float]],
    second: tuple[tuple[float, float], tuple[float, float]],
) -> bool:
    def orientation(
        a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]
    ) -> float:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    a, b = first
    c, d = second
    return (
        orientation(a, b, c) * orientation(a, b, d) < 0
        and orientation(c, d, a) * orientation(c, d, b) < 0
    )


def _score(
    positions: dict[str, list[float]],
    internal: tuple[GraphNode, ...],
    graph: ExplanationGraph,
    spec: GraphLayoutSpec,
) -> float:
    geometry = {node.node_id: node for node in internal}
    left = min(positions[node.node_id][0] - node.width / 2.0 for node in internal)
    right = max(positions[node.node_id][0] + node.width / 2.0 for node in internal)
    bottom = min(positions[node.node_id][1] - node.height / 2.0 for node in internal)
    top = max(positions[node.node_id][1] + node.height / 2.0 for node in internal)
    width, height = max(right - left, 1e-9), max(top - bottom, 1e-9)
    center_x, center_y = (left + right) / 2.0, (bottom + top) / 2.0
    total_area = sum(node.width * node.height for node in internal)
    upper_area = sum(
        node.width * node.height
        for node in internal
        if positions[node.node_id][1] >= center_y
    )
    right_area = sum(
        node.width * node.height
        for node in internal
        if positions[node.node_id][0] >= center_x
    )
    mass_imbalance = (
        abs(2.0 * upper_area - total_area) + abs(2.0 * right_area - total_area)
    ) / max(total_area, 1e-9)

    aspect_target = min(
        1.45,
        max(
            1.0,
            math.sqrt(
                sum(node.width / max(node.height, 1e-9) for node in internal)
                / len(internal)
            ),
        ),
    )
    aspect_penalty = abs(math.log((width / height) / aspect_target))
    pairs = internal_edge_pairs(graph)
    crossings = 0
    for index, first_pair in enumerate(pairs):
        first_nodes = set(first_pair)
        first_segment = (
            tuple(positions[first_pair[0]]),
            tuple(positions[first_pair[1]]),
        )
        for second_pair in pairs[index + 1 :]:
            if first_nodes.intersection(second_pair):
                continue
            second_segment = (
                tuple(positions[second_pair[0]]),
                tuple(positions[second_pair[1]]),
            )
            crossings += int(_segment_crosses(first_segment, second_segment))

    edge_lengths = [
        math.dist(positions[source], positions[target]) for source, target in pairs
    ]
    mean_edge = sum(edge_lengths) / max(len(edge_lengths), 1)
    edge_variation = (
        sum(abs(length - mean_edge) for length in edge_lengths)
        / max(len(edge_lengths) * mean_edge, 1e-9)
        if edge_lengths
        else 0.0
    )
    available_left, available_right, available_bottom, available_top = _internal_box(
        spec
    )
    occupied_fraction = (width * height) / max(
        (available_right - available_left) * (available_top - available_bottom), 1e-9
    )

    first_id = next(
        (
            edge.target_id
            for edge in graph.traversal_edges
            if edge.source_id == START_NODE_ID and edge.target_id in geometry
        ),
        None,
    )
    last_id = next(
        (
            edge.source_id
            for edge in graph.traversal_edges
            if edge.target_id == END_NODE_ID and edge.source_id in geometry
        ),
        None,
    )
    attachment_penalty = 0.0
    if first_id is not None:
        attachment_penalty += (positions[first_id][0] - left) / width
    if last_id is not None:
        attachment_penalty += (right - positions[last_id][0]) / width

    ordered_path = sorted(
        (node for node in internal if node.step_index is not None),
        key=lambda node: node.step_index,
    )
    vertical_deltas = [
        positions[right.node_id][1] - positions[left.node_id][1]
        for left, right in itertools.pairwise(ordered_path)
    ]
    switchbacks = sum(
        1 for left, right in itertools.pairwise(vertical_deltas) if left * right < 0.0
    )
    return (
        crossings * 3.5
        + max(0, switchbacks - 1) * 1.8
        + switchbacks * 0.15
        + mass_imbalance * 2.4
        + aspect_penalty * 1.6
        + edge_variation * 0.45
        + attachment_penalty * 0.25
        + max(0.0, 0.42 - occupied_fraction) * 1.5
    )


def _add_terminals(
    positions: dict[str, list[float]],
    internal: tuple[GraphNode, ...],
    spec: GraphLayoutSpec,
) -> None:
    internal_left = min(
        positions[node.node_id][0] - node.width / 2.0 for node in internal
    )
    internal_right = max(
        positions[node.node_id][0] + node.width / 2.0 for node in internal
    )
    terminal_y = spec.viewport.height_pt / 2.0
    start_radius = spec.geometry(START_NODE_ID).width_pt / 2.0
    end_radius = spec.geometry(END_NODE_ID).width_pt / 2.0
    positions[START_NODE_ID] = [
        internal_left - spec.terminal_gap_pt - start_radius,
        terminal_y,
    ]
    positions[END_NODE_ID] = [
        internal_right + spec.terminal_gap_pt + end_radius,
        terminal_y,
    ]
    occupied_left = positions[START_NODE_ID][0] - start_radius
    occupied_right = positions[END_NODE_ID][0] + end_radius
    desired_center = spec.viewport.width_pt / 2.0
    shift_x = desired_center - (occupied_left + occupied_right) / 2.0
    for point in positions.values():
        point[0] += shift_x

    horizontal_inset = max(spec.viewport.padding_pt, 15.0)
    available_width = spec.viewport.width_pt - (2.0 * horizontal_inset)
    occupied_width = occupied_right - occupied_left
    horizontal_scale = min(1.22, available_width / max(occupied_width, 1e-9))
    if horizontal_scale > 1.0:
        for point in positions.values():
            point[0] = desired_center + ((point[0] - desired_center) * horizontal_scale)


@dataclass(frozen=True, slots=True)
class NetworkXIsotropicLayout:
    name: str = "networkx_isotropic_candidates"
    seed_count: int = 3
    orientation_count: int = 12
    arf_max_iter: int = 1200
    forceatlas_max_iter: int = 650

    def layout(
        self,
        graph: ExplanationGraph,
        spec: GraphLayoutSpec,
        *,
        seed: int,
    ) -> tuple[GraphPosition, ...]:
        internal = tuple(node for node in graph.nodes if not node.terminal)
        if not internal:
            raise ValueError("graph layout requires at least one non-terminal node")
        candidates = generate_layout_candidates(
            graph,
            internal,
            seed=seed,
            seed_count=self.seed_count,
            arf_max_iter=self.arf_max_iter,
            forceatlas_max_iter=self.forceatlas_max_iter,
        )
        resolved: list[tuple[float, str, int, int, bool, dict[str, list[float]]]] = []
        for candidate in candidates:
            for orientation in range(self.orientation_count):
                angle = orientation * math.pi / self.orientation_count
                for reflected in (False, True):
                    raw = _rotate(candidate.positions, angle, reflected)
                    positions = _fit_isotropic(raw, internal, spec)
                    if not _symmetric_collision_relaxation(positions, internal, spec):
                        continue
                    score = _score(positions, internal, graph, spec)
                    resolved.append(
                        (
                            score,
                            candidate.engine,
                            candidate.seed,
                            orientation,
                            reflected,
                            positions,
                        )
                    )
        if not resolved:
            raise ValueError("graph layout could not resolve measured node collisions")
        resolved.sort(key=operator.itemgetter(slice(5)))
        positions: dict[str, list[float]] | None = None
        from .presentation import select_context_edges
        from .routing import route_graph_scene

        display_context_edges = select_context_edges(graph, max_per_source=2)

        for _, engine_name, _candidate_seed, _, _, internal_positions in resolved:
            routed_positions = {
                node_id: list(point) for node_id, point in internal_positions.items()
            }
            _add_terminals(routed_positions, internal, spec)
            provisional_positions = tuple(
                GraphPosition(
                    node.node_id,
                    routed_positions[node.node_id][0],
                    routed_positions[node.node_id][1],
                )
                for node in graph.nodes
            )
            provisional_scene = GraphScene(
                graph=graph,
                display_context_edges=display_context_edges,
                positions=provisional_positions,
                layout_spec=spec,
                engine=engine_name,
                seed=seed,
            )
            try:
                route_graph_scene(provisional_scene)
            except ValueError:
                continue
            positions = routed_positions
            break
        if positions is None:
            raise ValueError(
                "no isotropic graph-layout candidate supports collision-free routing"
            )
        return tuple(
            GraphPosition(
                node.node_id, positions[node.node_id][0], positions[node.node_id][1]
            )
            for node in graph.nodes
        )
