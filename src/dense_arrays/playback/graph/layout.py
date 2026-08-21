"""Deterministic point-space graph layout and terminal bracketing."""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from typing import Protocol

from ..models import PlaybackPlan
from .geometry import default_layout_spec
from .model import (
    END_NODE_ID,
    START_NODE_ID,
    ExplanationGraph,
    GraphLayoutSpec,
    GraphNode,
    GraphPosition,
    GraphScene,
)
from .presentation import select_context_edges
from .projection import project_explanation_graph


class GraphLayoutEngine(Protocol):
    name: str

    def layout(
        self, graph: ExplanationGraph, spec: GraphLayoutSpec, *, seed: int
    ) -> tuple[GraphPosition, ...]: ...


def _measured_graph(graph: ExplanationGraph, spec: GraphLayoutSpec) -> ExplanationGraph:
    return replace(
        graph,
        nodes=tuple(
            replace(
                node,
                width=spec.geometry(node.node_id).width_pt,
                height=spec.geometry(node.node_id).height_pt,
            )
            for node in graph.nodes
        ),
    )


def _overlap(
    left: GraphNode,
    left_pos: list[float],
    right: GraphNode,
    right_pos: list[float],
    clearance: float,
) -> tuple[float, float]:
    return (
        ((left.width + right.width) / 2.0)
        + clearance
        - abs(right_pos[0] - left_pos[0]),
        ((left.height + right.height) / 2.0)
        + clearance
        - abs(right_pos[1] - left_pos[1]),
    )


def _relax_internal_collisions(
    positions: dict[str, list[float]],
    nodes: tuple[GraphNode, ...],
    spec: GraphLayoutSpec,
) -> None:
    padding = spec.viewport.padding_pt
    min_x, max_x = padding, spec.viewport.width_pt - padding
    min_y, max_y = padding, spec.viewport.height_pt - padding
    ids = tuple(node.node_id for node in nodes)
    by_id = {node.node_id: node for node in nodes}
    for _ in range(240):
        moved = False
        for left_index, left_id in enumerate(ids):
            left = by_id[left_id]
            for right_id in ids[left_index + 1 :]:
                right = by_id[right_id]
                overlap_x, overlap_y = _overlap(
                    left,
                    positions[left_id],
                    right,
                    positions[right_id],
                    spec.node_clearance_pt,
                )
                if overlap_x <= 0.0 or overlap_y <= 0.0:
                    continue
                moved = True
                dx = positions[right_id][0] - positions[left_id][0]
                dy = positions[right_id][1] - positions[left_id][1]
                if overlap_x <= overlap_y:
                    direction = (
                        1.0 if dx > 0.0 or (dx == 0.0 and left_id < right_id) else -1.0
                    )
                    shift = overlap_x / 2.0 + 0.1
                    positions[left_id][0] -= direction * shift
                    positions[right_id][0] += direction * shift
                else:
                    direction = (
                        1.0 if dy > 0.0 or (dy == 0.0 and left_id < right_id) else -1.0
                    )
                    shift = overlap_y / 2.0 + 0.1
                    positions[left_id][1] -= direction * shift
                    positions[right_id][1] += direction * shift
                for node_id in (left_id, right_id):
                    node = by_id[node_id]
                    positions[node_id][0] = min(
                        max_x - node.width / 2.0,
                        max(min_x + node.width / 2.0, positions[node_id][0]),
                    )
                    positions[node_id][1] = min(
                        max_y - node.height / 2.0,
                        max(min_y + node.height / 2.0, positions[node_id][1]),
                    )
        if not moved:
            return
    for left_index, left_id in enumerate(ids):
        for right_id in ids[left_index + 1 :]:
            overlap_x, overlap_y = _overlap(
                by_id[left_id],
                positions[left_id],
                by_id[right_id],
                positions[right_id],
                spec.node_clearance_pt,
            )
            if overlap_x > 0.0 and overlap_y > 0.0:
                raise ValueError(
                    f"graph layout could not resolve node collision: {left_id!r}, {right_id!r}"
                )


def _raw_to_point_positions(
    raw: dict[str, tuple[float, float]],
    nodes: tuple[GraphNode, ...],
    spec: GraphLayoutSpec,
) -> dict[str, list[float]]:
    viewport = spec.viewport
    reserve_x = spec.terminal_gap_pt + spec.geometry(START_NODE_ID).width_pt
    available_w = viewport.width_pt - 2.0 * (viewport.padding_pt + reserve_x)
    available_h = viewport.height_pt - 2.0 * viewport.padding_pt
    if available_w <= 0.0 or available_h <= 0.0:
        raise ValueError("graph viewport is too small for terminal geometry")
    max_w = max(node.width for node in nodes)
    max_h = max(node.height for node in nodes)
    best = None
    for rotate in (False, True):
        candidate = {
            node_id: (-position[1], position[0]) if rotate else position
            for node_id, position in raw.items()
        }
        xs = [value[0] for value in candidate.values()]
        ys = [value[1] for value in candidate.values()]
        x_span = max(max(xs) - min(xs), 1e-9)
        y_span = max(max(ys) - min(ys), 1e-9)
        scale = min(
            max(1.0, available_w - max_w) / x_span,
            max(1.0, available_h - max_h) / y_span,
        )
        if best is None or scale > best[0]:
            best = (scale, candidate)
    if best is None:
        raise RuntimeError("graph layout did not produce a fitted candidate")
    scale, selected = best
    xs = [value[0] for value in selected.values()]
    ys = [value[1] for value in selected.values()]
    center_x = (min(xs) + max(xs)) / 2.0
    center_y = (min(ys) + max(ys)) / 2.0
    return {
        node_id: [
            viewport.width_pt / 2.0 + (value[0] - center_x) * scale,
            viewport.height_pt / 2.0 + (value[1] - center_y) * scale,
        ]
        for node_id, value in selected.items()
    }


@dataclass(frozen=True, slots=True)
class NetworkXForceAtlas2Layout:
    name: str = "networkx_forceatlas2"
    max_iter: int = 400
    scaling_ratio: float = 2.8
    gravity: float = 0.35

    def layout(
        self, graph: ExplanationGraph, spec: GraphLayoutSpec, *, seed: int
    ) -> tuple[GraphPosition, ...]:
        try:
            import networkx as nx
        except ImportError as exc:
            raise RuntimeError(
                "ForceAtlas2 playback requires the 'dense-arrays[playback]' extra"
            ) from exc
        internal = tuple(node for node in graph.nodes if not node.terminal)
        if not internal:
            raise ValueError("graph layout requires at least one non-terminal node")
        internal_ids = {node.node_id for node in internal}
        layout_graph = nx.Graph()
        layout_graph.add_nodes_from(internal_ids)
        for edge in graph.context_edges:
            if edge.source_id in internal_ids and edge.target_id in internal_ids:
                layout_graph.add_edge(edge.source_id, edge.target_id, weight=0.75)
        for edge in graph.traversal_edges:
            if edge.source_id in internal_ids and edge.target_id in internal_ids:
                existing = layout_graph.get_edge_data(
                    edge.source_id, edge.target_id, default={}
                )
                layout_graph.add_edge(
                    edge.source_id,
                    edge.target_id,
                    weight=max(1.15, float(existing.get("weight", 0.0))),
                )
        if len(internal) == 1:
            raw = {internal[0].node_id: (0.0, 0.0)}
        else:
            result = nx.forceatlas2_layout(
                layout_graph,
                max_iter=self.max_iter,
                scaling_ratio=self.scaling_ratio,
                gravity=self.gravity,
                node_size={
                    node.node_id: max(node.width, node.height) for node in internal
                },
                weight="weight",
                seed=seed,
            )
            raw = {
                node_id: (float(position[0]), float(position[1]))
                for node_id, position in result.items()
            }
        positions = _raw_to_point_positions(raw, internal, spec)
        _relax_internal_collisions(positions, internal, spec)
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
        available_left = spec.viewport.padding_pt
        available_right = spec.viewport.width_pt - spec.viewport.padding_pt
        if occupied_right - occupied_left > available_right - available_left + 1e-6:
            raise ValueError(
                "graph viewport cannot contain measured nodes and terminal gaps"
            )
        shift_x = (available_left + available_right) / 2.0 - (
            occupied_left + occupied_right
        ) / 2.0
        for position in positions.values():
            position[0] += shift_x
        return tuple(
            GraphPosition(
                node.node_id, positions[node.node_id][0], positions[node.node_id][1]
            )
            for node in graph.nodes
        )


from .isotropic_layout import NetworkXIsotropicLayout

_DEFAULT_LAYOUT_ENGINE = NetworkXIsotropicLayout()


@lru_cache(maxsize=64)
def _default_graph_scene(
    plan: PlaybackPlan,
    seed: int,
    spec: GraphLayoutSpec,
    max_context_edges_per_source: int,
) -> GraphScene:
    semantic_graph = project_explanation_graph(plan)
    measured_graph = _measured_graph(semantic_graph, spec)
    display_context = select_context_edges(
        measured_graph, max_per_source=max_context_edges_per_source
    )
    replace(measured_graph, context_edges=display_context)
    return GraphScene(
        measured_graph,
        display_context,
        _DEFAULT_LAYOUT_ENGINE.layout(measured_graph, spec, seed=seed),
        spec,
        _DEFAULT_LAYOUT_ENGINE.name,
        seed,
    )


def build_graph_scene(
    plan: PlaybackPlan,
    *,
    seed: int = 23,
    engine: GraphLayoutEngine | None = None,
    layout_spec: GraphLayoutSpec | None = None,
    max_context_edges_per_source: int = 2,
) -> GraphScene:
    semantic_graph = project_explanation_graph(plan)
    spec = layout_spec or default_layout_spec(semantic_graph)
    if engine is None:
        return _default_graph_scene(plan, seed, spec, max_context_edges_per_source)
    measured_graph = _measured_graph(semantic_graph, spec)
    display_context = select_context_edges(
        measured_graph, max_per_source=max_context_edges_per_source
    )
    display_graph = replace(measured_graph, context_edges=display_context)
    return GraphScene(
        measured_graph,
        display_context,
        engine.layout(display_graph, spec, seed=seed),
        spec,
        engine.name,
        seed,
    )
