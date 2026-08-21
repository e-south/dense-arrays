"""Topology-derived NetworkX layout candidates for compact playback graphs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .model import END_NODE_ID, START_NODE_ID, ExplanationGraph, GraphNode


@dataclass(frozen=True, slots=True)
class RawLayoutCandidate:
    engine: str
    seed: int
    positions: dict[str, tuple[float, float]]


def _layout_graph(nx: Any, graph: ExplanationGraph, internal_ids: set[str]) -> Any:
    layout_graph = nx.Graph()
    layout_graph.add_nodes_from(sorted(internal_ids))
    for edge in graph.context_edges:
        if edge.source_id in internal_ids and edge.target_id in internal_ids:
            layout_graph.add_edge(edge.source_id, edge.target_id, weight=1.0)
    for edge in graph.traversal_edges:
        if edge.source_id not in internal_ids or edge.target_id not in internal_ids:
            continue
        existing = layout_graph.get_edge_data(
            edge.source_id, edge.target_id, default={}
        )
        layout_graph.add_edge(
            edge.source_id,
            edge.target_id,
            weight=max(0.9, float(existing.get("weight", 0.0))),
        )
    return layout_graph


def _initial_positions(nx: Any, node_ids: tuple[str, ...], seed: int) -> dict[str, Any]:
    import math

    import numpy as np

    base = nx.circular_layout(node_ids, scale=1.0)
    phase = (seed % 360) * math.pi / 180.0
    cosine, sine = math.cos(phase), math.sin(phase)
    rng = np.random.default_rng(seed)
    positions: dict[str, Any] = {}
    for node_id in node_ids:
        x, y = base[node_id]
        positions[node_id] = np.asarray(
            (
                cosine * x - sine * y + rng.normal(0.0, 0.015),
                sine * x + cosine * y + rng.normal(0.0, 0.015),
            ),
            dtype=float,
        )
    return positions


def generate_layout_candidates(
    graph: ExplanationGraph,
    internal: tuple[GraphNode, ...],
    *,
    seed: int,
    seed_count: int,
    arf_max_iter: int,
    forceatlas_max_iter: int,
) -> tuple[RawLayoutCandidate, ...]:
    """Generate deterministic ARF and ForceAtlas2 equilibria."""
    try:
        import networkx as nx
    except ImportError as exc:
        raise RuntimeError(
            "isotropic playback requires the 'dense-arrays[playback]' extra"
        ) from exc

    internal_ids = {node.node_id for node in internal}
    node_ids = tuple(sorted(internal_ids))
    layout_graph = _layout_graph(nx, graph, internal_ids)
    if len(node_ids) == 1:
        return (RawLayoutCandidate("singleton", seed, {node_ids[0]: (0.0, 0.0)}),)

    median_size = sorted(max(node.width, node.height) for node in internal)[
        len(internal) // 2
    ]
    node_size = {
        node.node_id: max(node.width, node.height) / max(median_size, 1e-9)
        for node in internal
    }
    node_mass = {
        node_id: 1.0 + 0.12 * float(layout_graph.degree[node_id])
        for node_id in node_ids
    }
    candidates: list[RawLayoutCandidate] = []
    for index in range(seed_count):
        candidate_seed = seed + index * 37
        initial = _initial_positions(nx, node_ids, candidate_seed)
        try:
            result = nx.arf_layout(
                layout_graph,
                pos=initial,
                scaling=1.0,
                a=1.35,
                etol=1e-7,
                dt=0.0025,
                max_iter=arf_max_iter,
                seed=candidate_seed,
            )
            candidates.append(
                RawLayoutCandidate(
                    "networkx_arf",
                    candidate_seed,
                    {
                        node_id: (float(position[0]), float(position[1]))
                        for node_id, position in result.items()
                    },
                )
            )
        except (ArithmeticError, ValueError):
            pass

        result = nx.forceatlas2_layout(
            layout_graph,
            pos=initial,
            max_iter=forceatlas_max_iter,
            jitter_tolerance=0.2,
            scaling_ratio=2.2,
            gravity=0.8,
            distributed_action=True,
            strong_gravity=False,
            node_mass=node_mass,
            node_size=node_size,
            weight="weight",
            linlog=False,
            seed=candidate_seed,
        )
        candidates.append(
            RawLayoutCandidate(
                "networkx_forceatlas2",
                candidate_seed,
                {
                    node_id: (float(position[0]), float(position[1]))
                    for node_id, position in result.items()
                },
            )
        )

    if not candidates:
        raise ValueError("no graph-layout candidate converged")
    return tuple(candidates)


def internal_edge_pairs(graph: ExplanationGraph) -> tuple[tuple[str, str], ...]:
    pairs: set[tuple[str, str]] = set()
    for edge in (*graph.context_edges, *graph.traversal_edges):
        if edge.source_id in {START_NODE_ID, END_NODE_ID}:
            continue
        if edge.target_id in {START_NODE_ID, END_NODE_ID}:
            continue
        pair = tuple(sorted((edge.source_id, edge.target_id)))
        if pair[0] != pair[1]:
            pairs.add(pair)
    return tuple(sorted(pairs))
