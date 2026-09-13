"""Assemble routed context and traversal relations for a measured scene.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from functools import lru_cache

from .curves import (
    quadratic_arc_length,
    quadratic_arc_t,
    quadratic_point,
    quadratic_segment,
)
from .edge_routing import route_edge
from .labels import (
    EDGE_LABEL_FONT_SIZE_PT,
    LabelObstacles,
    place_edge_aware_label,
    place_leader_label,
)
from .model import GraphEdge, GraphRoutes, GraphScene, QuadraticCurve, RoutedEdge
from .obstacles import inside_viewport

__all__ = (
    "EDGE_LABEL_FONT_SIZE_PT",
    "edge_curve",
    "quadratic_arc_length",
    "quadratic_arc_t",
    "quadratic_point",
    "quadratic_segment",
    "route_graph_scene",
)


@lru_cache(maxsize=128)
def route_graph_scene(scene: GraphScene) -> GraphRoutes:
    """Route context and traversal edges and assign readable cost labels."""
    traversal_curves = []
    for edge in scene.graph.traversal_edges:
        curve = route_edge(scene, edge, allow_masked=True)
        traversal_curves.append((edge, curve))
    context = []
    for edge in scene.display_context_edges:
        try:
            curve = route_edge(
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
            label_position, label_bounds, label_anchor, label_font_size = _cost_label(
                scene, edge, curve, occupied_labels, other_curves
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
    """Route one edge around node obstacles without masking."""
    return route_edge(scene, edge)


def _cost_label(
    scene: GraphScene,
    edge: GraphEdge,
    curve: QuadraticCurve,
    occupied_labels: list[tuple[float, float, float, float]],
    other_curves: tuple[QuadraticCurve, ...],
) -> tuple[
    tuple[float, float], tuple[float, float, float, float], tuple[float, float], float
]:
    occupied = LabelObstacles(tuple(occupied_labels), other_curves)
    for placement in (place_edge_aware_label, place_leader_label):
        for font_size in (EDGE_LABEL_FONT_SIZE_PT, 11.0, 10.0, 9.0):
            try:
                position, bounds, anchor = placement(
                    scene, edge, curve, occupied, font_size
                )
            except ValueError:
                continue
            if placement is place_edge_aware_label and not inside_viewport(
                scene, bounds
            ):
                continue
            return position, bounds, anchor, font_size
    msg = (
        f"no collision-free cost label for edge {edge.source_id!r} "
        f"-> {edge.target_id!r}"
    )
    raise ValueError(msg)
