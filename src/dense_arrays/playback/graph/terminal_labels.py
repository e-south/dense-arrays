"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/graph/terminal_labels.py

Place centered terminal captions clear of measured nodes and visible edges.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .model import END_NODE_ID, START_NODE_ID
from .obstacles import curve_intersects_bounds, intersects, node_bounds

if TYPE_CHECKING:
    from collections.abc import Mapping

    from .model import GraphRoutes, GraphScene, GraphViewport, RoutedEdge


@dataclass(frozen=True, slots=True)
class TerminalLabel:
    """One immutable caption center and its reserved point-space bounds."""

    node_id: str
    text: str
    x: float
    y: float
    bounds: tuple[float, float, float, float]


def _caption_fits(
    bounds: tuple[float, float, float, float],
    viewport: GraphViewport,
    occupied: list[tuple[float, float, float, float]],
    curves: tuple[RoutedEdge, ...],
) -> bool:
    if (
        min(bounds[0], bounds[1]) < 0
        or bounds[2] > viewport.width_pt
        or bounds[3] > viewport.height_pt
    ):
        return False
    if any(intersects(bounds, box) for box in occupied):
        return False
    expanded = (bounds[0] - 2.0, bounds[1] - 2.0, bounds[2] + 2.0, bounds[3] + 2.0)
    return not any(curve_intersects_bounds(edge.curve, expanded) for edge in curves)


def plan_terminal_labels(
    scene: GraphScene,
    routes: GraphRoutes,
    extents: Mapping[str, tuple[float, float]],
    *,
    include_context: bool = True,
    include_edge_labels: bool = True,
) -> tuple[TerminalLabel, ...]:
    """Prefer above, then below, without moving labels sideways or shrinking them."""
    if not scene.graph.traversal_edges:
        return ()
    curves = (*routes.traversal, *(routes.context if include_context else ()))
    occupied = [
        node_bounds(scene, node.node_id, inflate=1.5) for node in scene.graph.nodes
    ]
    if include_edge_labels:
        occupied.extend(
            edge.label_bounds for edge in curves if edge.label_bounds is not None
        )
    viewport = scene.layout_spec.viewport
    result = []
    for text, node_id in (("Start", START_NODE_ID), ("End", END_NODE_ID)):
        x, y = scene.position(node_id)
        radius = scene.geometry(node_id).height_pt / 2.0
        width, height = extents[text]
        selected = None
        for step in range(3):
            gap = 3.5 + step * (height + 3.5)
            for side in (1, -1):
                center_y = y + side * (radius + gap + height / 2.0)
                bounds = (
                    x - width / 2,
                    center_y - height / 2,
                    x + width / 2,
                    center_y + height / 2,
                )
                if not _caption_fits(bounds, viewport, occupied, curves):
                    continue
                selected = TerminalLabel(node_id, text, x, center_y, bounds)
                break
            if selected is not None:
                break
        if selected is None:
            msg = (
                f"no collision-free terminal label position for {text!r}; "
                "enlarge the graph viewport"
            )
            raise ValueError(msg)
        result.append(selected)
        occupied.append(selected.bounds)
    return tuple(result)
