"""Terminal captions remain centered and clear of visible graph edges.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from dense_arrays.playback.graph.model import (
    END_NODE_ID,
    START_NODE_ID,
    ExplanationGraph,
    GraphEdge,
    GraphLayoutSpec,
    GraphNode,
    GraphPosition,
    GraphRoutes,
    GraphScene,
    GraphViewport,
    NodeGeometry,
    QuadraticCurve,
    RoutedEdge,
)
from dense_arrays.playback.graph.obstacles import intersects, node_bounds
from dense_arrays.playback.graph.terminal_labels import plan_terminal_labels


def _scene() -> GraphScene:
    nodes = (
        GraphNode(START_NODE_ID, None, terminal=True),
        GraphNode("part", 0),
        GraphNode(END_NODE_ID, None, terminal=True),
    )
    edges = (
        GraphEdge(START_NODE_ID, "part", None, 0, "traversal", 0),
        GraphEdge("part", END_NODE_ID, None, 0, "traversal", 1),
    )
    return GraphScene(
        ExplanationGraph(nodes, (), edges),
        (),
        (
            GraphPosition(START_NODE_ID, 30, 50),
            GraphPosition("part", 100, 50),
            GraphPosition(END_NODE_ID, 170, 50),
        ),
        GraphLayoutSpec(
            GraphViewport(200, 100),
            tuple(
                NodeGeometry(
                    n.node_id, 11 if n.terminal else 30, 11 if n.terminal else 14
                )
                for n in nodes
            ),
        ),
        "test",
        0,
    )


def _route(
    edge: GraphEdge,
    start: tuple[float, float],
    control: tuple[float, float],
    end: tuple[float, float],
) -> RoutedEdge:
    return RoutedEdge(edge, QuadraticCurve(start, control, end, start, end))


def test_terminal_labels_choose_above_or_below_without_horizontal_drift() -> None:
    scene = _scene()
    routes = GraphRoutes(
        (),
        (
            _route(scene.graph.traversal_edges[0], (30, 55.5), (30, 90), (85, 50)),
            _route(scene.graph.traversal_edges[1], (115, 50), (170, 10), (170, 44.5)),
        ),
    )
    labels = plan_terminal_labels(scene, routes, {"Start": (24, 12), "End": (20, 12)})
    assert labels[0].x == 30
    assert labels[0].y < 50
    assert labels[1].x == 170
    assert labels[1].y > 50


def test_terminal_choice_is_stable_and_hidden_context_does_not_block_it() -> None:
    scene = _scene()
    routes = GraphRoutes(
        (_route(scene.graph.traversal_edges[0], (30, 55.5), (30, 90), (85, 50)),), ()
    )
    extents = {"Start": (24, 12), "End": (20, 12)}
    assert plan_terminal_labels(scene, routes, extents, include_context=False)[0].y > 50
    assert plan_terminal_labels(scene, routes, extents)[0].y < 50
    assert plan_terminal_labels(scene, routes, extents) == plan_terminal_labels(
        scene, routes, extents
    )


def test_terminal_labels_fail_if_both_vertical_sides_are_obstructed() -> None:
    scene = _scene()
    edge = scene.graph.traversal_edges[0]
    routes = GraphRoutes(
        (),
        (
            _route(edge, (30, 56), (30, 75), (30, 100)),
            _route(edge, (30, 44), (30, 25), (30, 0)),
        ),
    )
    with pytest.raises(ValueError, match=r"terminal label.*Start"):
        plan_terminal_labels(scene, routes, {"Start": (24, 12), "End": (20, 12)})


def test_terminal_label_bounds_respect_other_nodes_and_viewport() -> None:

    scene = _scene()
    scene = replace(
        scene,
        positions=(
            scene.positions[0],
            GraphPosition("part", 30, 65),
            scene.positions[2],
        ),
    )
    labels = plan_terminal_labels(
        scene, GraphRoutes((), ()), {"Start": (24, 12), "End": (20, 12)}
    )
    assert labels[0].y < 50
    assert not intersects(labels[0].bounds, node_bounds(scene, "part"))
    assert all(0 <= label.bounds[1] < label.bounds[3] <= 100 for label in labels)


def test_hidden_edge_cost_labels_do_not_displace_terminal_captions() -> None:
    scene = _scene()
    routed = _route(scene.graph.traversal_edges[0], (35.5, 50), (60, 50), (85, 50))
    routed = replace(routed, label_bounds=(15, 58, 45, 75))
    routes = GraphRoutes((), (routed,))
    extents = {"Start": (24, 12), "End": (20, 12)}
    assert plan_terminal_labels(scene, routes, extents)[0].y < 50
    assert (
        plan_terminal_labels(scene, routes, extents, include_edge_labels=False)[0].y
        > 50
    )
