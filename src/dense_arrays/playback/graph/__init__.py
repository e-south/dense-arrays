"""Measured explanation-graph geometry for dense-array playback."""

from ..positions import journey_path_positions, radial_path_positions
from .geometry import (
    KMER_FONT_FAMILY,
    KMER_FONT_SIZE_PT,
    KMER_FONT_WEIGHT,
    default_layout_spec,
    matplotlib_layout_spec,
    node_box_width,
)
from .layout import (
    GraphLayoutEngine,
    NetworkXForceAtlas2Layout,
    build_graph_scene,
)
from .model import (
    DEFAULT_STEP_COLORS,
    DOWNSTREAM_FIXED_COLOR,
    END_NODE_ID,
    NODE_HEIGHT,
    START_NODE_ID,
    TERMINAL_DIAMETER,
    UPSTREAM_FIXED_COLOR,
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
from .presentation import select_context_edges, step_color
from .projection import project_explanation_graph
from .routing import edge_curve, quadratic_point, route_graph_scene

__all__ = (
    "DEFAULT_STEP_COLORS",
    "DOWNSTREAM_FIXED_COLOR",
    "END_NODE_ID",
    "KMER_FONT_FAMILY",
    "KMER_FONT_SIZE_PT",
    "KMER_FONT_WEIGHT",
    "NODE_HEIGHT",
    "START_NODE_ID",
    "TERMINAL_DIAMETER",
    "UPSTREAM_FIXED_COLOR",
    "ExplanationGraph",
    "GraphEdge",
    "GraphLayoutEngine",
    "GraphLayoutSpec",
    "GraphNode",
    "GraphPosition",
    "GraphRoutes",
    "GraphScene",
    "GraphViewport",
    "NetworkXForceAtlas2Layout",
    "NodeGeometry",
    "QuadraticCurve",
    "RoutedEdge",
    "build_graph_scene",
    "default_layout_spec",
    "edge_curve",
    "journey_path_positions",
    "matplotlib_layout_spec",
    "node_box_width",
    "project_explanation_graph",
    "quadratic_point",
    "radial_path_positions",
    "route_graph_scene",
    "select_context_edges",
    "step_color",
)
