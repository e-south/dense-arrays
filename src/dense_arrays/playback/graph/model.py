"""Immutable topology and measured scene contracts for playback graphs.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import dataclass

from .. import theme

DEFAULT_STEP_COLORS = theme.DEFAULT_STEP_COLORS
UPSTREAM_FIXED_COLOR = "#7D86D1"
DOWNSTREAM_FIXED_COLOR = "#C886D1"
START_NODE_ID = "__dense_arrays_start__"
END_NODE_ID = "__dense_arrays_end__"
NODE_HEIGHT = 13.0
TERMINAL_DIAMETER = 11.0


@dataclass(frozen=True, slots=True)
class GraphNode:
    """One semantic node; measured scenes populate compatibility extents."""

    node_id: str
    step_index: int | None
    width: float = 0.0
    height: float = 0.0
    terminal: bool = False
    sequence: str = ""


@dataclass(frozen=True, slots=True)
class GraphEdge:
    """One truthful directed relation in the explanation graph."""

    source_id: str
    target_id: str
    added_bases: int | None
    overlap_bases: int
    relation_kind: str
    traversal_index: int | None = None


@dataclass(frozen=True, slots=True)
class ExplanationGraph:
    """Store semantic nodes, context relations, and permitted traversal edges."""

    nodes: tuple[GraphNode, ...]
    context_edges: tuple[GraphEdge, ...]
    traversal_edges: tuple[GraphEdge, ...]


@dataclass(frozen=True, slots=True)
class NodeGeometry:
    """Store measured node dimensions in points."""

    node_id: str
    width_pt: float
    height_pt: float


@dataclass(frozen=True, slots=True)
class GraphViewport:
    """Store the drawing extent and its padding in points."""

    width_pt: float
    height_pt: float
    padding_pt: float = 7.0


@dataclass(frozen=True, slots=True)
class GraphLayoutSpec:
    """All visual inputs that can change a frozen layout."""

    viewport: GraphViewport
    node_geometries: tuple[NodeGeometry, ...]
    terminal_gap_pt: float = 10.0
    node_clearance_pt: float = 12.0
    route_clearance_pt: float = 3.0

    def geometry(self, node_id: str) -> NodeGeometry:
        """Return measured geometry for a known node identity."""
        for geometry in self.node_geometries:
            if geometry.node_id == node_id:
                return geometry
        msg = f"missing node geometry: {node_id!r}"
        raise KeyError(msg)


@dataclass(frozen=True, slots=True)
class GraphPosition:
    """Bind a node identity to its point-space coordinates."""

    node_id: str
    x: float
    y: float


@dataclass(frozen=True, slots=True)
class GraphScene:
    """Combine graph semantics, displayed context, and measured layout."""

    graph: ExplanationGraph
    display_context_edges: tuple[GraphEdge, ...]
    positions: tuple[GraphPosition, ...]
    layout_spec: GraphLayoutSpec
    engine: str
    seed: int

    def node(self, node_id: str) -> GraphNode:
        """Return a graph node by its stable identity."""
        for node in self.graph.nodes:
            if node.node_id == node_id:
                return node
        msg = f"unknown graph node: {node_id!r}"
        raise KeyError(msg)

    def geometry(self, node_id: str) -> NodeGeometry:
        """Return measured geometry for a known node identity."""
        return self.layout_spec.geometry(node_id)

    def position(self, node_id: str) -> tuple[float, float]:
        """Return point-space coordinates for a known node identity."""
        for position in self.positions:
            if position.node_id == node_id:
                return (position.x, position.y)
        msg = f"missing graph position: {node_id!r}"
        raise KeyError(msg)


@dataclass(frozen=True, slots=True)
class QuadraticCurve:
    """Keep canonical and endpoint-clipped geometry for one quadratic edge."""

    visible_start: tuple[float, float]
    control: tuple[float, float]
    visible_end: tuple[float, float]
    motion_start: tuple[float, float]
    motion_end: tuple[float, float]
    visible_control: tuple[float, float] | None = None
    visible_t_start: float = 0.0
    visible_t_end: float = 1.0


@dataclass(frozen=True, slots=True)
class RoutedEdge:
    """Pair one semantic edge with its curve and optional cost label."""

    edge: GraphEdge
    curve: QuadraticCurve
    label_position: tuple[float, float] | None = None
    label_bounds: tuple[float, float, float, float] | None = None
    label_font_size: float | None = None
    label_anchor: tuple[float, float] | None = None


@dataclass(frozen=True, slots=True)
class GraphRoutes:
    """Collect routed context and traversal edges for drawing."""

    context: tuple[RoutedEdge, ...]
    traversal: tuple[RoutedEdge, ...]
