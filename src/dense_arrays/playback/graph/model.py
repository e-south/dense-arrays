"""Immutable topology and measured scene contracts for playback graphs."""

from __future__ import annotations

from dataclasses import dataclass

DEFAULT_STEP_COLORS = (
    "#67BFA5",
    "#D883A4",
    "#7BA4D9",
    "#C08A56",
    "#5DA79F",
    "#D1B06C",
    "#74C0CB",
    "#86A5D8",
)
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
    nodes: tuple[GraphNode, ...]
    context_edges: tuple[GraphEdge, ...]
    traversal_edges: tuple[GraphEdge, ...]


@dataclass(frozen=True, slots=True)
class NodeGeometry:
    node_id: str
    width_pt: float
    height_pt: float


@dataclass(frozen=True, slots=True)
class GraphViewport:
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
        for geometry in self.node_geometries:
            if geometry.node_id == node_id:
                return geometry
        raise KeyError(f"missing node geometry: {node_id!r}")


@dataclass(frozen=True, slots=True)
class GraphPosition:
    node_id: str
    x: float
    y: float


@dataclass(frozen=True, slots=True)
class GraphScene:
    graph: ExplanationGraph
    display_context_edges: tuple[GraphEdge, ...]
    positions: tuple[GraphPosition, ...]
    layout_spec: GraphLayoutSpec
    engine: str
    seed: int

    def node(self, node_id: str) -> GraphNode:
        for node in self.graph.nodes:
            if node.node_id == node_id:
                return node
        raise KeyError(f"unknown graph node: {node_id!r}")

    def geometry(self, node_id: str) -> NodeGeometry:
        return self.layout_spec.geometry(node_id)

    def position(self, node_id: str) -> tuple[float, float]:
        for position in self.positions:
            if position.node_id == node_id:
                return (position.x, position.y)
        raise KeyError(f"missing graph position: {node_id!r}")


@dataclass(frozen=True, slots=True)
class QuadraticCurve:
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
    edge: GraphEdge
    curve: QuadraticCurve
    label_position: tuple[float, float] | None = None
    label_bounds: tuple[float, float, float, float] | None = None
    label_font_size: float | None = None
    label_anchor: tuple[float, float] | None = None


@dataclass(frozen=True, slots=True)
class GraphRoutes:
    context: tuple[RoutedEdge, ...]
    traversal: tuple[RoutedEdge, ...]
