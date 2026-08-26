"""Matplotlib-faithful point-space node measurement."""

from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

from ..typography import PUBLICATION_NUCLEOTIDE_TYPOGRAPHY
from .model import (
    TERMINAL_DIAMETER,
    ExplanationGraph,
    GraphLayoutSpec,
    GraphViewport,
    NodeGeometry,
)

KMER_FONT_FAMILY = PUBLICATION_NUCLEOTIDE_TYPOGRAPHY.family
KMER_FONT_SIZE_PT = PUBLICATION_NUCLEOTIDE_TYPOGRAPHY.graph_font_size_pt
KMER_FONT_WEIGHT = PUBLICATION_NUCLEOTIDE_TYPOGRAPHY.weight
KMER_PADDING_X_PT = 3.8
KMER_PADDING_Y_PT = 2.8


def _text_path_extent(text: str) -> tuple[float, float]:
    prop = FontProperties(
        family=KMER_FONT_FAMILY, size=KMER_FONT_SIZE_PT, weight=KMER_FONT_WEIGHT
    )
    bounds = TextPath((0, 0), text or " ", prop=prop).get_extents()
    return float(bounds.width), float(bounds.height)


def node_box_width(sequence: str) -> float:
    width, _height = _text_path_extent(sequence)
    return width + (2.0 * KMER_PADDING_X_PT)


def matplotlib_layout_spec(
    axis,
    graph: ExplanationGraph,
    *,
    kmer_font_size_pt: float = KMER_FONT_SIZE_PT,
) -> GraphLayoutSpec:
    renderer = axis.figure.canvas.get_renderer()
    bbox = axis.get_window_extent(renderer=renderer)
    dpi = float(axis.figure.dpi)
    px_to_pt = 72.0 / dpi
    prop = FontProperties(
        family=KMER_FONT_FAMILY, size=kmer_font_size_pt, weight=KMER_FONT_WEIGHT
    )
    geometries: list[NodeGeometry] = []
    for node in graph.nodes:
        if node.terminal:
            geometries.append(
                NodeGeometry(node.node_id, TERMINAL_DIAMETER, TERMINAL_DIAMETER)
            )
            continue
        width_px, height_px, _descent = renderer.get_text_width_height_descent(
            node.sequence or " ", prop, ismath=False
        )
        geometries.append(
            NodeGeometry(
                node.node_id,
                (float(width_px) * px_to_pt) + (2.0 * KMER_PADDING_X_PT),
                (float(height_px) * px_to_pt) + (2.0 * KMER_PADDING_Y_PT),
            )
        )
    return GraphLayoutSpec(
        GraphViewport(float(bbox.width) * px_to_pt, float(bbox.height) * px_to_pt),
        tuple(geometries),
    )


def default_layout_spec(graph: ExplanationGraph) -> GraphLayoutSpec:
    geometries: list[NodeGeometry] = []
    for node in graph.nodes:
        if node.terminal:
            geometries.append(
                NodeGeometry(node.node_id, TERMINAL_DIAMETER, TERMINAL_DIAMETER)
            )
            continue
        width, height = _text_path_extent(node.sequence)
        geometries.append(
            NodeGeometry(
                node.node_id,
                width + (2.0 * KMER_PADDING_X_PT),
                height + (2.0 * KMER_PADDING_Y_PT),
            )
        )
    return GraphLayoutSpec(GraphViewport(388.0, 166.0), tuple(geometries))
