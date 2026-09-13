"""Draw measured placement relations and nucleotide nodes.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .graph.geometry import (
    KMER_FONT_FAMILY,
    KMER_FONT_SIZE_PT,
    KMER_FONT_WEIGHT,
    matplotlib_layout_spec,
)
from .graph.layout import build_graph_scene
from .graph.model import END_NODE_ID, START_NODE_ID
from .graph.projection import project_explanation_graph
from .graph.routing import (
    EDGE_LABEL_FONT_SIZE_PT,
    quadratic_arc_t,
    quadratic_point,
    quadratic_segment,
    route_graph_scene,
)
from .theme import constraint_relation_color

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.path import Path

    from .graph.model import GraphRoutes, GraphScene, QuadraticCurve
    from .presentation import PlaybackDocument

_PAPER = "#ffffff"
_INK = "#4b5563"
_LINE = "#9cc9c1"
_TRAVERSED = "#50635f"
_ACTIVE = "#167a70"
_GRAPH_TEXT = "#1F2423"
_TERMINAL_FONT_SIZE_PT = 12.4
_SETTLED_PROGRESS = 0.999


def _matplotlib_path(curve: QuadraticCurve) -> Path:
    from matplotlib.path import Path

    return Path(
        (
            curve.visible_start,
            curve.visible_control or curve.control,
            curve.visible_end,
        ),
        (Path.MOVETO, Path.CURVE3, Path.CURVE3),
    )


def _matplotlib_segment_path(points: tuple[tuple[float, float], ...]) -> Path:
    from matplotlib.path import Path

    return Path(points, (Path.MOVETO, Path.CURVE3, Path.CURVE3))


def draw_graph(
    axis: Axes,
    document: PlaybackDocument,
    transition_index: int,
    progress: float,
    *,
    kmer_font_size_pt: float = KMER_FONT_SIZE_PT,
) -> None:
    """Draw context relations, permitted traversal, and sequence nodes."""
    semantic_graph = project_explanation_graph(document.plan)
    layout_spec = matplotlib_layout_spec(
        axis,
        semantic_graph,
        kmer_font_size_pt=kmer_font_size_pt,
    )
    scene = build_graph_scene(document.plan, layout_spec=layout_spec)
    routes = route_graph_scene(scene)
    axis.set_xlim(0, layout_spec.viewport.width_pt)
    axis.set_ylim(0, layout_spec.viewport.height_pt)
    axis.set_aspect("equal", adjustable="box")
    axis.axis("off")

    _draw_context_edges(axis, document, routes)
    active_geometry = _draw_traversal_edges(
        axis, document, routes, transition_index, progress
    )
    _draw_nodes(axis, document, scene, kmer_font_size_pt=kmer_font_size_pt)
    _draw_terminals(axis, scene)
    _draw_active_marker(axis, active_geometry)


def _draw_context_edges(
    axis: Axes, document: PlaybackDocument, routes: GraphRoutes
) -> None:
    from matplotlib.patches import FancyArrowPatch

    visible_context = tuple(
        routed
        for routed in routes.context
        if document.presentation.graph_detail == "full"
    )
    for routed in visible_context:
        curve = routed.curve
        declared_constraint = routed.edge.relation_kind == "declared_constraint"
        axis.add_patch(
            FancyArrowPatch(
                path=_matplotlib_path(curve),
                arrowstyle="-|>",
                mutation_scale=5.4,
                linewidth=1.7 if declared_constraint else 0.85,
                color=(
                    constraint_relation_color(document.presentation.color_profile)
                    if declared_constraint
                    else _LINE
                ),
                alpha=0.82 if declared_constraint else 0.34,
                zorder=0,
            )
        )


def _draw_traversal_edges(
    axis: Axes,
    document: PlaybackDocument,
    routes: GraphRoutes,
    transition_index: int,
    progress: float,
) -> tuple[QuadraticCurve, float] | None:
    from matplotlib.patches import FancyArrowPatch

    active_geometry = None
    for routed in routes.traversal:
        edge = routed.edge
        edge_index = edge.traversal_index
        if edge_index is None:
            msg = "traversal edge is missing its timeline index"
            raise ValueError(msg)
        completed = edge_index < transition_index or (
            edge_index == transition_index and progress >= _SETTLED_PROGRESS
        )
        active = edge_index == transition_index and progress < _SETTLED_PROGRESS
        color = _TRAVERSED if completed else _LINE
        width = 2.4 if completed else 1.45
        curve = routed.curve
        axis.add_patch(
            FancyArrowPatch(
                path=_matplotlib_path(curve),
                arrowstyle="-|>",
                mutation_scale=7.8,
                linewidth=width,
                color=color,
                alpha=1.0 if completed or active else 0.72,
                zorder=1,
            )
        )
        active_t = None
        if active:
            active_t = quadratic_arc_t(curve, progress)
            if active_t > curve.visible_t_start:
                partial_end = min(active_t, curve.visible_t_end)
                if partial_end > curve.visible_t_start:
                    axis.add_patch(
                        FancyArrowPatch(
                            path=_matplotlib_segment_path(
                                quadratic_segment(
                                    curve,
                                    curve.visible_t_start,
                                    partial_end,
                                )
                            ),
                            arrowstyle="-",
                            linewidth=2.4,
                            color=_ACTIVE,
                            alpha=1.0,
                            zorder=1.2,
                        )
                    )
        if (
            document.presentation.show_edge_costs
            and edge.added_bases is not None
            and routed.label_position is not None
        ):
            axis.text(
                routed.label_position[0],
                routed.label_position[1],
                str(edge.added_bases),
                ha="center",
                va="center",
                color=_GRAPH_TEXT,
                fontsize=routed.label_font_size or EDGE_LABEL_FONT_SIZE_PT,
                family=KMER_FONT_FAMILY,
                bbox={
                    "boxstyle": "round,pad=0.10",
                    "facecolor": _PAPER,
                    "edgecolor": _LINE,
                    "linewidth": 0.7,
                    "alpha": 0.97,
                },
                zorder=2.0,
            )
        if active:
            active_geometry = (curve, active_t)

    return active_geometry


def _draw_nodes(
    axis: Axes,
    document: PlaybackDocument,
    scene: GraphScene,
    *,
    kmer_font_size_pt: float,
) -> None:
    from matplotlib.patches import FancyBboxPatch

    steps = document.plan.steps
    for index, step in enumerate(steps):
        geometry = scene.geometry(step.placement_id)
        x, y = scene.position(step.placement_id)
        color = document.step_color(index)
        axis.add_patch(
            FancyBboxPatch(
                (x - geometry.width_pt / 2, y - geometry.height_pt / 2),
                geometry.width_pt,
                geometry.height_pt,
                boxstyle="round,pad=0.0,rounding_size=1.5",
                facecolor=color,
                edgecolor="none",
                linewidth=0.0,
                alpha=1.0,
                zorder=3,
            )
        )
        axis.text(
            x,
            y,
            step.placement_sequence,
            ha="center",
            va="center",
            color="#FFFFFF",
            fontsize=kmer_font_size_pt,
            family=KMER_FONT_FAMILY,
            fontweight=KMER_FONT_WEIGHT,
            zorder=4,
        )


def _draw_terminals(axis: Axes, scene: GraphScene) -> None:
    from matplotlib.patches import Circle

    terminals = (
        (("Start", START_NODE_ID), ("End", END_NODE_ID))
        if scene.graph.traversal_edges
        else ()
    )
    for label, node_id in terminals:
        x, y = scene.position(node_id)
        geometry = scene.geometry(node_id)
        radius = geometry.width_pt / 2.0
        axis.add_patch(
            Circle(
                (x, y),
                radius=radius,
                facecolor=_PAPER,
                edgecolor=_INK,
                linewidth=1.25,
                zorder=4,
            )
        )
        axis.text(
            x,
            y + radius + 3.5,
            label,
            ha="center",
            va="bottom",
            color=_GRAPH_TEXT,
            fontsize=_TERMINAL_FONT_SIZE_PT,
            family=KMER_FONT_FAMILY,
            zorder=4,
        )


def _draw_active_marker(
    axis: Axes, active_geometry: tuple[QuadraticCurve, float] | None
) -> None:
    if active_geometry is not None:
        curve, active_t = active_geometry
        point = quadratic_point(
            curve.motion_start,
            curve.control,
            curve.motion_end,
            active_t,
        )
        axis.scatter(
            (point[0],),
            (point[1],),
            s=26,
            facecolor=_ACTIVE,
            edgecolor="none",
            linewidth=0.0,
            zorder=1.5,
        )
