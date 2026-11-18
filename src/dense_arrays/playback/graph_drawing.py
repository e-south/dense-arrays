"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/graph_drawing.py

Draw measured placement relations and nucleotide nodes.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from .frame_schedule import PlaybackFrame
from .graph.arrowheads import draw_arrow
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
from .graph.terminal_labels import TerminalLabel, plan_terminal_labels
from .theme import (
    RESTING_COLOR,
    RESTING_TEXT_COLOR,
    blend_color,
    constraint_relation_color,
)
from .timeline import placement_progress

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.path import Path

    from .graph.model import GraphRoutes, GraphScene, QuadraticCurve
    from .graph_style import GraphArrowDimensions
    from .presentation import PlaybackDocument

_PAPER = "#ffffff"
_INK = "#4b5563"
_LINE = RESTING_COLOR
_ACTIVE = "#167a70"
_GRAPH_TEXT = "#1F2423"
_TERMINAL_FONT_SIZE_PT = 12.4
_SETTLED_PROGRESS = 0.999


def _transition_complete(index: int, frame: PlaybackFrame) -> bool:
    return index < frame.transition_index or (
        index == frame.transition_index and frame.progress >= _SETTLED_PROGRESS
    )


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
    terminal_labels = plan_terminal_labels(
        scene,
        routes,
        _terminal_label_extents(axis),
        include_context=document.presentation.graph_detail == "full",
        include_edge_labels=document.presentation.show_edge_costs,
    )
    axis.set_xlim(0, layout_spec.viewport.width_pt)
    axis.set_ylim(0, layout_spec.viewport.height_pt)
    axis.set_aspect("equal", adjustable="box")
    axis.axis("off")

    frame = PlaybackFrame(transition_index, progress)
    arrow_dimensions = document.presentation.graph_arrows.resolve(kmer_font_size_pt)
    _draw_context_edges(axis, document, routes, frame, arrow_dimensions)
    active_geometry = _draw_traversal_edges(
        axis, document, routes, frame, arrow_dimensions
    )
    _draw_nodes(
        axis,
        document,
        scene,
        PlaybackFrame(transition_index, progress),
        kmer_font_size_pt=kmer_font_size_pt,
    )
    _draw_terminals(axis, scene, terminal_labels, transition_index, progress)
    _draw_active_marker(axis, active_geometry)


def _draw_context_edges(
    axis: Axes,
    document: PlaybackDocument,
    routes: GraphRoutes,
    frame: PlaybackFrame,
    arrow_dimensions: GraphArrowDimensions,
) -> None:

    transition_index, progress = frame.transition_index, frame.progress
    visible_context = tuple(
        routed
        for routed in routes.context
        if document.presentation.graph_detail == "full"
    )
    for routed in visible_context:
        curve = routed.curve
        declared_constraint = routed.edge.relation_kind == "declared_constraint"
        emphasis = min(
            placement_progress(index, transition_index, progress)
            for index, step in enumerate(document.plan.steps)
            if step.placement_id in (routed.edge.source_id, routed.edge.target_id)
        )
        color = (
            blend_color(
                document.presentation.graph_arrows.resting_color,
                constraint_relation_color(document.presentation.color_profile)
                if document.presentation.color_profile == "constraints"
                else document.presentation.graph_arrows.completed_color,
                emphasis,
            )
            if declared_constraint
            else document.presentation.graph_arrows.resting_color
        )
        draw_arrow(
            axis,
            _matplotlib_path(curve),
            arrow_dimensions,
            color=color,
            zorder=0,
        )


def _draw_traversal_edges(
    axis: Axes,
    document: PlaybackDocument,
    routes: GraphRoutes,
    frame: PlaybackFrame,
    arrow_dimensions: GraphArrowDimensions,
) -> tuple[QuadraticCurve, float] | None:
    from matplotlib.patches import FancyArrowPatch

    transition_index, progress = frame.transition_index, frame.progress
    active_geometry = None
    for routed in routes.traversal:
        edge = routed.edge
        edge_index = edge.traversal_index
        if edge_index is None:
            msg = "traversal edge is missing its timeline index"
            raise ValueError(msg)
        terminal = edge.target_id == END_NODE_ID
        timeline_index = edge_index - 1 if terminal else edge_index
        completed = _transition_complete(timeline_index, frame)
        active = (
            not terminal
            and edge_index == transition_index
            and 0 < progress < _SETTLED_PROGRESS
        )
        color = (
            document.presentation.graph_arrows.completed_color
            if completed
            else document.presentation.graph_arrows.resting_color
        )
        curve = routed.curve
        draw_arrow(
            axis,
            _matplotlib_path(curve),
            arrow_dimensions,
            color=color,
            zorder=1,
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
                            linewidth=arrow_dimensions.shaft_width_pt,
                            capstyle="butt",
                            gid="graph-arrow-shaft",
                            color=document.presentation.graph_arrows.active_color,
                            alpha=document.presentation.graph_arrows.opacity,
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
                color=blend_color(
                    RESTING_TEXT_COLOR,
                    _GRAPH_TEXT,
                    placement_progress(edge_index, transition_index, progress),
                ),
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
    frame: PlaybackFrame,
    *,
    kmer_font_size_pt: float,
) -> None:
    from matplotlib.patches import FancyBboxPatch

    steps = document.plan.steps
    for index, step in enumerate(steps):
        geometry = scene.geometry(step.placement_id)
        x, y = scene.position(step.placement_id)
        emphasis = placement_progress(index, frame.transition_index, frame.progress)
        color = blend_color(RESTING_COLOR, document.step_color(index), emphasis)
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
            color=blend_color(RESTING_TEXT_COLOR, "#FFFFFF", emphasis),
            fontsize=kmer_font_size_pt,
            family=KMER_FONT_FAMILY,
            fontweight=KMER_FONT_WEIGHT,
            zorder=4,
        )


def _terminal_label_extents(axis: Axes) -> dict[str, tuple[float, float]]:
    from matplotlib.font_manager import FontProperties

    prop = FontProperties(family=KMER_FONT_FAMILY, size=_TERMINAL_FONT_SIZE_PT)
    renderer = axis.figure.canvas.get_renderer()
    scale = 72.0 / axis.figure.dpi
    return {
        text: (
            renderer.get_text_width_height_descent(text, prop, ismath=False)[0] * scale,
            _TERMINAL_FONT_SIZE_PT * 1.25,
        )
        for text in ("Start", "End")
    }


def _draw_terminals(
    axis: Axes,
    scene: GraphScene,
    labels: tuple[TerminalLabel, ...],
    transition_index: int,
    progress: float,
) -> None:
    from matplotlib.patches import Circle

    for label in labels:
        node_id = label.node_id
        index = 0 if node_id == START_NODE_ID else len(scene.graph.traversal_edges) - 2
        emphasis = (
            placement_progress(index, transition_index, progress)
            if node_id == START_NODE_ID
            else float(
                _transition_complete(index, PlaybackFrame(transition_index, progress))
            )
        )
        x, y = scene.position(node_id)
        geometry = scene.geometry(node_id)
        radius = geometry.width_pt / 2.0
        axis.add_patch(
            Circle(
                (x, y),
                radius=radius,
                facecolor=_PAPER,
                edgecolor=blend_color(RESTING_COLOR, _INK, emphasis),
                linewidth=1.25,
                zorder=4,
            )
        )
        axis.text(
            label.x,
            label.y,
            label.text,
            ha="center",
            va="center",
            color=blend_color(RESTING_TEXT_COLOR, _GRAPH_TEXT, emphasis),
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
