"""Truthful projection from playback placements to graph topology."""

from __future__ import annotations

from ..models import PlaybackPlan, PlaybackStep
from .model import END_NODE_ID, START_NODE_ID, ExplanationGraph, GraphEdge, GraphNode


def _added_bases(step: PlaybackStep) -> int:
    return sum(span.end - span.start for span in step.added_spans)


def _coordinate_overlap(source: PlaybackStep, target: PlaybackStep) -> int:
    if target.start < source.start:
        return 0
    overlap_start = target.start
    overlap_end = min(source.end, target.end)
    if overlap_end <= overlap_start:
        return 0
    source_offset = overlap_start - source.start
    length = overlap_end - overlap_start
    source_bases = source.placement_sequence[source_offset : source_offset + length]
    target_bases = target.placement_sequence[:length]
    return length if source_bases == target_bases else 0


def project_explanation_graph(plan: PlaybackPlan) -> ExplanationGraph:
    """Project all coordinate-compatible relations without display pruning."""
    placement_ids = tuple(step.placement_id for step in plan.steps)
    if len(placement_ids) != len(set(placement_ids)):
        raise ValueError("playback placement ids must be unique for graph projection")
    reserved = {START_NODE_ID, END_NODE_ID}.intersection(placement_ids)
    if reserved:
        raise ValueError(
            f"playback placement ids use reserved graph ids: {sorted(reserved)}"
        )

    nodes = [GraphNode(START_NODE_ID, None, terminal=True)]
    nodes.extend(
        GraphNode(step.placement_id, index, sequence=step.placement_sequence)
        for index, step in enumerate(plan.steps)
    )
    nodes.append(GraphNode(END_NODE_ID, None, terminal=True))

    traversal_edges: list[GraphEdge] = []
    previous_id = START_NODE_ID
    for index, step in enumerate(plan.steps):
        traversal_edges.append(
            GraphEdge(
                previous_id,
                step.placement_id,
                _added_bases(step),
                0,
                "realized_traversal",
                index,
            )
        )
        previous_id = step.placement_id
    traversal_edges.append(
        GraphEdge(
            previous_id, END_NODE_ID, None, 0, "realized_traversal", len(plan.steps)
        )
    )
    traversal_pairs = {(edge.source_id, edge.target_id) for edge in traversal_edges}

    context_edges: list[GraphEdge] = []
    for result in plan.constraint_results:
        if (
            result.upstream_placement_id not in placement_ids
            or result.downstream_placement_id not in placement_ids
        ):
            continue
        context_edges.append(
            GraphEdge(
                result.upstream_placement_id,
                result.downstream_placement_id,
                None,
                0,
                "declared_constraint",
            )
        )
    for source in plan.steps:
        for target in plan.steps:
            pair = (source.placement_id, target.placement_id)
            if source.placement_id == target.placement_id or pair in traversal_pairs:
                continue
            overlap = _coordinate_overlap(source, target)
            if overlap <= 0:
                continue
            context_edges.append(
                GraphEdge(
                    source.placement_id,
                    target.placement_id,
                    max(0, target.end - source.end),
                    overlap,
                    "realized_coordinate_overlap",
                )
            )
    context_edges.sort(
        key=lambda edge: (
            edge.source_id,
            0 if edge.relation_kind == "declared_constraint" else 1,
            -edge.overlap_bases,
            edge.added_bases or 0,
            edge.target_id,
        )
    )
    return ExplanationGraph(tuple(nodes), tuple(context_edges), tuple(traversal_edges))
