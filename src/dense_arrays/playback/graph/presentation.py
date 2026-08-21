"""Explicit display choices applied after semantic graph projection."""

from ..theme import step_color as resolve_step_color
from .model import ExplanationGraph, GraphEdge


def select_context_edges(
    graph: ExplanationGraph, *, max_per_source: int = 2
) -> tuple[GraphEdge, ...]:
    if max_per_source < 0:
        raise ValueError("max_per_source must be non-negative")
    grouped: dict[str, list[GraphEdge]] = {}
    for edge in graph.context_edges:
        grouped.setdefault(edge.source_id, []).append(edge)
    selected: list[GraphEdge] = []
    for source_id in sorted(grouped):
        candidates = sorted(
            grouped[source_id],
            key=lambda edge: (
                0 if edge.relation_kind == "declared_constraint" else 1,
                -edge.overlap_bases,
                edge.added_bases or 0,
                edge.target_id,
            ),
        )
        selected.extend(candidates[:max_per_source])
    return tuple(selected)


def step_color(step, index: int, profile: str = "categorical") -> str:
    return resolve_step_color(step, index, profile)
