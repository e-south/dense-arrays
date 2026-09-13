"""Validate geometry and supported evidence semantics for playback v1.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from dense_arrays._record_validation import validate_placement_sequence

if TYPE_CHECKING:
    from collections.abc import Sequence

    from .models import PlaybackPlan


class _Interval(Protocol):
    placement_id: str
    start: int

    @property
    def end(self) -> int: ...


def coordinate_order[T: _Interval](items: Sequence[T]) -> list[T]:
    """Return the v1 deterministic coordinate order."""
    return sorted(
        items, key=lambda item: (item.start, item.end - item.start, item.placement_id)
    )


def inferred_ordering(items: Sequence[_Interval]) -> str:
    """Infer ordering strength from actual placement intervals."""
    ordered = coordinate_order(items)
    cursor = ordered[0].end
    ambiguous = False
    layout_gap = False
    previous_start = ordered[0].start
    for item in ordered[1:]:
        layout_gap |= item.start > cursor
        ambiguous |= item.start == previous_start or item.end <= cursor
        cursor = max(cursor, item.end)
        previous_start = item.start
    if layout_gap:
        return "layout_only"
    if ambiguous:
        return "ambiguous"
    return "unique"


def reveal_intervals(
    *, start: int, end: int, revealed: list[bool]
) -> tuple[tuple[int, int], ...]:
    """Reveal an interval and return the maximal runs of newly covered bases."""
    spans: list[tuple[int, int]] = []
    run_start: int | None = None
    for coordinate in range(start, end):
        if not revealed[coordinate] and run_start is None:
            run_start = coordinate
        if revealed[coordinate] and run_start is not None:
            spans.append((run_start, coordinate))
            run_start = None
    if run_start is not None:
        spans.append((run_start, end))
    for coordinate in range(start, end):
        revealed[coordinate] = True
    return tuple(spans)


def _validate_step_layout(plan: PlaybackPlan) -> None:
    if list(plan.steps) != coordinate_order(plan.steps):
        msg = "playback steps must follow deterministic v1 coordinate order"
        raise ValueError(msg)
    revealed = [False] * len(plan.realized_sequence)
    predecessor_id = None
    for step in plan.steps:
        validate_placement_sequence(
            placement_id=step.placement_id,
            start=step.start,
            end=step.end,
            sequence=step.placement_sequence,
            realized_sequence=plan.realized_sequence,
        )
        if step.predecessor_placement_id != predecessor_id:
            msg = (
                f"step {step.index} predecessor must be the preceding "
                "coordinate-ordered placement"
            )
            raise ValueError(msg)
        expected = reveal_intervals(start=step.start, end=step.end, revealed=revealed)
        observed = tuple((span.start, span.end) for span in step.added_spans)
        if observed != expected:
            msg = (
                f"step {step.index} added_spans must exactly match "
                "newly covered placement bases"
            )
            raise ValueError(msg)
        predecessor_id = step.placement_id


def _validate_constraint_layout(plan: PlaybackPlan) -> None:
    by_id = {step.placement_id: step for step in plan.steps}
    result_ids = [item.constraint_id for item in plan.constraint_results]
    if len(result_ids) != len(set(result_ids)):
        msg = "constraint_id values must be unique within a playback plan"
        raise ValueError(msg)
    for result in plan.constraint_results:
        missing = {
            result.upstream_placement_id,
            result.downstream_placement_id,
        } - by_id.keys()
        if missing:
            msg = (
                f"constraint {result.constraint_id!r} references "
                f"unknown placements: {sorted(missing)}"
            )
            raise ValueError(msg)
        upstream = by_id[result.upstream_placement_id]
        downstream = by_id[result.downstream_placement_id]
        if result.actual_distance_bp != downstream.start - upstream.end:
            msg = (
                f"constraint {result.constraint_id!r} actual_distance_bp "
                "must match placement coordinates"
            )
            raise ValueError(msg)


def validate_plan(plan: PlaybackPlan) -> None:
    """Check all cross-record claims after local record validation."""
    if plan.authority != "placement_reconstructed":
        msg = (
            "playback v1 supports only placement_reconstructed authority; "
            "solver_selected requires a future trace contract"
        )
        raise ValueError(msg)
    placement_ids = [step.placement_id for step in plan.steps]
    if len(placement_ids) != len(set(placement_ids)):
        msg = "placement_id values must be unique within a playback plan"
        raise ValueError(msg)
    _validate_step_layout(plan)
    if plan.ordering_status != inferred_ordering(plan.steps):
        msg = "ordering_status must match the realized placement intervals"
        raise ValueError(msg)
    _validate_constraint_layout(plan)
    for notice in plan.notices:
        if notice.code == "solver_selected" or (
            notice.code in {"ambiguous_order", "layout_only"}
            and notice.code
            != {"ambiguous": "ambiguous_order", "layout_only": "layout_only"}.get(
                plan.ordering_status
            )
        ):
            msg = f"notice code {notice.code!r} contradicts playback v1 evidence"
            raise ValueError(msg)
