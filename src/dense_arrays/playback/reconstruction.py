"""Compile persisted placements into a truthful semantic playback plan.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import hashlib
import json

from dense_arrays._record_validation import records
from dense_arrays.realized import RealizedArray

from .models import (
    ConstraintResult,
    CoordinateSpan,
    NoticeLevel,
    OrderingStatus,
    PlaybackAuthority,
    PlaybackNotice,
    PlaybackPlan,
    PlaybackStep,
)
from .validation import coordinate_order, inferred_ordering, reveal_intervals


def _realization_digest(realized: RealizedArray) -> str:
    payload = {
        "coordinate_space": realized.coordinate_space,
        "constraints": [
            {
                "constraint_id": item.constraint_id,
                "downstream_placement_id": item.downstream_placement_id,
                "max_distance_bp": item.max_distance_bp,
                "min_distance_bp": item.min_distance_bp,
                "upstream_placement_id": item.upstream_placement_id,
            }
            for item in realized.constraints
        ],
        "placements": [
            {
                "feature_id": item.feature_id,
                "kind": item.kind.value,
                "orientation": item.orientation.value,
                "placement_id": item.placement_id,
                "sequence": item.sequence,
                "start": item.start,
            }
            for item in realized.placements
        ],
        "sequence": realized.sequence,
        "source_id": realized.source_id,
    }
    encoded = json.dumps(
        payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def reconstruct_playback(
    realized: RealizedArray, *, notices: tuple[PlaybackNotice, ...] = ()
) -> PlaybackPlan:
    """Build deterministic placement playback without claiming solver path authority.

    Parameters
    ----------
    realized
        Validated persisted sequence, placements, and declared constraints.
    notices
        Explicit caller-owned evidence qualifications to preserve in the plan.
        Metadata names do not imply recovery methods or other process claims.

    Returns
    -------
    PlaybackPlan
        A validated placement-reconstructed plan with evaluated constraints.

    Raises
    ------
    TypeError
        If the input or notice records have invalid types.
    ValueError
        If caller notices contradict reserved v1 evidence codes.
    """
    if not isinstance(realized, RealizedArray):
        msg = "realized must be a RealizedArray record"
        raise TypeError(msg)
    caller_notices = records(notices, PlaybackNotice, field_name="notices")
    ordered = coordinate_order(realized.placements)
    ordering_status = OrderingStatus(inferred_ordering(realized.placements))
    revealed = [False] * len(realized.sequence)
    steps: list[PlaybackStep] = []
    predecessor_id: str | None = None
    for index, placement in enumerate(ordered):
        steps.append(
            PlaybackStep(
                index=index,
                placement_id=placement.placement_id,
                feature_id=placement.feature_id,
                start=placement.start,
                end=placement.end,
                placement_kind=placement.kind.value,
                orientation=placement.orientation.value,
                placement_sequence=placement.sequence,
                added_spans=tuple(
                    CoordinateSpan(start=start, end=end)
                    for start, end in reveal_intervals(
                        start=placement.start,
                        end=placement.end,
                        revealed=revealed,
                    )
                ),
                predecessor_placement_id=predecessor_id,
                label=placement.label,
            )
        )
        predecessor_id = placement.placement_id

    placement_by_id = {
        placement.placement_id: placement for placement in realized.placements
    }
    constraint_results: list[ConstraintResult] = []
    for constraint in realized.constraints:
        upstream = placement_by_id[constraint.upstream_placement_id]
        downstream = placement_by_id[constraint.downstream_placement_id]
        actual_distance = downstream.start - upstream.end
        constraint_results.append(
            ConstraintResult(
                constraint_id=constraint.constraint_id,
                upstream_placement_id=constraint.upstream_placement_id,
                downstream_placement_id=constraint.downstream_placement_id,
                actual_distance_bp=actual_distance,
                min_distance_bp=constraint.min_distance_bp,
                max_distance_bp=constraint.max_distance_bp,
                passed=constraint.min_distance_bp
                <= actual_distance
                <= constraint.max_distance_bp,
                label=constraint.label,
            )
        )

    plan_notices = [
        PlaybackNotice(
            code="placement_reconstructed",
            message=(
                "Playback order is reconstructed from persisted placements; "
                "relations are not recorded solver-selected edges."
            ),
        )
    ]
    if ordering_status is OrderingStatus.AMBIGUOUS:
        plan_notices.append(
            PlaybackNotice(
                code="ambiguous_order",
                message=(
                    "Multiple placements share a start or one contains another; "
                    "the display uses the documented deterministic tie-break."
                ),
                level=NoticeLevel.WARNING,
            )
        )
    elif ordering_status is OrderingStatus.LAYOUT_ONLY:
        plan_notices.append(
            PlaybackNotice(
                code="layout_only",
                message=(
                    "The realized placements contain an internal uncovered span; "
                    "render the layout without an active-edge claim."
                ),
                level=NoticeLevel.WARNING,
            )
        )
    unrevealed_count = revealed.count(False)
    if unrevealed_count:
        plan_notices.append(
            PlaybackNotice(
                code="unannotated_sequence",
                message=(
                    f"{unrevealed_count} realized bases are not covered "
                    "by persisted placements."
                ),
            )
        )

    return PlaybackPlan(
        source_id=realized.source_id,
        source_digest=realized.source_digest,
        realization_digest=_realization_digest(realized),
        realized_sequence=realized.sequence,
        authority=PlaybackAuthority.PLACEMENT_RECONSTRUCTED,
        ordering_status=ordering_status,
        steps=tuple(steps),
        constraint_results=tuple(constraint_results),
        notices=(*plan_notices, *caller_notices),
    )
