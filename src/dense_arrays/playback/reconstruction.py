"""Compile persisted placements into a truthful semantic playback plan."""

from __future__ import annotations

import hashlib
import json

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


def _validate_realization(realized: RealizedArray) -> None:
    placement_ids = {placement.placement_id for placement in realized.placements}
    for placement in realized.placements:
        if placement.end > len(realized.sequence):
            msg = (
                f"placement {placement.placement_id!r} ends at {placement.end}, "
                f"beyond sequence length {len(realized.sequence)}"
            )
            raise ValueError(msg)
        observed = realized.sequence[placement.start : placement.end]
        if observed != placement.sequence:
            msg = (
                f"placement {placement.placement_id!r} is sequence-inconsistent: "
                f"expected {placement.sequence!r}, observed {observed!r}"
            )
            raise ValueError(msg)
    for constraint in realized.constraints:
        missing = {
            constraint.upstream_placement_id,
            constraint.downstream_placement_id,
        } - placement_ids
        if missing:
            msg = f"constraint {constraint.constraint_id!r} references unknown placements: {sorted(missing)}"
            raise ValueError(msg)


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


def _ordering_status(realized: RealizedArray) -> OrderingStatus:
    ordered = sorted(
        realized.placements,
        key=lambda item: (item.start, len(item.sequence), item.placement_id),
    )
    cursor = ordered[0].end
    ambiguous = False
    layout_gap = False
    previous_start = ordered[0].start
    for placement in ordered[1:]:
        if placement.start > cursor:
            layout_gap = True
        if placement.start == previous_start or placement.end <= cursor:
            ambiguous = True
        cursor = max(cursor, placement.end)
        previous_start = placement.start
    if layout_gap:
        return OrderingStatus.LAYOUT_ONLY
    if ambiguous:
        return OrderingStatus.AMBIGUOUS
    return OrderingStatus.UNIQUE


def _added_spans(
    *, start: int, end: int, revealed: list[bool]
) -> tuple[CoordinateSpan, ...]:
    spans: list[CoordinateSpan] = []
    run_start: int | None = None
    for coordinate in range(start, end):
        if not revealed[coordinate] and run_start is None:
            run_start = coordinate
        if revealed[coordinate] and run_start is not None:
            spans.append(CoordinateSpan(start=run_start, end=coordinate))
            run_start = None
    if run_start is not None:
        spans.append(CoordinateSpan(start=run_start, end=end))
    for coordinate in range(start, end):
        revealed[coordinate] = True
    return tuple(spans)


def reconstruct_playback(realized: RealizedArray) -> PlaybackPlan:
    """Build deterministic placement playback without claiming solver path authority."""
    _validate_realization(realized)
    ordered = sorted(
        realized.placements,
        key=lambda item: (item.start, len(item.sequence), item.placement_id),
    )
    ordering_status = _ordering_status(realized)
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
                added_spans=_added_spans(
                    start=placement.start,
                    end=placement.end,
                    revealed=revealed,
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

    notices = [
        PlaybackNotice(
            code="placement_reconstructed",
            message=(
                "Playback order is reconstructed from persisted placements; "
                "relations are not recorded solver-selected edges."
            ),
        )
    ]
    recovered_coordinates = [
        placement
        for placement in realized.placements
        if placement.metadata.get("coordinate_source") == "offset_raw"
    ]
    if recovered_coordinates:
        notices.append(
            PlaybackNotice(
                code="coordinate_recovered",
                message=(
                    f"{len(recovered_coordinates)} placement coordinate(s) were recovered "
                    "from offset_raw by exact realized-sequence agreement."
                ),
                level=NoticeLevel.WARNING,
            )
        )
    if ordering_status is OrderingStatus.AMBIGUOUS:
        notices.append(
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
        notices.append(
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
        notices.append(
            PlaybackNotice(
                code="unannotated_sequence",
                message=f"{unrevealed_count} realized bases are not covered by persisted placements.",
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
        notices=tuple(notices),
    )
