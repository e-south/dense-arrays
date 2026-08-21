"""Strict JSON serialization for realized arrays and playback plans."""

from __future__ import annotations

import json
from collections.abc import Mapping

from dense_arrays.realized import (
    REALIZED_ARRAY_SCHEMA_VERSION,
    DeclaredConstraint,
    Orientation,
    Placement,
    PlacementKind,
    RealizedArray,
)

from .models import (
    PLAYBACK_PLAN_SCHEMA_VERSION,
    ConstraintResult,
    CoordinateSpan,
    NoticeLevel,
    OrderingStatus,
    PlaybackAuthority,
    PlaybackNotice,
    PlaybackPlan,
    PlaybackStep,
)


def _object(value: object, *, context: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        msg = f"{context} must be a JSON object"
        raise TypeError(msg)
    return value


def _exact_keys(
    value: Mapping[str, object], *, expected: set[str], context: str
) -> None:
    missing = sorted(expected - set(value))
    unknown = sorted(set(value) - expected)
    if missing or unknown:
        msg = f"{context} has missing keys {missing} and unknown keys {unknown}"
        raise ValueError(msg)


def _list(value: object, *, context: str) -> list[object]:
    if not isinstance(value, list):
        msg = f"{context} must be a JSON array"
        raise TypeError(msg)
    return value


def realized_array_to_dict(realized: RealizedArray) -> dict[str, object]:
    """Return the canonical JSON-compatible realized-array representation."""
    return {
        "schema_version": realized.schema_version,
        "source_id": realized.source_id,
        "source_digest": realized.source_digest,
        "coordinate_space": realized.coordinate_space,
        "sequence": realized.sequence,
        "placements": [
            {
                "placement_id": item.placement_id,
                "feature_id": item.feature_id,
                "kind": item.kind.value,
                "sequence": item.sequence,
                "start": item.start,
                "end": item.end,
                "orientation": item.orientation.value,
                "label": item.label,
                "metadata": dict(item.metadata),
            }
            for item in realized.placements
        ],
        "constraints": [
            {
                "constraint_id": item.constraint_id,
                "upstream_placement_id": item.upstream_placement_id,
                "downstream_placement_id": item.downstream_placement_id,
                "min_distance_bp": item.min_distance_bp,
                "max_distance_bp": item.max_distance_bp,
                "label": item.label,
                "metadata": dict(item.metadata),
            }
            for item in realized.constraints
        ],
        "provenance": dict(realized.provenance),
    }


def realized_array_from_dict(value: Mapping[str, object]) -> RealizedArray:
    """Parse a realized array and reject missing or unknown fields."""
    expected = {
        "schema_version",
        "source_id",
        "source_digest",
        "coordinate_space",
        "sequence",
        "placements",
        "constraints",
        "provenance",
    }
    _exact_keys(value, expected=expected, context="realized_array")
    if value["schema_version"] != REALIZED_ARRAY_SCHEMA_VERSION:
        msg = f"unsupported realized-array schema: {value['schema_version']!r}"
        raise ValueError(msg)
    placements: list[Placement] = []
    placement_keys = {
        "placement_id",
        "feature_id",
        "kind",
        "sequence",
        "start",
        "end",
        "orientation",
        "label",
        "metadata",
    }
    for index, raw in enumerate(_list(value["placements"], context="placements")):
        item = _object(raw, context=f"placements[{index}]")
        _exact_keys(item, expected=placement_keys, context=f"placements[{index}]")
        placement = Placement(
            placement_id=str(item["placement_id"]),
            feature_id=str(item["feature_id"]),
            kind=PlacementKind(str(item["kind"])),
            sequence=str(item["sequence"]),
            start=int(item["start"]),
            orientation=Orientation(str(item["orientation"])),
            label=None if item["label"] is None else str(item["label"]),
            metadata=_object(item["metadata"], context=f"placements[{index}].metadata"),
        )
        if item["end"] != placement.end:
            msg = f"placements[{index}].end does not match start + sequence length"
            raise ValueError(msg)
        placements.append(placement)
    constraints: list[DeclaredConstraint] = []
    constraint_keys = {
        "constraint_id",
        "upstream_placement_id",
        "downstream_placement_id",
        "min_distance_bp",
        "max_distance_bp",
        "label",
        "metadata",
    }
    for index, raw in enumerate(_list(value["constraints"], context="constraints")):
        item = _object(raw, context=f"constraints[{index}]")
        _exact_keys(item, expected=constraint_keys, context=f"constraints[{index}]")
        constraints.append(
            DeclaredConstraint(
                constraint_id=str(item["constraint_id"]),
                upstream_placement_id=str(item["upstream_placement_id"]),
                downstream_placement_id=str(item["downstream_placement_id"]),
                min_distance_bp=int(item["min_distance_bp"]),
                max_distance_bp=int(item["max_distance_bp"]),
                label=None if item["label"] is None else str(item["label"]),
                metadata=_object(
                    item["metadata"], context=f"constraints[{index}].metadata"
                ),
            )
        )
    source_digest = value["source_digest"]
    return RealizedArray(
        source_id=str(value["source_id"]),
        source_digest=None if source_digest is None else str(source_digest),
        coordinate_space=str(value["coordinate_space"]),
        sequence=str(value["sequence"]),
        placements=tuple(placements),
        constraints=tuple(constraints),
        provenance=_object(value["provenance"], context="provenance"),
    )


def playback_plan_to_dict(plan: PlaybackPlan) -> dict[str, object]:
    """Return the canonical JSON-compatible playback-plan representation."""
    return {
        "schema_version": plan.schema_version,
        "source_id": plan.source_id,
        "source_digest": plan.source_digest,
        "realization_digest": plan.realization_digest,
        "realized_sequence": plan.realized_sequence,
        "authority": plan.authority.value,
        "ordering_status": plan.ordering_status.value,
        "steps": [
            {
                "index": item.index,
                "placement_id": item.placement_id,
                "feature_id": item.feature_id,
                "start": item.start,
                "end": item.end,
                "placement_kind": item.placement_kind,
                "orientation": item.orientation,
                "placement_sequence": item.placement_sequence,
                "added_spans": [
                    {"start": span.start, "end": span.end} for span in item.added_spans
                ],
                "predecessor_placement_id": item.predecessor_placement_id,
                "relation_kind": item.relation_kind,
                "label": item.label,
            }
            for item in plan.steps
        ],
        "constraint_results": [
            {
                "constraint_id": item.constraint_id,
                "upstream_placement_id": item.upstream_placement_id,
                "downstream_placement_id": item.downstream_placement_id,
                "actual_distance_bp": item.actual_distance_bp,
                "min_distance_bp": item.min_distance_bp,
                "max_distance_bp": item.max_distance_bp,
                "passed": item.passed,
                "label": item.label,
            }
            for item in plan.constraint_results
        ],
        "notices": [
            {"code": item.code, "message": item.message, "level": item.level.value}
            for item in plan.notices
        ],
    }


def playback_plan_from_dict(value: Mapping[str, object]) -> PlaybackPlan:
    """Parse a playback plan and reject missing or unknown fields."""
    expected = {
        "schema_version",
        "source_id",
        "source_digest",
        "realization_digest",
        "realized_sequence",
        "authority",
        "ordering_status",
        "steps",
        "constraint_results",
        "notices",
    }
    _exact_keys(value, expected=expected, context="playback_plan")
    if value["schema_version"] != PLAYBACK_PLAN_SCHEMA_VERSION:
        msg = f"unsupported playback-plan schema: {value['schema_version']!r}"
        raise ValueError(msg)
    steps: list[PlaybackStep] = []
    step_keys = {
        "index",
        "placement_id",
        "feature_id",
        "start",
        "end",
        "placement_kind",
        "orientation",
        "placement_sequence",
        "added_spans",
        "predecessor_placement_id",
        "relation_kind",
        "label",
    }
    for index, raw in enumerate(_list(value["steps"], context="steps")):
        item = _object(raw, context=f"steps[{index}]")
        _exact_keys(item, expected=step_keys, context=f"steps[{index}]")
        spans = tuple(
            CoordinateSpan(
                start=int(_object(span, context="added_span")["start"]),
                end=int(_object(span, context="added_span")["end"]),
            )
            for span in _list(
                item["added_spans"], context=f"steps[{index}].added_spans"
            )
        )
        predecessor = item["predecessor_placement_id"]
        steps.append(
            PlaybackStep(
                index=int(item["index"]),
                placement_id=str(item["placement_id"]),
                feature_id=str(item["feature_id"]),
                start=int(item["start"]),
                end=int(item["end"]),
                placement_kind=str(item["placement_kind"]),
                orientation=str(item["orientation"]),
                placement_sequence=str(item["placement_sequence"]),
                added_spans=spans,
                predecessor_placement_id=None
                if predecessor is None
                else str(predecessor),
                relation_kind=str(item["relation_kind"]),
                label=None if item["label"] is None else str(item["label"]),
            )
        )
    results: list[ConstraintResult] = []
    result_keys = {
        "constraint_id",
        "upstream_placement_id",
        "downstream_placement_id",
        "actual_distance_bp",
        "min_distance_bp",
        "max_distance_bp",
        "passed",
        "label",
    }
    for index, raw in enumerate(
        _list(value["constraint_results"], context="constraint_results")
    ):
        item = _object(raw, context=f"constraint_results[{index}]")
        _exact_keys(item, expected=result_keys, context=f"constraint_results[{index}]")
        results.append(
            ConstraintResult(
                constraint_id=str(item["constraint_id"]),
                upstream_placement_id=str(item["upstream_placement_id"]),
                downstream_placement_id=str(item["downstream_placement_id"]),
                actual_distance_bp=int(item["actual_distance_bp"]),
                min_distance_bp=int(item["min_distance_bp"]),
                max_distance_bp=int(item["max_distance_bp"]),
                passed=bool(item["passed"]),
                label=None if item["label"] is None else str(item["label"]),
            )
        )
    notices: list[PlaybackNotice] = []
    notice_keys = {"code", "message", "level"}
    for index, raw in enumerate(_list(value["notices"], context="notices")):
        item = _object(raw, context=f"notices[{index}]")
        _exact_keys(item, expected=notice_keys, context=f"notices[{index}]")
        notices.append(
            PlaybackNotice(
                code=str(item["code"]),
                message=str(item["message"]),
                level=NoticeLevel(str(item["level"])),
            )
        )
    source_digest = value["source_digest"]
    return PlaybackPlan(
        source_id=str(value["source_id"]),
        source_digest=None if source_digest is None else str(source_digest),
        realization_digest=str(value["realization_digest"]),
        realized_sequence=str(value["realized_sequence"]),
        authority=PlaybackAuthority(str(value["authority"])),
        ordering_status=OrderingStatus(str(value["ordering_status"])),
        steps=tuple(steps),
        constraint_results=tuple(results),
        notices=tuple(notices),
    )


def dumps_realized_array(realized: RealizedArray) -> str:
    """Serialize a realized array deterministically."""
    return json.dumps(
        realized_array_to_dict(realized),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def loads_realized_array(payload: str) -> RealizedArray:
    """Deserialize a strict realized-array JSON document."""
    return realized_array_from_dict(
        _object(json.loads(payload), context="realized_array")
    )


def dumps_playback_plan(plan: PlaybackPlan) -> str:
    """Serialize a playback plan deterministically."""
    return json.dumps(
        playback_plan_to_dict(plan),
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def loads_playback_plan(payload: str) -> PlaybackPlan:
    """Deserialize a strict playback-plan JSON document."""
    return playback_plan_from_dict(
        _object(json.loads(payload), context="playback_plan")
    )
