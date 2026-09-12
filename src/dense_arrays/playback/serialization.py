"""Strict JSON serialization for realized arrays and playback plans.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import json
from collections.abc import Mapping

from dense_arrays._record_validation import integer, mutable_json
from dense_arrays.realized import (
    REALIZED_ARRAY_SCHEMA_VERSION,
    DeclaredConstraint,
    Placement,
    RealizedArray,
)

from .models import (
    PLAYBACK_PLAN_SCHEMA_VERSION,
    ConstraintResult,
    CoordinateSpan,
    PlaybackNotice,
    PlaybackPlan,
    PlaybackStep,
)


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result = {}
    for key, value in pairs:
        if key in result:
            msg = f"duplicate JSON object key: {key!r}"
            raise ValueError(msg)
        result[key] = value
    return result


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
                "metadata": mutable_json(item.metadata),
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
                "metadata": mutable_json(item.metadata),
            }
            for item in realized.constraints
        ],
        "provenance": mutable_json(realized.provenance),
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
    value = _object(value, context="realized_array")
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
            placement_id=item["placement_id"],
            feature_id=item["feature_id"],
            kind=item["kind"],
            sequence=item["sequence"],
            start=item["start"],
            orientation=item["orientation"],
            label=None if item["label"] is None else item["label"],
            metadata=_object(item["metadata"], context=f"placements[{index}].metadata"),
        )
        integer(item["end"], field_name=f"placements[{index}].end", minimum=0)
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
                constraint_id=item["constraint_id"],
                upstream_placement_id=item["upstream_placement_id"],
                downstream_placement_id=item["downstream_placement_id"],
                min_distance_bp=item["min_distance_bp"],
                max_distance_bp=item["max_distance_bp"],
                label=None if item["label"] is None else item["label"],
                metadata=_object(
                    item["metadata"], context=f"constraints[{index}].metadata"
                ),
            )
        )
    source_digest = value["source_digest"]
    return RealizedArray(
        source_id=value["source_id"],
        source_digest=None if source_digest is None else source_digest,
        coordinate_space=value["coordinate_space"],
        sequence=value["sequence"],
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
    value = _object(value, context="playback_plan")
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
        spans = []
        for span_index, raw_span in enumerate(
            _list(item["added_spans"], context=f"steps[{index}].added_spans")
        ):
            context = f"steps[{index}].added_spans[{span_index}]"
            span = _object(raw_span, context=context)
            _exact_keys(span, expected={"start", "end"}, context=context)
            spans.append(CoordinateSpan(start=span["start"], end=span["end"]))
        predecessor = item["predecessor_placement_id"]
        steps.append(
            PlaybackStep(
                index=item["index"],
                placement_id=item["placement_id"],
                feature_id=item["feature_id"],
                start=item["start"],
                end=item["end"],
                placement_kind=item["placement_kind"],
                orientation=item["orientation"],
                placement_sequence=item["placement_sequence"],
                added_spans=spans,
                predecessor_placement_id=None if predecessor is None else predecessor,
                relation_kind=item["relation_kind"],
                label=None if item["label"] is None else item["label"],
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
        passed = item["passed"]
        if not isinstance(passed, bool):
            msg = f"constraint_results[{index}].passed must be a JSON boolean"
            raise TypeError(msg)
        results.append(
            ConstraintResult(
                constraint_id=item["constraint_id"],
                upstream_placement_id=item["upstream_placement_id"],
                downstream_placement_id=item["downstream_placement_id"],
                actual_distance_bp=item["actual_distance_bp"],
                min_distance_bp=item["min_distance_bp"],
                max_distance_bp=item["max_distance_bp"],
                passed=passed,
                label=None if item["label"] is None else item["label"],
            )
        )
    notices: list[PlaybackNotice] = []
    notice_keys = {"code", "message", "level"}
    for index, raw in enumerate(_list(value["notices"], context="notices")):
        item = _object(raw, context=f"notices[{index}]")
        _exact_keys(item, expected=notice_keys, context=f"notices[{index}]")
        notices.append(
            PlaybackNotice(
                code=item["code"],
                message=item["message"],
                level=item["level"],
            )
        )
    source_digest = value["source_digest"]
    return PlaybackPlan(
        source_id=value["source_id"],
        source_digest=None if source_digest is None else source_digest,
        realization_digest=value["realization_digest"],
        realized_sequence=value["realized_sequence"],
        authority=value["authority"],
        ordering_status=value["ordering_status"],
        steps=tuple(steps),
        constraint_results=tuple(results),
        notices=tuple(notices),
    )


def dumps_realized_array(realized: RealizedArray) -> str:
    """Serialize a realized array deterministically."""
    return json.dumps(
        realized_array_to_dict(realized),
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def loads_realized_array(payload: str) -> RealizedArray:
    """Deserialize a strict realized-array JSON document."""
    return realized_array_from_dict(
        _object(
            json.loads(payload, object_pairs_hook=_unique_object),
            context="realized_array",
        )
    )


def dumps_playback_plan(plan: PlaybackPlan) -> str:
    """Serialize a playback plan deterministically."""
    return json.dumps(
        playback_plan_to_dict(plan),
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def loads_playback_plan(payload: str) -> PlaybackPlan:
    """Deserialize a strict playback-plan JSON document."""
    return playback_plan_from_dict(
        _object(
            json.loads(payload, object_pairs_hook=_unique_object),
            context="playback_plan",
        )
    )
