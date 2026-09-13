"""Exercise persisted playback and realized-layout trust boundaries.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import json
from dataclasses import replace
from typing import TYPE_CHECKING

import pytest

from dense_arrays.playback import (
    CoordinateSpan,
    PlaybackNotice,
    dumps_playback_plan,
    dumps_realized_array,
    loads_playback_plan,
    loads_realized_array,
    playback_plan_from_dict,
    playback_plan_to_dict,
    realized_array_from_dict,
    realized_array_to_dict,
    reconstruct_playback,
)
from dense_arrays.realized import (
    DeclaredConstraint,
    Placement,
    PlacementKind,
    RealizedArray,
)

if TYPE_CHECKING:
    from collections.abc import Callable


def _realized() -> RealizedArray:
    return RealizedArray(
        source_id="record:1",
        sequence="AAACCC",
        placements=(
            Placement("left", "feature:left", PlacementKind.TFBS, "AAA", 0),
            Placement("right", "feature:right", PlacementKind.OTHER, "CCC", 3),
        ),
        constraints=(DeclaredConstraint("distance", "left", "right", 0, 0),),
    )


@pytest.mark.parametrize(
    "field,value",
    [
        ("placement_id", 12),
        ("feature_id", None),
        ("kind", "unknown"),
        ("kind", 1),
        ("orientation", "unknown"),
        ("sequence", 123),
        ("label", 7),
        ("metadata", {1: "invalid-key"}),
    ],
)
def test_realized_fields_have_the_same_python_and_json_contract(
    field: str, value: object
):
    realized = _realized()
    with pytest.raises((TypeError, ValueError)):
        replace(realized.placements[0], **{field: value})
    payload = realized_array_to_dict(realized)
    payload["placements"][0][field] = value
    with pytest.raises((TypeError, ValueError)):
        realized_array_from_dict(payload)


@pytest.mark.parametrize("value", [True, 0.5, "0", None, float("inf")])
def test_json_coordinates_do_not_coerce_values(value: object):
    payload = realized_array_to_dict(_realized())
    payload["placements"][0]["start"] = value
    with pytest.raises((TypeError, ValueError)):
        realized_array_from_dict(payload)


def test_realized_rejects_sequence_disagreement_at_construction():
    with pytest.raises(ValueError, match="sequence-inconsistent"):
        replace(_realized(), sequence="TTTCCC")


def test_realized_rejects_unknown_constraint_references_at_construction():
    realized = _realized()
    constraint = replace(realized.constraints[0], downstream_placement_id="absent")
    with pytest.raises(ValueError, match="unknown placements"):
        replace(realized, constraints=(constraint,))


def test_metadata_is_an_immutable_json_snapshot():
    source = {"evidence": {"positions": [0, 1], "quality": 0.5}}
    placement = replace(_realized().placements[0], metadata=source)
    realized = replace(
        _realized(),
        provenance=source,
        placements=(placement, _realized().placements[1]),
    )
    before = dumps_realized_array(realized)
    source["evidence"]["positions"].append(9)
    source["evidence"]["quality"] = 1.0
    assert dumps_realized_array(realized) == before
    with pytest.raises(TypeError):
        realized.provenance["evidence"]["quality"] = 4
    with pytest.raises(TypeError):
        placement.metadata["evidence"]["positions"][0] = 9
    assert dumps_realized_array(loads_realized_array(before)) == before


@pytest.mark.parametrize("value", [float("nan"), float("inf"), object(), {1: "bad"}])
def test_non_json_provenance_is_rejected(value: object):
    with pytest.raises((TypeError, ValueError)):
        replace(_realized(), provenance={"nested": [value]})


def test_metadata_cycles_are_rejected():
    circular = {}
    circular["self"] = circular
    with pytest.raises(ValueError, match="cyclic"):
        replace(_realized(), provenance=circular)


@pytest.mark.parametrize(
    "path,value",
    [
        (("realized_sequence",), "TTTCCC"),
        (("steps", 1, "placement_id"), "left"),
        (("steps", 1, "predecessor_placement_id"), "absent"),
        (("steps", 0, "added_spans", 0, "end"), 100),
        (("steps", 0, "added_spans", 0, "extra"), "unexpected"),
        (("steps", 0, "start"), False),
        (("steps", 0, "start"), 0.5),
        (("constraint_results", 0, "passed"), False),
        (("constraint_results", 0, "downstream_placement_id"), "absent"),
    ],
)
def test_all_nine_audit_mutations_are_rejected(
    path: tuple[str | int, ...], value: object
):
    payload = playback_plan_to_dict(reconstruct_playback(_realized()))
    target = payload
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises((TypeError, ValueError)):
        loads_playback_plan(json.dumps(payload))


@pytest.mark.parametrize(
    "field,value",
    [
        ("index", True),
        ("start", True),
        ("end", 3.0),
        ("placement_id", None),
        ("feature_id", 8),
        ("placement_kind", "unknown"),
        ("orientation", "unknown"),
        ("placement_sequence", "ZZZ"),
        ("label", False),
        ("relation_kind", "solver_selected"),
        ("predecessor_placement_id", 1),
        ("added_spans", ["invalid"]),
    ],
)
def test_python_step_fields_reject_malformed_values(field: str, value: object):
    step = reconstruct_playback(_realized()).steps[0]
    with pytest.raises((TypeError, ValueError)):
        replace(step, **{field: value})


@pytest.mark.parametrize(
    "field,value",
    [
        ("source_id", 1),
        ("source_digest", 2),
        ("realization_digest", "invalid"),
        ("authority", "solver_selected"),
        ("authority", "unknown"),
        ("ordering_status", "ambiguous"),
        ("realized_sequence", "AAANNN"),
        ("notices", ["invalid"]),
        ("constraint_results", ["invalid"]),
        ("steps", ["invalid"]),
    ],
)
def test_python_and_json_plan_boundaries_agree(field: str, value: object):
    plan = reconstruct_playback(_realized())
    with pytest.raises((TypeError, ValueError)):
        replace(plan, **{field: value})
    payload = playback_plan_to_dict(plan)
    payload[field] = value
    with pytest.raises((TypeError, ValueError)):
        playback_plan_from_dict(payload)


def test_reveal_spans_cannot_repeat_previous_bases_or_reveal_unrelated_bases():
    plan = reconstruct_playback(_realized())
    for spans in (
        (CoordinateSpan(0, 3), CoordinateSpan(3, 6)),
        (CoordinateSpan(3, 6), CoordinateSpan(3, 6)),
        (CoordinateSpan(3, 5),),
        (),
    ):
        with pytest.raises(ValueError, match="added_spans"):
            replace(
                plan, steps=(plan.steps[0], replace(plan.steps[1], added_spans=spans))
            )


def test_coordinate_order_and_predecessor_must_match_v1_order():
    plan = reconstruct_playback(_realized())
    with pytest.raises(ValueError, match="order"):
        replace(
            plan,
            steps=(replace(plan.steps[1], index=0), replace(plan.steps[0], index=1)),
        )
    with pytest.raises(ValueError, match="predecessor"):
        replace(
            plan,
            steps=(
                plan.steps[0],
                replace(plan.steps[1], predecessor_placement_id=None),
            ),
        )


@pytest.mark.parametrize(
    "field,value",
    [
        ("constraint_id", 1),
        ("upstream_placement_id", "right"),
        ("min_distance_bp", False),
        ("min_distance_bp", -1),
        ("min_distance_bp", 1),
        ("max_distance_bp", 0.5),
        ("actual_distance_bp", 0.5),
        ("passed", 1),
        ("label", 1),
    ],
)
def test_constraint_result_local_contract(field: str, value: object):
    result = reconstruct_playback(_realized()).constraint_results[0]
    with pytest.raises((TypeError, ValueError)):
        replace(result, **{field: value})


def test_constraint_distance_must_match_actual_layout():
    plan = reconstruct_playback(_realized())
    # The evaluation is internally valid but contradicts the placement geometry.
    result = replace(plan.constraint_results[0], actual_distance_bp=1, passed=False)
    with pytest.raises(ValueError, match="actual_distance_bp"):
        replace(plan, constraint_results=(result,))


def test_plan_constraint_ids_are_unique():
    plan = reconstruct_playback(_realized())
    with pytest.raises(ValueError, match="constraint_id"):
        replace(plan, constraint_results=plan.constraint_results * 2)


def test_legitimate_failed_constraint_and_layout_only_plan_round_trip():
    realized = replace(
        _realized(),
        sequence="AAATTTCCC",
        placements=(
            _realized().placements[0],
            replace(_realized().placements[1], start=6),
        ),
    )
    plan = reconstruct_playback(realized)
    assert plan.ordering_status == "layout_only"
    assert plan.constraint_results[0].passed is False
    assert plan.constraint_results[0].actual_distance_bp == 3
    assert loads_playback_plan(dumps_playback_plan(plan)) == plan


def test_contained_placement_adds_no_revealed_bases():
    realized = RealizedArray(
        "contained",
        "AAAAC",
        (
            Placement("long", "long", PlacementKind.OTHER, "AAAAC", 0),
            Placement("short", "short", PlacementKind.OTHER, "AAA", 1),
        ),
    )
    plan = reconstruct_playback(realized)
    assert plan.ordering_status == "ambiguous"
    assert plan.steps[1].added_spans == ()
    assert loads_playback_plan(dumps_playback_plan(plan)) == plan


@pytest.mark.parametrize(
    "field,value", [("code", 1), ("message", None), ("level", "unknown")]
)
def test_notice_fields_are_validated(field: str, value: object):
    with pytest.raises((TypeError, ValueError)):
        replace(
            PlaybackNotice("caller_evidence", "Caller-supplied evidence."),
            **{field: value},
        )


def test_reconstruction_carries_explicit_evidence_without_inferring_recovery():
    realized = _realized()
    realized = replace(
        realized,
        placements=(
            replace(
                realized.placements[0], metadata={"coordinate_source": "offset_raw"}
            ),
            realized.placements[1],
        ),
    )
    plan = reconstruct_playback(realized)
    assert all(notice.code != "coordinate_recovered" for notice in plan.notices)
    evidence = PlaybackNotice(
        "adapter_evidence", "Adapter reports reviewed coordinate recovery."
    )
    plan = reconstruct_playback(realized, notices=(evidence,))
    assert evidence in plan.notices


@pytest.mark.parametrize("loader", [loads_realized_array, loads_playback_plan])
def test_strict_json_rejects_duplicate_object_keys(loader: Callable[[str], object]):
    payload = (
        dumps_realized_array(_realized())
        if loader is loads_realized_array
        else dumps_playback_plan(reconstruct_playback(_realized()))
    )
    with pytest.raises(ValueError, match="duplicate"):
        loader(
            payload.replace(
                '"source_id":"record:1"',
                '"source_id":"record:2","source_id":"record:1"',
            )
        )


@pytest.mark.parametrize("code", ["solver_selected", "ambiguous_order", "layout_only"])
def test_reserved_notices_cannot_contradict_actual_plan_evidence(code: str):
    plan = reconstruct_playback(_realized())
    with pytest.raises(ValueError, match="contradicts"):
        replace(plan, notices=(PlaybackNotice(code, "Claim."),))


def test_enum_strings_serialize_and_caller_text_identity_is_preserved():
    placement = Placement(" placement ", "feature", "tfbs", "aaa", 0, "fwd")
    realized = RealizedArray("record", "AAA", (placement,))
    assert realized.placements[0].placement_id == " placement "
    assert realized_array_to_dict(realized)["placements"][0]["kind"] == "tfbs"
    assert loads_realized_array(dumps_realized_array(realized)) == realized


@pytest.mark.parametrize("value", [[], None, 1, "record"])
@pytest.mark.parametrize("parser", [realized_array_from_dict, playback_plan_from_dict])
def test_dictionary_entrypoints_reject_non_objects(
    parser: Callable[[object], object], value: object
):
    with pytest.raises(TypeError, match="JSON object"):
        parser(value)
