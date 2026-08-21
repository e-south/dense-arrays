"""Contract tests for realized-placement playback."""

from __future__ import annotations

import json

import pytest

from dense_arrays.playback import (
    OrderingStatus,
    PlaybackAuthority,
    dumps_playback_plan,
    dumps_realized_array,
    loads_playback_plan,
    loads_realized_array,
    reconstruct_playback,
    render_playback_html,
)
from dense_arrays.playback.graph_layout import journey_path_positions
from dense_arrays.realized import (
    DeclaredConstraint,
    Orientation,
    Placement,
    PlacementKind,
    RealizedArray,
)


def _placement(
    placement_id: str,
    sequence: str,
    start: int,
    *,
    kind: PlacementKind = PlacementKind.TFBS,
    label: str | None = None,
) -> Placement:
    return Placement(
        placement_id=placement_id,
        feature_id=f"feature:{placement_id}",
        kind=kind,
        sequence=sequence,
        start=start,
        orientation=Orientation.FORWARD,
        label=label,
    )


def _realized_array() -> RealizedArray:
    placements = (
        _placement("p1", "AAA", 0, label="TF A"),
        _placement("p2", "CCC", 3, label="TF B"),
        _placement(
            "upstream",
            "GGG",
            6,
            kind=PlacementKind.FIXED_ELEMENT,
            label="upstream",
        ),
        _placement(
            "downstream",
            "TTT",
            9,
            kind=PlacementKind.FIXED_ELEMENT,
            label="downstream",
        ),
    )
    return RealizedArray(
        source_id="fixture#array-1",
        source_digest="a" * 64,
        sequence="AAACCCGGGTTT",
        placements=placements,
        constraints=(
            DeclaredConstraint(
                constraint_id="fixed-pair:0",
                upstream_placement_id="upstream",
                downstream_placement_id="downstream",
                min_distance_bp=0,
                max_distance_bp=0,
                label="fixed pair",
            ),
        ),
    )


def test_reconstruction_compiles_unique_truthful_plan() -> None:
    plan = reconstruct_playback(_realized_array())

    assert plan.authority is PlaybackAuthority.PLACEMENT_RECONSTRUCTED
    assert plan.ordering_status is OrderingStatus.UNIQUE
    assert [step.placement_id for step in plan.steps] == [
        "p1",
        "p2",
        "upstream",
        "downstream",
    ]
    assert plan.constraint_results[0].passed is True
    assert plan.constraint_results[0].actual_distance_bp == 0
    assert plan.notices[0].code == "placement_reconstructed"


def test_contract_json_round_trips_deterministically() -> None:
    realized = _realized_array()
    plan = reconstruct_playback(realized)

    realized_json = dumps_realized_array(realized)
    plan_json = dumps_playback_plan(plan)

    assert dumps_realized_array(loads_realized_array(realized_json)) == realized_json
    assert dumps_playback_plan(loads_playback_plan(plan_json)) == plan_json


def test_contract_json_rejects_unknown_fields() -> None:
    payload = json.loads(dumps_realized_array(_realized_array()))
    payload["optimizer"] = {"backend": "not-public"}

    with pytest.raises(ValueError, match="unknown keys"):
        loads_realized_array(json.dumps(payload))


def test_reconstruction_rejects_sequence_inconsistency() -> None:
    realized = RealizedArray(
        source_id="fixture#invalid",
        sequence="AAAAAA",
        placements=(_placement("p1", "AAAC", 0),),
    )

    with pytest.raises(ValueError, match="sequence-inconsistent"):
        reconstruct_playback(realized)


def test_reconstruction_marks_equal_starts_ambiguous() -> None:
    realized = RealizedArray(
        source_id="fixture#ambiguous",
        sequence="AAACCC",
        placements=(
            _placement("short", "AAA", 0),
            _placement("long", "AAACCC", 0),
        ),
    )

    plan = reconstruct_playback(realized)

    assert plan.ordering_status is OrderingStatus.AMBIGUOUS
    assert [step.placement_id for step in plan.steps] == ["short", "long"]
    assert any(notice.code == "ambiguous_order" for notice in plan.notices)


def test_reconstruction_marks_internal_gaps_layout_only() -> None:
    realized = RealizedArray(
        source_id="fixture#layout-only",
        sequence="AAATTTCCC",
        placements=(
            _placement("left", "AAA", 0),
            _placement("right", "CCC", 6),
        ),
    )

    plan = reconstruct_playback(realized)

    assert plan.ordering_status is OrderingStatus.LAYOUT_ONLY
    assert any(notice.code == "layout_only" for notice in plan.notices)


def test_html_is_self_contained_and_preserves_authority() -> None:
    plan = reconstruct_playback(_realized_array())

    artifact = render_playback_html(
        plan,
        title="Packing explanation",
        label_overrides={"upstream": "fixed upstream element"},
    )

    assert "Packing explanation" in artifact
    assert "dense_arrays.playback_plan.v1" in artifact
    assert "placement_reconstructed" in artifact
    assert '<script id="playback-data" type="application/json">' in artifact
    assert "https://" not in artifact


def test_journey_layout_is_monotonic_and_slide_safe() -> None:
    positions = journey_path_positions(8)

    assert [x for x, _ in positions] == sorted(x for x, _ in positions)
    assert all(0.15 <= x <= 0.85 for x, _ in positions)
    assert all(0.15 <= y <= 0.85 for _, y in positions)
