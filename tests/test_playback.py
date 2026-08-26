"""Contract tests for realized-placement playback."""

from __future__ import annotations

import json
import re
from typing import TYPE_CHECKING, ClassVar, Self

import matplotlib.animation as mpl_animation
import matplotlib.pyplot as plt
import pytest
from matplotlib.patches import FancyArrowPatch

from dense_arrays.playback import (
    OrderingStatus,
    PlaybackAuthority,
    dumps_playback_plan,
    dumps_realized_array,
    loads_playback_plan,
    loads_realized_array,
    matplotlib_renderer,
    reconstruct_playback,
    render_playback_html,
)
from dense_arrays.playback.graph_layout import journey_path_positions
from dense_arrays.playback.html import PlaybackDocument
from dense_arrays.playback.models import ConstraintResult
from dense_arrays.playback.theme import PlaybackPresentation
from dense_arrays.playback.timeline import complement_sequence, revealed_indices
from dense_arrays.realized import (
    DeclaredConstraint,
    Orientation,
    Placement,
    PlacementKind,
    RealizedArray,
)

if TYPE_CHECKING:
    from pathlib import Path


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


@pytest.mark.parametrize("invalid", ["false", "true", 0, 1, None, {}, []])
def test_contract_json_rejects_non_boolean_constraint_result(invalid: object) -> None:
    payload = json.loads(dumps_playback_plan(reconstruct_playback(_realized_array())))
    payload["constraint_results"][0]["passed"] = invalid

    with pytest.raises(TypeError, match="must be a JSON boolean"):
        loads_playback_plan(json.dumps(payload))


def test_direct_constraint_result_rejects_non_boolean_passed() -> None:
    with pytest.raises(TypeError, match="passed must be a boolean"):
        ConstraintResult(
            constraint_id="fixture",
            upstream_placement_id="left",
            downstream_placement_id="right",
            actual_distance_bp=0,
            min_distance_bp=0,
            max_distance_bp=0,
            passed="false",  # type: ignore[arg-type]
        )


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
    assert "revealed_indices" in artifact
    assert "timeline_frames" in artifact


def test_full_graph_draws_declared_constraint_relation() -> None:
    document = PlaybackDocument(
        plan=reconstruct_playback(_realized_array()),
        title="fixture",
        presentation=PlaybackPresentation(
            graph_detail="full", color_profile="constraints"
        ),
    )
    figure, axis = plt.subplots()

    matplotlib_renderer._draw_graph(  # noqa: SLF001
        axis, document, transition_index=0, progress=0.0
    )

    assert any(
        isinstance(patch, FancyArrowPatch) and patch.get_linewidth() == 1.7
        for patch in axis.patches
    )
    plt.close(figure)


def test_revealed_positions_follow_added_spans_not_whole_placements() -> None:
    plan = reconstruct_playback(
        RealizedArray(
            source_id="fixture#gapped",
            sequence="AAATTTCCC",
            placements=(
                _placement("left", "AAA", 0),
                _placement("right", "CCC", 6),
            ),
        )
    )
    document = matplotlib_renderer.PlaybackDocument(plan=plan, title="fixture")

    assert revealed_indices(document.plan.steps, 1) == (0, 1, 2, 6, 7, 8)


def test_overlap_reveal_mask_preserves_complete_placement_bars() -> None:
    plan = reconstruct_playback(
        RealizedArray(
            source_id="fixture#overlap",
            sequence="AAAT",
            placements=(
                _placement("first", "AAA", 0),
                _placement("overlap", "AAT", 1),
            ),
        )
    )
    document = matplotlib_renderer.PlaybackDocument(plan=plan, title="fixture")
    figure, axis = plt.subplots()

    matplotlib_renderer._draw_fallback_duplex(axis, document, 1)  # noqa: SLF001
    bar_widths = [patch.get_width() for patch in axis.patches]
    bar_labels = {text.get_text() for text in axis.texts}
    artifact = render_playback_html(plan, title="Overlap fixture")

    assert revealed_indices(plan.steps, 1) == (0, 1, 2, 3)
    assert bar_widths == [3, 3]
    assert {"AAA", "AAT"} <= bar_labels
    assert "step.added_spans.map" not in artifact
    assert "esc(step.placement_sequence)" in artifact
    plt.close(figure)


def test_iupac_complement_is_complete() -> None:
    assert complement_sequence("ATCGRYSWKMBDHVN") == "TAGCYRSWMKVHDBN"


def test_html_payload_uses_python_owned_iupac_complement() -> None:
    alphabet = "ATCGRYSWKMBDHVN"
    plan = reconstruct_playback(
        RealizedArray(
            source_id="fixture#iupac",
            sequence=alphabet,
            placements=(_placement("iupac", alphabet, 0),),
        )
    )

    artifact = render_playback_html(plan, title="IUPAC fixture")
    match = re.search(
        r'<script id="playback-data" type="application/json">(.*?)</script>',
        artifact,
    )

    assert match is not None
    payload = json.loads(match.group(1))
    assert payload[0]["complement_sequence"] == "TAGCYRSWMKVHDBN"


class _FrameCountingWriter:
    instances: ClassVar[list[_FrameCountingWriter]] = []

    def __init__(self, **_kwargs: object) -> None:
        self.frame_count = 0
        self.instances.append(self)

    def saving(self, *_args: object, **_kwargs: object) -> _FrameCountingWriter:
        return self

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_args: object) -> None:
        return None

    def grab_frame(self, **_kwargs: object) -> None:
        self.frame_count += 1


@pytest.mark.parametrize(
    "renderer_name,writer_name,transition_seconds,expected_frames",
    [
        ("render_collection_gif", "PillowWriter", 0.0, 20),
        ("render_collection_gif", "PillowWriter", 1.0, 24),
        ("render_collection_mp4", "FFMpegWriter", 0.0, 20),
        ("render_collection_mp4", "FFMpegWriter", 1.0, 24),
    ],
)
def test_collection_renderers_honor_scene_transition_seconds(  # noqa: PLR0913, PLR0917
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    renderer_name: str,
    writer_name: str,
    transition_seconds: float,
    expected_frames: int,
) -> None:
    _FrameCountingWriter.instances.clear()
    monkeypatch.setattr(mpl_animation, writer_name, _FrameCountingWriter)
    monkeypatch.setattr(
        matplotlib_renderer, "_draw_document", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        matplotlib_renderer,
        "_transition_frame_counts",
        lambda *_args, **_kwargs: (3, 5),
    )
    plan = reconstruct_playback(_realized_array())
    document = matplotlib_renderer.PlaybackDocument(plan=plan, title="fixture")
    renderer = getattr(matplotlib_renderer, renderer_name)

    renderer(
        (document, document),
        tmp_path / f"playback.{renderer_name[-3:]}",
        fps=2,
        seconds_per_step=1.0,
        hold_seconds=0.0,
        lead_seconds=0.0,
        scene_transition_seconds=transition_seconds,
    )

    assert _FrameCountingWriter.instances[-1].frame_count == expected_frames


def test_journey_layout_is_monotonic_and_slide_safe() -> None:
    positions = journey_path_positions(8)

    assert [x for x, _ in positions] == sorted(x for x, _ in positions)
    assert all(0.15 <= x <= 0.85 for x, _ in positions)
    assert all(0.15 <= y <= 0.85 for _, y in positions)
