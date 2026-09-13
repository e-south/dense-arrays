"""Contract tests for realized-placement playback."""

from __future__ import annotations

import json
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
    export,
    gif_writer,
    loads_playback_plan,
    loads_realized_array,
    matplotlib_renderer,
    reconstruct_playback,
)
from dense_arrays.playback.duplex_drawing import draw_duplex
from dense_arrays.playback.graph_drawing import draw_graph
from dense_arrays.playback.graph_layout import journey_path_positions
from dense_arrays.playback.models import ConstraintResult
from dense_arrays.playback.presentation import PlaybackDocument
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
    with pytest.raises(ValueError, match="sequence-inconsistent"):
        RealizedArray(
            source_id="fixture#invalid",
            sequence="AAAAAA",
            placements=(_placement("p1", "AAAC", 0),),
        )


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


def test_full_graph_draws_declared_constraint_relation() -> None:
    document = PlaybackDocument(
        plan=reconstruct_playback(_realized_array()),
        title="fixture",
        presentation=PlaybackPresentation(
            graph_detail="full", color_profile="constraints"
        ),
    )
    figure, axis = plt.subplots()

    draw_graph(axis, document, transition_index=0, progress=0.0)

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

    draw_duplex(axis, document, 1)
    bar_widths = [patch.get_width() for patch in axis.patches]
    bar_labels = [
        "".join(
            text.get_text()
            for text in axis.texts
            if text.get_position()[1]
            == pytest.approx(patch.get_y() + patch.get_height() / 2)
        )
        for patch in axis.patches
    ]

    assert revealed_indices(plan.steps, 1) == (0, 1, 2, 3)
    assert bar_widths == [3, 3]
    assert bar_labels == ["AAA", "AAT"]
    plt.close(figure)


def test_iupac_complement_is_complete() -> None:
    assert complement_sequence("ATCGRYSWKMBDHVN") == "TAGCYRSWMKVHDBN"


class _FrameCountingWriter:
    instances: ClassVar[list[_FrameCountingWriter]] = []

    def __init__(self, **_kwargs: object) -> None:
        self.frame_count = 0
        self.instances.append(self)

    @classmethod
    def isAvailable(cls) -> bool:  # noqa: N802 - Matplotlib writer interface
        return True

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
        ("render_collection_gif", "EvidencePillowWriter", 0.0, 18),
        ("render_collection_gif", "EvidencePillowWriter", 1.0, 20),
        ("render_collection_mp4", "FFMpegWriter", 0.0, 18),
        ("render_collection_mp4", "FFMpegWriter", 1.0, 20),
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
    writer_module = (
        gif_writer if renderer_name == "render_collection_gif" else mpl_animation
    )
    monkeypatch.setattr(writer_module, writer_name, _FrameCountingWriter)
    monkeypatch.setattr(export, "draw_document", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        export,
        "transition_frame_counts",
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
