"""Renderer-independent playback-plan contracts.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

from dense_arrays._record_validation import (
    digest,
    enum_value,
    integer,
    normalized_dna,
    records,
    required_text,
)
from dense_arrays.realized import Orientation, PlacementKind

from .validation import validate_plan

PLAYBACK_PLAN_SCHEMA_VERSION = "dense_arrays.playback_plan.v1"


class PlaybackAuthority(StrEnum):
    """Truth level of ordering and relation claims in a playback plan."""

    SOLVER_SELECTED = "solver_selected"
    PLACEMENT_RECONSTRUCTED = "placement_reconstructed"


class OrderingStatus(StrEnum):
    """Strength of the order inferred from persisted placements."""

    UNIQUE = "unique"
    AMBIGUOUS = "ambiguous"
    LAYOUT_ONLY = "layout_only"


class NoticeLevel(StrEnum):
    """Severity of a renderer-visible playback notice."""

    INFO = "info"
    WARNING = "warning"


@dataclass(frozen=True, slots=True)
class CoordinateSpan:
    """Zero-based, half-open sequence span."""

    start: int
    end: int

    def __post_init__(self) -> None:
        """Require a nonempty interval of nonnegative integer coordinates."""
        for field_name in ("start", "end"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                msg = f"{field_name} must be an integer >= 0"
                raise ValueError(msg)
        if self.end <= self.start:
            msg = "span.end must be greater than span.start"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class PlaybackStep:
    """One semantic placement event in a playback plan."""

    index: int
    placement_id: str
    feature_id: str
    start: int
    end: int
    placement_kind: str
    orientation: str
    placement_sequence: str
    added_spans: tuple[CoordinateSpan, ...]
    predecessor_placement_id: str | None = None
    relation_kind: str = "coordinate_precedence"
    label: str | None = None

    def __post_init__(self) -> None:
        """Validate placement fields and each locally declared reveal span."""
        for field_name in ("index", "start", "end"):
            integer(
                getattr(self, field_name), field_name=f"step.{field_name}", minimum=0
            )
        for field_name in ("placement_id", "feature_id", "relation_kind"):
            required_text(getattr(self, field_name), field_name=f"step.{field_name}")
        for field_name in ("label", "predecessor_placement_id"):
            value = getattr(self, field_name)
            if value is not None:
                required_text(value, field_name=f"step.{field_name}")
        object.__setattr__(
            self,
            "placement_kind",
            enum_value(
                self.placement_kind, PlacementKind, field_name="step.placement_kind"
            ).value,
        )
        object.__setattr__(
            self,
            "orientation",
            enum_value(
                self.orientation, Orientation, field_name="step.orientation"
            ).value,
        )
        object.__setattr__(
            self,
            "placement_sequence",
            normalized_dna(
                self.placement_sequence, field_name="step.placement_sequence"
            ),
        )
        if self.relation_kind != "coordinate_precedence":
            msg = "playback v1 supports only coordinate_precedence relations"
            raise ValueError(msg)
        if self.end <= self.start:
            msg = "step.end must be greater than step.start"
            raise ValueError(msg)
        if len(self.placement_sequence) != self.end - self.start:
            msg = "step placement_sequence length must equal end - start"
            raise ValueError(msg)
        spans = records(self.added_spans, CoordinateSpan, field_name="step.added_spans")
        previous_end = self.start
        for span in spans:
            if span.start < previous_end or span.end > self.end:
                msg = (
                    "step.added_spans must be ordered, disjoint, "
                    "and inside the placement"
                )
                raise ValueError(msg)
            previous_end = span.end
        object.__setattr__(self, "added_spans", spans)


@dataclass(frozen=True, slots=True)
class ConstraintResult:
    """Evaluation of one declared constraint on the realized layout."""

    constraint_id: str
    upstream_placement_id: str
    downstream_placement_id: str
    actual_distance_bp: int
    min_distance_bp: int
    max_distance_bp: int
    passed: bool
    label: str | None = None

    def __post_init__(self) -> None:
        """Check that the reported result agrees with its allowed distance."""
        for field_name in (
            "constraint_id",
            "upstream_placement_id",
            "downstream_placement_id",
        ):
            required_text(getattr(self, field_name), field_name=field_name)
        integer(self.actual_distance_bp, field_name="actual_distance_bp")
        for field_name in ("min_distance_bp", "max_distance_bp"):
            integer(getattr(self, field_name), field_name=field_name, minimum=0)
        if self.min_distance_bp > self.max_distance_bp:
            msg = "min_distance_bp must be <= max_distance_bp"
            raise ValueError(msg)
        if self.upstream_placement_id == self.downstream_placement_id:
            msg = "a constraint must reference two different placements"
            raise ValueError(msg)
        if self.label is not None:
            required_text(self.label, field_name="label")
        if not isinstance(self.passed, bool):
            msg = "constraint result passed must be a boolean"
            raise TypeError(msg)
        if self.passed != (
            self.min_distance_bp <= self.actual_distance_bp <= self.max_distance_bp
        ):
            msg = "constraint result passed must match the declared distance range"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class PlaybackNotice:
    """Structured qualification that renderers must preserve."""

    code: str
    message: str
    level: NoticeLevel = NoticeLevel.INFO

    def __post_init__(self) -> None:
        """Require an explicit code, message, and supported notice level."""
        required_text(self.code, field_name="notice.code")
        required_text(self.message, field_name="notice.message")
        object.__setattr__(
            self,
            "level",
            enum_value(self.level, NoticeLevel, field_name="notice.level"),
        )


@dataclass(frozen=True, slots=True)
class PlaybackPlan:
    """Immutable semantic timeline consumed by all playback renderers."""

    source_id: str
    source_digest: str | None
    realization_digest: str
    realized_sequence: str
    authority: PlaybackAuthority
    ordering_status: OrderingStatus
    steps: tuple[PlaybackStep, ...]
    constraint_results: tuple[ConstraintResult, ...] = ()
    notices: tuple[PlaybackNotice, ...] = ()
    schema_version: str = field(default=PLAYBACK_PLAN_SCHEMA_VERSION, init=False)

    def __post_init__(self) -> None:
        """Freeze input records and validate their combined evidence claims."""
        required_text(self.source_id, field_name="source_id")
        if self.source_digest is not None:
            object.__setattr__(
                self,
                "source_digest",
                digest(self.source_digest, field_name="source_digest"),
            )
        object.__setattr__(
            self,
            "realization_digest",
            digest(self.realization_digest, field_name="realization_digest"),
        )
        object.__setattr__(
            self,
            "realized_sequence",
            normalized_dna(self.realized_sequence, field_name="realized_sequence"),
        )
        object.__setattr__(
            self,
            "authority",
            enum_value(self.authority, PlaybackAuthority, field_name="authority"),
        )
        object.__setattr__(
            self,
            "ordering_status",
            enum_value(
                self.ordering_status, OrderingStatus, field_name="ordering_status"
            ),
        )
        object.__setattr__(
            self, "steps", records(self.steps, PlaybackStep, field_name="steps")
        )
        object.__setattr__(
            self,
            "constraint_results",
            records(
                self.constraint_results,
                ConstraintResult,
                field_name="constraint_results",
            ),
        )
        object.__setattr__(
            self, "notices", records(self.notices, PlaybackNotice, field_name="notices")
        )
        if not self.steps:
            msg = "a playback plan must contain at least one step"
            raise ValueError(msg)
        expected_indices = list(range(len(self.steps)))
        actual_indices = [step.index for step in self.steps]
        if actual_indices != expected_indices:
            msg = "playback step indices must be contiguous and zero-based"
            raise ValueError(msg)
        validate_plan(self)
