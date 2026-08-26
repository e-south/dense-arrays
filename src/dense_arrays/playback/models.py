"""Renderer-independent playback-plan contracts."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

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
        if (
            isinstance(self.index, bool)
            or not isinstance(self.index, int)
            or self.index < 0
        ):
            msg = "step.index must be an integer >= 0"
            raise ValueError(msg)
        if self.end <= self.start:
            msg = "step.end must be greater than step.start"
            raise ValueError(msg)
        if len(self.placement_sequence) != self.end - self.start:
            msg = "step placement_sequence length must equal end - start"
            raise ValueError(msg)
        object.__setattr__(self, "added_spans", tuple(self.added_spans))


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
        if not isinstance(self.passed, bool):
            msg = "constraint result passed must be a boolean"
            raise TypeError(msg)


@dataclass(frozen=True, slots=True)
class PlaybackNotice:
    """Structured qualification that renderers must preserve."""

    code: str
    message: str
    level: NoticeLevel = NoticeLevel.INFO


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
        object.__setattr__(self, "steps", tuple(self.steps))
        object.__setattr__(self, "constraint_results", tuple(self.constraint_results))
        object.__setattr__(self, "notices", tuple(self.notices))
        if not self.steps:
            msg = "a playback plan must contain at least one step"
            raise ValueError(msg)
        expected_indices = list(range(len(self.steps)))
        actual_indices = [step.index for step in self.steps]
        if actual_indices != expected_indices:
            msg = "playback step indices must be contiguous and zero-based"
            raise ValueError(msg)
