"""Renderer-neutral documents and visible placement evidence.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from .models import ConstraintResult, OrderingStatus, PlaybackPlan
from .theme import PlaybackPresentation, step_color, validate_color


def _freeze_overrides(
    values: Mapping[str, str], placement_ids: set[str], name: str
) -> Mapping[str, str]:
    if not isinstance(values, Mapping):
        msg = f"{name} must be a mapping keyed by placement ID"
        raise TypeError(msg)
    if set(values) - placement_ids:
        msg = f"{name} contains unknown placement IDs"
        raise ValueError(msg)
    if any(
        not isinstance(value, str) or not value.strip() for value in values.values()
    ):
        msg = f"{name} values must be nonempty text"
        raise ValueError(msg)
    return MappingProxyType(dict(values))


@dataclass(frozen=True, slots=True)
class PlaybackDocument:
    """One scene, caller-owned styles, and artifact metadata."""

    plan: PlaybackPlan
    title: str
    subtitle: str = ""
    label_overrides: Mapping[str, str] = field(default_factory=dict)
    presentation: PlaybackPresentation = field(default_factory=PlaybackPresentation)
    color_overrides: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate document fields and freeze caller-owned presentation maps."""
        if not isinstance(self.plan, PlaybackPlan) or not self.plan.steps:
            msg = "playback document requires a PlaybackPlan with steps"
            raise ValueError(msg)
        if not isinstance(self.title, str) or not self.title.strip():
            msg = "title must be nonempty text"
            raise ValueError(msg)
        if not isinstance(self.subtitle, str):
            msg = "subtitle must be text"
            raise TypeError(msg)
        if not isinstance(self.presentation, PlaybackPresentation):
            msg = "presentation must be PlaybackPresentation"
            raise TypeError(msg)
        placement_ids = {step.placement_id for step in self.plan.steps}
        for name in ("label_overrides", "color_overrides"):
            values = _freeze_overrides(getattr(self, name), placement_ids, name)
            object.__setattr__(self, name, values)
        for color in self.color_overrides.values():
            validate_color(color)

    def step_label(self, index: int) -> str:
        """Resolve the label for a placement from the caller's explicit map."""
        step = self.plan.steps[index]
        return self.label_overrides.get(
            step.placement_id, step.label or step.placement_id
        )

    def step_color(self, index: int) -> str:
        """Resolve a caller color, then the selected generic palette."""
        step = self.plan.steps[index]
        return self.color_overrides.get(
            step.placement_id, step_color(step, index, self.presentation.color_profile)
        )


@dataclass(frozen=True, slots=True)
class PlaybackEvidence:
    """Mandatory interpretation and optional provenance for one scene."""

    qualification: str
    constraints: tuple[str, ...]
    animate_chain: bool
    notices: tuple[str, ...]


def resolve_evidence(document: PlaybackDocument) -> PlaybackEvidence:
    """Resolve visible authority, ordering, and failed-constraint qualifications."""
    plan = document.plan
    order = {
        OrderingStatus.UNIQUE: "Unique coordinate order",
        OrderingStatus.AMBIGUOUS: "Ambiguous order; deterministic tie-break",
        OrderingStatus.LAYOUT_ONLY: "Layout only; gaps prevent a placement chain",
    }[plan.ordering_status]
    failures = tuple(
        _constraint_evidence(result)
        for result in plan.constraint_results
        if not result.passed
    )
    return PlaybackEvidence(
        qualification=f"Reconstructed from placements · {order}",
        constraints=failures,
        animate_chain=plan.ordering_status != OrderingStatus.LAYOUT_ONLY,
        notices=tuple(notice.message for notice in plan.notices)
        if document.presentation.show_authority_notice
        else (),
    )


def _constraint_evidence(result: ConstraintResult) -> str:
    status = "PASSED" if result.passed else "FAILED"
    return (
        f"{status} {result.constraint_id}: {result.actual_distance_bp} bp "
        f"(required {result.min_distance_bp}..{result.max_distance_bp} bp)"
    )


def evidence_metadata(document: PlaybackDocument) -> str:
    """Retain all distance results and selected notices beyond canvas summaries."""
    evidence = resolve_evidence(document)
    distances = tuple(
        _constraint_evidence(result) for result in document.plan.constraint_results
    )
    return "\n".join((evidence.qualification, *distances, *evidence.notices))


def collection_comment(documents: tuple[PlaybackDocument, ...]) -> str:
    """Describe each scene and its complete evidence in native media metadata."""
    return "\n\n".join(
        "\n".join(
            value
            for value in (
                document.title,
                document.subtitle,
                evidence_metadata(document),
            )
            if value
        )
        for document in documents
    )


@dataclass(frozen=True, slots=True)
class DistanceBracket:
    """A declared distance rendered in realized sequence coordinates."""

    start: int
    end: int
    label: str
    failed: bool


def resolve_distance_brackets(
    document: PlaybackDocument,
) -> tuple[DistanceBracket, ...]:
    """Resolve only declared constraints under the caller's bracket control."""
    if document.presentation.show_distance_bracket == "never":
        return ()
    steps = {step.placement_id: step for step in document.plan.steps}
    return tuple(
        DistanceBracket(
            steps[result.upstream_placement_id].end,
            steps[result.downstream_placement_id].start,
            f"{result.constraint_id}: {result.actual_distance_bp} bp "
            f"(required {result.min_distance_bp}..{result.max_distance_bp} bp)",
            not result.passed,
        )
        for result in document.plan.constraint_results
    )
