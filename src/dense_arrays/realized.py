"""Public contracts for a persisted, realized dense array.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from typing import TYPE_CHECKING

from ._record_validation import (
    digest,
    enum_value,
    immutable_json_mapping,
    normalized_dna,
    records,
    required_text,
    validate_placement_sequence,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

REALIZED_ARRAY_SCHEMA_VERSION = "dense_arrays.realized_array.v1"


class PlacementKind(StrEnum):
    """Semantic kind of a feature placed on a realized sequence."""

    TFBS = "tfbs"
    FIXED_ELEMENT = "fixed_element"
    OTHER = "other"


class Orientation(StrEnum):
    """Orientation of the already-oriented placement sequence."""

    FORWARD = "fwd"
    REVERSE = "rev"
    UNSPECIFIED = "unspecified"


@dataclass(frozen=True, slots=True)
class Placement:
    """One feature placement in zero-based, half-open coordinates."""

    placement_id: str
    feature_id: str
    kind: PlacementKind
    sequence: str
    start: int
    orientation: Orientation = Orientation.UNSPECIFIED
    label: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate placement identity and freeze caller metadata."""
        object.__setattr__(
            self,
            "placement_id",
            required_text(self.placement_id, field_name="placement_id"),
        )
        object.__setattr__(
            self,
            "feature_id",
            required_text(self.feature_id, field_name="feature_id"),
        )
        object.__setattr__(
            self,
            "sequence",
            normalized_dna(self.sequence, field_name="placement.sequence"),
        )
        object.__setattr__(
            self,
            "kind",
            enum_value(self.kind, PlacementKind, field_name="placement.kind"),
        )
        object.__setattr__(
            self,
            "orientation",
            enum_value(
                self.orientation, Orientation, field_name="placement.orientation"
            ),
        )
        if (
            isinstance(self.start, bool)
            or not isinstance(self.start, int)
            or self.start < 0
        ):
            msg = "placement.start must be an integer >= 0"
            raise ValueError(msg)
        if self.label is not None:
            object.__setattr__(
                self, "label", required_text(self.label, field_name="label")
            )
        object.__setattr__(self, "metadata", immutable_json_mapping(self.metadata))

    @property
    def end(self) -> int:
        """Return the exclusive placement end coordinate."""
        return self.start + len(self.sequence)


@dataclass(frozen=True, slots=True)
class DeclaredConstraint:
    """A declared distance constraint between two realized placements."""

    constraint_id: str
    upstream_placement_id: str
    downstream_placement_id: str
    min_distance_bp: int
    max_distance_bp: int
    label: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Validate the referenced pair and allowed distance interval."""
        for field_name in (
            "constraint_id",
            "upstream_placement_id",
            "downstream_placement_id",
        ):
            object.__setattr__(
                self,
                field_name,
                required_text(getattr(self, field_name), field_name=field_name),
            )
        for field_name in ("min_distance_bp", "max_distance_bp"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                msg = f"{field_name} must be an integer >= 0"
                raise ValueError(msg)
        if self.min_distance_bp > self.max_distance_bp:
            msg = "min_distance_bp must be <= max_distance_bp"
            raise ValueError(msg)
        if self.upstream_placement_id == self.downstream_placement_id:
            msg = "a constraint must reference two different placements"
            raise ValueError(msg)
        if self.label is not None:
            object.__setattr__(
                self, "label", required_text(self.label, field_name="label")
            )
        object.__setattr__(self, "metadata", immutable_json_mapping(self.metadata))


@dataclass(frozen=True, slots=True)
class RealizedArray:
    """A sequence and the persisted placements known to realize it."""

    source_id: str
    sequence: str
    placements: tuple[Placement, ...]
    constraints: tuple[DeclaredConstraint, ...] = ()
    source_digest: str | None = None
    coordinate_space: str = "realized_sequence"
    provenance: Mapping[str, object] = field(default_factory=dict)
    schema_version: str = field(default=REALIZED_ARRAY_SCHEMA_VERSION, init=False)

    def __post_init__(self) -> None:
        """Validate the sequence, placement alignment, and constraint references."""
        object.__setattr__(
            self, "source_id", required_text(self.source_id, field_name="source_id")
        )
        object.__setattr__(
            self,
            "sequence",
            normalized_dna(self.sequence, field_name="realized_array.sequence"),
        )
        object.__setattr__(
            self,
            "placements",
            records(self.placements, Placement, field_name="placements"),
        )
        object.__setattr__(
            self,
            "constraints",
            records(self.constraints, DeclaredConstraint, field_name="constraints"),
        )
        object.__setattr__(
            self,
            "coordinate_space",
            required_text(self.coordinate_space, field_name="coordinate_space"),
        )
        if not self.placements:
            msg = "a realized array must contain at least one placement"
            raise ValueError(msg)
        placement_ids = [placement.placement_id for placement in self.placements]
        if len(placement_ids) != len(set(placement_ids)):
            msg = "placement_id values must be unique within a realized array"
            raise ValueError(msg)
        constraint_ids = [constraint.constraint_id for constraint in self.constraints]
        if len(constraint_ids) != len(set(constraint_ids)):
            msg = "constraint_id values must be unique within a realized array"
            raise ValueError(msg)
        for placement in self.placements:
            validate_placement_sequence(
                placement_id=placement.placement_id,
                start=placement.start,
                end=placement.end,
                sequence=placement.sequence,
                realized_sequence=self.sequence,
            )
        for constraint in self.constraints:
            missing = {
                constraint.upstream_placement_id,
                constraint.downstream_placement_id,
            } - set(placement_ids)
            if missing:
                msg = (
                    f"constraint {constraint.constraint_id!r} references unknown "
                    f"placements: {sorted(missing)}"
                )
                raise ValueError(msg)
        if self.source_digest is not None:
            object.__setattr__(
                self,
                "source_digest",
                digest(self.source_digest, field_name="source_digest"),
            )
        object.__setattr__(self, "provenance", immutable_json_mapping(self.provenance))
