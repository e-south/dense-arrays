"""Public contracts for a persisted, realized dense array."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import StrEnum
from types import MappingProxyType

REALIZED_ARRAY_SCHEMA_VERSION = "dense_arrays.realized_array.v1"
_IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


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


def _required_text(value: str, *, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        msg = f"{field_name} must be a non-empty string"
        raise ValueError(msg)
    return text


def _normalized_dna(value: str, *, field_name: str) -> str:
    sequence = _required_text(value, field_name=field_name).upper()
    invalid = sorted(set(sequence) - _IUPAC_DNA)
    if invalid:
        msg = f"{field_name} contains non-IUPAC DNA symbols: {invalid}"
        raise ValueError(msg)
    return sequence


def _immutable_mapping(value: Mapping[str, object]) -> Mapping[str, object]:
    return MappingProxyType({str(key): item for key, item in value.items()})


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
        object.__setattr__(
            self,
            "placement_id",
            _required_text(self.placement_id, field_name="placement_id"),
        )
        object.__setattr__(
            self,
            "feature_id",
            _required_text(self.feature_id, field_name="feature_id"),
        )
        object.__setattr__(
            self,
            "sequence",
            _normalized_dna(self.sequence, field_name="placement.sequence"),
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
                self, "label", _required_text(self.label, field_name="label")
            )
        object.__setattr__(self, "metadata", _immutable_mapping(self.metadata))

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
        for field_name in (
            "constraint_id",
            "upstream_placement_id",
            "downstream_placement_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _required_text(getattr(self, field_name), field_name=field_name),
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
                self, "label", _required_text(self.label, field_name="label")
            )
        object.__setattr__(self, "metadata", _immutable_mapping(self.metadata))


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
        object.__setattr__(
            self, "source_id", _required_text(self.source_id, field_name="source_id")
        )
        object.__setattr__(
            self,
            "sequence",
            _normalized_dna(self.sequence, field_name="realized_array.sequence"),
        )
        object.__setattr__(self, "placements", tuple(self.placements))
        object.__setattr__(self, "constraints", tuple(self.constraints))
        object.__setattr__(
            self,
            "coordinate_space",
            _required_text(self.coordinate_space, field_name="coordinate_space"),
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
        if self.source_digest is not None:
            digest = self.source_digest.lower().strip()
            if len(digest) != 64 or any(
                char not in "0123456789abcdef" for char in digest
            ):
                msg = "source_digest must be a lowercase SHA-256 hex digest"
                raise ValueError(msg)
            object.__setattr__(self, "source_digest", digest)
        object.__setattr__(self, "provenance", _immutable_mapping(self.provenance))
