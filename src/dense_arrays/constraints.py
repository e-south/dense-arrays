"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/constraints.py

Validated promoter and regulator requirements for motif packing.

Module Author(s): Virgile Andreani, Eric J. South
Maintainer(s): Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from numbers import Real
from typing import Self

from ._record_validation import integer, required_text
from .problem import discrete_integer, motif_library

_INTERVAL_ENDPOINTS = 2


@dataclass(frozen=True)
class PromoterConstraint:
    """Promoter constraint (up/downstream elements, positions and spacing)."""

    upstream_index: int
    downstream_index: int
    upstream_pos: tuple[int | None, int | None]
    downstream_pos: tuple[int | None, int | None]
    spacer_length: tuple[int | None, int | None]

    def __init__(
        self: Self,
        *,
        upstream_index: int,
        downstream_index: int,
        upstream_pos: int | tuple[int | None, int | None] | None = None,
        downstream_pos: int | tuple[int | None, int | None] | None = None,
        spacer_length: int | tuple[int | None, int | None] | None = None,
    ) -> None:
        upstream_index = discrete_integer(upstream_index, "upstream_index", minimum=0)
        downstream_index = discrete_integer(
            downstream_index, "downstream_index", minimum=0
        )
        if upstream_index == downstream_index:
            msg = "Promoter indices must identify distinct library entries"
            raise ValueError(msg)
        object.__setattr__(self, "upstream_index", upstream_index)
        object.__setattr__(self, "downstream_index", downstream_index)
        object.__setattr__(
            self, "upstream_pos", _interval(upstream_pos, "upstream_pos", minimum=0)
        )
        object.__setattr__(
            self,
            "downstream_pos",
            _interval(downstream_pos, "downstream_pos", minimum=0),
        )
        object.__setattr__(
            self, "spacer_length", _interval(spacer_length, "spacer_length")
        )


@dataclass(frozen=True)
class RegulatorRequirements:
    """A caller-independent snapshot of entry labels and coverage requirements."""

    mapping: tuple[tuple[int, str], ...]
    min_counts: tuple[tuple[str, int], ...]
    min_required: int | None


@dataclass(frozen=True)
class CountConstraint:
    """Inclusive occurrence bounds on a validated set of library indices."""

    indices: tuple[int, ...]
    minimum: int | None
    maximum: int | None


@dataclass(frozen=True)
class CoverageConstraint:
    """Minimum represented groups, each containing validated library indices."""

    groups: tuple[tuple[int, ...], ...]
    minimum: int


@dataclass(frozen=True)
class FixedOccurrence:
    """One supplied identity, a declared strand and an inclusive start window."""

    index: int
    orientation: str
    start: tuple[int | None, int | None]
    origin: str = "start"


@dataclass(frozen=True)
class SpacingConstraint:
    """Signed downstream-start minus upstream-end bounds for fixed identities."""

    upstream: int
    downstream: int
    interval: tuple[int, int]


def occurrence_indices(indices: list[int], available: int) -> tuple[int, ...]:
    """Validate an explicit nonempty set of supplied occurrence indices."""
    if not isinstance(indices, (list, tuple)) or not indices:
        msg = "occurrence indices must be a nonempty list or tuple"
        raise ValueError(msg)
    result = tuple(discrete_integer(i, "index", minimum=0) for i in indices)
    if max(result) >= available or len(set(result)) != len(result):
        msg = "occurrence indices must be unique and within the available library"
        raise ValueError(msg)
    return result


def count_bounds(
    minimum: int | None, maximum: int | None, available: int
) -> tuple[int | None, int | None]:
    """Validate inclusive integer count bounds without conflating omission/zero."""
    if minimum is None and maximum is None:
        msg = "at least one count bound is required"
        raise ValueError(msg)
    low = None if minimum is None else discrete_integer(minimum, "minimum", minimum=0)
    high = None if maximum is None else discrete_integer(maximum, "maximum", minimum=0)
    if low is not None and low > available:
        msg = f"minimum exceeds {available} available occurrences"
        raise ValueError(msg)
    if low is not None and high is not None and low > high:
        msg = "minimum must not exceed maximum"
        raise ValueError(msg)
    return low, high


def _interval(
    value: object, name: str, *, minimum: int | None = None
) -> tuple[int | None, int | None]:
    if isinstance(value, tuple):
        if len(value) != _INTERVAL_ENDPOINTS:
            msg = f"{name} must be an integer, None, or a two-item tuple"
            raise ValueError(msg)
        lower, upper = value
    else:
        lower = upper = value
    lower = None if lower is None else discrete_integer(lower, name, minimum=minimum)
    upper = None if upper is None else discrete_integer(upper, name, minimum=minimum)
    if lower is not None and upper is not None and lower > upper:
        msg = f"{name} minimum must not exceed maximum"
        raise ValueError(msg)
    return lower, upper


def _normalize_regulator_mapping(
    nb_motifs: int,
    regulator_by_index: list[str] | dict[int, str],
) -> dict[int, str]:
    if isinstance(regulator_by_index, list):
        if len(regulator_by_index) != nb_motifs:
            msg = "regulator_by_index list length must match number of motifs"
            raise ValueError(msg)
        mapping = dict(enumerate(regulator_by_index))
    elif isinstance(regulator_by_index, dict):
        for index in regulator_by_index:
            discrete_integer(index, "regulator_by_index index", minimum=0)
        if set(regulator_by_index.keys()) != set(range(nb_motifs)):
            msg = "regulator_by_index dict must cover all motif indices"
            raise ValueError(msg)
        mapping = dict(regulator_by_index)
    else:
        msg = "regulator_by_index must be a list or dict"
        raise TypeError(msg)
    if any(
        not isinstance(label, str) or not label or label.strip() != label
        for label in mapping.values()
    ):
        msg = "regulator_by_index labels must be non-empty strings"
        raise ValueError(msg)
    return mapping


def _normalize_min_counts(
    mapping: dict[int, str],
    required: set[str],
    min_count_by_regulator: dict[str, int] | None,
) -> dict[str, int]:
    min_counts = {
        key: discrete_integer(value, "min_count_by_regulator")
        for key, value in (min_count_by_regulator or {}).items()
    }
    for regulator, count in min_counts.items():
        if count <= 0:
            msg = f"min_count_by_regulator must be > 0 (got {regulator}={count})"
            raise ValueError(msg)
        if regulator not in set(mapping.values()):
            msg = f"min_count_by_regulator regulator not in mapping: {regulator}"
            raise ValueError(msg)
    for regulator in required:
        min_counts[regulator] = max(1, min_counts.get(regulator, 1))

    counts = Counter(mapping.values())
    for regulator, count in min_counts.items():
        if counts[regulator] < count:
            msg = (
                f"Regulator '{regulator}' has only {counts[regulator]} motifs, "
                f"cannot satisfy min_count={count}."
            )
            raise ValueError(msg)
    return min_counts


def _normalize_min_required(
    min_required_regulators: int | None,
    available: set[str],
) -> int | None:
    if min_required_regulators is None:
        return None
    min_required_regulators = discrete_integer(
        min_required_regulators, "min_required_regulators"
    )
    if min_required_regulators <= 0:
        msg = "min_required_regulators must be > 0 (use None to disable)."
        raise ValueError(msg)
    if min_required_regulators > len(available):
        msg = "min_required_regulators exceeds available regulators"
        raise ValueError(msg)
    return min_required_regulators


@dataclass(frozen=True)
class Length:
    """A packing maximum or explicit final-length requirement, never both."""

    maximum: int | None = None
    exact: int | None = None

    def __post_init__(self) -> None:
        """Require exactly one positive integer length."""
        if (self.maximum is None) == (self.exact is None):
            msg = "length requires exactly one of maximum or exact"
            raise ValueError(msg)
        for name in ("maximum", "exact"):
            if (value := getattr(self, name)) is not None:
                integer(value, field_name=f"length.{name}", minimum=1)


@dataclass(frozen=True)
class Avoid:
    """Exclude literal DNA matches except those wholly within named fixed intervals."""

    id: str
    patterns: tuple[str, ...]
    strands: str = "both"
    except_placements: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Freeze strict ACGT patterns and explicit interval exception identities."""
        required_text(self.id, field_name="requirement.id")
        object.__setattr__(self, "patterns", tuple(motif_library(self.patterns)))
        if len(set(self.patterns)) != len(self.patterns):
            msg = "avoid patterns must be unique"
            raise ValueError(msg)
        if self.strands not in {"forward", "both"}:
            msg = "avoid strands must be forward or both"
            raise ValueError(msg)
        if not isinstance(self.except_placements, (list, tuple)):
            msg = "except_placements must be a sequence of fixed part IDs"
            raise TypeError(msg)
        for name in self.except_placements:
            required_text(name, field_name="except_placements")
        if len(set(self.except_placements)) != len(self.except_placements):
            msg = "except_placements must contain unique fixed part IDs"
            raise ValueError(msg)
        object.__setattr__(self, "except_placements", tuple(self.except_placements))


@dataclass(frozen=True)
class GC:
    """Inclusive GC fraction bounds on the final sequence or added padding."""

    id: str
    scope: str
    min: float
    max: float

    def __post_init__(self) -> None:
        """Require finite fractions and a declared counting region."""
        required_text(self.id, field_name="requirement.id")
        if self.scope not in {"sequence", "padding"}:
            msg = "gc scope must be sequence or padding"
            raise ValueError(msg)
        for name in ("min", "max"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real):
                msg = f"gc.{name} must be a finite fraction"
                raise TypeError(msg)
            if not math.isfinite(value) or not 0 <= value <= 1:
                msg = f"gc.{name} must be a finite fraction in [0, 1]"
                raise ValueError(msg)
        if self.min > self.max:
            msg = "gc minimum exceeds maximum"
            raise ValueError(msg)
