"""Validated promoter and regulator requirements for motif packing.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Self

from .problem import discrete_integer

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
