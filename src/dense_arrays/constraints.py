"""
--------------------------------------------------------------------------------
<dense-array project>

Constraint helpers for dense-arrays.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Self


@dataclass
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
        self.upstream_index = upstream_index
        self.downstream_index = downstream_index
        self.upstream_pos = (
            upstream_pos
            if isinstance(upstream_pos, tuple)
            else (upstream_pos, upstream_pos)
        )
        self.downstream_pos = (
            downstream_pos
            if isinstance(downstream_pos, tuple)
            else (downstream_pos, downstream_pos)
        )
        self.spacer_length = (
            spacer_length
            if isinstance(spacer_length, tuple)
            else (spacer_length, spacer_length)
        )


def _normalize_regulator_mapping(
    nb_motifs: int,
    regulator_by_index: list[str] | dict[int, str],
) -> dict[int, str]:
    if isinstance(regulator_by_index, list):
        if len(regulator_by_index) != nb_motifs:
            msg = "regulator_by_index list length must match number of motifs"
            raise ValueError(msg)
        mapping = {i: str(label).strip() for i, label in enumerate(regulator_by_index)}
    elif isinstance(regulator_by_index, dict):
        if set(regulator_by_index.keys()) != set(range(nb_motifs)):
            msg = "regulator_by_index dict must cover all motif indices"
            raise ValueError(msg)
        mapping = {
            int(i): str(label).strip() for i, label in regulator_by_index.items()
        }
    else:
        msg = "regulator_by_index must be a list or dict"
        raise TypeError(msg)
    if any(not label for label in mapping.values()):
        msg = "regulator_by_index labels must be non-empty strings"
        raise ValueError(msg)
    return mapping


def _normalize_min_counts(
    mapping: dict[int, str],
    required: set[str],
    min_count_by_regulator: dict[str, int] | None,
) -> dict[str, int]:
    min_counts = {str(k): int(v) for k, v in (min_count_by_regulator or {}).items()}
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
    if min_required_regulators <= 0:
        msg = "min_required_regulators must be > 0 (use None to disable)."
        raise ValueError(msg)
    if min_required_regulators > len(available):
        msg = "min_required_regulators exceeds available regulators"
        raise ValueError(msg)
    return min_required_regulators
