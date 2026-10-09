"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/retention/selection.py

Deterministic candidate equivalence, representatives and ordered retention.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

from .bands import ScoreBands, assign_score_bands

if TYPE_CHECKING:
    from dense_arrays.parts.candidates import Candidate
    from dense_arrays.parts.motifs import Motif
    from dense_arrays.parts.preparation import Retention
    from dense_arrays.parts.sampling import Uniqueness


def select_candidates(
    candidates: tuple[Candidate, ...],
    uniqueness: Uniqueness,
    retention: Retention,
    *,
    motif: Motif | None = None,
    score_bands: ScoreBands | None = None,
) -> tuple[Candidate, ...]:
    """Choose best-score representatives, then retain by the declared policy."""
    selected = _select_candidates(candidates, uniqueness, retention, motif=motif)
    return (
        selected if score_bands is None else assign_score_bands(selected, score_bands)
    )


def _select_candidates(
    candidates: tuple[Candidate, ...],
    uniqueness: Uniqueness,
    retention: Retention,
    *,
    motif: Motif | None,
) -> tuple[Candidate, ...]:
    groups: dict[str, list[Candidate]] = {}
    for candidate in candidates:
        key = eligible_key(candidate, uniqueness)
        if key is None:
            continue
        groups.setdefault(key, []).append(candidate)
    representatives = {}
    for members in groups.values():
        if len({c.part.group for c in members}) > 1:
            msg = "equivalent candidates belong to different groups"
            raise ValueError(msg)
        chosen = min(members, key=_score_order)
        for candidate in members:
            representatives[candidate.index] = chosen.index
    selected = [c for c in candidates if representatives.get(c.index) == c.index]
    if retention.policy == "mmr":
        from .mmr import select_mmr  # noqa: PLC0415 - policy-specific dependency

        if motif is None:
            msg = "MMR requires a resolved motif"
            raise ValueError(msg)
        choices = select_mmr(
            [replace(c, representative=c.index) for c in selected], retention, motif
        )
        return tuple(
            choices[c.index]
            if c.index in choices
            else replace(c, representative=representatives.get(c.index))
            for c in candidates
        )
    selected.sort(
        key=_score_order if retention.policy == "top_score" else lambda c: c.index
    )
    ranks = {candidate.index: i for i, candidate in enumerate(selected, 1)}
    return tuple(
        replace(
            c,
            representative=representatives.get(c.index),
            rank=ranks.get(c.index),
            retained=c.index in ranks and ranks[c.index] <= retention.count,
        )
        for c in candidates
    )


def eligible_key(candidate: Candidate, uniqueness: Uniqueness) -> str | None:
    """Use the same admitted equivalence for mining progress and final retention."""
    if candidate.reasons or candidate.error:
        return None
    key = (
        candidate.part.sequence
        if uniqueness.key == "sequence"
        else candidate.part.core_sequence
    )
    if key is None:
        msg = "core uniqueness requires a scored core for every eligible candidate"
        raise ValueError(msg)
    return key


def _score_order(candidate: Candidate) -> tuple[float, int]:
    score = candidate.score
    return -(score.raw if score else 0), candidate.index
