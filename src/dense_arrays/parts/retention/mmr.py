"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/retention/mmr.py

Greedy MMR with linear retained state and explicit core-distance semantics.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

import numpy as np

from .models import MMRDecision

if TYPE_CHECKING:
    from dense_arrays.parts.candidates import Candidate
    from dense_arrays.parts.motifs import Motif
    from dense_arrays.parts.preparation import Retention


ZERO_WEIGHT_TOLERANCE = 1e-6


def _weights(motif: Motif) -> np.ndarray:
    probabilities = np.array(motif.probabilities, dtype=float)
    background = np.array(motif.background, dtype=float)
    background /= background.sum()
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    log_ratio = np.zeros_like(probabilities)
    positive = probabilities > 0
    np.log2(probabilities, out=log_ratio, where=positive)
    log_ratio -= np.log2(background)
    information = (probabilities * log_ratio).sum(axis=1)
    maximum = -np.log2(background.min())
    weights = 1 - np.clip(information / maximum, 0, 1)
    return np.ones(motif.width) if weights.sum() <= ZERO_WEIGHT_TOLERANCE else weights


def _percentiles(scores: np.ndarray) -> np.ndarray:
    if len(scores) == 1 or scores.min() == scores.max():
        return np.ones(len(scores))
    unique, inverse, counts = np.unique(scores, return_inverse=True, return_counts=True)
    del unique
    ends = np.cumsum(counts) - 1
    starts = ends - counts + 1
    ranks = (starts + ends) / (2 * (len(scores) - 1))
    return ranks[inverse]


def select_mmr(
    candidates: list[Candidate], retention: Retention, motif: Motif
) -> dict[int, Candidate]:
    """Rank admitted representatives and record every pool/choice decision."""
    policy = retention.mmr
    pool_limit = policy.pool_limit(retention.count)
    if any(c.score is None for c in candidates):
        msg = "MMR requires scored candidates"
        raise ValueError(msg)
    ranked = sorted(
        candidates,
        key=lambda c: (-c.score.raw, c.part.core_sequence, c.part.sequence, c.index),
    )
    decisions = {}
    pool = []
    for candidate in ranked:
        hit = candidate.score
        if hit is None or (
            hit.fraction_of_max is None
            and (
                policy.score_scaling == "fraction_of_max_clipped"
                or policy.minimum_fraction_of_max is not None
            )
        ):
            msg = "MMR requires scored candidates with a positive theoretical maximum"
            raise ValueError(msg)
        if len(hit.core) != motif.width:
            msg = "MMR cores must have the motif's uniform core length"
            raise ValueError(msg)
        if (
            policy.minimum_fraction_of_max is not None
            and max(0, hit.fraction_of_max) < policy.minimum_fraction_of_max
        ):
            status = "below_score"
        elif len(pool) >= pool_limit:
            status = "beyond_limit"
        else:
            pool.append(candidate)
            continue
        decisions[candidate.index] = replace(candidate, selection=MMRDecision(status))
    if not pool:
        return decisions
    scores = np.array([c.score.raw for c in pool])
    relevance = (
        _percentiles(scores)
        if policy.score_scaling == "score_percentile"
        else np.array([max(0, c.score.fraction_of_max) for c in pool])
    )
    for i, candidate in enumerate(pool):
        decisions[candidate.index] = replace(
            candidate, selection=MMRDecision("included", float(relevance[i]))
        )
    cores = np.array([list(c.part.core_sequence) for c in pool])
    weights = _weights(motif)
    selected = np.zeros(len(pool), dtype=bool)
    minimum_distance = np.full(len(pool), np.inf)
    for rank in range(1, min(retention.count, len(pool)) + 1):
        similarity = np.zeros(len(pool)) if rank == 1 else 1 / (1 + minimum_distance)
        utility = (
            policy.relevance_weight * relevance
            - (1 - policy.relevance_weight) * similarity
        )
        utility[selected] = -np.inf
        # The pre-ranked pool implements raw-score/core/sequence/index tie breaking.
        chosen = int(np.argmax(utility))
        candidate = pool[chosen]
        decisions[candidate.index] = replace(
            candidate,
            rank=rank,
            retained=True,
            selection=MMRDecision(
                "included",
                float(relevance[chosen]),
                float(utility[chosen]),
                None if rank == 1 else float(minimum_distance[chosen]),
                None if rank == 1 else float(similarity[chosen]),
            ),
        )
        selected[chosen] = True
        distances = ((cores != cores[chosen]) * weights).sum(axis=1)
        np.minimum(minimum_distance, distances, out=minimum_distance)
    return decisions
