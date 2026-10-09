"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/eligibility.py

Candidate eligibility from recorded scores and shared sequence requirements.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from .screening import PWMExclusion
from .screening.sequence import passes_sequence

if TYPE_CHECKING:
    from dense_arrays.constraints import GC, Avoid

    from .candidates import Candidate
    from .sampling import Eligibility


def rejection_reasons(
    candidate: Candidate,
    eligibility: Eligibility,
    screening: tuple[Avoid | GC | PWMExclusion, ...],
    *,
    requires_hit: bool,
) -> tuple[str, ...]:
    """Evaluate every applicable screen and keep absent hits separate from scores."""
    reasons = []
    if requires_hit:
        score = candidate.score
        if score is None:
            reasons.append("no_qualifying_hit")
        elif (
            eligibility.best_hit_score_min_exclusive is not None
            and score.raw <= eligibility.best_hit_score_min_exclusive
        ):
            reasons.append("best_hit_score")
    for rule in screening:
        if isinstance(rule, PWMExclusion):
            observed = tuple(
                item for item in candidate.screening if item.rule_id == rule.id
            )
            if len(observed) != len(rule.motifs):
                msg = "PWM exclusion requires complete recorded observations"
                raise ValueError(msg)
            # Evaluate all observations so an undefined score ratio cannot hide.
            rejected = [rule.rejects(item.hit) for item in observed]
            if any(rejected):
                reasons.append(rule.id)
        elif not passes_sequence(rule, candidate.part.sequence):
            reasons.append(rule.id)
    return tuple(reasons)
