"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/retention/bands.py

Empirical score bands on eligible representatives, with intact boundary ties.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import pairwise
from math import ceil
from typing import TYPE_CHECKING

from dense_arrays._record_validation import object_fields
from dense_arrays.parts.scoring.configuration import finite

if TYPE_CHECKING:
    from dense_arrays.parts.candidates import Candidate

BAND_POLICY = "upper_rank_include_ties.v1"


@dataclass(frozen=True)
class ScoreBands:
    """Cumulative upper fractions of score-ranked eligible representatives.

    Fractions are strictly increasing in (0, 1). The final remainder is implicit.
    A boundary uses rank ceil(fraction * total) and includes every equal score.
    These bands describe one recipe's empirical scores, not biological activity.
    """

    upper_fractions: tuple[float, ...]

    def __post_init__(self) -> None:
        """Require explicit cumulative fractions without inferred tier defaults."""
        if (
            not isinstance(self.upper_fractions, (list, tuple))
            or not self.upper_fractions
        ):
            msg = "score_bands.upper_fractions requires a nonempty ordered array"
            raise ValueError(msg)
        values = tuple(
            finite(v, "score_bands.upper_fractions") for v in self.upper_fractions
        )
        if any(not 0 < v < 1 for v in values) or any(
            a >= b for a, b in pairwise(values)
        ):
            msg = "score band upper fractions must be strictly increasing in (0, 1)"
            raise ValueError(msg)
        object.__setattr__(self, "upper_fractions", values)

    def to_dict(self) -> dict[str, object]:
        """Encode cumulative fractions; the plan binds the classification version."""
        return {"upper_fractions": list(self.upper_fractions)}

    @classmethod
    def from_dict(cls, value: object) -> ScoreBands:
        """Parse only the declared score-band configuration."""
        return cls(**object_fields(value, {"upper_fractions"}, "score bands"))


def score_cutoffs(
    candidates: tuple[Candidate, ...], policy: ScoreBands
) -> tuple[float | None, ...]:
    """Resolve boundaries from every eligible representative, before pool caps."""
    scores = []
    for candidate in candidates:
        if candidate.representative != candidate.index:
            continue
        score = candidate.score
        if score is None:
            msg = "score bands require scored eligible representatives"
            raise ValueError(msg)
        scores.append(score.raw)
    scores.sort(reverse=True)
    return tuple(
        scores[ceil(Fraction(str(fraction)) * len(scores)) - 1] if scores else None
        for fraction in policy.upper_fractions
    )


def assign_score_bands(
    candidates: tuple[Candidate, ...], policy: ScoreBands
) -> tuple[Candidate, ...]:
    """Annotate representatives without changing eligibility or retained identity."""
    cutoffs = score_cutoffs(candidates, policy)
    thresholds = [-value for value in cutoffs if value is not None]
    return tuple(
        replace(c, score_band=bisect_left(thresholds, -c.score.raw) + 1)
        if c.representative == c.index
        else c
        for c in candidates
    )
