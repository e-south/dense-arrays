"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/pools/filters.py

Typed selection of recorded preparation outcomes.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    integer,
    object_fields,
    records,
    required_text,
)
from dense_arrays.artifacts.preparation.sets import SetAccounting

if TYPE_CHECKING:
    from dense_arrays.artifacts.pool_records import PoolSummary
    from dense_arrays.artifacts.preparation.candidates import PoolCandidate

CANDIDATE_FILTER_SCHEMA = "dense_arrays.candidate-filter.v1"
OUTCOMES = {
    "retained",
    "not_selected",
    "duplicate_discarded",
    "eligibility_rejected",
    "execution_error",
}


@dataclass(frozen=True)
class CandidateFilter:
    """OR within selectors; AND across fields, with recipe-local score bands."""

    indices: tuple[int, ...] = ()
    outcomes: tuple[str, ...] = ()
    reasons: tuple[str, ...] = ()
    recipes: tuple[str, ...] = ()
    score_bands: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        """Freeze explicit selectors and reject ambiguous or unknown values."""
        for name in ("indices", "outcomes", "reasons", "recipes", "score_bands"):
            numeric = name in {"indices", "score_bands"}
            values = records(
                getattr(self, name), int if numeric else str, field_name=name
            )
            for value in values:
                if numeric:
                    integer(value, field_name=name, minimum=1)
                else:
                    required_text(value, field_name=name)
            if len(set(values)) != len(values):
                msg = f"{name} must not repeat values"
                raise ValueError(msg)
            object.__setattr__(self, name, values)
        if set(self.outcomes) - OUTCOMES:
            msg = "unknown candidate outcomes"
            raise ValueError(msg)

    @property
    def identities(self) -> int:
        """Count explicit selector values held for the query."""
        return (
            len(self.indices)
            + len(self.outcomes)
            + len(self.reasons)
            + len(self.recipes)
            + len(self.score_bands)
        )

    def validate(self, summary: PoolSummary) -> None:
        """Resolve indices and observed rejection reasons in the immutable pool."""
        if summary.preparation is None:
            msg = "pool has no saved candidate evidence"
            raise ValueError(msg)
        if any(index > summary.source_parts for index in self.indices):
            msg = "unknown candidate indices in this pool"
            raise ValueError(msg)
        if set(self.reasons) - set(summary.preparation.rejections):
            msg = "unknown rejection reasons in this pool"
            raise ValueError(msg)
        if self.recipes and (
            not isinstance(summary.preparation, SetAccounting)
            or set(self.recipes) - set(summary.preparation.recipes)
        ):
            msg = "unknown preparation recipe IDs in this pool"
            raise ValueError(msg)
        self._validate_score_bands(summary)

    def _validate_score_bands(self, summary: PoolSummary) -> None:
        if not self.score_bands:
            return
        accounting = summary.preparation
        if isinstance(accounting, SetAccounting):
            if len(self.recipes) != 1:
                msg = (
                    "score band queries require exactly one recipe in a preparation set"
                )
                raise ValueError(msg)
            accounting = accounting.recipes[self.recipes[0]]
        if accounting.score_bands is None or any(
            band > len(accounting.score_bands["bands"]) for band in self.score_bands
        ):
            msg = "unknown score bands for this preparation recipe"
            raise ValueError(msg)

    def matches(self, record: PoolCandidate) -> bool:
        """Evaluate saved decisions without sampling or scoring."""
        return (
            (not self.indices or record.candidate.index in self.indices)
            and (not self.recipes or record.candidate.recipe_id in self.recipes)
            and (not self.outcomes or record.outcome in self.outcomes)
            and (
                not self.score_bands or record.candidate.score_band in self.score_bands
            )
            and (
                not self.reasons
                or any(r in record.candidate.reasons for r in self.reasons)
            )
        )

    def to_dict(self) -> dict[str, object]:
        """Serialize a reusable predicate for both interfaces."""
        return {
            "schema": CANDIDATE_FILTER_SCHEMA,
            "indices": list(self.indices),
            "outcomes": list(self.outcomes),
            "reasons": list(self.reasons),
            **({"recipes": list(self.recipes)} if self.recipes else {}),
            **({"score_bands": list(self.score_bands)} if self.score_bands else {}),
        }

    @classmethod
    def from_dict(cls, value: object) -> CandidateFilter:
        """Read only the declared candidate-filter schema."""
        data = object_fields(
            value,
            {"schema", "indices", "outcomes", "reasons", "recipes", "score_bands"},
            "candidate filter",
        )
        if data.pop("schema", None) != CANDIDATE_FILTER_SCHEMA:
            msg = "unsupported candidate-filter schema"
            raise ValueError(msg)
        return cls(**data)
