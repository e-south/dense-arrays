"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/sampling.py

Candidate mining, eligibility and uniqueness contracts, independent of execution.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass, field
from fractions import Fraction
from math import ceil

from dense_arrays._record_validation import integer, required_text
from dense_arrays.constraints import Length
from dense_arrays.parts.background.contracts import ConditionalLimits
from dense_arrays.parts.mining import POLICIES, sample_length
from dense_arrays.parts.motifs.models import numeric_row
from dense_arrays.parts.scoring.configuration import finite

_DEFAULT_BATCH_BASES = 1_000_000
_DEFAULT_TOTAL_BASES = 100_000_000


@dataclass(frozen=True)
class Background:
    """Declared independent-base DNA distribution and output group."""

    base_probabilities: tuple[float, ...] = (0.25, 0.25, 0.25, 0.25)
    group: str = "background"

    def __post_init__(self) -> None:
        """Require a normalized ACGT distribution and a usable group label."""
        object.__setattr__(
            self,
            "base_probabilities",
            numeric_row(self.base_probabilities, probability=True),
        )
        required_text(self.group, field_name="background.group")


@dataclass(frozen=True)
class LengthRange:
    """Inclusive positive integer bounds for a uniform length prior."""

    minimum: int
    maximum: int

    def __post_init__(self) -> None:
        """Reject reversed, fractional, boolean or unbounded length intervals."""
        integer(self.minimum, field_name="length.minimum", minimum=1)
        integer(self.maximum, field_name="length.maximum", minimum=self.minimum)


@dataclass(frozen=True)
class Sampling:
    """Full-sequence proposal strategy, length and optional sampling background."""

    length: Length | LengthRange
    strategy: str = "stochastic"
    base_probabilities: tuple[float, ...] | None = None
    limits: ConditionalLimits | None = None

    def __post_init__(self) -> None:
        """Reject implicit lengths or unsupported proposal strategies."""
        if not isinstance(self.length, (Length, LengthRange)) or (
            isinstance(self.length, Length) and self.length.exact is None
        ):
            msg = "sampling requires Length(exact=...) or LengthRange(minimum, maximum)"
            raise ValueError(msg)
        if self.strategy not in POLICIES:
            msg = (
                "supported sampling strategies: stochastic, consensus, "
                "background, conditional"
            )
            raise ValueError(msg)
        if self.strategy == "conditional":
            if self.limits is None:
                object.__setattr__(self, "limits", ConditionalLimits())
            elif not isinstance(self.limits, ConditionalLimits):
                msg = "conditional sampling requires ConditionalLimits"
                raise TypeError(msg)
        elif self.limits is not None:
            msg = "sampling limits apply only to conditional sampling"
            raise ValueError(msg)
        if self.base_probabilities is not None:
            object.__setattr__(
                self,
                "base_probabilities",
                numeric_row(self.base_probabilities, probability=True),
            )

    @property
    def minimum_length(self) -> int:
        """Smallest admitted full-sequence length."""
        return (
            self.length.minimum
            if isinstance(self.length, LengthRange)
            else self.length.exact
        )

    @property
    def maximum_length(self) -> int:
        """Largest admitted full-sequence length, used for work admission."""
        return (
            self.length.maximum
            if isinstance(self.length, LengthRange)
            else self.length.exact
        )

    def candidate_length(self, *, seed: int, index: int) -> int:
        """Resolve one independent length draw; exact lengths consume no entropy."""
        if self.strategy == "conditional":
            msg = "conditional lengths must be drawn jointly with their sequence"
            raise ValueError(msg)
        if isinstance(self.length, LengthRange):
            return sample_length(
                minimum=self.length.minimum,
                maximum=self.length.maximum,
                seed=seed,
                index=index,
            )
        return self.length.exact

    @property
    def policy(self) -> str:
        """Bind the versioned random stream and proposal interpretation."""
        return POLICIES[self.strategy]

    @property
    def records_proposal(self) -> bool:
        """Preserve earlier stochastic records while binding explicit new proposals."""
        return self.strategy != "stochastic" or self.base_probabilities is not None


@dataclass(frozen=True)
class MiningTarget:
    """Desired eligible supply, separate from effort caps and retained-part counts."""

    eligible_unique: int | None = None
    minimum_candidates: int = 0
    max_retained_fraction: float | None = None

    def __post_init__(self) -> None:
        """Require an explicit positive supply target and nonnegative effort floor."""
        if (self.eligible_unique is None) == (self.max_retained_fraction is None):
            msg = (
                "mining target requires exactly one of eligible_unique "
                "or max_retained_fraction"
            )
            raise ValueError(msg)
        if self.eligible_unique is not None:
            integer(
                self.eligible_unique,
                field_name="mining_target.eligible_unique",
                minimum=1,
            )
        if self.max_retained_fraction is not None:
            value = finite(
                self.max_retained_fraction, "mining_target.max_retained_fraction"
            )
            if not 0 < value <= 1:
                msg = "mining_target.max_retained_fraction must be in (0, 1]"
                raise ValueError(msg)
            object.__setattr__(self, "max_retained_fraction", float(value))
        integer(
            self.minimum_candidates,
            field_name="mining_target.minimum_candidates",
            minimum=0,
        )

    def resolve(self, retained_count: int) -> dict[str, int]:
        """Expose the stopping criteria independently of candidate outcomes."""
        integer(retained_count, field_name="retain.count", minimum=0)
        return {
            "eligible_unique": self.eligible_unique
            if self.eligible_unique is not None
            else ceil(
                Fraction(retained_count) / Fraction(str(self.max_retained_fraction))
            ),
            "minimum_candidates": self.minimum_candidates,
        }


@dataclass(frozen=True)
class CandidateBudget:
    """Maximum candidate effort, independent of the requested retained count."""

    candidates: int
    seconds: float | None = None
    batch_size: int = 1000
    batch_bases: int = field(default=_DEFAULT_BATCH_BASES, kw_only=True)
    total_bases: int = field(default=_DEFAULT_TOTAL_BASES, kw_only=True)

    def __post_init__(self) -> None:
        """Require finite effort and bounded scoring batches."""
        integer(self.candidates, field_name="budget.candidates", minimum=1)
        integer(self.batch_size, field_name="budget.batch_size", minimum=1)
        integer(self.batch_bases, field_name="budget.batch_bases", minimum=1)
        integer(self.total_bases, field_name="budget.total_bases", minimum=1)
        if self.seconds is not None:
            if finite(self.seconds, "budget.seconds") <= 0:
                msg = "budget.seconds must be positive"
                raise ValueError(msg)
            object.__setattr__(self, "seconds", float(self.seconds))

    def admit(self, maximum_length: int) -> None:
        """Admit declared candidate bases without equating bases with RAM usage."""
        integer(maximum_length, field_name="sampling maximum length", minimum=1)
        batch = min(self.candidates, self.batch_size) * maximum_length
        if batch > self.batch_bases:
            msg = (
                f"candidate batch bound {batch} bases exceeds "
                f"budget.batch_bases={self.batch_bases}; reduce budget.batch_size "
                "or sampling length, or explicitly raise budget.batch_bases"
            )
            raise ValueError(msg)
        total = self.candidates * maximum_length
        if total > self.total_bases:
            msg = (
                f"candidate total bound {total} bases exceeds "
                f"budget.total_bases={self.total_bases}; reduce budget.candidates "
                "or sampling length, or explicitly raise budget.total_bases"
            )
            raise ValueError(msg)

    def to_dict(self) -> dict[str, int | float | None]:
        """Encode explicit cap changes while preserving earlier default requests."""
        return {
            "candidates": self.candidates,
            "seconds": self.seconds,
            "batch_size": self.batch_size,
            **(
                {"batch_bases": self.batch_bases}
                if self.batch_bases != _DEFAULT_BATCH_BASES
                else {}
            ),
            **(
                {"total_bases": self.total_bases}
                if self.total_bases != _DEFAULT_TOTAL_BASES
                else {}
            ),
        }


@dataclass(frozen=True)
class Eligibility:
    """An optional strict score cutoff applied after a qualifying motif hit."""

    best_hit_score_min_exclusive: float | None = None

    def __post_init__(self) -> None:
        """Keep absent cutoffs distinct from a declared zero threshold."""
        if self.best_hit_score_min_exclusive is not None:
            object.__setattr__(
                self,
                "best_hit_score_min_exclusive",
                finite(
                    self.best_hit_score_min_exclusive,
                    "eligibility.best_hit_score_min_exclusive",
                ),
            )


@dataclass(frozen=True)
class Uniqueness:
    """Choose sequence or oriented core equivalence within one preparation recipe."""

    key: str = "sequence"

    def __post_init__(self) -> None:
        """Reject unsupported per-recipe equivalence keys."""
        if self.key not in {"sequence", "core"}:
            msg = "uniqueness.key must be sequence or core"
            raise ValueError(msg)
