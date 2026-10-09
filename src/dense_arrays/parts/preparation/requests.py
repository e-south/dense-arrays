"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/preparation/requests.py

Individual preparation separates candidate effort, eligibility and retention.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field, replace

from dense_arrays._record_validation import integer
from dense_arrays.constraints import GC, Avoid
from dense_arrays.parts.filters import PartFilter
from dense_arrays.parts.models import PartTable
from dense_arrays.parts.motifs.artifacts import PWMArtifact
from dense_arrays.parts.retention import MMR, ScoreBands
from dense_arrays.parts.sampling import (
    Background,
    CandidateBudget,
    Eligibility,
    MiningTarget,
    Sampling,
    Uniqueness,
)
from dense_arrays.parts.scoring.configuration import FimoScoring
from dense_arrays.parts.screening import PWMExclusion


@dataclass(frozen=True)
class Retention:
    """Curated filtering or explicit sampled-pool count and selection policy."""

    select: PartFilter | None = None
    count: int | None = None
    policy: str | None = None
    rank_by: str | None = None
    mmr: MMR | None = None

    def __post_init__(self) -> None:
        """Reject mixed curated filtering and sampled selection contracts."""
        if (self.policy == "mmr") != (self.mmr is not None):
            msg = "retain.mmr is required only for the mmr policy"
            raise ValueError(msg)
        if self.mmr is not None and not isinstance(self.mmr, MMR):
            msg = "retain.mmr must be MMR"
            raise TypeError(msg)
        if self.select is not None and not isinstance(self.select, PartFilter):
            msg = "retention.select must be PartFilter"
            raise TypeError(msg)
        if self.count is None:
            if self.policy is not None or self.rank_by is not None:
                msg = "sampled retention requires count and policy"
                raise ValueError(msg)
            return
        integer(self.count, field_name="retain.count", minimum=0)
        if self.select is not None:
            msg = "sampled retention cannot include a curated part filter"
            raise ValueError(msg)
        if self.policy not in {"first_eligible", "top_score", "mmr"}:
            msg = "supported retention policies: first_eligible, top_score, mmr"
            raise ValueError(msg)
        if self.rank_by != (
            "best_hit_score" if self.policy in {"top_score", "mmr"} else None
        ):
            msg = "retention rank_by must match its selection policy"
            raise ValueError(msg)


@dataclass(frozen=True)
class PreparationSpec:
    """Curated import or bounded candidate preparation, through one operation."""

    source: PartTable | PWMArtifact | Background
    retain: Retention = field(default_factory=Retention)
    sampling: Sampling | None = None
    budget: CandidateBudget | None = None
    scoring: FimoScoring | None = None
    eligibility: Eligibility = field(default_factory=Eligibility)
    uniqueness: Uniqueness = field(default_factory=Uniqueness)
    screening: tuple[Avoid | GC | PWMExclusion, ...] = ()
    seed: int = 0
    mining_target: MiningTarget | None = None
    score_bands: ScoreBands | None = None

    def __post_init__(self) -> None:
        """Validate applicable policies before opening sources or invoking tools."""
        if not isinstance(self.source, (PartTable, PWMArtifact, Background)):
            msg = "preparation requires PartTable, PWMArtifact or Background"
            raise TypeError(msg)
        for name, kind in (
            ("retain", Retention),
            ("eligibility", Eligibility),
            ("uniqueness", Uniqueness),
        ):
            if not isinstance(getattr(self, name), kind):
                msg = f"preparation.{name} must be {kind.__name__}"
                raise TypeError(msg)
        integer(self.seed, field_name="preparation.seed", minimum=0)
        if self.mining_target is not None and not isinstance(
            self.mining_target, MiningTarget
        ):
            msg = "preparation.mining_target must be MiningTarget"
            raise TypeError(msg)
        self._validate_screens()
        self._validate_score_bands()
        self._validate_source_policies()

    def _validate_score_bands(self) -> None:
        if self.score_bands is not None:
            if not isinstance(self.score_bands, ScoreBands):
                msg = "preparation.score_bands must be ScoreBands"
                raise TypeError(msg)
            if not isinstance(self.source, PWMArtifact):
                msg = "score bands require PWM preparation with explicit scoring"
                raise ValueError(msg)

    def _validate_screens(self) -> None:
        if not isinstance(self.screening, (tuple, list)):
            msg = "screening must be an ordered list of sequence requirements"
            raise TypeError(msg)
        object.__setattr__(self, "screening", tuple(self.screening))
        for rule in self.screening:
            if (
                not isinstance(rule, (Avoid, GC, PWMExclusion))
                or (isinstance(rule, GC) and rule.scope != "sequence")
                or (isinstance(rule, Avoid) and rule.except_placements)
            ):
                msg = (
                    "preparation screening supports sequence GC and avoid "
                    "without placement exceptions"
                )
                raise ValueError(msg)
        if {r.id for r in self.screening} & {"no_qualifying_hit", "best_hit_score"}:
            msg = "screening ID is reserved for built-in score eligibility"
            raise ValueError(msg)
        if len({r.id for r in self.screening}) != len(self.screening):
            msg = "preparation screening IDs must be unique"
            raise ValueError(msg)

    def _validate_source_policies(self) -> None:
        if isinstance(self.source, PartTable):
            if (
                any(
                    v is not None
                    for v in (
                        self.sampling,
                        self.budget,
                        self.scoring,
                        self.retain.count,
                        self.mining_target,
                    )
                )
                or self.screening
                or self.seed != 0
                or self.eligibility != Eligibility()
                or self.uniqueness != Uniqueness()
            ):
                msg = "curated preparation accepts only its part filter"
                raise ValueError(msg)
            return
        if (
            not isinstance(self.sampling, Sampling)
            or not isinstance(self.budget, CandidateBudget)
            or self.retain.count is None
        ):
            msg = "sampled preparation requires sampling, budget and explicit retention"
            raise ValueError(msg)
        if self.sampling.strategy == "conditional" and not isinstance(
            self.source, Background
        ):
            msg = "conditional sampling requires a Background source"
            raise ValueError(msg)
        if isinstance(self.source, Background):
            if (
                self.sampling.strategy not in {"stochastic", "conditional"}
                or self.sampling.base_probabilities is not None
            ):
                msg = (
                    "Background owns its base distribution and uses "
                    "stochastic or conditional sampling"
                )
                raise ValueError(msg)
            if (
                self.scoring is not None
                or self.eligibility != Eligibility()
                or self.uniqueness.key != "sequence"
                or self.retain.policy != "first_eligible"
            ):
                msg = (
                    "background preparation requires sequence uniqueness and "
                    "first_eligible retention without a score cutoff"
                )
                raise ValueError(msg)
        elif not isinstance(self.scoring, FimoScoring):
            msg = "PWM preparation requires explicit FimoScoring"
            raise TypeError(msg)

    def with_changes(self, **changes: object) -> PreparationSpec:
        """Validate an immutable edit without opening inputs or preparing parts."""
        return replace(self, **changes)
