"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/preparation/sampled.py

Resolved sampled preparation binds source evidence without generating candidates.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from functools import cached_property
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import mutable_json, object_fields, semantic_digest
from dense_arrays.constraints import GC, Avoid
from dense_arrays.parts import Background, LengthRange, PreparationSpec, PWMArtifact
from dense_arrays.parts.background import model_identity
from dense_arrays.parts.mining import LENGTH_POLICY
from dense_arrays.parts.motifs.windows import WINDOW_POLICY
from dense_arrays.parts.retention.bands import BAND_POLICY
from dense_arrays.parts.retention.pool import POOL_SIZE_POLICY, PoolSize

from .motifs import MotifSource, resolve_motif
from .requests import preparation_from_dict, preparation_to_dict
from .screening import (
    BoundExclusion,
    resolve_screens,
    screening_preview,
    validate_screens,
)

if TYPE_CHECKING:
    from pathlib import Path

SAMPLED_PLAN_SCHEMA = "dense_arrays.preparation_plan.v2"
SAMPLED_POLICIES = MappingProxyType(
    {
        "sampling": "candidate_shake256.v1",
        "representative": "best_score_then_first.v1",
        "retention": "explicit_rank.v1",
        "evidence": "full_candidates.v1",
    }
)


@dataclass(frozen=True, repr=False)
class SampledPreparation:
    """Complete sampled request and resolved sources with unknown retained yield."""

    request: PreparationSpec
    source: Background | MotifSource
    screens: tuple[BoundExclusion, ...] = ()

    def __post_init__(self) -> None:
        """Reject contradictory source families or request/scorer bindings."""
        if not isinstance(self.request, PreparationSpec):
            msg = "sampled plan requires PreparationSpec"
            raise TypeError(msg)
        if isinstance(self.source, Background):
            if self.request.source != self.source:
                msg = "sampled background source disagrees with request"
                raise ValueError(msg)
        elif isinstance(self.source, MotifSource):
            request = self.request
            if (
                not isinstance(request.source, PWMArtifact)
                or request.source.path != self.source.input.path
                or request.scoring != self.source.scoring.settings
                or request.source.window != self.source.window
                or request.source.format != self.source.input.format
            ):
                msg = "sampled PWM source disagrees with request"
                raise ValueError(msg)
            if request.sampling.minimum_length < self.source.motif.width:
                msg = (
                    f"minimum candidate length {request.sampling.minimum_length} "
                    f"is shorter than selected motif width {self.source.motif.width}; "
                    "increase sampling.length or explicitly select a shorter motif "
                    "with source.window (PWMArtifact.window in Python)"
                )
                raise ValueError(msg)
            batch = min(request.budget.candidates, request.budget.batch_size)
            windows = self.source.window_bound(request.sampling.maximum_length, batch)
            if windows > request.scoring.limits.windows:
                msg = (
                    "candidate batch exceeds the scoring window limit, "
                    "including calibration"
                )
                raise ValueError(msg)
            if request.source.motif_ids and request.source.motif_ids != (
                self.source.input.motif.motif_id,
            ):
                msg = "sampled motif selector disagrees with resolved input"
                raise ValueError(msg)
        else:
            msg = "unsupported resolved preparation source"
            raise TypeError(msg)
        if not isinstance(self.screens, (tuple, list)) or any(
            not isinstance(screen, BoundExclusion) for screen in self.screens
        ):
            msg = "sampled screens require ordered BoundExclusion records"
            raise TypeError(msg)
        object.__setattr__(self, "screens", tuple(self.screens))
        validate_screens(self.request, self.screens)

    @property
    def policies(self) -> dict[str, str]:
        """Bind the chosen retention algorithm without changing earlier plans."""
        return (
            dict(SAMPLED_POLICIES)
            | {"sampling": self.request.sampling.policy}
            | (
                {"length": self.length_policy}
                if isinstance(self.request.sampling.length, LengthRange)
                else {}
            )
            | (
                {"motif_window": WINDOW_POLICY}
                if (isinstance(self.source, MotifSource) and self.source.window)
                or any(m.window for screen in self.screens for m in screen.motifs)
                else {}
            )
            | ({"screening": "pwm_exclusion.v1"} if self.screens else {})
            | (
                {"pool_size": POOL_SIZE_POLICY}
                if self.request.retain.mmr is not None
                and isinstance(self.request.retain.mmr.pool_size, PoolSize)
                else {}
            )
            | (
                {"score_bands": BAND_POLICY}
                if self.request.score_bands is not None
                else {}
            )
            | (
                {"mining": "eligible_unique_batch_stop.v1"}
                if self.request.mining_target is not None
                else {}
            )
            | (
                {"retention": self.request.retain.mmr.algorithm}
                if self.request.retain.mmr is not None
                else {}
            )
        )

    @property
    def base_probabilities(self) -> tuple[float, ...]:
        """Resolve sampling probabilities independently of scoring background."""
        if isinstance(self.source, Background):
            return self.source.base_probabilities
        return (
            self.request.sampling.base_probabilities
            or self.source.input.motif.background
        )

    @property
    def length_policy(self) -> str:
        """Conditional length is part of the joint sequence draw."""
        return (
            self.request.sampling.policy
            if self.request.sampling.strategy == "conditional"
            else LENGTH_POLICY
        )

    @property
    def compiled_rules(self) -> tuple[Avoid | GC, ...]:
        """Share native sequence constraints; external scorers stay post-proposal."""
        return tuple(r for r in self.request.screening if isinstance(r, (Avoid, GC)))

    @cached_property
    def conditional_model_id(self) -> str | None:
        """Bind a conditional model without building its counting tables."""
        if self.request.sampling.strategy != "conditional":
            return None
        return model_identity(
            minimum=self.request.sampling.minimum_length,
            maximum=self.request.sampling.maximum_length,
            probabilities=self.base_probabilities,
            screening=self.compiled_rules,
        )

    @property
    def proposal_preview(self) -> dict[str, object]:
        """Disclose explicit proposal geometry and the source of background draws."""
        sampling = self.request.sampling
        if not sampling.records_proposal:
            return {}
        return {
            "proposal": {
                "strategy": sampling.strategy,
                "policy": sampling.policy,
                "base_probabilities": self.base_probabilities,
                "background_source": "sampling"
                if sampling.base_probabilities is not None
                else "background"
                if isinstance(self.source, Background)
                else "motif_artifact",
                "motif_placement": "none"
                if sampling.strategy in {"background", "conditional"}
                else "embedded_forward",
                **(
                    {
                        "model_id": self.conditional_model_id,
                        "constraint_ids": tuple(r.id for r in self.compiled_rules),
                        "construction_limits": asdict(sampling.limits),
                        "distribution": (
                            "declared_background_conditioned_on_constraints"
                        ),
                        "valid_mass": None,
                    }
                    if sampling.strategy == "conditional"
                    else {}
                ),
            }
        }

    @property
    def identity_count(self) -> int:
        """Count embedded model positions and declared sequence screens."""
        return (
            sum(
                m.input.motif.width + 1 + (m.motif.width + 1 if m.window else 0)
                for screen in self.screens
                for m in screen.motifs
            )
            + len(self.request.screening)
            + (
                len(self.request.score_bands.upper_fractions)
                if self.request.score_bands is not None
                else 0
            )
            + (
                self.source.input.motif.width
                + 1
                + (self.source.motif.width + 1 if self.source.window else 0)
                if isinstance(self.source, MotifSource)
                else 1
            )
        )

    @property
    def plan_id(self) -> str:
        """Bind source/scoring semantics independently of local file locations."""
        request = preparation_to_dict(self.request)
        if isinstance(self.source, MotifSource):
            request["source"].pop("path")
            request["scoring"].pop("executable")
            request["scoring"].pop("background")
        for rule in request["screening"]:
            if rule.get("kind") == "pwm_exclusion":
                for motif in rule["motifs"]:
                    motif.pop("path")
                rule["scoring"].pop("executable")
                rule["scoring"].pop("background")
        source = (
            {
                "input_sha256": self.source.input.sha256,
                "scoring_id": self.source.scoring.binding_id,
            }
            if isinstance(self.source, MotifSource)
            else {"kind": "background"}
        )
        return semantic_digest(
            {
                "schema": SAMPLED_PLAN_SCHEMA,
                "request": request,
                "source": source,
                **(
                    {
                        "screens": [
                            {
                                "rule_id": screen.rule_id,
                                "motifs": [
                                    {
                                        "input_sha256": m.input.sha256,
                                        "scoring_id": m.scoring.binding_id,
                                    }
                                    for m in screen.motifs
                                ],
                            }
                            for screen in self.screens
                        ]
                    }
                    if self.screens
                    else {}
                ),
                "policies": self.policies,
            }
        )

    @property
    def mining_target(self) -> dict[str, int] | None:
        """Resolve desired eligible supply separately from the candidate effort cap."""
        target = self.request.mining_target
        return None if target is None else target.resolve(self.request.retain.count)

    @property
    def preview(self) -> MappingProxyType:
        """Report effort and requested retention without predicting eligibility."""
        return MappingProxyType(
            {
                "source_parts": 0,
                "source_motifs": int(isinstance(self.source, MotifSource)),
                "retained_parts": None,
                "requested_retention": self.request.retain.count,
                "retained_count_status": "unknown",
                "candidate_budget": self.request.budget.candidates,
                **(
                    {
                        "score_bands": {
                            **self.request.score_bands.to_dict(),
                            "policy": BAND_POLICY,
                            "population": "eligible_unique",
                            "scoring_id": self.source.scoring.binding_id,
                            "counts": None,
                        }
                    }
                    if self.request.score_bands is not None
                    else {}
                ),
                **(
                    {"mining_target": self.mining_target}
                    if self.mining_target is not None
                    else {}
                ),
                "candidate_bases_bound": self.request.budget.candidates
                * self.request.sampling.maximum_length,
                **(
                    {
                        "retention": {
                            **self.request.retain.mmr.to_dict(),
                            **(
                                {
                                    "pool_limit": self.request.retain.mmr.pool_limit(
                                        self.request.retain.count
                                    )
                                }
                                if isinstance(
                                    self.request.retain.mmr.pool_size, PoolSize
                                )
                                else {}
                            ),
                            "distance_work_bound": min(
                                self.request.budget.candidates,
                                self.request.retain.mmr.pool_limit(
                                    self.request.retain.count
                                ),
                            )
                            * min(
                                self.request.retain.count,
                                self.request.retain.mmr.pool_limit(
                                    self.request.retain.count
                                ),
                                self.request.budget.candidates,
                            )
                            * self.source.motif.width,
                        }
                    }
                    if self.request.retain.mmr is not None
                    else {}
                ),
                **(
                    {
                        "sampled_length": {
                            "minimum": self.request.sampling.minimum_length,
                            "maximum": self.request.sampling.maximum_length,
                            "distribution": "uniform_prior_conditioned_on_constraints"
                            if self.request.sampling.strategy == "conditional"
                            else "uniform_integer",
                            "policy": self.length_policy,
                        }
                    }
                    if isinstance(self.request.sampling.length, LengthRange)
                    else {}
                ),
                **(
                    {"motif_window": self.source.selection.to_dict()}
                    if isinstance(self.source, MotifSource) and self.source.selection
                    else {}
                ),
                **self.proposal_preview,
                **(
                    {
                        "motif_import": {
                            "format": self.source.input.format,
                            "score_matrix": "unavailable"
                            if self.source.motif.log_odds is None
                            else "supplied",
                        }
                    }
                    if isinstance(self.source, MotifSource)
                    and self.source.input.format != "artifact"
                    else {}
                ),
                **screening_preview(self.request, self.screens),
                "required_tools": ("fimo",)
                if isinstance(self.source, MotifSource) or self.screens
                else (),
                "screening_stages": (
                    "sampling",
                    "scoring",
                    "eligibility",
                    "uniqueness",
                    "retention",
                )
                if isinstance(self.source, MotifSource) or self.screens
                else ("sampling", "eligibility", "uniqueness", "retention"),
            }
        )

    def verify_inputs(self) -> None:
        """Recheck exact source and executable bytes before mining."""
        if isinstance(self.source, MotifSource):
            self.source.verify()
        for screen in self.screens:
            for motif in screen.motifs:
                motif.verify()

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Save all effective policies and bound source evidence."""
        source = {"kind": "background"}
        if isinstance(self.source, MotifSource):
            source = self.source.to_dict(base=base)
        return {
            "schema": SAMPLED_PLAN_SCHEMA,
            "plan_id": self.plan_id,
            "request": preparation_to_dict(self.request, base=base),
            "source": source,
            **(
                {"screens": [screen.to_dict(base=base) for screen in self.screens]}
                if self.screens
                else {}
            ),
            "policies": self.policies,
            "preview": mutable_json(dict(self.preview)),
        }

    @classmethod
    def from_dict(
        cls, value: object, *, base: Path | None = None
    ) -> SampledPreparation:
        """Restore complete evidence without reading files or querying a tool."""
        keys = {"schema", "plan_id", "request", "source", "policies", "preview"}
        data = object_fields(value, keys | {"screens"}, "sampled preparation plan")
        encoded_screens = data.pop("screens", None)
        if encoded_screens is not None and not isinstance(encoded_screens, list):
            msg = "sampled screens must be an array"
            raise TypeError(msg)
        if set(data) != keys or data["schema"] != SAMPLED_PLAN_SCHEMA:
            msg = "unsupported or incomplete sampled preparation plan"
            raise ValueError(msg)
        request = preparation_from_dict(data["request"], base=base)
        if isinstance(request.source, Background):
            source = request.source
        else:
            source = MotifSource.from_dict(data["source"], base=base)
        screens = tuple(
            BoundExclusion.from_dict(x, base=base) for x in (encoded_screens or [])
        )
        result = cls(request, source, screens)
        if encoded_screens is not None:
            data["screens"] = encoded_screens
        if result.to_dict(base=base) != data:
            msg = "sampled preparation digest or resolved fields disagree"
            raise ValueError(msg)
        return result


def resolve_sampled(request: PreparationSpec) -> SampledPreparation:
    """Validate an input model and scoring availability without sampling."""
    request.budget.admit(request.sampling.maximum_length)
    request, screens = resolve_screens(request)
    if isinstance(request.source, Background):
        return SampledPreparation(request, request.source, screens)
    source, bound = resolve_motif(request.source, request.scoring)
    resolved = replace(request, source=source, scoring=bound.scoring.settings)
    return SampledPreparation(resolved, bound, screens)
