"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/preparation.py

Compose bounded candidate generation, scoring, screening and native publication.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays.artifacts.preparation.records import recount
from dense_arrays.artifacts.preparation.sets import (
    SetAccounting,
    qualify_candidate,
)
from dense_arrays.artifacts.preparation.storage import pool_destination, publish
from dense_arrays.parts import Part
from dense_arrays.parts.background import compile_background
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.eligibility import rejection_reasons
from dense_arrays.parts.mining import Proposal, propose_sequence, sample_sequence
from dense_arrays.parts.retention.collisions import validate_collisions
from dense_arrays.parts.retention.pool import PoolSize
from dense_arrays.parts.retention.selection import eligible_key, select_candidates
from dense_arrays.parts.scoring import ScoringError, scan_fimo
from dense_arrays.parts.screening import PWMExclusion, ScreenObservation
from dense_arrays.planning.preparation.sampled import MotifSource
from dense_arrays.planning.preparation.sets import SetPreparation

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.preparation.records import PoolAccounting
    from dense_arrays.parts.background import CountedBackground
    from dense_arrays.parts.pools import PoolHandle
    from dense_arrays.parts.scoring import FimoBinding, FimoHit, FimoResult
    from dense_arrays.planning.preparation import PreparationPlan
    from dense_arrays.planning.preparation.sampled import SampledPreparation


def _batch(
    plan: SampledPreparation,
    start: int,
    stop: int,
    deadline: float | None,
    sampler: CountedBackground | None = None,
) -> tuple[Candidate, ...]:
    request = plan.request
    motif = plan.source.motif if isinstance(plan.source, MotifSource) else None
    proposals = tuple(
        Proposal(sampler.draw(seed=request.seed, index=i), None, None)
        if sampler is not None
        else propose_sequence(
            length=request.sampling.candidate_length(seed=request.seed, index=i),
            seed=request.seed,
            index=i,
            probabilities=plan.base_probabilities,
            motif=motif,
            strategy=request.sampling.strategy,
        )
        if request.sampling.records_proposal
        else Proposal(
            sample_sequence(
                length=request.sampling.candidate_length(seed=request.seed, index=i),
                seed=request.seed,
                index=i,
                probabilities=plan.base_probabilities,
                motif=motif,
            ),
            None,
            None,
        )
        for i in range(start, stop)
    )
    sequences = tuple(p.sequence for p in proposals)
    scores = (None,) * len(sequences)
    error = None
    observations = [[] for _ in sequences]
    if (
        deadline is not None
        and (isinstance(plan.source, MotifSource) or plan.screens)
        and time.monotonic() >= deadline
    ):
        return ()
    try:
        if isinstance(plan.source, MotifSource):
            scores = _scan(plan.source.scoring, sequences, deadline).hits
        rules = {
            rule.id: rule
            for rule in request.screening
            if isinstance(rule, PWMExclusion)
        }
        for screen in plan.screens:
            rule = rules[screen.rule_id]
            for bound in screen.motifs:
                hits = _screen_hits(bound.scoring, sequences, deadline)
                for row, hit in zip(observations, hits, strict=True):
                    row.append(
                        ScreenObservation(screen.rule_id, bound.scoring.binding_id, hit)
                    )
                _check_exclusion_scale(rule, hits)
    except ScoringError as err:
        error = str(err)
    candidates = []
    group = motif.motif_id if motif is not None else plan.source.group
    for index, proposal, score, observed in zip(
        range(start, stop), proposals, scores, observations, strict=True
    ):
        part = Part(
            f"candidate_{index}",
            proposal.sequence,
            group=group,
            source="pwm_artifact" if motif else "background",
            core_start=score.start if score else None,
            core_end=score.end if score else None,
            core_orientation=score.strand if score else None,
            metadata={
                **({"score": score.to_dict()} if score else {}),
                **(
                    {
                        "proposal": {
                            **proposal.evidence(request.sampling.strategy),
                            **(
                                {"model_id": sampler.model_id}
                                if sampler is not None
                                else {}
                            ),
                        }
                    }
                    if request.sampling.records_proposal
                    else {}
                ),
            },
        )
        candidate = Candidate(index, part, error=error, screening=tuple(observed))
        candidates.append(
            candidate
            if error
            else replace(
                candidate,
                reasons=rejection_reasons(
                    candidate,
                    request.eligibility,
                    request.screening,
                    requires_hit=motif is not None,
                ),
            )
        )
    return tuple(candidates)


def _screen_hits(
    binding: FimoBinding, sequences: tuple[str, ...], deadline: float | None
) -> tuple[FimoHit | None, ...]:
    """Sequences shorter than an exclusion motif have no full-width windows."""
    positions = tuple(
        i for i, seq in enumerate(sequences) if len(seq) >= binding.motif.width
    )
    hits = [None] * len(sequences)
    if positions:
        scored = _scan(binding, tuple(sequences[i] for i in positions), deadline).hits
        for position, hit in zip(positions, scored, strict=True):
            hits[position] = hit
    return tuple(hits)


def _check_exclusion_scale(
    rule: PWMExclusion, hits: tuple[FimoHit | None, ...]
) -> None:
    """Keep undefined score ratios distinct from candidate rejection."""
    if rule.score_field == "fraction_of_max" and any(
        hit is not None and hit.fraction_of_max is None for hit in hits
    ):
        reason = "undefined_score"
        raise ScoringError(
            reason, "exclusion ratio requires a positive theoretical maximum"
        )


def _scan(
    binding: FimoBinding, sequences: tuple[str, ...], deadline: float | None
) -> FimoResult:
    """Give every source/screen call the remaining preparation time budget."""
    if deadline is not None:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            reason = "timeout"
            raise ScoringError(reason, "preparation time budget exhausted")
        binding = replace(
            binding,
            settings=replace(
                binding.settings,
                limits=replace(
                    binding.settings.limits,
                    seconds=min(binding.settings.limits.seconds, remaining),
                ),
            ),
        )
    return scan_fimo(binding, sequences)


def execute_preparation(plan: PreparationPlan, out: Path) -> PoolHandle:
    """Own the output before mining; publish all evidence in one transaction."""
    sources = (
        plan.resolved.recipes.values()
        if isinstance(plan.resolved, SetPreparation)
        else (plan.resolved,)
    )
    for source in sources:
        source.request.budget.admit(source.request.sampling.maximum_length)
    plan.verify_inputs()
    with pool_destination(out) as connection:
        if isinstance(plan.resolved, SetPreparation):
            candidates, records = [], {}
            for name, source in plan.resolved.recipes.items():
                decided, accounting = _mine(source)
                offset = len(candidates)
                candidates.extend(qualify_candidate(c, name, offset) for c in decided)
                records[name] = accounting
            candidates, accounting = tuple(candidates), SetAccounting(records)
            validate_collisions(
                candidates,
                sequences=plan.resolved.sequence_collisions,
                cores=plan.resolved.core_collisions,
            )
        else:
            candidates, accounting = _mine(plan.resolved)
        plan.verify_inputs()
        return publish(connection, out, plan, candidates, accounting)


def _mine(source: SampledPreparation) -> tuple[tuple[Candidate, ...], PoolAccounting]:
    """Execute one independent budget and retention policy without publication."""
    request = source.request
    deadline = (
        None
        if request.budget.seconds is None
        else time.monotonic() + request.budget.seconds
    )
    candidates = []
    reason = "candidate_budget"
    target = source.mining_target
    unique = set()
    sampler, construction = None, None
    for start in range(1, request.budget.candidates + 1, request.budget.batch_size):
        if (
            target is not None
            and target["eligible_unique"] == 0
            and target["minimum_candidates"] == 0
        ):
            reason = "mining_target"
            break
        if deadline is not None and time.monotonic() >= deadline:
            reason = "time_budget"
            break
        if request.sampling.strategy == "conditional" and construction is None:
            compiled = compile_background(
                minimum=request.sampling.minimum_length,
                maximum=request.sampling.maximum_length,
                probabilities=source.base_probabilities,
                screening=source.compiled_rules,
                limits=request.sampling.limits,
                deadline=deadline,
            )
            sampler, construction = compiled.sampler, compiled.report
            if sampler is None:
                reason = f"construction_{construction.status}"
                break
        stop = min(start + request.budget.batch_size, request.budget.candidates + 1)
        batch = _batch(source, start, stop, deadline, sampler)
        if not batch:
            reason = "time_budget"
            break
        candidates.extend(batch)
        if any(c.error for c in batch):
            reason = "execution_error"
            break
        if target is not None:
            unique.update(
                key
                for c in batch
                if (key := eligible_key(c, request.uniqueness)) is not None
            )
            if (
                len(unique) >= target["eligible_unique"]
                and len(candidates) >= target["minimum_candidates"]
            ):
                reason = "mining_target"
                break
    decided = select_candidates(
        tuple(candidates),
        request.uniqueness,
        request.retain,
        motif=source.source.motif if isinstance(source.source, MotifSource) else None,
        score_bands=request.score_bands,
    )
    accounting = recount(
        decided,
        target=request.retain.count,
        budget=request.budget.candidates,
        stop_reason=reason,
        mmr=request.retain.mmr is not None,
        pool_sizing=request.retain.mmr.pool_size
        if request.retain.mmr is not None
        and isinstance(request.retain.mmr.pool_size, PoolSize)
        else None,
        mining_target=target,
        score_bands=request.score_bands,
        scoring_id=source.source.scoring.binding_id
        if isinstance(source.source, MotifSource)
        else None,
        construction=construction,
    )
    return decided, accounting
