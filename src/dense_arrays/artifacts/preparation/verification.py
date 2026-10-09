"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/verification.py

Verify sampled pools from persisted evidence without sampling or scoring.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, object_fields
from dense_arrays.artifacts.preparation.records import PoolAccounting, recount
from dense_arrays.artifacts.preparation.sets import (
    SetAccounting,
    local_candidate,
)
from dense_arrays.artifacts.preparation.storage import read_candidates, result_identity
from dense_arrays.artifacts.reading import ReadBudget, Verification
from dense_arrays.parts.eligibility import rejection_reasons
from dense_arrays.parts.mining import consensus
from dense_arrays.parts.retention.collisions import validate_collisions
from dense_arrays.parts.retention.pool import PoolSize
from dense_arrays.parts.retention.selection import eligible_key, select_candidates
from dense_arrays.parts.screening.sequence import passes_sequence
from dense_arrays.planning.preparation.sampled import MotifSource
from dense_arrays.planning.preparation.screening import verify_observations
from dense_arrays.planning.preparation.sets import SetPreparation

if TYPE_CHECKING:
    import sqlite3
    from pathlib import Path

    from dense_arrays.artifacts.pool_records import PoolSummary
    from dense_arrays.parts.candidates import Candidate
    from dense_arrays.planning.preparation import PreparationPlan
    from dense_arrays.planning.preparation.sampled import SampledPreparation


def verify_sampled(
    connection: sqlite3.Connection,
    path: Path,
    plan: PreparationPlan,
    summary: PoolSummary,
    budget: ReadBudget,
) -> tuple[Verification, tuple[Candidate, ...]]:
    """Recount stage decisions, score geometry and retained-record joins."""
    from dense_arrays.artifacts.pools import (  # noqa: PLC0415 - resolve owner after module initialization
        iter_parts,
    )

    if summary.preparation is None or not plan.sampled:
        msg = "sampled pool requires a matching plan and stage accounting"
        raise ValueError(msg)
    budget.retain(plan.identity_count)
    candidates = read_candidates(connection, budget)
    if isinstance(plan.resolved, SetPreparation):
        if not isinstance(summary.preparation, SetAccounting) or tuple(
            summary.preparation.recipes
        ) != tuple(plan.resolved.recipes):
            msg = "preparation set requires matching per-recipe accounting"
            raise ValueError(msg)
        offset = 0
        for name, source in plan.resolved.recipes.items():
            recorded = summary.preparation.recipes[name]
            stop = offset + recorded.counts["processed"]
            local = tuple(
                local_candidate(c, name, offset) for c in candidates[offset:stop]
            )
            verify_decisions(local, source, recorded, budget=budget)
            offset = stop
        if offset != len(candidates):
            msg = "recipe candidate counts do not cover the pool"
            raise ValueError(msg)
        validate_collisions(
            candidates,
            sequences=plan.resolved.sequence_collisions,
            cores=plan.resolved.core_collisions,
        )
    else:
        if not isinstance(summary.preparation, PoolAccounting):
            msg = "single-recipe pool requires single-recipe accounting"
            raise TypeError(msg)
        verify_decisions(candidates, plan.resolved, summary.preparation, budget=budget)
    if (
        result_identity(plan.plan_id, candidates, summary.preparation)
        != summary.pool_id
        or summary.plan_id != plan.plan_id
    ):
        msg = "sampled pool identity disagrees with its preparation"
        raise ValueError(msg)
    order = (
        {name: index for index, name in enumerate(summary.preparation.recipes)}
        if isinstance(summary.preparation, SetAccounting)
        else {}
    )
    retained = sorted(
        (c for c in candidates if c.retained),
        key=lambda c: (order.get(c.recipe_id, 0), c.rank),
    )
    records = iter_parts(path, pool_id=summary.pool_id, budget=budget)
    try:
        for ordinal, (candidate, record) in enumerate(
            zip(retained, records, strict=True), 1
        ):
            if record.ordinal != ordinal or record.part != candidate.part:
                msg = "retained part does not match its sampled candidate"
                raise ValueError(msg)
    finally:
        records.close()
    return (
        Verification(
            "sampled_preparation_and_retained_parts",
            budget.examined,
            budget.bytes_checked,
        ),
        candidates,
    )


def verify_decisions(
    candidates: tuple[Candidate, ...],
    source: SampledPreparation,
    recorded: PoolAccounting,
    *,
    budget: ReadBudget | None = None,
) -> None:
    """Check one recipe's local evidence without invoking sampling or scoring."""
    request = source.request
    model = source.source
    motif = model.motif if isinstance(model, MotifSource) else None
    group = motif.motif_id if motif else model.group
    _verify_construction(source, recorded)
    for candidate in candidates:
        _verify_proposal(candidate, source)
        if candidate.error and motif is None and not source.screens:
            msg = "candidate scoring error requires a declared scoring operation"
            raise ValueError(msg)
        verify_observations(candidate, source.screens)
        expected_reasons = (
            ()
            if candidate.error
            else rejection_reasons(
                candidate,
                request.eligibility,
                request.screening,
                requires_hit=motif is not None,
            )
        )
        if (
            candidate.reasons != expected_reasons
            or not request.sampling.minimum_length
            <= len(candidate.part.sequence)
            <= request.sampling.maximum_length
            or candidate.part.group != group
            or candidate.part.part_id != f"candidate_{candidate.index}"
            or candidate.recipe_id is not None
            or candidate.recipe_index is not None
        ):
            msg = "sampled candidate meaning disagrees with its preparation recipe"
            raise ValueError(msg)
        if isinstance(model, MotifSource):
            model.scoring.verify_hit(candidate.score)
        elif candidate.score is not None or candidate.part.core_start is not None:
            msg = "background candidates cannot carry a primary score or motif core"
            raise ValueError(msg)
    _admit_mmr(candidates, source, budget or ReadBudget())
    reset = tuple(
        replace(
            c,
            representative=None,
            rank=None,
            retained=False,
            selection=None,
            score_band=None,
        )
        for c in candidates
    )
    selected = select_candidates(
        reset,
        request.uniqueness,
        request.retain,
        motif=motif,
        score_bands=request.score_bands,
    )
    _verify_mining_stop(reset, source, recorded.stop_reason)
    accounting = recount(
        selected,
        target=request.retain.count,
        budget=request.budget.candidates,
        stop_reason=recorded.stop_reason,
        mmr=request.retain.mmr is not None,
        pool_sizing=request.retain.mmr.pool_size
        if request.retain.mmr is not None
        and isinstance(request.retain.mmr.pool_size, PoolSize)
        else None,
        mining_target=source.mining_target,
        score_bands=request.score_bands,
        scoring_id=model.scoring.binding_id if isinstance(model, MotifSource) else None,
        construction=recorded.construction,
    )
    if selected != candidates or accounting != recorded:
        msg = "sampled selection or stage accounting disagrees"
        raise ValueError(msg)


def _admit_mmr(
    candidates: tuple[Candidate, ...], source: SampledPreparation, budget: ReadBudget
) -> None:
    """Bound selection replay from supplied candidates, not claimed admission counts."""
    request = source.request
    policy = request.retain.mmr
    if policy is None:
        return
    keys = set()
    for candidate in candidates:
        key = eligible_key(candidate, request.uniqueness)
        if key is not None and key not in keys:
            budget.retain()
            keys.add(key)
    # Score cutoffs may shrink this bound. Each retained rank updates the whole pool.
    pool = min(len(keys), policy.pool_limit(request.retain.count))
    budget.compare(pool * min(request.retain.count, pool))


def _verify_construction(source: SampledPreparation, recorded: PoolAccounting) -> None:
    """Check recorded counting scope and allowances without rerunning generation."""
    conditional = source.request.sampling.strategy == "conditional"
    report = recorded.construction
    if not conditional:
        if report is not None:
            msg = "construction evidence is not declared by this preparation"
            raise ValueError(msg)
        return
    if report is None:
        if recorded.counts["processed"] or recorded.stop_reason not in {
            "mining_target",
            "time_budget",
        }:
            msg = "conditional candidates require construction evidence"
            raise ValueError(msg)
        return
    limits = source.request.sampling.limits
    if report.model_id != source.conditional_model_id or any(
        getattr(report, key) > getattr(limits, key)
        for key in ("states", "automaton_states", "mass_bits")
    ):
        msg = "conditional construction disagrees with its model or work limits"
        raise ValueError(msg)


def _verify_mining_stop(
    candidates: tuple[Candidate, ...], source: SampledPreparation, stop_reason: str
) -> None:
    """Require the first attained complete batch, using saved eligible equivalence."""
    target = source.mining_target
    if target is None:
        return
    if target["eligible_unique"] == 0 and target["minimum_candidates"] == 0:
        if candidates or stop_reason != "mining_target":
            msg = "zero mining target requires stopping before the first batch"
            raise ValueError(msg)
        return
    unique = set()
    request = source.request
    for start in range(0, len(candidates), request.budget.batch_size):
        expected_end = min(start + request.budget.batch_size, request.budget.candidates)
        batch = candidates[start:expected_end]
        if any(c.error for c in batch):
            return
        unique.update(
            key
            for c in batch
            if (key := eligible_key(c, request.uniqueness)) is not None
        )
        processed = start + len(batch)
        if (
            len(unique) >= target["eligible_unique"]
            and processed >= target["minimum_candidates"]
        ):
            if (
                processed != expected_end
                or processed != len(candidates)
                or stop_reason != "mining_target"
            ):
                msg = "mining target must stop at the first attained complete batch"
                raise ValueError(msg)
            return


def _verify_proposal(candidate: Candidate, plan: SampledPreparation) -> None:
    """Check declared construction geometry and base support without random replay."""
    sampling = plan.request.sampling
    evidence = candidate.part.metadata.get("proposal")
    if not sampling.records_proposal:
        if evidence is not None:
            msg = "proposal evidence is not declared by this preparation plan"
            raise ValueError(msg)
        if not isinstance(plan.source, MotifSource):
            _verify_support(candidate.part.sequence, plan, start=None, end=None)
        return
    fields = {"schema", "strategy", "policy", "start", "end"}
    if sampling.strategy == "conditional":
        fields.add("model_id")
    data = object_fields(evidence, fields, "proposal evidence")
    if (
        set(data) != fields
        or data["schema"] != "dense_arrays.proposal.v1"
        or data["strategy"] != sampling.strategy
        or data["policy"] != sampling.policy
    ):
        msg = "proposal evidence disagrees with the declared sampling policy"
        raise ValueError(msg)
    _verify_geometry_and_support(data, candidate.part.sequence, plan)


def _verify_geometry_and_support(
    data: dict[str, object], sequence: str, plan: SampledPreparation
) -> None:
    """Validate construction intervals, conditioned rules and positive base support."""
    sampling = plan.request.sampling
    motif = plan.source.motif if isinstance(plan.source, MotifSource) else None
    start, end = data["start"], data["end"]
    if sampling.strategy in {"background", "conditional"}:
        if start is not None or end is not None:
            msg = "background proposal cannot declare an intended motif interval"
            raise ValueError(msg)
        if sampling.strategy == "conditional":
            _verify_conditional_sequence(data["model_id"], sequence, plan)
    else:
        integer(start, field_name="proposal.start", minimum=0)
        integer(end, field_name="proposal.end", minimum=1)
        if end != start + motif.width or end > len(sequence):
            msg = "proposal interval does not fit the motif or candidate sequence"
            raise ValueError(msg)
        if sampling.strategy == "consensus" and sequence[start:end] != consensus(motif):
            msg = "proposal motif interval does not contain the declared consensus"
            raise ValueError(msg)
    _verify_support(sequence, plan, start=start, end=end)


def _verify_support(
    sequence: str, plan: SampledPreparation, *, start: int | None, end: int | None
) -> None:
    """Reject bases with zero probability under the declared proposal distribution."""
    motif = plan.source.motif if isinstance(plan.source, MotifSource) else None
    for index, base in enumerate(sequence):
        weights = (
            motif.probabilities[index - start]
            if start is not None and start <= index < end
            else plan.base_probabilities
        )
        if weights["ACGT".index(base)] == 0:
            msg = "proposal contains a base outside its declared sampling support"
            raise ValueError(msg)


def _verify_conditional_sequence(
    model_id: str, sequence: str, plan: SampledPreparation
) -> None:
    """Independently apply native constraints even when later scoring failed."""
    if model_id != plan.conditional_model_id or any(
        not passes_sequence(rule, sequence) for rule in plan.compiled_rules
    ):
        msg = "conditional proposal violates its compiled constraints or model"
        raise ValueError(msg)
