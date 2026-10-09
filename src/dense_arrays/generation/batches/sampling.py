"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/batches/sampling.py

Finite identity-based sampling with independent named random streams.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

import math
from collections import Counter, defaultdict
from heapq import heapify, heappop, heappush

from dense_arrays._record_validation import semantic_digest
from dense_arrays.planning import (
    BatchSampling,
    CandidateBatch,
    FeedbackSnapshot,
    Fixed,
    GenerationPlan,
)
from dense_arrays.planning.batches.models import CONSTRAINED_SAMPLING_POLICY
from dense_arrays.planning.batches.validation import validate_eligible


def sample_batch(
    plan: GenerationPlan,
    policy: BatchSampling,
    *,
    stream: str,
    feedback: FeedbackSnapshot | None = None,
) -> CandidateBatch:
    """Select once without replacement, retaining each required fixed occurrence."""
    parts = plan.request.parts
    fixed = {r.part_id for r in plan.request.requirements if isinstance(r, Fixed)}
    if not len(fixed) <= policy.size <= len(parts):
        msg = "batch size must cover fixed occurrences and not exceed eligible parts"
        raise ValueError(msg)
    validate_eligible(parts, policy)
    context = {"policy": policy.policy_id, "seed": policy.seed, "stream": stream}
    weights = None if feedback is None else feedback.weights(parts)
    ordered = sorted(
        parts,
        key=lambda p: _priority(p.part_id, context, weights),
    )
    selected = [p for p in ordered if p.part_id in fixed]
    if policy.policy_id == CONSTRAINED_SAMPLING_POLICY:
        from .constrained import select_constrained  # noqa: PLC0415

        selected = select_constrained(ordered, selected, policy)
        return CandidateBatch(
            tuple(p.part_id for p in selected),
            plan.collection_id,
            policy,
            stream,
            feedback=feedback,
        )
    remaining = [p for p in ordered if p.part_id not in fixed]
    if policy.strategy == "uniform":
        selected.extend(remaining[: policy.size - len(selected)])
    else:
        by_group = defaultdict(list)
        for part in reversed(remaining):
            by_group[part.group].append(part)
        counts = Counter(p.group for p in selected)
        queue = [
            (counts[g], semantic_digest({**context, "group": g}), g) for g in by_group
        ]
        heapify(queue)
        while len(selected) < policy.size:
            count, priority, group = heappop(queue)
            selected.append(by_group[group].pop())
            if by_group[group]:
                heappush(queue, (count + 1, priority, group))
    return CandidateBatch(
        tuple(p.part_id for p in selected),
        plan.collection_id,
        policy,
        stream,
        feedback=feedback,
    )


def _priority(part_id: str, context: dict, weights: dict[str, float] | None) -> tuple:
    """Use weighted exponential priorities; equal weights retain hash ordering."""
    key = semantic_digest({**context, "part_id": part_id})
    if weights is None:
        return (key, part_id)
    uniform = 1 - (int(key[:13], 16) + 0.5) / 2**52
    rank = math.log(-math.log(uniform)) - math.log(weights[part_id])
    return (rank, key, part_id)
