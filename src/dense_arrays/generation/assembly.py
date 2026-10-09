"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/assembly.py

Assemble exact-length candidates while retaining a reversible coordinate transform.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

import time
from dataclasses import dataclass, replace

from dense_arrays.generation.acceptance import evaluate
from dense_arrays.generation.randomness import PADDING_POLICY, padding_dna
from dense_arrays.planning import GC, Avoid, GenerationPlan
from dense_arrays.realized import RealizedArray


@dataclass(frozen=True)
class CandidateResult:
    """Final candidate or last rejection, with effort separate from solver attempts."""

    realized: RealizedArray | None
    requirements: tuple[dict[str, object], ...]
    trials: int
    reason: str
    last_rejected: RealizedArray | None = None


def finalize(
    packed: RealizedArray, plan: GenerationPlan, *, attempt: int, deadline: float
) -> CandidateResult:
    """Bound final assembly checks without claiming padding infeasibility."""
    padding = None if plan.request.assembly is None else plan.request.assembly.padding
    count = (
        0
        if plan.request.length.exact is None
        else plan.request.length.exact - len(packed.sequence)
    )
    maximum = padding.max_trials if padding is not None and count > 0 else 1
    screen_ids = {r.id for r in plan.request.requirements if isinstance(r, (Avoid, GC))}
    results = ()
    last_rejected = None
    for trial in range(1, maximum + 1):
        if time.monotonic() >= deadline:
            return CandidateResult(
                None, results, trial - 1, "active_time_limit", last_rejected
            )
        realized = assemble(packed, plan, attempt=attempt, trial=trial)
        results = evaluate(realized, plan)
        if any(not r["passed"] and r["id"] not in screen_ids for r in results):
            msg = "solver result failed independent requirement recount"
            raise ValueError(msg)
        if all(r["passed"] for r in results):
            return CandidateResult(realized, results, trial, "accepted")
        last_rejected = realized
    reason = (
        "padding_trials_exhausted"
        if padding is not None and count
        else "screening_rejection"
    )
    return CandidateResult(None, results, maximum, reason, last_rejected)


def assemble(
    packed: RealizedArray, plan: GenerationPlan, *, attempt: int, trial: int
) -> RealizedArray:
    """Generate one deterministic proposal with an explicit coordinate transform."""
    if plan.request.assembly is None:
        return packed
    padding = plan.request.assembly.padding
    count = plan.request.length.exact - len(packed.sequence)
    if count < 0 or (count and padding is None):
        msg = "packing length cannot satisfy the declared assembly policy"
        raise ValueError(msg)
    dna, stream_id = padding_dna(
        count, seed=plan.request.seed, attempt=attempt, trial=trial
    )
    side = None if padding is None else padding.side
    shift = count if side == "left" else 0
    sequence = dna + packed.sequence if side == "left" else packed.sequence + dna
    provenance = dict(packed.provenance)
    provenance["assembly"] = {
        "policy": PADDING_POLICY,
        "side": side,
        "padding_length": count,
        "packed_start": shift,
        "packed_length": len(packed.sequence),
        "trial": trial,
        "stream_id": stream_id,
    }
    return replace(
        packed,
        sequence=sequence,
        placements=tuple(replace(p, start=p.start + shift) for p in packed.placements),
        provenance=provenance,
    )
