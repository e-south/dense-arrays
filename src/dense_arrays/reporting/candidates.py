"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/candidates.py

Verify recorded candidate geometry and screening without regenerating it.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import replace

from dense_arrays._record_validation import mutable_json, semantic_digest
from dense_arrays.artifacts.batches import BatchDecision
from dense_arrays.artifacts.records import Attempt
from dense_arrays.generation.acceptance import evaluate
from dense_arrays.generation.packing import restore_packing
from dense_arrays.planning import GC, Avoid, GenerationPlan
from dense_arrays.planning.batches.bindings import offered_batch
from dense_arrays.playback.serialization import realized_array_to_dict
from dense_arrays.realized import RealizedArray


def verify_candidate(
    attempt: Attempt,
    plan: GenerationPlan,
    *,
    run_id: str,
    decision: BatchDecision | None = None,
) -> tuple[str, str] | None:
    """Return final sequence/record digests, or explicit absent final evidence."""
    plan = _candidate_plan(attempt, plan, run_id, decision)
    candidate = attempt.candidate
    if candidate is None:
        return None
    packed, final = candidate.packed, candidate.final
    expected_id = f"{run_id}/{attempt.cell_id}/d{attempt.attempt_id:08d}"
    if packed.source_id != expected_id:
        msg = "candidate identity does not match its run and attempt"
        raise ValueError(msg)
    restore_packing(packed, plan, batch_id=attempt.evidence.get("batch_id"))
    assembly = plan.request.assembly
    padding = None if assembly is None else assembly.padding
    padded = padding is not None and len(packed.sequence) < plan.request.length.exact
    maximum = padding.max_trials if padded else 1
    trials = attempt.evidence.get("assembly_trials")
    if trials is None or not 0 <= trials <= maximum:
        msg = "candidate assembly trials exceed the declared effort"
        raise ValueError(msg)
    if final is None:
        if attempt.evidence.get("code") != "active_time_limit" or attempt.evidence.get(
            "requirements"
        ):
            msg = "unevaluated candidate must record an active-time stop without checks"
            raise ValueError(msg)
        return None
    evidence = evaluate(final, plan, batch_id=attempt.evidence.get("batch_id"))
    if (
        tuple(mutable_json(r) for r in attempt.evidence.get("requirements", ()))
        != evidence
    ):
        msg = "candidate requirement evidence does not match independent recount"
        raise ValueError(msg)
    _verify_transform(packed, final, trials=trials, assembled=assembly is not None)
    failed = {r["id"] for r in evidence if not r["passed"]}
    screen_ids = {r.id for r in plan.request.requirements if isinstance(r, (Avoid, GC))}
    if failed - screen_ids:
        msg = "candidate violates a packing requirement"
        raise ValueError(msg)
    code = attempt.evidence.get("code")
    valid = {
        "accepted": not failed,
        "duplicate": not failed,
        "rejected": bool(failed)
        and trials == maximum
        and code == ("padding_trials_exhausted" if padded else "screening_rejection"),
        "no_candidate": bool(failed)
        and trials < maximum
        and code == "active_time_limit",
    }.get(attempt.outcome, False)
    if not valid:
        msg = "candidate outcome disagrees with screening or assembly effort"
        raise ValueError(msg)
    return (
        semantic_digest(
            {"schema": "dense_arrays.sequence.v1", "sequence": final.sequence}
        ),
        semantic_digest(realized_array_to_dict(final)),
    )


def _verify_transform(
    packed: RealizedArray, final: RealizedArray, *, trials: int, assembled: bool
) -> None:
    if not assembled:
        matches = final == packed and trials == 1
    else:
        transform = final.provenance["assembly"]
        shift = transform["packed_start"]
        matches = (
            final.sequence[shift : shift + len(packed.sequence)] == packed.sequence
            and final.placements
            == tuple(replace(p, start=p.start + shift) for p in packed.placements)
            and trials == transform["trial"]
        )
    if not matches:
        msg = "final candidate does not preserve its packing and assembly trial"
        raise ValueError(msg)


def _candidate_plan(
    attempt: Attempt, plan: GenerationPlan, run_id: str, decision: BatchDecision | None
) -> GenerationPlan:
    """Resolve a candidate's saved offered geometry and batch coordinates."""
    if decision is not None:
        if (
            decision.cell_id != attempt.cell_id
            or decision.run_id != run_id
            or decision.index != attempt.evidence.get("batch_index")
        ):
            msg = "candidate does not match its runtime batch decision"
            raise ValueError(msg)
        plan = decision.bind(plan)
    batch = offered_batch(
        plan.request,
        attempt.evidence.get("batch_id"),
        index=attempt.evidence.get("batch_index") if decision is None else None,
    )
    if attempt.evidence.get("batch_id") != (None if batch is None else batch.batch_id):
        msg = "attempt batch identity does not match its offered plan"
        raise ValueError(msg)
    return plan
