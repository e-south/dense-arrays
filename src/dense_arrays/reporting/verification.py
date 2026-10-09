"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/verification.py

Independently verify native evidence at one immutable revision.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter
from contextlib import closing
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, mutable_json, semantic_digest
from dense_arrays.artifacts.reading import ReadBudget, ReadLimits, Verification
from dense_arrays.artifacts.run_plans import cell_plans, validate_run_binding
from dense_arrays.artifacts.store import reader, stored_plan
from dense_arrays.generation.acceptance import evaluate
from dense_arrays.planning.batches.bindings import membership_size
from dense_arrays.playback.serialization import realized_array_to_dict
from dense_arrays.reporting.batch_accounting import BatchAccounting
from dense_arrays.reporting.candidates import verify_candidate
from dense_arrays.reporting.objectives import ObjectiveHistory
from dense_arrays.reporting.readers import RecordView, read_records
from dense_arrays.reporting.runtime_batches import read_runtime_batches
from dense_arrays.reporting.search import SearchHistory

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.planning import GenerationPlan
    from dense_arrays.reporting.runtime_batches import RuntimeBatches
    from dense_arrays.reporting.summary import RunSummary


def verify_run(
    path: Path,
    summary: RunSummary,
    limits: ReadLimits | None = None,
    *,
    budget: ReadBudget | None = None,
) -> Verification:
    """Recount outcomes and validate placement, requirement and provenance joins."""
    budget = budget if budget is not None else ReadBudget(limits or ReadLimits())
    with reader(path) as connection:
        payload = connection.execute("SELECT payload FROM plan WHERE id=1").fetchone()[
            0
        ]
        budget.examine(payload)
        plan = stored_plan(connection, max_identities=budget.limits.identities)
        plans = cell_plans(plan)
        budget.retain(
            sum(
                len(p.request.parts)
                + len(p.request.requirements)
                + membership_size(p.request)
                for p in plans.values()
            )
        )
    validate_run_binding(plan, summary.to_dict())
    bound = max(summary.counts["started"] + 1, summary.accepted + 1)
    excluded = {
        (cell, e.sequence_id): e.design_ref
        for cell, p in plans.items()
        for e in p.exclusions
    }
    budget.retain(len(excluded))
    runtime = read_runtime_batches(path, summary, plans, budget)
    accepted = _verify_attempts(
        path, summary, plans, budget=budget, excluded=excluded, runtime=runtime
    )
    runtime.finish()
    _verify_designs(
        path,
        summary,
        bound,
        plans,
        accepted,
        budget=budget,
        excluded=excluded,
        runtime=runtime,
    )
    return Verification(
        "committed_plan_attempts_designs", budget.examined, budget.bytes_checked
    )


def _verify_attempts(  # noqa: PLR0913 - one snapshot and shared evidence budget
    path: Path,
    summary: RunSummary,
    plans: dict[str, GenerationPlan],
    *,
    budget: ReadBudget,
    excluded: dict[tuple[str, str], str],
    runtime: RuntimeBatches,
) -> dict[int, tuple[str, str | None, str | None]]:
    counts = Counter()
    cells = {name: Counter() for name in plans}
    batches = {
        name: BatchAccounting(p.request.schedule or p.request.resampling)
        for name, p in plans.items()
    }
    accepted = {}
    sequences = set()
    objectives = ObjectiveHistory(plans, budget)
    searches = SearchHistory(plans, budget)
    query = RecordView(
        path,
        summary.revision,
        "attempts",
        summary.counts["started"] + 1,
        run_id=summary.run_id,
    )
    with closing(read_records(query, budget)) as records:
        for record in records:
            if record.cell_id not in plans:
                msg = "unsupported attempt schema or cell"
                raise ValueError(msg)
            counts[record.outcome] += 1
            counts["started"] += 1
            cells[record.cell_id][record.outcome] += 1
            cells[record.cell_id]["started"] += 1
            _check_cell_ordinal(
                record.evidence.get("cell_attempt"),
                cells[record.cell_id]["started"],
                required=bool(summary.cells),
            )
            batches[record.cell_id].observe(record)
            decision = runtime.observe(record, plans[record.cell_id])
            candidate = verify_candidate(
                record,
                plans[record.cell_id],
                run_id=summary.run_id,
                decision=decision,
            )
            objectives.observe(record, plans[record.cell_id], decision)
            searches.observe(record, plans[record.cell_id])
            if candidate is not None and record.outcome == "duplicate":
                if (record.cell_id, candidate[0]) not in sequences and (
                    record.cell_id,
                    candidate[0],
                ) not in excluded:
                    msg = (
                        "duplicate candidate has no earlier accepted "
                        "or excluded sequence"
                    )
                    raise ValueError(msg)
                if (
                    record.evidence.get("code") == "parent_duplicate"
                    and candidate[0] != record.evidence["candidate_sequence_id"]
                ):
                    msg = "parent duplicate sequence does not match its candidate"
                    raise ValueError(msg)
            if (
                record.evidence.get("code") == "parent_duplicate"
                and excluded.get(
                    (record.cell_id, record.evidence["candidate_sequence_id"])
                )
                != record.evidence["matched_design_ref"]
            ):
                msg = "parent duplicate does not match the frozen exclusion evidence"
                raise ValueError(msg)
            if record.outcome == "accepted":
                budget.retain()
                accepted[record.attempt_id] = (
                    record.design_ref,
                    None if candidate is None else candidate[1],
                    record.evidence.get("batch_id")
                    if plans[record.cell_id].request.schedule is not None
                    or plans[record.cell_id].request.resampling is not None
                    else None,
                )
                if candidate is not None:
                    budget.retain()
                    sequences.add((record.cell_id, candidate[0]))
    if any(counts[name] != value for name, value in summary.counts.items()) or set(
        counts
    ) - set(summary.counts):
        msg = "stored attempts do not reconcile with committed counters"
        raise ValueError(msg)
    _reconcile_cells(summary, cells)
    return accepted


def _reconcile_cells(summary: RunSummary, cells: dict[str, Counter]) -> None:
    """Check each cell independently of the aggregate attempt total."""
    for name, cell in summary.cells.items():
        if any(cells[name][key] != count for key, count in cell.counts.items()):
            msg = "stored attempts do not reconcile with cell counters"
            raise ValueError(msg)


def _check_cell_ordinal(value: object, expected: int, *, required: bool) -> None:
    """Preserve logical cell work independently of interleaved run ordinals."""
    if not required:
        return
    integer(value, field_name="cell_attempt", minimum=1)
    if value != expected:
        msg = "cell attempt ordinals are not contiguous"
        raise ValueError(msg)


def _verify_designs(  # noqa: PLR0913 - one verification scope and shared budget
    path: Path,
    summary: RunSummary,
    bound: int,
    plans: dict[str, GenerationPlan],
    accepted: dict[int, tuple[str, str | None, str | None]],
    *,
    budget: ReadBudget,
    excluded: dict[tuple[str, str], str],
    runtime: RuntimeBatches,
) -> None:
    seen_sequences = set()
    query = RecordView(path, summary.revision, "designs", bound, run_id=summary.run_id)
    with closing(read_records(query, budget)) as records:
        for design in records:
            if (
                design.run_id != summary.run_id
                or design.cell_id not in plans
                or design.plan_id != plans[design.cell_id].plan_id
            ):
                msg = "design provenance does not match the run and plan"
                raise ValueError(msg)
            attempt = accepted.pop(design.attempt_id, None)
            if (
                attempt is None
                or attempt[0] != design.reference
                or design.realized.source_id != design.reference
                or attempt[2] != design.batch_id
            ):
                msg = "design and accepted-attempt identities do not join"
                raise ValueError(msg)
            root = runtime.design_plan(
                design.cell_id, design.batch_id, plans[design.cell_id]
            )
            evidence = evaluate(design.realized, root, batch_id=design.batch_id)
            if tuple(
                mutable_json(r) for r in design.requirements
            ) != evidence or not all(r["passed"] for r in evidence):
                msg = "design requirements do not match independent recount"
                raise ValueError(msg)
            if attempt[1] is not None and attempt[1] != semantic_digest(
                realized_array_to_dict(design.realized)
            ):
                msg = "accepted design differs from its recorded final candidate"
                raise ValueError(msg)
            if (design.cell_id, design.sequence_id) in seen_sequences:
                msg = "duplicate accepted final sequence in the same cell"
                raise ValueError(msg)
            if (design.cell_id, design.sequence_id) in excluded:
                msg = "accepted child design duplicates a parent or ancestor sequence"
                raise ValueError(msg)
            budget.retain()
            seen_sequences.add((design.cell_id, design.sequence_id))
    if accepted or len(seen_sequences) != summary.accepted:
        msg = "accepted designs are missing or counters disagree"
        raise ValueError(msg)
