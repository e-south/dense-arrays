"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/execution.py

Bounded local execution coordinating pure generation and native commits.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from contextlib import closing
from typing import TYPE_CHECKING

from dense_arrays.artifacts.candidates import CandidateEvidence
from dense_arrays.artifacts.records import Design, RunHandle
from dense_arrays.artifacts.run_plans import cell_plans
from dense_arrays.artifacts.search import batch_search_finished
from dense_arrays.artifacts.store import RunWriter, create_run
from dense_arrays.generation.acceptance import realize
from dense_arrays.generation.assembly import finalize
from dense_arrays.generation.heuristic import HeuristicResult
from dense_arrays.generation.packing import build_optimizer, restore_packing
from dense_arrays.planning import MatrixPlan
from dense_arrays.solver import SolveReport, SolveStatus

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

    from dense_arrays.artifacts.batches import BatchDecision
    from dense_arrays.artifacts.run_plans import RunPlan
    from dense_arrays.generation.batches.progress import BatchWork
    from dense_arrays.generation.packing import PackingEngine
    from dense_arrays.planning import GenerationPlan
    from dense_arrays.realized import RealizedArray
    from dense_arrays.reporting.readers import RecordView


class RunExecutionError(RuntimeError):
    """Execution failed after creating a run; inspect ``run`` for committed evidence."""

    def __init__(self, message: str, *, run: RunHandle) -> None:
        super().__init__(message)
        self.run = run

    @property
    def artifact(self) -> Path:
        """Locate the committed run without interpreting the exception message."""
        return self.run.path


def execute(plan: RunPlan, out: Path) -> RunHandle:
    """Execute an already resolved request in one explicitly new destination."""
    if out.exists():
        msg = f"output destination already exists: {out}"
        raise FileExistsError(msg)
    started = time.monotonic()
    plan.admit_work()
    plan.verify_inputs()
    with create_run(plan, out) as writer:
        continue_run(plan, writer, started)
        return writer.handle


def continue_run(
    plan: RunPlan,
    writer: RunWriter,
    started: float,
    *,
    packings: Generator[RealizedArray, None, None] | None = None,
    replay: RecordView | None = None,
) -> None:
    """Share generation and measured interruption handling for new/resumed runs."""
    try:
        if packings is not None or replay is not None:
            writer.begin_resume(active_seconds=time.monotonic() - started)
        if any(
            p.request.schedule is not None or p.request.resampling is not None
            for p in cell_plans(plan).values()
        ):
            from dense_arrays.workflow import schedules  # noqa: PLC0415

            schedules.generate_schedules(plan, writer, started, replay=replay)
        elif isinstance(plan, MatrixPlan):
            from dense_arrays.workflow.matrices import generate_matrix  # noqa: PLC0415

            generate_matrix(plan, writer, started, replay=replay)
        else:
            _generate(plan, writer, started, packings=packings)
    except (Exception, KeyboardInterrupt) as err:
        writer.refresh()
        elapsed = time.monotonic() - started
        interrupted = isinstance(err, KeyboardInterrupt)
        if writer.state["counts"]["in_progress"]:
            writer.publish(
                writer.state["counts"]["started"],
                "interrupted_unresolved" if interrupted else "error",
                {
                    "code": "interrupted" if interrupted else "execution_error",
                    "detail": str(err),
                },
                active_seconds=elapsed,
            )
        if writer.state["counts"]["accepted"] == writer.state["target"]:
            writer.finish("completed", "target_attained", active_seconds=elapsed)
        else:
            writer.finish(
                "stopped" if interrupted else "failed",
                "interrupted" if interrupted else "execution_error",
                active_seconds=elapsed,
            )
        if interrupted:
            err.add_note(
                f"Committed results remain inspectable at {writer.handle.path}"
            )
            raise
        raise RunExecutionError(str(err), run=writer.handle) from err


def _generate(
    plan: GenerationPlan,
    writer: RunWriter,
    started: float,
    *,
    packings: Generator[RealizedArray, None, None] | None = None,
) -> None:
    limits = plan.request.limits
    optimizer = None
    while writer.state["counts"]["accepted"] < plan.request.target.count:
        elapsed = time.monotonic() - started
        if elapsed >= limits.active_seconds:
            writer.finish("stopped", "active_time_limit", active_seconds=elapsed)
            return
        if writer.state["counts"]["started"] >= limits.attempts:
            writer.finish("stopped", "attempt_limit", active_seconds=elapsed)
            return
        seconds = min(limits.solver_seconds, limits.active_seconds - elapsed)
        if optimizer is None:
            optimizer = _restore_optimizer(
                plan,
                writer,
                seconds=seconds,
                deadline=started + limits.active_seconds,
                packings=packings,
            )
            continue  # Account for model construction before reserving search.
        terminal = search_attempt(plan, writer, optimizer, started, seconds=seconds)
        if terminal is not None:
            writer.finish(*terminal, active_seconds=time.monotonic() - started)
            return
    writer.finish(
        "completed", "target_attained", active_seconds=time.monotonic() - started
    )


def search_attempt(  # noqa: PLR0913 - exact attempt and scheduling context
    plan: GenerationPlan,
    writer: RunWriter,
    optimizer: PackingEngine,
    started: float,
    *,
    seconds: float,
    cell_id: str = "default",
    batch_work: BatchWork | None = None,
    batch_decision: BatchDecision | None = None,
) -> tuple[str, str] | None:
    """Persist one cell's bounded search and final acceptance as a coherent attempt."""
    elapsed = time.monotonic() - started
    batch = plan.request.batch
    attempt = writer.reserve(
        active_seconds=elapsed,
        cell_id=cell_id,
        batch_decision=batch_decision,
        **({"batch_id": batch.batch_id} if batch is not None else {}),
        **({} if batch_work is None else batch_work.to_dict()),
    )
    cell_attempt = writer.state.get("cells", {}).get(cell_id, writer.state)["counts"][
        "started"
    ]
    report = optimizer.solve_report(time_limit_seconds=seconds)
    evidence = {
        **(
            {"packing_objective": optimizer.packing_objective}
            if plan.request.packing_preference
            else {}
        ),
        **({"batch_id": batch.batch_id} if batch is not None else {}),
        **({} if batch_work is None else batch_work.to_dict()),
        **_search_evidence(report),
    }
    elapsed = time.monotonic() - started
    if report.solution is None:
        counts = writer.state.get("cells", {}).get(cell_id, writer.state)["counts"]
        enumerated = sum(counts[k] for k in ("accepted", "rejected", "duplicate"))
        terminal = _search_stop(
            report,
            enumerated=bool(enumerated)
            if batch_work is None
            else batch_work.enumerated,
        )
        writer.publish(
            attempt,
            "error" if terminal[0] == "failed" else "no_candidate",
            evidence,
            active_seconds=elapsed,
            terminal=None
            if batch_work is not None
            and batch_search_finished(evidence, on_unproven=batch_work.on_unproven)
            else terminal,
        )
        return terminal
    design_id = f"d{attempt:08d}"
    realized = realize(
        report.solution,
        plan,
        source_id=f"{writer.handle.run_id}/{cell_id}/{design_id}",
    )
    candidate = finalize(
        realized,
        plan,
        attempt=cell_attempt,
        deadline=started + writer.limits.active_seconds,
    )
    evidence.update(
        assembly_trials=candidate.trials,
        candidate=CandidateEvidence(
            realized, candidate.realized or candidate.last_rejected
        ).to_dict(),
        requirements=list(candidate.requirements),
    )
    elapsed = time.monotonic() - started
    if candidate.realized is None:
        evidence.update(
            code=candidate.reason, requirements=list(candidate.requirements)
        )
        writer.publish(
            attempt,
            "rejected" if candidate.reason != "active_time_limit" else "no_candidate",
            evidence,
            active_seconds=elapsed,
            terminal=("stopped", candidate.reason)
            if candidate.reason == "active_time_limit"
            else None,
        )
        if candidate.reason == "active_time_limit":
            return ("stopped", candidate.reason)
        optimizer.forbid(report.solution)
        return None
    design = Design(
        writer.handle.run_id,
        cell_id,
        design_id,
        writer.state.get("cells", {}).get(cell_id, writer.state)["plan_id"],
        attempt,
        candidate.realized,
        candidate.requirements,
        batch_id=batch.batch_id if batch_work is not None else None,
    )
    writer.publish(attempt, "accepted", evidence, active_seconds=elapsed, design=design)
    progress = writer.state.get("cells", {}).get(cell_id, writer.state)
    if progress["counts"]["accepted"] < plan.request.target.count:
        optimizer.forbid(report.solution)
    return None


def _search_evidence(report: SolveReport | HeuristicResult) -> dict[str, object]:
    """Keep heuristic observations distinct from exact-backend proof reports."""
    if isinstance(report, HeuristicResult):
        return {
            "heuristic": report.evidence.to_dict(),
            "proof_scope": None,
            "termination_reason": f"heuristic_{report.evidence.status}",
        }
    return {
        "solver_status": report.status.value,
        "backend_status": report.backend_status,
        "proof_scope": report.proof_scope,
        "termination_reason": report.termination_reason,
        "detail": report.detail,
    }


def _search_stop(
    report: SolveReport | HeuristicResult, *, enumerated: bool
) -> tuple[str, str]:
    if isinstance(report, HeuristicResult):
        return "stopped", f"heuristic_{report.evidence.status}"
    if report.status is SolveStatus.INFEASIBLE:
        return "stopped", "batch_exhausted" if enumerated else "batch_infeasible"
    failed = report.status in {SolveStatus.BACKEND_ERROR, SolveStatus.INVALID_RESULT}
    return "failed" if failed else "stopped", f"solver_{report.status.value}"


def _restore_optimizer(
    plan: GenerationPlan,
    writer: RunWriter,
    *,
    seconds: float,
    deadline: float,
    packings: Generator[RealizedArray, None, None] | None,
) -> PackingEngine:
    """Rebuild the offered model and its exclusions within the active allowance."""
    optimizer = build_optimizer(plan, seconds=seconds)
    if optimizer.solver_identity is not None:
        writer.bind_solver(optimizer.solver_identity)
    if packings is not None:
        with closing(packings):
            for packed in packings:
                if time.monotonic() >= deadline:
                    break
                optimizer.forbid(restore_packing(packed, plan))
    return optimizer
