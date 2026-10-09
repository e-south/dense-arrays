"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/recovery.py

Resume an unchanged measured interruption from its verified committed prefix.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays.artifacts.errors import integrity_boundary
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.reading import ReadLimits
from dense_arrays.artifacts.records import RunHandle
from dense_arrays.artifacts.recovery import RecoveryError, own_run
from dense_arrays.artifacts.run_plans import cell_plans, excluded_sequences, run_limits
from dense_arrays.artifacts.search import batch_search_finished
from dense_arrays.artifacts.store import RunWriter, latest, stored_plan
from dense_arrays.planning import MatrixPlan
from dense_arrays.reporting.plans.reading import plan_identities
from dense_arrays.reporting.readers import RecordView
from dense_arrays.reporting.summary import RunSummary
from dense_arrays.reporting.verification import verify_run
from dense_arrays.workflow.execution import continue_run

if TYPE_CHECKING:
    from collections.abc import Generator
    from pathlib import Path

    from dense_arrays.artifacts.run_plans import RunPlan
    from dense_arrays.realized import RealizedArray


def resume_run(path: Path) -> RunHandle:
    """Acquire ownership, verify all saved evidence, then admit bounded work."""
    started = time.monotonic()
    with own_run(path) as connection:
        with integrity_boundary(path):
            state = latest(connection)
            summary = RunSummary.from_manifest(state)
            plan = stored_plan(connection)

        with integrity_boundary(path):
            if summary.counts["started"] > run_limits(plan).attempts:
                msg = "recorded attempts exceed the immutable plan budget"
                raise ValueError(msg)
            limits = ReadLimits(
                records=1
                + summary.counts["started"]
                + summary.accepted
                + (summary.batch_count or 0),
                identities=(
                    plan_identities(plan)
                    + 3 * summary.accepted
                    + 1
                    + (summary.batch_count or 0)
                    * max(
                        (
                            3 * len(p.request.parts) + 2
                            for p in cell_plans(plan).values()
                        ),
                        default=0,
                    )
                    + sum(
                        2 * len(p.request.parts)
                        for p in cell_plans(plan).values()
                        if p.request.resampling is not None
                    )
                    + sum(
                        3 * len(p.request.parts)
                        for p in cell_plans(plan).values()
                        if p.request.packing_preference
                    )
                    + sum(
                        p.request.search == "greedy" for p in cell_plans(plan).values()
                    )
                ),
            )
            verify_run(path, summary, limits)
        plan.verify_inputs()
        handle = RunHandle(path, summary.run_id)
        if summary.state == "completed":
            return handle
        _admit(path, plan, summary)
        plan.admit_work()
        candidates = _candidates(path, plan, summary, limits)
        writer = RunWriter(
            connection,
            handle,
            state,
            excluded_sequences(plan),
            run_limits(plan),
        )
        started -= summary.active_seconds
        if isinstance(plan, MatrixPlan) or any(
            p.request.schedule is not None or p.request.resampling is not None
            for p in cell_plans(plan).values()
        ):
            continue_run(plan, writer, started, replay=candidates)
        else:
            continue_run(plan, writer, started, packings=_packings(candidates))
        return handle


def _packings(view: RecordView) -> Generator[RealizedArray, None, None]:
    """Own the replay reader, including early closure when time expires."""
    with view.records() as records:
        for attempt in records:
            if attempt.candidate is not None:
                yield attempt.candidate.packed


def _admit(path: Path, plan: RunPlan, summary: RunSummary) -> None:
    if summary.state in {"created", "running"}:
        code = "active_time_unknown"
        message = "abandoned execution has unknown active time; create a new linked run"
        raise RecoveryError(code, message, artifact=path)
    if summary.counts["started"] >= run_limits(plan).attempts:
        code = "attempt_budget_exhausted"
        message = "original attempt budget is exhausted; create a new linked run"
        raise RecoveryError(code, message, artifact=path)
    if summary.active_seconds >= run_limits(plan).active_seconds:
        code = "active_budget_exhausted"
        message = "original active-time budget is exhausted; create a new linked run"
        raise RecoveryError(code, message, artifact=path)
    if (
        summary.state != "stopped"
        or summary.termination_reason != "interrupted"
        or not summary.resumable
    ):
        code = "not_resumable"
        message = (
            "recorded termination is not recoverable; "
            "inspect outcomes and create a new linked run"
        )
        raise RecoveryError(code, message, artifact=path)
    if (
        summary.producer is None
        or replace(summary.producer, solver=None) != Producer.capture()
    ):
        code = "producer_changed"
        message = (
            "execution environment differs from the recorded producer; "
            "use the original environment or a new linked run"
        )
        raise RecoveryError(code, message, artifact=path)


def _candidates(
    path: Path, plan: RunPlan, summary: RunSummary, limits: ReadLimits
) -> RecordView:
    scheduled = {
        cell: recipe.request.schedule or recipe.request.resampling
        for cell, recipe in cell_plans(plan).items()
        if recipe.request.schedule is not None or recipe.request.resampling is not None
    }
    view = RecordView(
        path,
        summary.revision,
        "attempts",
        None,
        run_id=summary.run_id,
        read_limits=limits,
    )
    with view.records() as records:
        for attempt in records:
            if (
                summary.cells
                and summary.cells[attempt.cell_id].termination_reason != "interrupted"
            ):
                continue
            if attempt.outcome in {"accepted", "duplicate", "rejected"}:
                if attempt.candidate is None:
                    code = "missing_candidate"
                    message = (
                        "saved packing evidence is unavailable; create a new linked run"
                    )
                    raise RecoveryError(code, message, artifact=path)
            elif (
                attempt.cell_id in scheduled
                and attempt.outcome == "no_candidate"
                and batch_search_finished(
                    attempt.evidence, on_unproven=scheduled[attempt.cell_id].on_unproven
                )
            ):
                # Verification already checked this batch's binding and frontier.
                # The declared boundary closes this batch; later batches may run.
                continue
            elif attempt.outcome != "interrupted_unresolved":
                code = "search_terminated"
                message = (
                    "recorded search termination cannot be resumed; "
                    "inspect its solver outcome"
                )
                raise RecoveryError(code, message, artifact=path)
    return view
