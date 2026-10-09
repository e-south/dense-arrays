"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/matrices.py

Local round-robin generation under one matrix run's shared effort budget.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from typing import TYPE_CHECKING

from dense_arrays.generation.packing import restore_packing
from dense_arrays.workflow.execution import _restore_optimizer, search_attempt

if TYPE_CHECKING:
    from dense_arrays.artifacts.store import RunWriter
    from dense_arrays.generation.packing import PackingEngine
    from dense_arrays.planning import MatrixPlan
    from dense_arrays.reporting.readers import RecordView


def generate_matrix(
    plan: MatrixPlan,
    writer: RunWriter,
    started: float,
    *,
    replay: RecordView | None = None,
) -> None:
    """Offer each active cell one attempt per round; never redistribute targets."""
    optimizers, after_cell = _restore_cells(plan, writer, started, replay)
    order = list(plan.cells)
    if after_cell is not None:
        pivot = (
            next(i for i, cell in enumerate(order) if cell.cell_id == after_cell) + 1
        )
        order = order[pivot:] + order[:pivot]
    while True:
        active = [
            c
            for c in order
            if writer.state["cells"][c.cell_id]["state"] in {"created", "running"}
        ]
        if not active:
            attained = writer.state["counts"]["accepted"] == plan.total
            failed = any(c["state"] == "failed" for c in writer.state["cells"].values())
            writer.finish(
                "completed" if attained else "failed" if failed else "stopped",
                "target_attained"
                if attained
                else "cell_failure"
                if failed
                else "cells_exhausted",
                active_seconds=time.monotonic() - started,
            )
            return
        for cell in active:
            elapsed = time.monotonic() - started
            reason = (
                "active_time_limit"
                if elapsed >= writer.limits.active_seconds
                else "attempt_limit"
                if writer.state["counts"]["started"] >= writer.limits.attempts
                else None
            )
            if reason:
                writer.finish("stopped", reason, active_seconds=elapsed)
                return
            seconds = min(
                writer.limits.solver_seconds, writer.limits.active_seconds - elapsed
            )
            if cell.cell_id not in optimizers:
                optimizers[cell.cell_id] = _restore_optimizer(
                    cell.plan,
                    writer,
                    seconds=seconds,
                    deadline=started + writer.limits.active_seconds,
                    packings=None,
                )
                # Recheck the shared clock before admitting an attempt.
                elapsed = time.monotonic() - started
                if elapsed >= writer.limits.active_seconds:
                    writer.finish(
                        "stopped", "active_time_limit", active_seconds=elapsed
                    )
                    return
                seconds = min(
                    writer.limits.solver_seconds, writer.limits.active_seconds - elapsed
                )
            search_attempt(
                cell.plan,
                writer,
                optimizers[cell.cell_id],
                started,
                seconds=seconds,
                cell_id=cell.cell_id,
            )
            if writer.state["cells"][cell.cell_id]["state"] not in {
                "created",
                "running",
            }:
                optimizers.pop(cell.cell_id, None)


def _restore_cells(
    plan: MatrixPlan, writer: RunWriter, started: float, replay: RecordView | None
) -> tuple[dict[str, PackingEngine], str | None]:
    """Scan the committed attempt history once, restoring only open cell models."""
    optimizers = {}
    after_cell = None
    if replay is None:
        return optimizers, after_cell
    plans = {cell.cell_id: cell.plan for cell in plan.cells}
    deadline = started + writer.limits.active_seconds
    with replay.records() as records:
        for attempt in records:
            after_cell = attempt.cell_id
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            if (
                attempt.candidate is None
                or writer.state["cells"][attempt.cell_id]["state"] != "running"
            ):
                continue
            if attempt.cell_id not in optimizers:
                optimizers[attempt.cell_id] = _restore_optimizer(
                    plans[attempt.cell_id],
                    writer,
                    seconds=min(writer.limits.solver_seconds, remaining),
                    deadline=deadline,
                    packings=None,
                )
            if time.monotonic() >= deadline:
                break
            optimizers[attempt.cell_id].forbid(
                restore_packing(attempt.candidate.packed, plans[attempt.cell_id])
            )
    return optimizers, after_cell
