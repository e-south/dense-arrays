"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/schedules.py

Round-robin execution and replay of bounded prepared batch schedules.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from typing import TYPE_CHECKING

from dense_arrays.artifacts.batches import BatchDecision
from dense_arrays.artifacts.records import Attempt
from dense_arrays.artifacts.run_plans import cell_plans
from dense_arrays.artifacts.search import batch_search_finished
from dense_arrays.artifacts.store import checked_payload
from dense_arrays.generation.packing import restore_packing
from dense_arrays.workflow.batch_cursor import BatchCursor
from dense_arrays.workflow.execution import _restore_optimizer, search_attempt

if TYPE_CHECKING:
    from dense_arrays.artifacts.run_plans import RunPlan
    from dense_arrays.artifacts.store import RunWriter
    from dense_arrays.reporting.readers import RecordView


def _progress(writer: RunWriter, cell: str) -> dict[str, object]:
    return writer.state.get("cells", {}).get(cell, writer.state)


def _open(writer: RunWriter, cell: str) -> bool:
    progress = _progress(writer, cell)
    return (
        progress["state"] in {"created", "running"}
        and progress["counts"]["accepted"] < progress["target"]
    )


def _stop(writer: RunWriter, cell: str, reason: str, started: float) -> None:
    elapsed = time.monotonic() - started
    if "cells" in writer.state:
        if _open(writer, cell):
            writer.stop_cell(cell, "stopped", reason, active_seconds=elapsed)
    else:
        writer.finish("stopped", reason, active_seconds=elapsed)


def generate_schedules(
    plan: RunPlan, writer: RunWriter, started: float, *, replay: RecordView | None
) -> None:
    """Advance selections only at declared batch boundaries, sharing global effort."""
    cursors = {name: BatchCursor(p) for name, p in cell_plans(plan).items()}
    after = _restore(cursors, writer, started, replay)
    order = list(cursors)
    if after is not None:
        pivot = order.index(after) + 1
        order = order[pivot:] + order[:pivot]
    while True:
        if writer.state["counts"]["accepted"] == writer.state["target"]:
            writer.finish(
                "completed",
                "target_attained",
                active_seconds=time.monotonic() - started,
            )
            return
        active = [name for name in order if _open(writer, name)]
        if not active:
            if "cells" in writer.state:
                failed = any(
                    c["state"] == "failed" for c in writer.state["cells"].values()
                )
                writer.finish(
                    "failed" if failed else "stopped",
                    "cell_failure" if failed else "cells_exhausted",
                    active_seconds=time.monotonic() - started,
                )
            return
        for cell in active:
            cursor = cursors[cell]
            if not cursor.advance():
                _stop(
                    writer,
                    cell,
                    "batch_limit"
                    if cursor.original.request.resampling is not None
                    else "batch_schedule_exhausted",
                    started,
                )
                continue
            reason = _step(cursor, cell, writer, started)
            if reason is not None:
                writer.finish(
                    "stopped", reason, active_seconds=time.monotonic() - started
                )
                return


def _step(
    cursor: BatchCursor, cell: str, writer: RunWriter, started: float
) -> str | None:
    elapsed = time.monotonic() - started
    if elapsed >= writer.limits.active_seconds:
        return "active_time_limit"
    if writer.state["counts"]["started"] >= writer.limits.attempts:
        return "attempt_limit"
    remaining = writer.limits.active_seconds - elapsed
    if cursor.optimizer is None:
        cursor.select(cell)
        cursor.optimizer = _restore_optimizer(
            cursor.plan,
            writer,
            seconds=min(writer.limits.solver_seconds, remaining),
            deadline=started + writer.limits.active_seconds,
            packings=None,
        )
    remaining = writer.limits.active_seconds - (time.monotonic() - started)
    if remaining <= 0:
        return "active_time_limit"
    terminal = search_attempt(
        cursor.plan,
        writer,
        cursor.optimizer,
        started,
        seconds=min(writer.limits.solver_seconds, remaining),
        cell_id=cell,
        batch_work=cursor.work,
        batch_decision=cursor.decision(
            writer.handle.run_id, cell, writer.state["counts"]["started"]
        ),
    )
    record = checked_payload(
        writer.connection.execute(
            "SELECT payload,digest FROM attempts WHERE attempt=? "
            "ORDER BY revision DESC LIMIT 1",
            (writer.state["counts"]["started"],),
        ).fetchone()
    )
    cursor.observe(Attempt.from_dict(record))
    if terminal is not None and not batch_search_finished(
        record["evidence"],
        on_unproven=cursor.policy.on_unproven if cursor.policy else "stop",
    ):
        if "cells" not in writer.state:
            writer.finish(*terminal, active_seconds=time.monotonic() - started)
    elif terminal is not None and cursor.policy is None:
        _stop(writer, cell, terminal[1], started)
    return None


def _restore(
    cursors: dict[str, BatchCursor],
    writer: RunWriter,
    started: float,
    replay: RecordView | None,
) -> str | None:
    """Recount frontiers, then restore current batches' committed exclusions."""
    if replay is None:
        return None
    after = None
    deadline = started + writer.limits.active_seconds
    with replay.records() as records:
        for attempt in records:
            if time.monotonic() >= deadline:
                return after
            after = attempt.cell_id
            _observe_saved(cursors[after], attempt, writer)
    for cursor in cursors.values():
        cursor.advance()
    with replay.records() as records:
        for attempt in records:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            cursor = cursors[attempt.cell_id]
            if (
                not _open(writer, attempt.cell_id)
                or cursor.ended
                or attempt.candidate is None
                or attempt.evidence.get("batch_index", 1) != cursor.index
            ):
                continue
            if cursor.optimizer is None:
                cursor.optimizer = _restore_optimizer(
                    cursor.plan,
                    writer,
                    seconds=min(writer.limits.solver_seconds, remaining),
                    deadline=deadline,
                    packings=None,
                )
            if time.monotonic() < deadline:
                cursor.optimizer.forbid(
                    restore_packing(attempt.candidate.packed, cursor.plan)
                )
    return after


def _observe_saved(cursor: BatchCursor, attempt: Attempt, writer: RunWriter) -> None:
    """Restore each recorded runtime membership once before observing outcomes."""
    if cursor.original.request.resampling is not None and (
        cursor.batch is None or cursor.index != attempt.evidence["batch_index"]
    ):
        saved = checked_payload(
            writer.connection.execute(
                "SELECT payload,digest FROM batches WHERE cell_id=? AND batch_index=?",
                (attempt.cell_id, attempt.evidence["batch_index"]),
            ).fetchone()
        )
        cursor.restore(BatchDecision.from_dict(saved))
    cursor.observe(attempt)
