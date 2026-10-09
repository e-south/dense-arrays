"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/bundles/batches.py

Bound portable runtime memberships; feedback histories remain unavailable.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.batches import BatchDecision
from dense_arrays.artifacts.run_state import origin_bindings
from dense_arrays.artifacts.store import checked_payload
from dense_arrays.planning.batches.bindings import encoded_batch_size

if TYPE_CHECKING:
    import sqlite3
    from collections.abc import Iterator, Mapping

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.planning import PlanEvidence
    from dense_arrays.reporting.bundles.views import BundleView


def load_batches(
    connection: sqlite3.Connection,
    manifest: Mapping[str, object],
    plans: dict[str, PlanEvidence],
    budget: ReadBudget,
) -> dict[tuple[str, str, str], BatchDecision]:
    """Check membership, policy, origin and count without drawing or replaying."""
    result = {}
    count = manifest.get("batches", 0)
    if "batches" not in manifest:
        return result
    origins = set(origin_bindings(manifest["source_runs"]))
    prefixes = {}
    for source in manifest["source_runs"]:
        prefixes[source["run_id"]] = max(
            prefixes.get(source["run_id"], 0), source["counts"]["started"]
        )
    for row in connection.execute(
        "SELECT run_id,cell_id,batch_index,batch_id,payload,digest FROM batches"
    ):
        budget.examine(row[-2])
        decision = BatchDecision.from_dict(checked_payload(row[-2:]))
        key = (decision.run_id, decision.cell_id, decision.batch.batch_id)
        if (row[0], row[1], row[2], row[3]) != (
            decision.run_id,
            decision.cell_id,
            decision.index,
            decision.batch.batch_id,
        ) or (decision.run_id, decision.cell_id, decision.plan_id) not in origins:
            msg = "portable batch decision origin or index mismatch"
            raise ValueError(msg)
        if decision.after_attempt >= prefixes[decision.run_id]:
            msg = "portable batch decision exceeds the original attempt prefix"
            raise ValueError(msg)
        decision.bind(plans[decision.plan_id])
        if key in result:
            msg = "portable runtime membership repeats an identity"
            raise ValueError(msg)
        feedback = decision.batch.feedback
        budget.retain(
            1
            + len(decision.batch.part_ids)
            + (0 if feedback is None else len(feedback.used) + len(feedback.failed))
        )
        result[key] = decision
    if len(result) != count:
        msg = "portable batch count does not match its manifest"
        raise ValueError(msg)
    return result


def read_batch_page(
    connection: sqlite3.Connection, query: BundleView, budget: ReadBudget
) -> Iterator[BatchDecision]:
    """Page immutable saved decisions by the contained insertion ordinal."""
    if "batches" not in query.summary.manifest:
        return
    start = 0 if query.after is None else query.after.ordinal
    count = 0
    for row in connection.execute(
        "SELECT rowid,payload,digest FROM batches WHERE rowid>? ORDER BY rowid LIMIT ?",
        (start, -1 if query.limit is None else query.limit),
    ):
        budget.examine(row[1])
        value = checked_payload(row[1:])
        size = 2 + encoded_batch_size(value.get("batch"))
        budget.retain(size)
        budget.identities -= size
        record = BatchDecision.from_dict(value)
        budget.position = row[0]
        count += 1
        yield record
    if query.limit is None and start + count != query.summary.manifest["batches"]:
        msg = "portable batch records do not reconcile with the manifest"
        raise ValueError(msg)
