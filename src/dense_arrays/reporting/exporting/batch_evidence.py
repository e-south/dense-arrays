"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/exporting/batch_evidence.py

Copy runtime membership evidence once per origin and logical selection.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.batches import BatchDecision
from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE
from dense_arrays.artifacts.store import checked_payload, reader
from dense_arrays.reporting.bundles import BundleView
from dense_arrays.reporting.summary import RunSummary

if TYPE_CHECKING:
    import sqlite3

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.reporting.collections.sources import SourceView


def copy_batches(
    target: sqlite3.Connection, sources: tuple[SourceView, ...], budget: ReadBudget
) -> int:
    """Preserve recorded membership and feedback without claiming attempt replay."""
    for source in sources:
        portable = isinstance(source, BundleView)
        with reader(
            source.path, **({"filename": BUNDLE_DATABASE} if portable else {})
        ) as connection:
            if portable:
                count = source.summary.manifest.get("batches", 0)
                query, parameters = "SELECT payload,digest FROM batches", ()
            else:
                row = connection.execute(
                    "SELECT payload,digest FROM commits WHERE revision=?",
                    (source.revision,),
                ).fetchone()
                summary = RunSummary.from_manifest(checked_payload(row))
                if summary.run_id != source.run_id:
                    msg = "runtime batch source identity changed"
                    raise ValueError(msg)
                count = summary.batch_count or 0
                query, parameters = (
                    (
                        "SELECT payload,digest FROM batches "
                        "WHERE revision<=? ORDER BY ordinal"
                    ),
                    (source.revision,),
                )
            if not count:
                continue
            observed = 0
            for row in connection.execute(query, parameters):
                budget.examine(row[0])
                payload = checked_payload(row)
                decision = BatchDecision.from_dict(payload)
                if not portable and decision.run_id != source.run_id:
                    msg = "batch decision does not belong to its source"
                    raise ValueError(msg)
                observed += 1
                key = (decision.run_id, decision.cell_id, decision.index)
                existing = target.execute(
                    "SELECT payload,digest FROM batches "
                    "WHERE run_id=? AND cell_id=? AND batch_index=?",
                    key,
                ).fetchone()
                if existing is not None:
                    if checked_payload(existing) != payload:
                        msg = "conflicting runtime batch evidence"
                        raise ValueError(msg)
                    continue
                target.execute(
                    "INSERT INTO batches VALUES (?,?,?,?,?,?)",
                    (
                        *key,
                        decision.batch.batch_id,
                        canonical_json(payload),
                        semantic_digest(payload),
                    ),
                )
            if observed != count:
                msg = "runtime batch source count mismatch"
                raise ValueError(msg)
    return target.execute("SELECT count(*) FROM batches").fetchone()[0]
