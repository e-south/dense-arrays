"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/pools/reading.py

Bounded candidate reads in stable mining order.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.pool_records import PoolSummary
from dense_arrays.artifacts.preparation.candidates import PoolCandidate
from dense_arrays.artifacts.store import checked_payload, reader
from dense_arrays.parts.candidates import Candidate

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.reporting.readers import RecordView


def iter_candidates(view: RecordView, budget: ReadBudget) -> Iterator[PoolCandidate]:
    """Stream validated records; a page does not claim full-pool verification."""
    with reader(view.path, filename="pool.sqlite3") as connection:
        summary = PoolSummary.from_dict(
            checked_payload(
                connection.execute(
                    "SELECT payload,digest FROM manifest WHERE id=1"
                ).fetchone()
            )
        )
        if summary.pool_id != view.pool_id or summary.preparation is None:
            msg = "candidate source changed or has no saved candidate evidence"
            raise ValueError(msg)
        if view.select is not None:
            view.select.validate(summary)
        previous = view.start_ordinal
        returned = 0
        for ordinal, payload, checksum in connection.execute(
            "SELECT ordinal,payload,digest FROM candidates "
            "WHERE ordinal>? ORDER BY ordinal",
            (previous,),
        ):
            budget.examine(payload)
            candidate = Candidate.from_dict(checked_payload((payload, checksum)))
            if (
                ordinal != previous + 1
                or ordinal != candidate.index
                or ordinal > summary.source_parts
            ):
                msg = "candidate ordinal or count disagrees with the pool"
                raise ValueError(msg)
            previous = ordinal
            record = PoolCandidate(summary.pool_id, candidate)
            if view.select is not None and not view.select.matches(record):
                continue
            budget.position = ordinal
            yield record
            returned += 1
            if view.limit is not None and returned >= view.limit:
                return
        if previous < summary.source_parts:
            msg = "candidate count disagrees with the pool"
            raise ValueError(msg)
