"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/selections/materialization.py

Materialize selections from complete, revision-bound design queries.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from contextlib import closing
from typing import TYPE_CHECKING

from dense_arrays._record_validation import mutable_json, semantic_digest
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.artifacts.run_state import manifest_cell_ids
from dense_arrays.artifacts.store import checked_payload, reader
from dense_arrays.reporting.bundles import BundleView
from dense_arrays.reporting.bundles.reading import read_bundle
from dense_arrays.reporting.collections import LibraryView
from dense_arrays.reporting.collections.reading import read_library
from dense_arrays.reporting.readers import RecordView, read_records

from .allocation import Allocation
from .requests import resolve_quotas

if TYPE_CHECKING:
    from .requests import LibrarySelection
from .snapshots import SelectionSnapshot, SelectionSource


def source_binding(
    source: RecordView | BundleView, budget: ReadBudget
) -> SelectionSource:
    """Pin the canonical committed manifest as well as its revision number."""
    budget.retain()
    if isinstance(source, BundleView):
        budget.examine()
        manifest = mutable_json(source.summary.manifest)
        budget.retain(len(manifest["source_runs"]))
        return SelectionSource(
            "bundle",
            source.summary.bundle_id,
            0,
            semantic_digest(manifest),
            source.path,
            tuple(
                dict.fromkeys(
                    f"{s['run_id']}/{c}"
                    for s in manifest["source_runs"]
                    for c in manifest_cell_ids(s)
                )
            ),
        )
    with reader(source.path) as connection:
        row = connection.execute(
            "SELECT payload,digest FROM commits WHERE revision=?", (source.revision,)
        ).fetchone()
        budget.examine(None if row is None else row[0])
        manifest = checked_payload(row)
    if manifest["run_id"] != source.run_id or manifest["revision"] != source.revision:
        msg = "selection source identity changed"
        raise ValueError(msg)
    budget.retain()
    return SelectionSource(
        "run",
        source.run_id,
        source.revision,
        semantic_digest(manifest),
        source.path,
        tuple(f"{source.run_id}/{c}" for c in manifest_cell_ids(manifest)),
    )


def materialize(
    query: RecordView | LibraryView | BundleView, request: LibrarySelection
) -> SelectionSnapshot:
    """Select under one shared work budget, never using a display page as input."""
    if query.view != "designs" or query.limit is not None or query.after is not None:
        msg = "selection requires a complete design query without pagination"
        raise ValueError(msg)
    budget = ReadBudget(query.read_limits)
    raw_sources = query.inputs if isinstance(query, LibraryView) else (query,)
    sources = tuple(source_binding(s, budget) for s in raw_sources)
    allocation = Allocation(
        request.take,
        resolve_quotas(request.take, (c for s in sources for c in s.cells)),
        budget,
    )
    iterator = (
        read_library(query, budget)
        if isinstance(query, LibraryView)
        else read_bundle(query, budget)
        if isinstance(query, BundleView)
        else read_records(query, budget)
    )
    with closing(iterator):
        for ordinal, design in enumerate(iterator):
            allocation.observe(design, ordinal)
    members, counts = allocation.finish()
    return SelectionSnapshot(
        sources,
        request,
        members,
        counts,
        "sha256_priority.v1"
        if request.take and request.take.policy == "random"
        else "source_order.v1",
        query.read_limits,
    )
