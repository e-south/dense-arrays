"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/bundles/reading.py

Stream selected bundle evidence through canonical filters and projections.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE
from dense_arrays.artifacts.bundles.storage import read_evidence, read_summary
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.records import Design
from dense_arrays.artifacts.run_state import manifest_cell_ids
from dense_arrays.artifacts.store import checked_payload, reader
from dense_arrays.reporting.bundles.batches import read_batch_page
from dense_arrays.reporting.collections.filtering import _unambiguous
from dense_arrays.reporting.design_filters import DesignFilter
from dense_arrays.reporting.design_queries import bundle_design_aliases
from dense_arrays.reporting.plans.bundles import read_bundle_plans
from dense_arrays.reporting.projections import project
from dense_arrays.reporting.summary import RunSummary

if TYPE_CHECKING:
    import sqlite3
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.reporting.bundles.views import BundleView


def load_plans(
    connection: sqlite3.Connection, query: BundleView, budget: ReadBudget
) -> dict:
    """Load the bounded declared plan inventory for annotations and verification."""
    budget.retain(len(query.summary.manifest["plans"]))
    plans = {
        identity: read_evidence(connection, identity, budget)
        for identity in query.summary.manifest["plans"]
    }
    for origin in query.summary.manifest["source_runs"]:
        summary = RunSummary.from_manifest(origin)
        for cell in summary.cell_summaries.values():
            if (
                cell.plan_id not in plans
                or plans[cell.plan_id].request.target.count != cell.target
            ):
                msg = "bundle cell plan identity or target does not match its origin"
                raise ValueError(msg)
    return plans


def read_bundle(query: BundleView, budget: ReadBudget) -> Iterator:
    """Read one immutable bundle, allowing continuation within a design's placements."""
    current = read_summary(query.path, query.read_limits)
    if current.bundle_id != query.summary.bundle_id:
        msg = "bundle changed since the query was bound"
        raise ValueError(msg)
    with reader(query.path, filename=BUNDLE_DATABASE) as connection:
        if query.view == "batches":
            yield from read_batch_page(connection, query, budget)
            return
        if query.view == "plans":
            yield from read_bundle_plans(connection, query, budget)
            return
        needs_parts = query.view == "placements" or bool(
            query.select and (query.select.part_ids or query.select.groups)
        )
        plans = load_plans(connection, query, budget) if needs_parts else {}
        annotations = {
            identity: ({p.part_id: p for p in plan.request.parts}, plan.collection_id)
            for identity, plan in plans.items()
        }
        selected = resolve_filter(connection, query, plans, budget)
        start = query.after.ordinal if query.after else 0
        offset = query.after.offset if query.after else 0
        emitted = scanned = 0
        for row in connection.execute(
            "SELECT ordinal,design_ref,run_id,cell_id,local_id,plan_id,payload,digest "
            "FROM designs WHERE ordinal>? ORDER BY ordinal",
            (start - bool(offset),),
        ):
            ordinal, design = decode_design(row, budget)
            scanned += 1
            parts, collection = annotations.get(design.plan_id, ({}, ""))
            if selected is not None and not selected.matches(design, parts, collection):
                continue
            rows = project(design, query.view, parts, collection)
            for index, value in enumerate(rows):
                if ordinal == start and index < offset:
                    continue
                budget.position = ordinal
                budget.offset = (
                    index + 1
                    if query.view == "placements"
                    and index + 1 < len(design.realized.placements)
                    else 0
                )
                emitted += 1
                yield value
                if query.limit is not None and emitted >= query.limit:
                    return
        if not start and scanned != query.summary.designs:
            msg = "contained designs do not reconcile with the bundle manifest"
            raise ArtifactIntegrityError(msg, artifact=query.path)


def decode_design(row: tuple, budget: ReadBudget) -> tuple[int, Design]:
    """Check the data and its SQL identity columns before using either."""
    ordinal, ref, run_id, cell_id, local_id, plan_id, payload, checksum = row
    budget.examine(payload)
    design = Design.from_dict(checked_payload((payload, checksum)))
    if (ref, run_id, cell_id, local_id, plan_id) != (
        design.reference,
        design.run_id,
        design.cell_id,
        design.design_id,
        design.plan_id,
    ):
        msg = "bundle design identity disagrees with its stored index"
        raise ValueError(msg)
    return ordinal, design


def resolve_filter(
    connection: sqlite3.Connection, query: BundleView, plans: dict, budget: ReadBudget
) -> DesignFilter | None:
    """Reject missing or ambiguous aliases using the same collection semantics."""
    selected = query.select
    if selected is None:
        return None
    offered = {
        name: {v: set() for v in getattr(selected, name)}
        for name in ("design_ids", "cells", "part_ids", "groups")
    }
    budget.retain(selected.identities)

    def offer(name: str, local: str, full: str) -> None:
        for label in {local, full} & offered[name].keys():
            if full not in offered[name][label]:
                budget.retain()
                offered[name][label].add(full)

    for run in query.summary.manifest["source_runs"]:
        for cell in manifest_cell_ids(run):
            offer("cells", cell, f"{run['run_id']}/{cell}")
    if selected.design_ids:
        for local, full in bundle_design_aliases(connection, selected.design_ids):
            offer("design_ids", local, full)
    for name, local, full in _part_aliases(plans):
        offer(name, local, full)
    resolved = _unambiguous(offered)
    budget.identities -= selected.identities + sum(
        len(matches) for labels in offered.values() for matches in labels.values()
    )
    result = DesignFilter(**resolved, metrics=selected.metrics)
    budget.retain(result.identities)
    return result


def _part_aliases(plans: dict) -> Iterator[tuple[str, str, str]]:
    """Use one collection digest per bound plan when resolving part predicates."""
    for plan in plans.values():
        collection = plan.collection_id
        for part in plan.request.parts:
            yield "part_ids", part.part_id, f"{collection}/{part.part_id}"
            if part.group is not None:
                yield "groups", part.group, part.group
