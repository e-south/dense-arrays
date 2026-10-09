"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/plans/bundles.py

Read included plan inventories and explicitly chosen single-plan evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE
from dense_arrays.artifacts.bundles.storage import read_evidence
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.artifacts.store import reader

from .filters import PlanFilter

if TYPE_CHECKING:
    import sqlite3
    from collections.abc import Iterator
    from pathlib import Path

    from dense_arrays.artifacts.bundles.models import BundleSummary
    from dense_arrays.artifacts.reading import ReadLimits
    from dense_arrays.planning import PlanEvidence
    from dense_arrays.reporting.bundles import BundleView


def selected_plans(
    summary: BundleSummary, selected: PlanFilter | None
) -> tuple[str, ...]:
    """Keep inventory order and reject every unknown requested identity."""
    if selected is not None and not isinstance(selected, PlanFilter):
        msg = "plan views require PlanFilter"
        raise TypeError(msg)
    inventory = tuple(summary.manifest["plans"])
    if selected is None or not selected.plan_ids:
        return inventory
    requested = set(selected.plan_ids)
    if unknown := requested - set(inventory):
        msg = f"unknown plan IDs: {sorted(unknown)}"
        raise ValueError(msg)
    return tuple(identity for identity in inventory if identity in requested)


def read_bundle_plan(
    path: Path,
    summary: BundleSummary,
    limits: ReadLimits,
    selected: PlanFilter | None = None,
) -> PlanEvidence:
    """Require an unambiguous plan instead of silently choosing an inventory entry."""
    identities = selected_plans(summary, selected)
    if not identities:
        msg = "bundle contains no included plans"
        raise ValueError(msg)
    if len(identities) != 1:
        msg = (
            "bundle contains multiple plans; "
            "inspect view='plans' and select one plan ID"
        )
        raise ValueError(msg)
    budget = ReadBudget(limits)
    budget.retain(len(identities) + (selected.identities if selected else 0))
    with reader(path, filename=BUNDLE_DATABASE) as connection:
        return read_evidence(connection, identities[0], budget)


def read_bundle_plans(
    connection: sqlite3.Connection, query: BundleView, budget: ReadBudget
) -> Iterator[PlanEvidence]:
    """Page immutable evidence while releasing temporary per-plan identity state."""
    identities = selected_plans(query.summary, query.select)
    budget.retain(len(identities) + (query.select.identities if query.select else 0))
    start = query.after.ordinal if query.after else 0
    if start > len(identities):
        msg = "plan cursor is beyond the included inventory"
        raise ValueError(msg)
    returned = 0
    for ordinal, identity in enumerate(identities, 1):
        if ordinal <= start:
            continue
        before = budget.identities
        plan = read_evidence(connection, identity, budget)
        retained = budget.identities - before
        budget.position = ordinal
        try:
            yield plan
        finally:
            budget.identities -= retained
        del plan
        returned += 1
        if query.limit is not None and returned >= query.limit:
            return
