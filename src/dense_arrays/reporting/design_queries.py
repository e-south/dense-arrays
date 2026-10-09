"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/design_queries.py

Bind design predicates and placement annotations to their stored plan.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json
from dense_arrays.artifacts.errors import InvalidQueryError
from dense_arrays.artifacts.run_plans import cell_plans
from dense_arrays.artifacts.store import stored_plan

if TYPE_CHECKING:
    import sqlite3
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.parts.models import Part
    from dense_arrays.reporting.readers import RecordView


def design_context(
    connection: sqlite3.Connection,
    view: RecordView,
    budget: ReadBudget,
    cell_ids: tuple[str, ...],
) -> dict[str, tuple[dict[str, Part], str]]:
    """Load annotations only when the query needs them; validate explicit IDs."""
    selected = view.select
    annotations = {}
    part_ids, groups = set(), set()
    if view.view == "placements" or (
        selected is not None and (selected.part_ids or selected.groups)
    ):
        budget.examine()
        plan = stored_plan(
            connection, max_identities=budget.limits.identities - budget.identities
        )
        for child in cell_plans(plan).values():
            budget.retain(len(child.request.parts))
            parts = {p.part_id: p for p in child.request.parts}
            annotations[child.plan_id] = (parts, child.collection_id)
            part_ids.update(parts)
            part_ids.update(f"{child.collection_id}/{p}" for p in parts)
            groups.update(p.group for p in parts.values())
    if selected is None:
        return annotations
    budget.retain(selected.identities)
    for name, labels, available in (
        (
            "cells",
            selected.cells,
            set(cell_ids) | {f"{view.run_id}/{c}" for c in cell_ids},
        ),
        (
            "part IDs",
            selected.part_ids,
            part_ids,
        ),
        ("groups", selected.groups, groups),
    ):
        if missing := set(labels) - available:
            msg = f"unknown {name}: {sorted(missing)}"
            raise InvalidQueryError(msg)
    if selected.design_ids:
        # Resolve all requested aliases in one native pass, not one pass per ID.
        missing = set(selected.design_ids)
        for local, full in design_aliases(
            connection, view.revision, selected.design_ids
        ):
            missing.discard(local)
            missing.discard(full)
        if missing:
            msg = f"unknown design IDs: {sorted(missing)}"
            raise InvalidQueryError(msg)
    return annotations


def design_aliases(
    connection: sqlite3.Connection, revision: int, requested: tuple[str, ...]
) -> Iterator[tuple[str, str]]:
    """Resolve requested design labels in one native SQL pass."""
    query = """WITH requested AS (SELECT value FROM json_each(?)),
            offered AS (
                SELECT json_extract(payload,'$.design_id') AS local,
                    json_extract(payload,'$.run_id') || '/' ||
                    json_extract(payload,'$.cell_id') || '/' ||
                    json_extract(payload,'$.design_id') AS full
                FROM designs WHERE revision<=?)
            SELECT local,full FROM offered WHERE
                local IN requested OR full IN requested"""
    yield from connection.execute(query, (canonical_json(list(requested)), revision))


def bundle_design_aliases(
    connection: sqlite3.Connection, requested: tuple[str, ...]
) -> Iterator[tuple[str, str]]:
    """Resolve only requested labels through the bundle's indexed identities."""
    yield from connection.execute(
        "WITH requested AS (SELECT value FROM json_each(?)) "
        "SELECT local_id,design_ref FROM designs WHERE local_id IN requested "
        "OR design_ref IN requested",
        (canonical_json(list(requested)),),
    )
