"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/collections/reading.py

Bounded full-record deduplication before filtering and scalar projection.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from contextlib import closing
from typing import TYPE_CHECKING

from dense_arrays._record_validation import semantic_digest
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.reporting.projections import project

from .filtering import resolve_filter
from .sources import annotations, read_designs

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Design
    from dense_arrays.reporting.projections import PlacementRecord, SequenceRecord

    from .views import LibraryView


def read_library(
    view: LibraryView, budget: ReadBudget
) -> Iterator[Design | SequenceRecord | PlacementRecord]:
    """Read one source at a time; replay prefixes to recover exact union state."""
    budget.retain(len(view.inputs))
    selected = resolve_filter(view.inputs, view.select, budget)
    seen: dict[str, str] = {}
    position = 0
    returned = 0
    start = 0 if view.after is None else view.after.ordinal
    needed = view.view == "placements" or bool(
        selected and (selected.part_ids or selected.groups)
    )
    for source in view.inputs:
        with (
            annotations(source, budget, needed=needed) as lookup,
            closing(read_designs(source, budget)) as records,
        ):
            for design in records:
                if needed and design.plan_id not in lookup:
                    msg = "design plan is not bound by its source artifact"
                    raise ArtifactIntegrityError(msg, artifact=source.path)
                parts, collection = lookup.get(design.plan_id, ({}, ""))
                fingerprint = semantic_digest(design.to_dict())
                previous = seen.get(design.reference)
                if previous is not None:
                    if previous != fingerprint:
                        msg = f"conflicting content for design {design.reference}"
                        raise ArtifactIntegrityError(msg, artifact=source.path)
                    continue
                budget.retain()
                seen[design.reference] = fingerprint
                if selected is not None and not selected.matches(
                    design, parts, collection
                ):
                    continue
                for record in project(design, view.view, parts, collection):
                    position += 1
                    if position <= start:
                        continue
                    budget.position, budget.offset = position, 0
                    yield record
                    returned += 1
                    if view.limit is not None and returned >= view.limit:
                        return
