"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/selections/views.py

Lazy projections of saved membership with source and record integrity checks.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from contextlib import closing
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, semantic_digest
from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE
from dense_arrays.artifacts.bundles.storage import read_evidence
from dense_arrays.artifacts.cursors import Cursor
from dense_arrays.artifacts.reading import ReadBudget, ReadCost
from dense_arrays.artifacts.run_plans import cell_plans
from dense_arrays.artifacts.store import reader, stored_plan
from dense_arrays.reporting.bundles import BundleView
from dense_arrays.reporting.bundles.reading import read_bundle
from dense_arrays.reporting.collections import LibraryView
from dense_arrays.reporting.collections.reading import read_library
from dense_arrays.reporting.projections import project
from dense_arrays.reporting.readers import Records, RecordView, read_records

from .materialization import source_binding
from .membership import Membership
from .snapshots import SelectionSnapshot

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadLimits
    from dense_arrays.artifacts.records import Design
    from dense_arrays.parts.models import Part
    from dense_arrays.reporting.projections import PlacementRecord, SequenceRecord


@dataclass(frozen=True, repr=False)
class SelectionView:
    """A saved ordered selection over pinned evidence; no sampling on iteration."""

    query: RecordView | LibraryView | BundleView
    snapshot: SelectionSnapshot
    view: str = "designs"
    limit: int | None = 100
    after: Cursor | None = None

    def __post_init__(self) -> None:
        """Keep the complete evidence query separate from output pagination."""
        if (
            not isinstance(self.query, (RecordView, LibraryView, BundleView))
            or self.query.view != "designs"
            or self.query.limit is not None
            or self.query.select is not None
            or self.query.after is not None
        ):
            msg = "saved selections require a complete unfiltered design query"
            raise ValueError(msg)
        if not isinstance(self.snapshot, SelectionSnapshot):
            msg = "saved selection views require SelectionSnapshot"
            raise TypeError(msg)
        if self.view not in {"designs", "sequences", "placements"}:
            msg = "saved selections support designs, sequences and placements"
            raise ValueError(msg)
        if self.limit is not None:
            integer(self.limit, field_name="limit", minimum=1)
        if self.after is not None and (
            not isinstance(self.after, Cursor)
            or self.after.offset
            or self.after.advance(0) != self.cursor
        ):
            msg = "cursor does not match the saved selection and projection"
            raise ValueError(msg)

    @property
    def inputs(self) -> tuple[RecordView | BundleView, ...]:
        """Expose bound native evidence for report and portable bundle writers."""
        return (
            self.query.inputs if isinstance(self.query, LibraryView) else (self.query,)
        )

    @property
    def read_limits(self) -> ReadLimits:
        """Reuse one work limit declaration for evidence and saved membership."""
        return self.query.read_limits

    @property
    def revision(self) -> int:
        """Selection revision is zero; native revisions are in the source bindings."""
        return 0

    @property
    def record_schema(self) -> str:
        """Keep scalar joins identical to projections of native evidence."""
        return replace(self.query, view=self.view).record_schema

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Expose identities and scope without leaking local source locations."""
        return tuple(
            {
                f"{s.kind}_id": s.source_id,
                "revision": s.revision,
                "manifest_digest": s.manifest_digest,
                "selection_id": self.snapshot.snapshot_id,
                "scope": "saved_selection",
                "record_schema": self.record_schema,
            }
            for s in self.snapshot.sources
        )

    @property
    def cursor(self) -> Cursor:
        """Bind saved membership and output projection independently of pagination."""
        return Cursor(
            self.snapshot.snapshot_id,
            0,
            semantic_digest(
                {
                    "schema": "dense_arrays.selection_query.v1",
                    "view": self.view,
                    "record_schema": self.record_schema,
                }
            ),
        )

    @property
    def cost(self) -> ReadCost:
        """Disclose source replay and conservative annotation-plan reads."""
        estimate = self.query.cost.records_estimate
        if estimate is not None:
            estimate += len(self.inputs)
            if self.view == "placements":
                estimate += (
                    self.snapshot.selected
                    if any(isinstance(source, BundleView) for source in self.inputs)
                    else len(self.inputs)
                )
        return ReadCost(
            self.snapshot.snapshot_id, 0, "scan", self.view, estimate, self.read_limits
        )

    def records(self) -> Records:
        """Open a fresh, bounded read without drawing a new sample."""
        budget = ReadBudget(self.read_limits)
        return Records(read_selection(self, budget), budget, self.cursor, self.limit)

    def __repr__(self) -> str:
        """Describe output scope without scanning evidence or expanding membership."""
        return (
            f"SelectionView({self.view}, selected={self.snapshot.selected}, "
            f"limit={self.limit})"
        )


def check_sources(query: SelectionView, budget: ReadBudget) -> None:
    """Require the original committed manifests, including for empty selections."""
    sources = tuple(source_binding(s, budget) for s in query.inputs)
    if tuple(s.content() for s in sources) != tuple(
        s.content() for s in query.snapshot.sources
    ):
        msg = "saved selection source revision or committed content changed"
        raise ValueError(msg)


def _annotations(
    design: Design, query: SelectionView, budget: ReadBudget
) -> tuple[dict[str, Part], str]:
    """Read one bound plan at a time to recover selected placement annotations."""
    retained = budget.identities
    for source in query.inputs:
        if isinstance(source, BundleView):
            if design.plan_id not in source.summary.manifest["plans"]:
                continue
            with reader(source.path, filename=BUNDLE_DATABASE) as connection:
                plan = read_evidence(connection, design.plan_id, budget)
            break
        if source.run_id == design.run_id:
            budget.examine()
            with reader(source.path) as connection:
                plan = stored_plan(
                    connection,
                    max_identities=budget.limits.identities - budget.identities,
                )
            plan = cell_plans(plan)[design.cell_id]
            break
    else:
        msg = "saved selection has no source for placement annotations"
        raise ValueError(msg)
    if plan.plan_id != design.plan_id:
        msg = "saved design plan disagrees with its source"
        raise ValueError(msg)
    budget.identities = retained
    budget.retain(len(plan.request.parts))
    return {p.part_id: p for p in plan.request.parts}, plan.collection_id


def read_selection(
    query: SelectionView, budget: ReadBudget
) -> Iterator[Design | SequenceRecord | PlacementRecord]:
    """Validate saved content and order before emitting selected projections."""
    check_sources(query, budget)
    membership = Membership(query.snapshot, budget)
    raw = query.query
    iterator = (
        read_library(raw, budget)
        if isinstance(raw, LibraryView)
        else read_bundle(raw, budget)
        if isinstance(raw, BundleView)
        else read_records(raw, budget)
    )
    position, returned = 0, 0
    parts, collection, plan_id = {}, "", None
    with closing(iterator):
        for design in iterator:
            if not membership.matches(design):
                continue
            if query.view == "placements" and design.plan_id != plan_id:
                budget.identities -= len(parts)
                parts, collection = _annotations(design, query, budget)
                plan_id = design.plan_id
            for record in project(design, query.view, parts, collection):
                position += 1
                if query.after and position <= query.after.ordinal:
                    continue
                budget.position = position
                yield record
                returned += 1
                if query.limit is not None and returned >= query.limit:
                    return
    membership.finish()
