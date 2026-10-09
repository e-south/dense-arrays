"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/collections/views.py

Bounded multi-source descriptors using the shared record iterator contract.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field, replace

from dense_arrays._record_validation import integer, semantic_digest
from dense_arrays.artifacts.cursors import Cursor
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits
from dense_arrays.reporting.bundles.views import BundleView
from dense_arrays.reporting.design_filters import DesignFilter
from dense_arrays.reporting.readers import Records, RecordView

from .sources import SourceView, annotation_cost, design_count


@dataclass(frozen=True, repr=False)
class LibraryView:
    """An ordered union, deduplicated by full design identity before projection."""

    inputs: tuple[SourceView, ...]
    view: str
    limit: int | None = 100
    select: DesignFilter | None = None
    read_limits: ReadLimits = field(default_factory=ReadLimits)
    after: Cursor | None = None

    def __post_init__(self) -> None:
        """Require explicit artifact snapshots and a matching continuation token."""
        if (
            not isinstance(self.inputs, (list, tuple))
            or not self.inputs
            or any(
                not isinstance(r, (RecordView, BundleView))
                or r.view != "designs"
                or r.select is not None
                or r.limit is not None
                or r.after is not None
                for r in self.inputs
            )
        ):
            msg = "library views require complete design snapshots"
            raise ValueError(msg)
        object.__setattr__(self, "inputs", tuple(self.inputs))
        if self.view not in {"designs", "sequences", "placements"}:
            msg = "combined sources support designs, sequences and placements"
            raise ValueError(msg)
        if self.limit is not None:
            integer(self.limit, field_name="limit", minimum=1)
        if not isinstance(self.read_limits, ReadLimits):
            msg = "library views require ReadLimits"
            raise TypeError(msg)
        if self.select is not None and not isinstance(self.select, DesignFilter):
            msg = "library views require DesignFilter"
            raise TypeError(msg)
        budget = ReadBudget(self.read_limits)
        budget.retain(len(self.inputs) + (self.select.identities if self.select else 0))
        if self.after is not None and (
            not isinstance(self.after, Cursor)
            or self.after.offset
            or self.after.advance(0) != self.cursor
        ):
            msg = "cursor does not match this library query or source order"
            raise ValueError(msg)
        if self.limit is not None:
            # Reject an unrepresentable token before streaming a page.
            self.cursor.token()

    @property
    def revision(self) -> int:
        """Union descriptors use revision zero; each run has its own pinned revision."""
        return 0

    @property
    def record_schema(self) -> str:
        """Use the same wire records for native and combined queries."""
        return replace(self.inputs[0], view=self.view).record_schema

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Expose every source snapshot in the supplied order, including repeats."""
        return tuple(
            replace(r, view=self.view, select=self.select).sources[0]
            for r in self.inputs
        )

    @property
    def cursor(self) -> Cursor:
        """Bind the entire ordered source set, predicate and projection."""
        bindings = tuple((r.cursor.source_id, r.cursor.revision) for r in self.inputs)
        return Cursor(
            semantic_digest(
                {"schema": "dense_arrays.library_sources.v1", "runs": bindings}
            ),
            self.revision,
            semantic_digest(
                {
                    "schema": "dense_arrays.library_query.v1",
                    "view": self.view,
                    "record_schema": self.record_schema,
                    "filter": None if self.select is None else self.select.to_dict(),
                    "order": "source_then_native_ordinal.v1",
                    "duplicates": "full_design_record.v1",
                }
            ),
            bindings=bindings,
        )

    @property
    def cost(self) -> ReadCost:
        """Include replayed prefixes; union pagination reconstructs its identity set."""
        annotations = self.view == "placements" or bool(
            self.select and (self.select.part_ids or self.select.groups)
        )
        filter_plans = bool(
            self.select and (self.select.part_ids or self.select.groups)
        )
        count = (
            None
            if any(design_count(r) is None for r in self.inputs)
            else sum(
                design_count(r)
                + annotation_cost(r) * (int(annotations) + int(filter_plans))
                for r in self.inputs
            )
        )
        return ReadCost(
            self.cursor.source_id,
            self.revision,
            "scan",
            self.view,
            count,
            self.read_limits,
        )

    def records(self) -> Records:
        """Open a fresh, closeable union reader with one shared work budget."""
        from .reading import read_library  # noqa: PLC0415

        budget = ReadBudget(self.read_limits)
        return Records(read_library(self, budget), budget, self.cursor, self.limit)

    def __repr__(self) -> str:
        """Describe a union without opening any source files."""
        return (
            f"LibraryView({self.view}, sources={len(self.inputs)}, "
            f"limit={self.limit}, mode=scan)"
        )
