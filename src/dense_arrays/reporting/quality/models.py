"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/models.py

Lazy quality reports bind the same design scope as inspection and export.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    integer,
    mutable_json,
    semantic_digest,
)
from dense_arrays.artifacts.bundles.models import BundleSummary
from dense_arrays.artifacts.cursors import Cursor
from dense_arrays.artifacts.reading import ReadCost
from dense_arrays.artifacts.records import COMPOSITION_POLICY
from dense_arrays.reporting.bundles.views import BundleView
from dense_arrays.reporting.collections import LibraryView
from dense_arrays.reporting.collections.sources import (
    SourceView,
    annotation_cost,
    design_count,
)
from dense_arrays.reporting.readers import RecordView
from dense_arrays.reporting.selections.snapshots import SelectionSnapshot
from dense_arrays.reporting.selections.views import SelectionView

from .evidence import origins

if TYPE_CHECKING:
    from collections.abc import Mapping

    from dense_arrays.artifacts.reading import ReadLimits
    from dense_arrays.reporting.summary import RunSummary

QUALITY_POLICY = COMPOSITION_POLICY
QUALITY_SCHEMA = "dense_arrays.quality.v3"


@dataclass(frozen=True, repr=False)
class QualityReport:
    """Exact selected aggregates with source attainment and paginated usage tables."""

    query: RecordView | LibraryView | BundleView
    summaries: tuple[RunSummary | BundleSummary, ...]
    limit: int = 100
    after: Cursor | None = None
    snapshot: SelectionSnapshot | None = None

    def __post_init__(self) -> None:
        """Validate snapshot bindings and continuation before scanning evidence."""
        integer(self.limit, field_name="limit", minimum=1)
        if not isinstance(self.query, (RecordView, LibraryView, BundleView)) or (
            self.query.view != "designs"
            or self.query.limit is not None
            or self.query.after is not None
        ):
            msg = "quality requires a complete design query"
            raise ValueError(msg)
        object.__setattr__(self, "summaries", tuple(self.summaries))
        if self.snapshot is not None and (
            not isinstance(self.snapshot, SelectionSnapshot)
            or self.query.select is not None
        ):
            msg = "quality accepts one saved snapshot or a design filter"
            raise ValueError(msg)
        if len(self.summaries) != len(self.inputs):
            msg = "quality requires a summary for every source snapshot"
            raise ValueError(msg)
        _ = self.origins
        if self.after is not None and (
            not isinstance(self.after, Cursor)
            or self.after.offset
            or self.after.advance(0) != self.cursor
        ):
            msg = "cursor does not match this quality query or its snapshots"
            raise ValueError(msg)
        self.cursor.token()

    @property
    def inputs(self) -> tuple[SourceView, ...]:
        """Expose supplied artifact snapshots in their declared order."""
        return (
            self.query.inputs if isinstance(self.query, LibraryView) else (self.query,)
        )

    @cached_property
    def origins(self) -> Mapping[str, RunSummary]:
        """Original run summaries, deduplicated with one consistent revision each."""
        from types import MappingProxyType  # noqa: PLC0415

        return MappingProxyType(origins(self.inputs, self.summaries, self.read_limits))

    @property
    def read_limits(self) -> ReadLimits:
        """Use one set of work limits for source evidence and report state."""
        return self.query.read_limits

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Retain ordered source identities and the declared predicate."""
        return (
            SelectionView(self.query, self.snapshot).sources
            if self.snapshot
            else self.query.sources
        )

    @property
    def cursor(self) -> Cursor:
        """Bind metric policy and usage ordering to the complete design query."""
        base = self.query.cursor
        return Cursor(
            base.source_id,
            base.revision,
            semantic_digest(
                {
                    "schema": "dense_arrays.quality_query.v1",
                    "policy": QUALITY_POLICY,
                    "design_query": base.query_id,
                    **(
                        {"selection_id": self.snapshot.snapshot_id}
                        if self.snapshot
                        else {}
                    ),
                    "order": "occurrences_descending_source_order.v1",
                }
            ),
            bindings=base.bindings,
        )

    @property
    def cost(self) -> ReadCost:
        """Disclose all source plans, designs and checked search histories."""
        predicate_plans = (
            sum(annotation_cost(s) for s in self.inputs)
            if self.query.select
            and (self.query.select.part_ids or self.query.select.groups)
            else 0
        )
        plan_reads = sum(
            len(s.manifest["source_runs"]) if isinstance(s, BundleSummary) else 1
            for s in self.summaries
        )
        native_attempts = sum(
            s.counts["started"]
            for s in self.summaries
            if not isinstance(s, BundleSummary)
        )
        count = (
            None
            if any(design_count(s) is None for s in self.inputs)
            else plan_reads
            + predicate_plans
            + sum(design_count(s) for s in self.inputs)
            + native_attempts
            + (len(self.inputs) if self.snapshot else 0)
        )
        return ReadCost(
            self.cursor.source_id,
            self.cursor.revision,
            "scan",
            "quality",
            count,
            self.read_limits,
        )

    @cached_property
    def _data(self) -> Mapping[str, object]:
        from .collection import collect  # noqa: PLC0415

        return immutable_json_mapping(collect(self))

    def to_dict(self) -> dict[str, object]:
        """Page usage tables without changing aggregate populations or source status."""
        value = mutable_json(self._data)
        start = 0 if self.after is None else self.after.ordinal
        stop = start + self.limit
        total = len(value["part_usage"])
        if start > total:
            msg = "quality cursor is beyond the selected population"
            raise ValueError(msg)
        for metrics in [value, *value["cells"]]:
            for name in ("part_usage", "group_usage"):
                metrics[name] = metrics[name][start:stop]
        value["next_cursor"] = (
            self.cursor.advance(stop).token() if stop < total else None
        )
        value["usage_order"] = "occurrences_descending_source_order.v1"
        value["usage_offset"] = start
        value["cost"] = self.cost.to_dict()
        return value

    def __repr__(self) -> str:
        """Describe report scope without computing aggregates."""
        return f"QualityReport(sources={len(self.inputs)}, limit={self.limit})"
