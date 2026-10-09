"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/bundles/views.py

Lazy bundle record descriptors with the common lifetime and cursor contract.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, semantic_digest
from dense_arrays.artifacts.cursors import Cursor
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits
from dense_arrays.reporting.design_filters import DesignFilter
from dense_arrays.reporting.plans.bundles import selected_plans
from dense_arrays.reporting.plans.filters import PlanFilter
from dense_arrays.reporting.readers import Records

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.bundles.models import BundleSummary


@dataclass(frozen=True, repr=False)
class BundleView:
    """Read one immutable selected collection without opening its original sources."""

    path: Path
    summary: BundleSummary
    view: str
    limit: int | None = 100
    select: DesignFilter | PlanFilter | None = None
    read_limits: ReadLimits = field(default_factory=ReadLimits)
    after: Cursor | None = None

    def __post_init__(self) -> None:
        """Reject unsupported projections and mismatched continuation queries."""
        if self.view not in {"designs", "sequences", "placements", "plans", "batches"}:
            msg = "unsupported bundle record view"
            raise ValueError(msg)
        if self.limit is not None:
            integer(self.limit, field_name="limit", minimum=1)
        if self.view == "batches" and self.select is not None:
            msg = "batch record views do not accept design or plan filters"
            raise ValueError(msg)
        expected = PlanFilter if self.view == "plans" else DesignFilter
        if self.select is not None and not isinstance(self.select, expected):
            msg = f"bundle {self.view} queries require {expected.__name__}"
            raise TypeError(msg)
        if self.after is not None and self.after.advance(0) != self.cursor:
            msg = "cursor does not match the bundle identity or query"
            raise ValueError(msg)
        if self.after is not None and self.after.offset and self.view != "placements":
            msg = "cursor offset is only valid for placement rows"
            raise ValueError(msg)

    @property
    def revision(self) -> int:
        """Immutable bundles use revision zero across all record projections."""
        return 0

    @property
    def record_schema(self) -> str:
        """Preserve the same record encodings as native run queries."""
        return {
            "designs": "dense_arrays.design_record.v1",
            "sequences": "dense_arrays.sequence_record.v1",
            "placements": "dense_arrays.placement_record.v1",
            "plans": "dense_arrays.plan_evidence.v1",
            "batches": "dense_arrays.batch_decision.v1",
        }[self.view]

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Identify this immutable collection and its selected source context."""
        return (
            {
                "bundle_id": self.summary.bundle_id,
                "revision": 0,
                "filter": None if self.select is None else self.select.to_dict(),
                "scope": "all_matching_records",
                "record_schema": self.record_schema,
            },
        )

    @property
    def cursor(self) -> Cursor:
        """Bind a continuation to the manifest, projection and predicate."""
        return Cursor(
            self.summary.bundle_id,
            0,
            semantic_digest(
                {
                    "schema": "dense_arrays.bundle_query.v1",
                    "view": self.view,
                    "filter": None if self.select is None else self.select.to_dict(),
                    "order": "bundle_ordinal.v1",
                }
            ),
        )

    @property
    def cost(self) -> ReadCost:
        """Include plan annotations and filter resolution in the scan bound."""
        if self.view == "batches":
            return ReadCost(
                self.summary.bundle_id,
                0,
                "indexed",
                self.view,
                self.summary.manifest.get("batches", 0),
                self.read_limits,
            )
        if self.view == "plans":
            return ReadCost(
                self.summary.bundle_id,
                0,
                "scan",
                self.view,
                len(selected_plans(self.summary, self.select)),
                self.read_limits,
            )
        return ReadCost(
            self.summary.bundle_id,
            0,
            "scan",
            self.view,
            self.summary.designs + 2 * len(self.summary.manifest["plans"]),
            self.read_limits,
        )

    def records(self) -> Records:
        """Open an independent closeable reader only when iteration begins."""
        from .reading import read_bundle  # noqa: PLC0415

        budget = ReadBudget(self.read_limits)
        return Records(read_bundle(self, budget), budget, self.cursor, self.limit)

    def __repr__(self) -> str:
        """Describe the bound query without exposing manifests or reading records."""
        return (
            f"BundleView({self.view}, bundle={self.summary.bundle_id[:12]}, "
            f"limit={self.limit}, mode=scan)"
        )
