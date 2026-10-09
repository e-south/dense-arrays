"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/comparison.py

Compare saved-library quality with one shared read budget and explicit scope.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

from dense_arrays._record_validation import (
    immutable_json_mapping,
    mutable_json,
    semantic_digest,
)
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits

from .collection import collect
from .differences import MetricDifference, metric_differences
from .models import QUALITY_POLICY, QualityReport
from .snapshots import QualitySnapshot


@dataclass(frozen=True, repr=False)
class QualityComparison:
    """Descriptive aggregate differences; no causal or significance inference."""

    before: QualityReport | QualitySnapshot
    after: QualityReport | QualitySnapshot
    read_limits: ReadLimits = field(default_factory=ReadLimits)

    def __post_init__(self) -> None:
        """Keep each side's independently declared population and source revisions."""
        if not isinstance(
            self.before, (QualityReport, QualitySnapshot)
        ) or not isinstance(self.after, (QualityReport, QualitySnapshot)):
            msg = "quality comparison requires two reports or saved quality snapshots"
            raise TypeError(msg)
        if not isinstance(self.read_limits, ReadLimits):
            msg = "quality comparison requires ReadLimits"
            raise TypeError(msg)

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Retain ordered source bindings from both sides without merging targets."""
        return (*self.before.sources, *self.after.sources)

    @property
    def cost(self) -> ReadCost:
        """Estimate the two complete evidence scans independently of display pages."""
        a, b = self.before.cost, self.after.cost
        count = (
            None
            if a.records_estimate is None or b.records_estimate is None
            else a.records_estimate + b.records_estimate
        )
        return ReadCost(
            semantic_digest(
                {"before": _identity(self.before), "after": _identity(self.after)}
            ),
            0,
            "scan",
            "quality_comparison",
            count,
            self.read_limits,
        )

    @cached_property
    def _reports(self) -> tuple[dict, dict]:
        budget = ReadBudget(self.read_limits)
        reports = []
        for report in (self.before, self.after):
            value = _read_report(report, budget)
            reports.append(immutable_json_mapping(value))
        return tuple(reports)

    @cached_property
    def metrics(self) -> tuple[MetricDifference, ...]:
        """Expose typed after-minus-before observations with separate denominators."""
        return metric_differences(*self._reports, policy=QUALITY_POLICY)

    def to_dict(self) -> dict[str, object]:
        """Preserve population and source attainment alongside every metric delta."""
        before, after = self._reports
        scopes = tuple(
            {
                **_scope(value),
                "evidence": "saved_report"
                if isinstance(source, QualitySnapshot)
                else "artifact_records",
            }
            for value, source in zip(
                (before, after), (self.before, self.after), strict=True
            )
        )
        same = scopes[0]["selection"] == scopes[1]["selection"]
        different = before["selection"]["designs"] != after["selection"]["designs"]
        return {
            "schema": "dense_arrays.quality_comparison.v1",
            "mode": "descriptive",
            "metric_scope": "aggregate_metrics",
            "before": scopes[0],
            "after": scopes[1],
            "population_changed": False if same else True if different else None,
            "scope_changed": not same,
            "metrics": [metric.to_dict() for metric in self.metrics],
            "cost": self.cost.to_dict(),
            "examined": before["examined"] + after["examined"],
        }

    def __repr__(self) -> str:
        """Describe scope without scanning evidence or constructing differences."""
        return (
            f"QualityComparison(before={type(self.before).__name__}, "
            f"after={type(self.after).__name__})"
        )


def _identity(report: QualityReport | QualitySnapshot) -> dict[str, object]:
    if isinstance(report, QualitySnapshot):
        return {"snapshot_id": report.snapshot_id}
    return {
        "source_id": report.cursor.source_id,
        "revision": report.cursor.revision,
        "query_id": report.cursor.query_id,
    }


def _read_report(report: QualityReport | QualitySnapshot, budget: ReadBudget) -> dict:
    """Preserve per-report limits inside the comparison's aggregate work cap."""
    total, local = budget.limits, report.read_limits
    budget.limits = ReadLimits(
        records=min(total.records, budget.examined + local.records),
        identities=min(total.identities, budget.identities + local.identities),
        pairs=min(total.pairs, local.pairs),
    )
    try:
        if isinstance(report, QualitySnapshot):
            budget.examine()
            budget.retain(report.identities)
            return {**report.data, "examined": 1}
        return collect(report, budget=budget)
    finally:
        budget.limits = total


def _scope(value: dict) -> dict[str, object]:
    selection = value["selection"]
    return {
        "policy": value["policy"],
        "population": value["population"],
        "selection": mutable_json(selection),
        "selected_designs": selection["designs"],
        "distinct_sequences": selection["distinct_sequences"],
        "occurrences": value["concentration"]["occurrence_denominator"],
        "source_runs": mutable_json(value["source_runs"]),
        "search_availability": value["search"]["availability"],
    }
