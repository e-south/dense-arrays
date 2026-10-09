"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/snapshots.py

Immutable recorded metrics for comparison after source artifacts are unavailable.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    mutable_json,
    semantic_digest,
)
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits

if TYPE_CHECKING:
    from collections.abc import Mapping

from .models import QUALITY_POLICY, QUALITY_SCHEMA, QualityReport
from .validation import report_entries, validate_quality


@dataclass(frozen=True, repr=False)
class QualitySnapshot:
    """Recorded quality values; validates internal consistency, not source truth."""

    data: Mapping[str, object]
    read_limits: ReadLimits = field(
        default_factory=ReadLimits, repr=False, compare=False
    )
    snapshot_id: str = field(init=False)
    identities: int = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Bound and validate recorded aggregate evidence before freezing it."""
        if not isinstance(self.read_limits, ReadLimits):
            msg = "quality snapshots require ReadLimits"
            raise TypeError(msg)
        entries = report_entries(self.data)
        ReadBudget(self.read_limits).retain(entries)
        value = validate_quality(
            self.data, schema=QUALITY_SCHEMA, policy=QUALITY_POLICY
        )
        object.__setattr__(self, "data", immutable_json_mapping(value))
        object.__setattr__(self, "snapshot_id", semantic_digest(value))
        object.__setattr__(self, "identities", entries)

    @classmethod
    def from_report(cls, report: QualityReport) -> QualitySnapshot:
        """Capture one complete report's aggregates and declared usage page."""
        if not isinstance(report, QualityReport):
            msg = "from_report requires QualityReport"
            raise TypeError(msg)
        return cls(report.to_dict(), report.read_limits)

    @classmethod
    def from_dict(
        cls, value: object, *, read_limits: ReadLimits | None = None
    ) -> QualitySnapshot:
        """Read the native quality schema; unknown metric policies remain explicit."""
        return cls(value, read_limits or ReadLimits())

    def to_dict(self) -> dict[str, object]:
        """Preserve the original report's native schema and population declarations."""
        return mutable_json(self.data)

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Expose recorded origin bindings without requiring their locations."""
        return tuple(mutable_json(v) for v in self.data["selection"]["sources"])

    @property
    def cost(self) -> ReadCost:
        """Saved metrics require a single document read and no source scan."""
        return ReadCost(self.snapshot_id, 0, "manifest", "quality", 1, self.read_limits)

    def __repr__(self) -> str:
        """Keep notebook display bounded even for large stored usage tables."""
        return (
            f"QualitySnapshot({self.snapshot_id[:12]}, "
            f"designs={self.data['selection']['designs']})"
        )
