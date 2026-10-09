"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/pools/quality.py

Sampled pool quality from verified committed candidate evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays.artifacts.pools import pool_summary
from dense_arrays.artifacts.preparation.reading import read_verified_sampled
from dense_arrays.artifacts.reading import ReadCost, ReadLimits
from dense_arrays.reporting.pools.diversity import summarize_diversity

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.pool_records import PoolSummary


@dataclass(frozen=True, repr=False)
class PoolQualityReport:
    """An exact sampled-pool report with visible verification work limits."""

    path: Path
    summary: PoolSummary
    read_limits: ReadLimits = field(default_factory=ReadLimits)

    @property
    def cost(self) -> ReadCost:
        """Describe the full evidence scan needed to recount the quality report."""
        return ReadCost(
            self.summary.pool_id,
            0,
            "scan",
            "quality",
            1 + self.summary.source_parts + self.summary.retained_parts,
            self.read_limits,
        )

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Bind the immutable source pool independently of its file location."""
        return (
            {
                "pool_id": self.summary.pool_id,
                "plan_id": self.summary.plan_id,
                "revision": 0,
            },
        )

    def to_dict(self) -> dict[str, object]:
        """Verify and recount saved decisions without invoking preparation."""
        if pool_summary(self.path) != self.summary:
            msg = "pool summary changed since binding the quality report"
            raise ValueError(msg)
        evidence = read_verified_sampled(self.path, self.summary, self.read_limits)
        diversity = summarize_diversity(evidence)
        return {
            **self.summary.preparation.to_dict(),
            "schema": "dense_arrays.pool_quality.v1",
            "pool_id": self.summary.pool_id,
            "plan_id": self.summary.plan_id,
            "state": self.summary.state,
            **({"diversity": diversity} if diversity else {}),
        }

    def __repr__(self) -> str:
        """Avoid reading artifacts while displaying a notebook value."""
        return (
            f"PoolQualityReport({self.summary.pool_id[:12]}, "
            f"records={self.cost.records_estimate})"
        )
