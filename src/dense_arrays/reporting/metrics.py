"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/metrics.py

Compute composition and interval-union metrics from persisted final placements.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer

if TYPE_CHECKING:
    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.realized import RealizedArray

COMPOSITION_METRICS = (
    "length",
    "gc_fraction",
    "placement_count",
    "density",
    "compression",
    "packed_span",
    "padding_length",
)


def covered_intervals(realized: RealizedArray) -> tuple[tuple[int, int], ...]:
    """Merge overlapping/touching half-open placements, counting overlap once."""
    merged: list[tuple[int, int]] = []
    for placement in sorted(realized.placements, key=lambda p: (p.start, p.end)):
        if merged and placement.start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], placement.end))
        else:
            merged.append((placement.start, placement.end))
    return tuple(merged)


def design_metrics(realized: RealizedArray) -> dict[str, int | float]:
    """Separate final coverage density from packed-span compression and padding."""
    length = len(realized.sequence)
    span = max(p.end for p in realized.placements) - min(
        p.start for p in realized.placements
    )
    covered = sum(end - start for start, end in covered_intervals(realized))
    return {
        "length": length,
        "gc_fraction": (realized.sequence.count("G") + realized.sequence.count("C"))
        / length,
        "placement_count": len(realized.placements),
        "density": covered / length,
        "compression": sum(len(p.sequence) for p in realized.placements) / span,
        "packed_span": span,
        "padding_length": length - span,
    }


def validate_metric_value(name: str, value: float) -> None:
    """Validate a finite histogram value under the native composition definitions."""
    if name in {"length", "packed_span", "placement_count", "padding_length"}:
        integer(value, field_name=name, minimum=0 if name == "padding_length" else 1)
    elif name in {"gc_fraction", "density"} and not 0 <= value <= 1:
        msg = f"{name} must be between 0 and 1"
        raise ValueError(msg)
    elif name == "compression" and value <= 0:
        msg = "compression must be positive"
        raise ValueError(msg)


@dataclass
class Distribution:
    """An exact frequency table with bounded lookup state and an explicit population."""

    counts: Counter = field(default_factory=Counter)

    def observe(self, value: float, budget: ReadBudget) -> None:
        """Charge each new value before retaining its bin."""
        if value not in self.counts:
            budget.retain()
        self.counts[value] += 1

    def merge(self, other: Distribution, budget: ReadBudget) -> None:
        """Merge exact frequency tables without retaining individual observations."""
        for value, count in other.counts.items():
            if value not in self.counts:
                budget.retain()
            self.counts[value] += count

    def to_dict(self, *, denominator: str = "accepted_designs") -> dict[str, object]:
        """Report null values for an empty denominator, never an invented zero."""
        count = self.counts.total()
        return {
            "count": count,
            "status": "exact",
            "denominator": denominator,
            "min": min(self.counts) if count else None,
            "max": max(self.counts) if count else None,
            "mean": math.fsum(value * n for value, n in self.counts.items()) / count
            if count
            else None,
            "histogram": [
                {"value": value, "count": n} for value, n in sorted(self.counts.items())
            ],
            "reason": None if count else "empty_population",
        }
