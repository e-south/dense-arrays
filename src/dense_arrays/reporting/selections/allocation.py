"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/selections/allocation.py

Resolve exact quotas and retain a bounded sample in source order.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import heapq
from typing import TYPE_CHECKING

from dense_arrays._record_validation import semantic_digest

from .snapshots import SelectedDesign

if TYPE_CHECKING:
    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Design

    from .requests import Take


class Allocation:
    """Retain only requested membership; inspect all candidates for exact counts."""

    def __init__(
        self, take: Take | None, resolved: dict[str, int | None], budget: ReadBudget
    ) -> None:
        """Charge allocation counters before creating independent reservoirs."""
        budget.retain(len(resolved))
        self.take, self.quotas, self.budget = take, resolved, budget
        self.available = dict.fromkeys(resolved, 0)
        self.samples = {key: [] for key in resolved}

    def observe(self, design: Design, ordinal: int) -> None:
        """Use SHA-256 priorities, not global or runtime-version-dependent RNG state."""
        per_cell = self.take is not None and self.take.per_cell is not None
        key = f"{design.run_id}/{design.cell_id}" if per_cell else "total"
        if key not in self.quotas:
            return
        self.available[key] += 1
        count = self.quotas[key]
        if count == 0:
            return
        sample = self.samples[key]
        random = self.take is not None and self.take.policy == "random"
        if not random and count is not None and len(sample) >= count:
            return
        priority = (
            int(
                semantic_digest(
                    {
                        "schema": "dense_arrays.selection_priority.v1",
                        "seed": self.take.seed,
                        "reference": design.reference,
                    }
                ),
                16,
            )
            if random
            else ordinal
        )
        if (
            count is not None
            and len(sample) >= count
            and (-priority, -ordinal) <= sample[0][:2]
        ):
            return
        member = SelectedDesign(design.reference, semantic_digest(design.to_dict()))
        entry = (-priority, -ordinal, member)
        if count is None or len(sample) < count:
            self.budget.retain()
            heapq.heappush(sample, entry)
        else:
            heapq.heapreplace(sample, entry)

    def finish(self) -> tuple[tuple[SelectedDesign, ...], dict[str, dict[str, int]]]:
        """Publish chosen records in encounter order, independently of priorities."""
        members = tuple(
            entry[2]
            for entry in sorted(
                (entry for sample in self.samples.values() for entry in sample),
                key=lambda e: -e[1],
            )
        )
        counts = {}
        for key, count in self.quotas.items():
            available, selected = self.available[key], len(self.samples[key])
            requested = available if count is None else count
            counts[key] = {
                "requested": requested,
                "available": available,
                "selected": selected,
                "shortfall": requested - selected,
            }
        return members, counts
