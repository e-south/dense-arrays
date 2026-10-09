"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/accumulation.py

Accumulate exact composition and usage from persisted designs.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter
from itertools import pairwise
from typing import TYPE_CHECKING

from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.reporting.metrics import (
    COMPOSITION_METRICS,
    Distribution,
    covered_intervals,
    design_metrics,
)

if TYPE_CHECKING:
    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Design
    from dense_arrays.planning import GenerationPlan, PlanEvidence


class QualityAccumulator:
    """Accumulate exact single-cell metrics with bounded lookup state."""

    def __init__(
        self, plan: GenerationPlan | PlanEvidence | None, budget: ReadBudget
    ) -> None:
        self.plan, self.budget = plan, budget
        self.collection_id = None if plan is None else plan.collection_id
        parts = () if plan is None else plan.request.parts
        requirements = () if plan is None else plan.request.requirements
        budget.retain(len(parts) + len(requirements))
        self.parts = {p.part_id: p for p in parts}
        self.collections = {}
        self.groups = dict.fromkeys(p.group for p in parts if p.group is not None)
        budget.retain(len(self.groups))
        self.occurrences, self.containing = Counter(), Counter()
        self.group_occurrences, self.group_containing = Counter(), Counter()
        self.sequences: set[str] = set()
        self.composition: dict[str, Distribution] = {}
        self.padding_sides = Counter()
        self.coverage, self.lengths = Counter(), Counter()
        self.requirements = {
            r.id: {
                "id": r.id,
                "passed": 0,
                "failed": 0,
                "missing": 0,
                "not_applicable": 0,
            }
            for r in requirements
        }
        self.accepted = 0

    def observe(self, design: Design) -> None:
        """Read supplied identities and requirement evidence without screening."""
        self.accepted += 1
        if design.sequence_id not in self.sequences:
            self.budget.retain()
            self.sequences.add(design.sequence_id)
        chosen = [p.feature_id for p in design.realized.placements]
        if set(chosen) - self.parts.keys() or len(chosen) != len(set(chosen)):
            msg = "quality placements contain unknown or repeated part identities"
            raise ArtifactIntegrityError(msg)
        self.occurrences.update(chosen)
        self.containing.update(set(chosen))
        groups = [
            self.parts[part_id].group
            for part_id in chosen
            if self.parts[part_id].group is not None
        ]
        self.group_occurrences.update(groups)
        self.group_containing.update(set(groups))
        for name, value in design_metrics(design.realized).items():
            self.composition.setdefault(name, Distribution()).observe(
                value, self.budget
            )
        assembly = design.realized.provenance.get("assembly")
        side = (
            "none"
            if not assembly or not assembly["padding_length"]
            else assembly["side"]
        )
        self.padding_sides[side] += 1
        for start, end in covered_intervals(design.realized):
            self._event(self.coverage, start, 1)
            self._event(self.coverage, end, -1)
        self._event(self.lengths, len(design.realized.sequence), 1)
        self._requirements(design)

    def _event(self, counter: Counter, position: int, count: int) -> None:
        if position not in counter:
            self.budget.retain()
        counter[position] += count

    def _requirements(self, design: Design) -> None:
        seen = set()
        for result in design.requirements:
            identity = result["id"]
            if (
                identity not in self.requirements
                or identity in seen
                or not isinstance(result["passed"], bool)
            ):
                msg = "invalid or unknown quality requirement evidence"
                raise ArtifactIntegrityError(msg)
            seen.add(identity)
            count = self.requirements[identity]
            count["passed" if result["passed"] else "failed"] += 1
            count["not_applicable"] += result.get("status") == "not_applicable"
        for identity in self.requirements.keys() - seen:
            self.requirements[identity]["missing"] += 1

    def usage(self, *, groups: bool = False) -> list[dict[str, object]]:
        """Rank source identities by use; ties retain declared source order."""
        source = self.groups if groups else self.parts
        counts = self.group_occurrences if groups else self.occurrences
        containing = self.group_containing if groups else self.containing
        denominator = self.occurrences.total()
        rows = []
        for identity in sorted(source, key=lambda key: -counts[key]):
            row = {
                "group" if groups else "part_id": identity,
                "occurrences": counts[identity],
                "designs": containing[identity],
                "occurrence_denominator": denominator,
                "design_denominator": self.accepted,
                "occurrence_share": counts[identity] / denominator
                if denominator
                else None,
                "design_fraction": containing[identity] / self.accepted
                if self.accepted
                else None,
                "reason": None if self.accepted else "empty_population",
            }
            if not groups:
                row["part_id"] = self.parts[identity].part_id
                row["group"] = self.parts[identity].group
                if identity in self.collections:
                    row["collection_id"] = self.collections[identity]
                    row["part_ref"] = identity
            rows.append(row)
        return rows

    def occupancy(self) -> list[dict[str, object]]:
        """Use designs reaching each position as the occupancy denominator."""
        boundaries = sorted({0, *self.coverage, *self.lengths})
        covered, available, rows = 0, self.accepted, []
        for start, end in pairwise(boundaries):
            covered += self.coverage[start]
            available -= self.lengths[start]
            if available:
                row = {
                    "start": start,
                    "end": end,
                    "designs": covered,
                    "denominator": available,
                    "fraction": covered / available,
                }
                if (
                    rows
                    and rows[-1]["end"] == start
                    and rows[-1]["designs"] == covered
                    and rows[-1]["denominator"] == available
                ):
                    rows[-1]["end"] = end
                else:
                    rows.append(row)
        return rows

    def to_dict(self, *, selected: bool = False) -> dict[str, object]:
        """Keep the eligible reference pool and all denominators in the report."""
        denominator = self.occurrences.total()
        return {
            "supply": {
                "eligible_parts": len(self.parts),
                "eligible_groups": len(self.groups),
                "unused_parts": len(self.parts)
                - sum(n > 0 for n in self.containing.values()),
                "unused_groups": len(self.groups)
                - sum(n > 0 for n in self.group_containing.values()),
            },
            "concentration": {
                "highest_part_occurrence_share": max(self.occurrences.values())
                / denominator
                if denominator
                else None,
                "occurrence_denominator": denominator,
                "reference": "eligible_parts_in_source_plans"
                if self.plan is None
                else "eligible_parts_in_plan",
                "reason": None if denominator else "empty_population",
            },
            "composition": {
                name: self.composition.get(name, Distribution()).to_dict(
                    denominator="selected_designs" if selected else "accepted_designs"
                )
                for name in COMPOSITION_METRICS
            },
            "padding_sides": dict(self.padding_sides),
            "occupancy": self.occupancy(),
            "part_usage": self.usage(),
            "group_usage": self.usage(groups=True),
            "requirements": list(self.requirements.values()),
        }

    def merge(self, source: QualityAccumulator, *, cell_ref: str) -> None:
        """Combine a source cell without collapsing part or requirement namespaces."""
        self.accepted += source.accepted
        self._merge_usage(source)
        for sequence in source.sequences - self.sequences:
            self.budget.retain()
            self.sequences.add(sequence)
        for name, distribution in source.composition.items():
            self.composition.setdefault(name, Distribution()).merge(
                distribution, self.budget
            )
        self.padding_sides.update(source.padding_sides)
        for attribute in ("coverage", "lengths"):
            for position, count in getattr(source, attribute).items():
                self._event(getattr(self, attribute), position, count)
        for identity, counts in source.requirements.items():
            self.budget.retain()
            self.requirements[f"{cell_ref}/{identity}"] = {
                **counts,
                "cell_ref": cell_ref,
            }

    def _merge_usage(self, source: QualityAccumulator) -> None:
        """Merge shared eligible collections and literal supplied group labels."""
        collection = source.collection_id
        for identity, part in source.parts.items():
            reference = f"{collection}/{identity}"
            if reference not in self.parts:
                self.budget.retain()
                self.parts[reference] = part
                self.collections[reference] = collection
            elif self.parts[reference] != part:
                msg = f"conflicting part evidence for {reference}"
                raise ArtifactIntegrityError(msg)
            self.occurrences[reference] += source.occurrences[identity]
            self.containing[reference] += source.containing[identity]
        for group in source.groups:
            if group not in self.groups:
                self.budget.retain()
                self.groups[group] = None
            self.group_occurrences[group] += source.group_occurrences[group]
            self.group_containing[group] += source.group_containing[group]
