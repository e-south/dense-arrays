"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/batches/packing.py

Pack only offered parts while preserving eligible-collection placement identities.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import replace

from dense_arrays.generation.heuristic import HeuristicResult
from dense_arrays.parts import PartSelector
from dense_arrays.planning import (
    Fixed,
    GenerationPlan,
    GroupCoverage,
    Occurrences,
    Spacing,
)
from dense_arrays.planning.requirements import Requirement
from dense_arrays.solution import DenseArray
from dense_arrays.solver import SolveReport, SolverIdentity


def _lower(plan: GenerationPlan) -> GenerationPlan:
    batch = plan.request.batch
    by_id = {p.part_id: p for p in plan.request.parts}
    offered = tuple(by_id[key] for key in batch.part_ids)
    membership = set(batch.part_ids)
    groups = {p.group for p in offered}
    rules = []
    impossible = False
    for rule in plan.request.requirements:
        lowered, unavailable = _lower_rule(rule, plan, membership, groups)
        impossible |= unavailable
        if lowered is not None:
            rules.append(lowered)
    if impossible:
        # Express a real contradiction through supported constraints, so CBC owns
        # its infeasibility report just as it does for geometric contradictions.
        rules = [
            Occurrences("batch_required", PartSelector(part_ids=batch.part_ids), min=1),
            Occurrences(
                "batch_unavailable", PartSelector(part_ids=batch.part_ids), max=0
            ),
        ]
    return GenerationPlan(
        plan.request.with_changes(
            parts=offered,
            requirements=tuple(rules),
            batch=None,
            lineage=None,
            exclude=None,
        )
    )


def _lower_rule(
    rule: Requirement,
    plan: GenerationPlan,
    membership: set[str],
    groups: set[str | None],
) -> tuple[Requirement | None, bool]:
    """Intersect selector capacity; distinguish absent requirements from screens."""
    if isinstance(rule, Fixed):
        return (rule, False) if rule.part_id in membership else (None, True)
    if isinstance(rule, Spacing):
        return (
            (rule, False)
            if {rule.upstream, rule.downstream} <= membership
            else (None, True)
        )
    if isinstance(rule, Occurrences):
        eligible = {
            plan.request.parts[i].part_id
            for i in rule.select.indices(plan.request.parts)
        }
        selected = tuple(key for key in plan.request.batch.part_ids if key in eligible)
        if (rule.min or 0) > len(selected):
            return None, True
        return (
            (replace(rule, select=PartSelector(part_ids=selected)), False)
            if selected
            else (None, False)
        )
    if isinstance(rule, GroupCoverage):
        selected = tuple(group for group in rule.groups if group in groups)
        return (
            (None, True)
            if rule.min > len(selected)
            else (replace(rule, groups=selected), False)
        )
    # Final-sequence screens stay with the full request, outside this model.
    return None, False


class BatchOptimizer:
    """Translate compact-model offsets at the workflow boundary, never by DNA text."""

    def __init__(self, plan: GenerationPlan, *, seconds: float) -> None:
        """Build adjacency only for the declared offered subset."""
        from dense_arrays.generation.packing import build_optimizer  # noqa: PLC0415

        lowered = _lower(plan)
        self._optimizer = build_optimizer(lowered, seconds=seconds)
        self._library = [p.sequence for p in plan.request.parts]
        by_id = {p.part_id: i for i, p in enumerate(plan.request.parts)}
        self._indices = tuple(by_id[p.part_id] for p in lowered.request.parts)
        self._length = plan.request.length.maximum or plan.request.length.exact

    @property
    def solver_identity(self) -> SolverIdentity | None:
        """Report the backend observed from the compact model."""
        return self._optimizer.solver_identity

    @property
    def packing_objective(self) -> dict[str, object]:
        """Keep objective keys in supplied part identities, not compact indices."""
        return self._optimizer.packing_objective

    def solve_report(
        self, *, time_limit_seconds: float
    ) -> SolveReport | HeuristicResult:
        """Lift placements into original collection indices after the actual solve."""
        report = self._optimizer.solve_report(time_limit_seconds=time_limit_seconds)
        if report.solution is None:
            return report
        forward, reverse = [None] * len(self._library), [None] * len(self._library)
        for local, original in enumerate(self._indices):
            forward[original] = report.solution.offsets_fwd[local]
            reverse[original] = report.solution.offsets_rev[local]
        return replace(
            report, solution=DenseArray(self._library, self._length, forward, reverse)
        )

    def forbid(self, solution: DenseArray) -> None:
        """Reapply a recorded exclusion to the same compact model on replay."""
        self._optimizer.forbid(
            DenseArray(
                [self._library[i] for i in self._indices],
                self._length,
                [solution.offsets_fwd[i] for i in self._indices],
                [solution.offsets_rev[i] for i in self._indices],
            )
        )
