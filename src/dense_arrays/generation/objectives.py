"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/objectives.py

Apply replayable part-usage preference at the exact packing boundary.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.objectives import PackingObjective
from dense_arrays.model import part_usage_weights

if TYPE_CHECKING:
    from dense_arrays.optimizer import Optimizer
    from dense_arrays.solution import DenseArray
    from dense_arrays.solver import SolveReport, SolverIdentity


class UsageBalancedOptimizer:
    """Count every excluded proposal once, independently of acceptance or strand."""

    def __init__(self, optimizer: Optimizer, part_ids: tuple[str, ...]) -> None:
        """Start a fresh offered-batch history; replay uses the same exclusion seam."""
        self._optimizer = optimizer
        self._usage = dict.fromkeys(part_ids, 0)
        self._proposed = 0
        self._weights = (1.0,) * len(part_ids)

    @property
    def solver_identity(self) -> SolverIdentity:
        """Preserve the actual backend identity."""
        return self._optimizer.solver_identity

    @property
    def packing_objective(self) -> dict[str, object]:
        """Return detached pre-solve evidence without changing search state."""
        return PackingObjective(self._usage, self._proposed).to_dict()

    def solve_report(self, *, time_limit_seconds: float) -> SolveReport:
        """Keep the exact solver's proof and termination outcomes intact."""
        return self._optimizer.solve_report(time_limit_seconds=time_limit_seconds)

    def forbid(self, solution: DenseArray) -> None:
        """Exclude the path, then update the next attempt's secondary objective."""
        self._optimizer.forbid(solution)
        for part_id, forward, reverse in zip(
            self._usage, solution.offsets_fwd, solution.offsets_rev, strict=True
        ):
            if forward is not None or reverse is not None:
                self._usage[part_id] += 1
        self._proposed += 1
        weights = part_usage_weights(tuple(self._usage.values()))
        for index, (before, after) in enumerate(
            zip(self._weights, weights, strict=True)
        ):
            if before != after:
                self._optimizer.set_motif_weight(index, after)
        self._weights = weights
