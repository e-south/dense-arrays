"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/heuristic.py

One bounded greedy proposal per offered batch, with explicit replay state.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays.artifacts.search import HeuristicEvidence
from dense_arrays.errors import InfeasibleError
from dense_arrays.greedy import realize_greedy
from dense_arrays.problem import PackingProblem

if TYPE_CHECKING:
    from dense_arrays.planning import GenerationPlan
    from dense_arrays.solution import DenseArray


@dataclass(frozen=True)
class HeuristicResult:
    """A candidate and its method evidence; intentionally no solver status."""

    evidence: HeuristicEvidence
    solution: DenseArray | None = None


class GreedySearch:
    """Reuse the packing primitive and consume at most one proposal per batch."""

    solver_identity = None

    def __init__(self, plan: GenerationPlan) -> None:
        """Allocate input validation without building a mathematical solver."""
        request = plan.request
        self._problem = PackingProblem.create(
            [p.sequence for p in request.parts],
            request.length.maximum or request.length.exact,
            request.strands,
        )
        self._consumed = False

    def solve_report(self, *, time_limit_seconds: float) -> HeuristicResult:
        """Return only a completed greedy proposal; never infer global feasibility."""
        if self._consumed:
            return HeuristicResult(HeuristicEvidence("exhausted"))
        try:
            solution = realize_greedy(
                self._problem, deadline=time.monotonic() + time_limit_seconds
            )
        except InfeasibleError:
            return HeuristicResult(HeuristicEvidence("exhausted"))
        except TimeoutError:
            return HeuristicResult(HeuristicEvidence("time_limit"))
        return HeuristicResult(HeuristicEvidence("candidate"), solution)

    def forbid(self, _solution: DenseArray) -> None:
        """Consume the sole proposal when a packing is committed or replayed."""
        self._consumed = True
