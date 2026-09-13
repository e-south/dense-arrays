"""Distinct outcomes of optimization and solver result validation.

Module Author(s): Eric J. South
"""


class OptimizationError(RuntimeError):
    """Base exception for unsuccessful optimization execution."""


class InfeasibleError(ValueError):
    """The solver proved that no feasible path exists."""


class UnprovenSolutionError(OptimizationError):
    """A feasible incumbent exists, but optimality was not proved."""


class SolverBackendError(OptimizationError):
    """The backend could not establish an optimization result."""


class InvalidSolverResultError(OptimizationError):
    """A reported optimal result violates the path or array contract."""
