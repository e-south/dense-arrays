"""Pack DNA motif libraries and describe their realized placements.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .errors import (
        InfeasibleError,
        InvalidSolverResultError,
        OptimizationError,
        SolverBackendError,
        UnprovenSolutionError,
    )
    from .optimizer import Optimizer
    from .solution import DenseArray

__all__ = [
    "DenseArray",
    "InfeasibleError",
    "InvalidSolverResultError",
    "OptimizationError",
    "Optimizer",
    "SolverBackendError",
    "UnprovenSolutionError",
]


def __getattr__(name: str) -> object:
    """Load public classes on demand so playback does not import a solver.

    Returns
    -------
    object
        The requested public class.

    Raises
    ------
    AttributeError
        If the name is not part of the public package interface.
    """
    if name not in __all__:
        msg = f"module {__name__!r} has no attribute {name!r}"
        raise AttributeError(msg)
    module = {"Optimizer": "optimizer", "DenseArray": "solution"}.get(name, "errors")
    return getattr(import_module(f"{__name__}.{module}"), name)
