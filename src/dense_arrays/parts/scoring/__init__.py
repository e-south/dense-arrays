"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/scoring/__init__.py

Optional motif scoring with bound inputs and finite execution limits.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .binding import FimoBinding, bind_fimo
from .configuration import FimoScoring, ScoringLimits
from .fimo import scan_fimo
from .process import ScoringError
from .records import FimoHit, FimoResult

__all__ = [
    "FimoBinding",
    "FimoHit",
    "FimoResult",
    "FimoScoring",
    "ScoringError",
    "ScoringLimits",
    "bind_fimo",
    "scan_fimo",
]
