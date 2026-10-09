"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/matrices/__init__.py

Bounded matrix expansion and explicit cell targets.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .allocation import Allocation
from .requests import MatrixSpec, Variant
from .resolution import MatrixPlan

__all__ = ["Allocation", "MatrixPlan", "MatrixSpec", "Variant"]
