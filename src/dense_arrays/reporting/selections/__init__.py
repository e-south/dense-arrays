"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/selections/__init__.py

Reproducible allocation of an inspected design population.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .requests import LibrarySelection, Take
from .snapshots import SelectionShortfall, SelectionSnapshot

__all__ = ["LibrarySelection", "SelectionShortfall", "SelectionSnapshot", "Take"]
