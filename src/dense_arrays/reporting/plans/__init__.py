"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/plans/__init__.py

Inspect resolved plan evidence and compare design semantics.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .comparison import PlanChange, PlanComparison
from .filters import PlanFilter
from .reading import check_plan_limits, editable_request, read_plan

__all__ = [
    "PlanChange",
    "PlanComparison",
    "PlanFilter",
    "check_plan_limits",
    "editable_request",
    "read_plan",
]
