"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/__init__.py

Exact selected-library reports over native run snapshots.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .comparison import QualityComparison
from .differences import MetricDifference
from .models import QualityReport
from .snapshots import QualitySnapshot

__all__ = ["MetricDifference", "QualityComparison", "QualityReport", "QualitySnapshot"]
