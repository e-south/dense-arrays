"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/__init__.py

Typed, read-only views of committed workflow evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dense_arrays.artifacts.bundles import BundleSummary
from dense_arrays.artifacts.reading import ReadCost, ReadLimitError, ReadLimits
from dense_arrays.parts.filters import Range

from .bundles import BundleView
from .collections import LibraryView
from .design_filters import DesignFilter
from .diagnostics import Diagnostic, DiagnosticReport
from .filters import AttemptFilter
from .plans import PlanChange, PlanComparison, PlanFilter
from .plans.requests import RequestReport
from .pools import CandidateFilter, PoolQualityReport, PoolQualitySnapshot
from .projections import PlacementRecord, SequenceRecord
from .quality import MetricDifference, QualityComparison, QualityReport, QualitySnapshot
from .readers import RecordView
from .selections import LibrarySelection, SelectionShortfall, SelectionSnapshot, Take
from .selections.views import SelectionView
from .summary import RunSummary

__all__ = [
    "AttemptFilter",
    "BundleSummary",
    "BundleView",
    "CandidateFilter",
    "DesignFilter",
    "Diagnostic",
    "DiagnosticReport",
    "LibrarySelection",
    "LibraryView",
    "MetricDifference",
    "PlacementRecord",
    "PlanChange",
    "PlanComparison",
    "PlanFilter",
    "PoolQualityReport",
    "PoolQualitySnapshot",
    "QualityComparison",
    "QualityReport",
    "QualitySnapshot",
    "Range",
    "ReadCost",
    "ReadLimitError",
    "ReadLimits",
    "RecordView",
    "RequestReport",
    "RunSummary",
    "SelectionShortfall",
    "SelectionSnapshot",
    "SelectionView",
    "SequenceRecord",
    "Take",
]
