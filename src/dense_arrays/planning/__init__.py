"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/__init__.py

Typed design requests and immutable generation previews.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .batches import (
    BatchSampling,
    BatchSchedule,
    CandidateBatch,
    FeedbackPolicy,
    FeedbackSnapshot,
    Resampling,
)
from .evidence import PlanEvidence
from .extension import ExtensionSpec, ParentRun
from .libraries import LibraryExclusion
from .lineage import Lineage, RunReference
from .matrices import Allocation, MatrixPlan, MatrixSpec, Variant
from .models import Assembly, DesignSpec, Length, Limits, Padding, Target
from .preparation import PreparationPlan
from .requirements import (
    GC,
    Avoid,
    Fixed,
    GroupCoverage,
    Occurrences,
    Spacing,
    StartWindow,
)
from .resolution import GenerationPlan

__all__ = [
    "GC",
    "Allocation",
    "Assembly",
    "Avoid",
    "BatchSampling",
    "BatchSchedule",
    "CandidateBatch",
    "DesignSpec",
    "ExtensionSpec",
    "FeedbackPolicy",
    "FeedbackSnapshot",
    "Fixed",
    "GenerationPlan",
    "GroupCoverage",
    "Length",
    "LibraryExclusion",
    "Limits",
    "Lineage",
    "MatrixPlan",
    "MatrixSpec",
    "Occurrences",
    "Padding",
    "ParentRun",
    "PlanEvidence",
    "PreparationPlan",
    "Resampling",
    "RunReference",
    "Spacing",
    "StartWindow",
    "Target",
    "Variant",
]
