"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/__init__.py

Part identities, selectors and explicit preparation requests.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .background.contracts import ConditionalLimits
from .bound import BoundParts
from .filters import PartFilter
from .models import Normalization, Part, PartSelector, PartTable
from .motifs.artifacts import PWMArtifact
from .motifs.windows import MotifWindow
from .pools import PoolHandle, PoolSource
from .preparation import PreparationSet, PreparationSpec, Retention
from .provenance import ImportReport
from .retention import MMR, PoolSize, ScoreBands
from .sampling import (
    Background,
    CandidateBudget,
    Eligibility,
    LengthRange,
    MiningTarget,
    Sampling,
    Uniqueness,
)
from .scoring.configuration import FimoScoring, ScoringLimits
from .screening import PWMExclusion

__all__ = [
    "MMR",
    "Background",
    "BoundParts",
    "CandidateBudget",
    "ConditionalLimits",
    "Eligibility",
    "FimoScoring",
    "ImportReport",
    "LengthRange",
    "MiningTarget",
    "MotifWindow",
    "Normalization",
    "PWMArtifact",
    "PWMExclusion",
    "Part",
    "PartFilter",
    "PartSelector",
    "PartTable",
    "PoolHandle",
    "PoolSize",
    "PoolSource",
    "PreparationSet",
    "PreparationSpec",
    "Retention",
    "Sampling",
    "ScoreBands",
    "ScoringLimits",
    "Uniqueness",
]
