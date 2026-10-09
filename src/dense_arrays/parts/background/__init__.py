"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/background/__init__.py

Exact conditional background generation from native sequence requirements.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .compiler import BackgroundCompilation, compile_background, model_identity
from .contracts import (
    CONDITIONAL_POLICY,
    ConditionalLimits,
    ConstructionLimitError,
    ConstructionReport,
)
from .counting import CountedBackground

__all__ = [
    "CONDITIONAL_POLICY",
    "BackgroundCompilation",
    "ConditionalLimits",
    "ConstructionLimitError",
    "ConstructionReport",
    "CountedBackground",
    "compile_background",
    "model_identity",
]
