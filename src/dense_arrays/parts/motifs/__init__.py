"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/__init__.py

Validated motif inputs and explicit supplied-matrix score semantics.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .artifacts import MotifInput, PWMArtifact, read_artifact
from .models import Motif
from .scoring import MotifHit, MotifScore, best_hit, score_core

__all__ = [
    "Motif",
    "MotifHit",
    "MotifInput",
    "MotifScore",
    "PWMArtifact",
    "best_hit",
    "read_artifact",
    "score_core",
]
