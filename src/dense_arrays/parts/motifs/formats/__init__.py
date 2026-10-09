"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/formats/__init__.py

Strict file adapters return probability models without inventing score matrices.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .jaspar import read_jaspar
from .meme import read_meme

__all__ = ["read_jaspar", "read_meme"]
