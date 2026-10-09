"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/__init__.py

Portable collections of supplied sequences, parts and realized placements.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from .models import ArrayCollection, ArrayFilter, ArrayRecord, CollectionSummary
from .reading import CollectionView

__all__ = [
    "ArrayCollection",
    "ArrayFilter",
    "ArrayRecord",
    "CollectionSummary",
    "CollectionView",
]
