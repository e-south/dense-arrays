"""
--------------------------------------------------------------------------------
<dense-array project>

Dense Arrays allows to create densely packed string arrays from a library of motifs.
It models the String Packing Problem as an Integer Linear Programming problem.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
--------------------------------------------------------------------------------
"""

__all__ = ["DenseArray", "Optimizer"]

from .optimizer import Optimizer
from .solution import DenseArray
