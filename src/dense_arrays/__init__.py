"""Dense Arrays allows to create densely packed string arrays from a library of motifs.

It models the String Packing Problem as an Integer Linear Programming problem.
"""

__all__ = ["DenseArray", "Optimizer"]

from .optimize import DenseArray, Optimizer
