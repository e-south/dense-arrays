"""Tests of the optimize.py module."""

import itertools as it
from collections import Counter

import dense_arrays as da
import pytest


@pytest.mark.parametrize(
    "motifa, motifb, shift",
    [
        ("ATGCATTA", "CATTATG", 3),
        ("ATGCATTA", "TATGA", 6),
        ("ATGTTAACT", "TTAA", 8),
        ("AGC", "CAG", 2),
        ("CAG", "AGC", 1),
        ("CGT", "CAG", 3),
        ("", "CGT", 0),
    ],
)
def test_shift_metric(motifa: str, motifb: str, shift: int):
    assert da.optimize.shift_metric(motifa, motifb) == shift


@pytest.mark.parametrize("strands", ["single", "double"])
def test_simple_optimal(strands: str):
    opt = da.Optimizer(
        ["ATGC", "CGT", "ATTA", "TTATTA"], sequence_length=8, strands=strands
    )
    best = opt.optimal()
    assert best.nb_motifs == 3
    assert best.sequence == "CGTTATTA"


@pytest.mark.parametrize(
    "strands, diverse, ns",
    [
        ("single", False, [1, 4, 9, 1]),
        ("single", True, [1, 4, 9, 1]),
        ("double", False, [1, 8, 34, 8]),
        ("double", True, [1, 8, 34, 8]),
    ],
)
def test_simple_solutions(strands: str, diverse: bool, ns: list[int]):
    opt = da.Optimizer(
        ["ATGC", "CGT", "ATTA", "TTATTA"], sequence_length=8, strands=strands
    )
    iterator = opt.solutions_diverse() if diverse else opt.solutions()
    solutions = list(it.islice(iterator, sum(ns) + 1))
    # There is the correct number of solutions
    assert len(solutions) == sum(ns)
    # There is the correct number of solution of each score
    sizes = Counter(sol.nb_motifs for sol in solutions)
    for i, n in enumerate(ns):
        assert sizes[i] == n
    # Solutions are ordered by score
    for sola, solb in it.pairwise(solutions):
        assert sola.nb_motifs >= solb.nb_motifs
