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
        ("ATGC", "ATGC", 4),
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


@pytest.mark.parametrize(
    "strands, sequence_length, noprom, prom",
    [
        ("single", 10, 4, 3),
        ("double", 8, 4, 3),
    ],
)
def test_promoter_constraints(
    strands: str, sequence_length: int, noprom: int, prom: int
):
    opt = da.Optimizer(
        ["GCA", "CCC", "ATGC", "CATT"], sequence_length=sequence_length, strands=strands
    )
    sol_noprom = opt.optimal()
    assert sol_noprom.nb_motifs == noprom
    opt.add_promoter_constraints(
        upstream="ATGC", downstream="CCC", upstream_pos=(0, 2), spacer_length=(0, 3)
    )
    sol_prom = opt.optimal()
    assert sol_prom.nb_motifs == prom


def test_side_bias():
    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="double")
    opt.add_side_biases(left=["AAA"], right=["CCC"])
    sol_left = opt.optimal()
    norm_offset_indices_left = [
        (offset, index % opt.nb_motifs)
        for offset, index in sol_left.offset_indices_in_order()
    ]
    assert norm_offset_indices_left == [(0, 0), (3, 1)]

    opt.add_side_biases(left=["CCC"], right=["AAA"])
    sol_right = opt.optimal()
    norm_offset_indices_right = [
        (offset, index % opt.nb_motifs)
        for offset, index in sol_right.offset_indices_in_order()
    ]
    assert norm_offset_indices_right == [(0, 1), (3, 0)]


def test_side_bias_with_same_promoter():
    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="double")
    opt.add_promoter_constraints(
        upstream="AAA", downstream="CCC", upstream_pos=(0, 2), spacer_length=(0, 2)
    )

    opt.add_side_biases(left=["AAA"], right=["CCC"])
    sol_left = opt.optimal()
    norm_offset_indices_left = [
        (offset, index % opt.nb_motifs)
        for offset, index in sol_left.offset_indices_in_order()
    ]
    assert norm_offset_indices_left == [(0, 0), (3, 1)]

    opt.add_side_biases(left=["CCC"], right=["AAA"])
    sol_right = opt.optimal()
    norm_offset_indices_right = [
        (offset, index % opt.nb_motifs)
        for offset, index in sol_right.offset_indices_in_order()
    ]
    assert norm_offset_indices_right == [(0, 0), (3, 1)]


def test_side_bias_with_other_promoter():
    library = ["GGGT", "CTTC", "TAGG", "AATC", "TCTA"]
    opt = da.Optimizer(library, sequence_length=14, strands="double")

    sol = opt.optimal()
    assert sol.nb_motifs == 5
    assert sol.offset_indices_in_order() == [(0, 3), (3, 1), (5, 4), (7, 2), (9, 0)]

    opt.add_promoter_constraints(
        upstream="GGGT", downstream="CTTC", upstream_pos=(0, 3), spacer_length=(1, 4)
    )

    opt.add_side_biases(left=library[::2], right=library[1::2])
    sol_left = opt.optimal()
    assert sol_left.nb_motifs == 4
    assert sol_left.offset_indices_in_order() == [(0, 2), (2, 0), (5, 4), (9, 1)]

    opt.add_side_biases(left=library[1::2], right=library[::2])
    sol_right = opt.optimal()
    assert sol_right.nb_motifs == 4
    assert sol_right.offset_indices_in_order() == [(0, 0), (4, 3), (7, 1), (10, 7)]
