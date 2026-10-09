"""Packing length and coordinate origin remain explicit before assembly.

Author: Eric J. South.
"""

import pytest

from dense_arrays import InfeasibleError, Optimizer


def test_exact_packing_length_does_not_invent_padding():
    assert Optimizer(["AAA"], 4, "single").optimal().sequence == "AAA"
    optimizer = Optimizer(["AAA"], 4, "single", length_mode="exact")
    with pytest.raises(InfeasibleError):
        optimizer.optimal()
    solution = Optimizer(["AAA", "CCC"], 6, "single", length_mode="exact").optimal()
    assert len(solution.sequence) == 6


def test_end_relative_position_supports_final_left_padding_coordinates():
    optimizer = Optimizer(["AACC", "CCGT", "TT"], 10, "single")
    optimizer.add_fixed_occurrence(0, orientation="forward", origin="end", start=-6)
    solution = optimizer.optimal()
    assert solution.offsets_fwd[0] - len(solution.sequence) == -6
    pad_left = 10 - len(solution.sequence)
    assert solution.offsets_fwd[0] + pad_left == 4


def test_greedy_cannot_ignore_exact_packing_length():
    optimizer = Optimizer(["AAA"], 4, length_mode="exact")
    with pytest.raises(ValueError, match="support"):
        optimizer.approximate()
