"""Fixed geometry names occurrence indices and supports either orientation.

Author: Eric J. South.
"""

import pytest

from dense_arrays import Optimizer


def test_repeated_string_does_not_change_the_named_fixed_occurrence():
    optimizer = Optimizer(["AAA", "AAA", "AAC"], 9, "double")
    optimizer.add_fixed_occurrence(1, orientation="reverse", start=0)
    optimizer.add_count_constraint([0], maximum=0)
    solution = optimizer.optimal()
    assert solution.offsets_rev[1] == 0
    assert solution.offsets_fwd[1] is None
    assert solution.offsets_fwd[0] is solution.offsets_rev[0] is None
    assert solution.sequence.startswith("TTT")


def test_spacing_uses_downstream_start_minus_upstream_end():
    optimizer = Optimizer(["AACC", "CCGT", "GTTA"], 8, "single")
    optimizer.add_fixed_occurrence(0, orientation="forward", start=0)
    optimizer.add_fixed_occurrence(2, orientation="forward", start=(4, 5))
    optimizer.add_spacing_constraint(0, 2, minimum=0, maximum=0)
    solution = optimizer.optimal()
    assert solution.offsets_fwd == [0, 2, 4]
    assert solution.sequence == "AACCGTTA"


def test_negative_spacing_preserves_an_intentional_overlap():
    optimizer = Optimizer(["AACC", "CCGT"], 6, "single")
    optimizer.add_fixed_occurrence(0, orientation="forward")
    optimizer.add_fixed_occurrence(1, orientation="forward")
    optimizer.add_spacing_constraint(0, 1, minimum=-2, maximum=-2)
    solution = optimizer.optimal()
    assert solution.offsets_fwd == [0, 2]


def test_invalid_fixed_geometry_fails_before_model_build():
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    with pytest.raises(ValueError, match="reverse"):
        optimizer.add_fixed_occurrence(0, orientation="reverse")
    with pytest.raises(ValueError, match="integer"):
        optimizer.add_fixed_occurrence(0, orientation="forward", start=True)
    with pytest.raises(ValueError, match="fixed"):
        optimizer.add_spacing_constraint(0, 1, minimum=0, maximum=1)
    assert optimizer.model is None


def test_greedy_cannot_drop_fixed_requirements():
    optimizer = Optimizer(["AAA", "CCC"], 6)
    optimizer.add_fixed_occurrence(1, orientation="forward")
    with pytest.raises(ValueError, match="support"):
        optimizer.approximate()
