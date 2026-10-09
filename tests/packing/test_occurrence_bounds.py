"""Count and coverage constraints use supplied identities, not string matches.

Author: Eric J. South.
"""

import pytest

from dense_arrays import Optimizer


def test_exact_and_overlapping_bounds_are_enforced_by_cbc():
    optimizer = Optimizer(["AAA", "AAA", "CCC", "GGG"], 12, "single")
    optimizer.add_count_constraint([0, 1, 2], minimum=2, maximum=2)
    optimizer.add_count_constraint([1, 2, 3], minimum=2)
    optimizer.add_group_coverage([[0, 1, 2], [3]], minimum=2)
    solution = optimizer.optimal()
    selected = {
        i for i, offset in enumerate(solution.offsets_fwd) if offset is not None
    }
    assert len(selected & {0, 1, 2}) == 2
    assert len(selected & {1, 2, 3}) >= 2
    assert 3 in selected


def test_zero_maximum_excludes_only_named_occurrence():
    optimizer = Optimizer(["AAA", "AAA", "CCC"], 9, "double")
    optimizer.add_count_constraint([0], maximum=0)
    solution = optimizer.optimal()
    assert solution.offsets_fwd[0] is solution.offsets_rev[0] is None
    assert solution.nb_motifs == 2


@pytest.mark.parametrize(
    "bounds",
    [
        {},
        {"minimum": True},
        {"maximum": 1.5},
        {"minimum": 3},
        {"minimum": 2, "maximum": 1},
        {"maximum": -1},
    ],
)
def test_invalid_bounds_fail_before_model_build(bounds: dict[str, object]):
    optimizer = Optimizer(["AAA", "CCC"], 6)
    with pytest.raises(ValueError, match=r"bound|integer|available|maximum"):
        optimizer.add_count_constraint([0, 1], **bounds)
    assert optimizer.model is None


def test_constraints_cannot_be_added_to_built_model_or_ignored_by_greedy():
    optimizer = Optimizer(["AAA", "CCC"], 6)
    optimizer.add_count_constraint([0], minimum=1)
    with pytest.raises(ValueError, match="support"):
        optimizer.approximate()
    optimizer.build_model()
    with pytest.raises(RuntimeError, match="before"):
        optimizer.add_count_constraint([1], minimum=1)


def test_coverage_requires_distinct_nonempty_groups():
    optimizer = Optimizer(["AAA", "CCC"], 6)
    with pytest.raises(ValueError, match="group"):
        optimizer.add_group_coverage([[0], [0]], minimum=2)
