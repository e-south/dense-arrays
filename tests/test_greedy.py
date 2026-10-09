"""Greedy realization preserves path entries and distinct occurrences.

Author: Eric J. South.
"""

import json
from pathlib import Path

import pytest

from dense_arrays import Optimizer


def test_greedy_cooperative_deadline_stops_before_search():
    from dense_arrays.greedy import realize_greedy  # noqa: PLC0415
    from dense_arrays.problem import PackingProblem  # noqa: PLC0415

    problem = PackingProblem.create(["AAA", "CCC"], 6, "single")
    with pytest.raises(TimeoutError, match="greedy"):
        realize_greedy(problem, deadline=0)


def test_greedy_preserves_released_fixture_geometry():
    fixture = json.loads(
        (
            Path(__file__).parent / "fixtures/workflow/dense-arrays-greedy-v1.json"
        ).read_text()
    )
    for case in fixture["cases"]:
        result = Optimizer(
            case["library"], case["length"], case["strands"]
        ).approximate()
        assert {
            "sequence": result.sequence,
            "offsets_fwd": result.offsets_fwd,
            "offsets_rev": result.offsets_rev,
            "occurrences": result.nb_motifs,
        } == case["expected"]


@pytest.mark.parametrize(
    "constraint", ["regulator", "promoter", "left", "right", "weight", "forbid"]
)
def test_approximate_rejects_unsupported_configuration(constraint: str):
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    if constraint == "regulator":
        optimizer.add_regulator_constraints(["R1", "R2"], required={"R2"})
    elif constraint == "promoter":
        optimizer.add_promoter_constraints(
            upstream="CCC", downstream="AAA", upstream_pos=0, spacer_length=0
        )
    elif constraint in {"left", "right"}:
        optimizer.add_side_biases(**{constraint: ["CCC"]})
    else:
        optimizer.build_model("CBC")
        if constraint == "weight":
            optimizer.set_motif_weight(0, 2.0)
        else:
            optimizer.forbid(optimizer.solve())
    with pytest.raises(ValueError, match=r"approximate.*support"):
        optimizer.approximate()


@pytest.mark.parametrize("strands", ["single", "double"])
@pytest.mark.parametrize("length, count", [(3, 1), (6, 2)])
def test_duplicate_entries_need_distinct_occurrences(
    strands: str, length: int, count: int
):
    solution = Optimizer(["AAA", "AAA"], length, strands).approximate()
    assert solution.nb_motifs == count
    offsets = [offset for offset, _ in solution.offset_indices_in_order()]
    assert offsets == [0] if count == 1 else offsets == [0, 3]


def test_overlong_entry_does_not_leave_gaps():
    solution = Optimizer(["AAAAAC", "ACG", "CGT"], 5, "single").approximate()
    assert solution.sequence == "ACGT"
    assert solution.offsets_fwd == [None, 0, 1]
    assert solution.nb_motifs == 2


def test_contained_motif_is_not_incidental_coverage():
    solution = Optimizer(["ATGTTAACT", "TTAA"], 9, "single").approximate()
    assert solution.nb_motifs == 1


def test_prefix_expansion_is_a_valid_path_entry():
    solution = Optimizer(["AAA", "AA"], 3, "single").approximate()
    assert solution.sequence == "AAA"
    assert solution.offsets_fwd == [0, 0]
    assert solution.nb_motifs == 2


def test_greedy_does_not_substitute_a_solver(monkeypatch: pytest.MonkeyPatch):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("approximate() invoked model construction")

    monkeypatch.setattr(Optimizer, "build_model", forbidden)
    assert Optimizer(["ACG", "CGT"], 4, "double").approximate().nb_motifs == 2
