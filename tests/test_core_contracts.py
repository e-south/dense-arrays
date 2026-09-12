"""Input and state contracts for the optimizer and direct results.

Author: Eric J. South.
"""

import pytest

from dense_arrays import DenseArray, Optimizer
from dense_arrays.constraints import PromoterConstraint


@pytest.mark.parametrize("length", [True, 3.5, "3"])
@pytest.mark.parametrize("owner", [Optimizer, DenseArray])
def test_lengths_are_discrete(length: object, owner: type):
    args = (["AAA"], length, [0], [None]) if owner is DenseArray else (["AAA"], length)
    with pytest.raises(ValueError, match="sequence_length"):
        owner(*args)


@pytest.mark.parametrize("count", [True, 1.9, "1"])
@pytest.mark.parametrize("field", ["min_count_by_regulator", "min_required_regulators"])
def test_regulator_counts_are_discrete(count: object, field: str):
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    kwargs = {field: {"R": count} if field == "min_count_by_regulator" else count}
    with pytest.raises(ValueError, match=field):
        optimizer.add_regulator_constraints(["R", "R"], **kwargs)
    optimizer.add_regulator_constraints(["R", "R"], required={"R"})


@pytest.mark.parametrize(
    "interval", [(1,), (1, 2, 3), (3, 1), (True, 3), (1.2, 3), [1, 3]]
)
def test_invalid_intervals_fail_at_configuration(interval: object):
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    with pytest.raises(ValueError, match="upstream_pos"):
        optimizer.add_promoter_constraints(
            upstream="AAA", downstream="CCC", upstream_pos=interval
        )
    assert not optimizer.promoters
    assert optimizer.model is None


@pytest.mark.parametrize("indices", [(True, 1), (-1, 1), (0, 0), (0, 1.5)])
def test_direct_promoter_indices(indices: tuple[object, object]):
    with pytest.raises(ValueError, match=r"index|indices"):
        PromoterConstraint(upstream_index=indices[0], downstream_index=indices[1])


def test_negative_spacer_supports_intentional_overlap():
    optimizer = Optimizer(["AAA", "AAT"], 4, "single")
    optimizer.add_promoter_constraints(
        upstream="AAA", downstream="AAT", upstream_pos=0, spacer_length=-2
    )
    assert optimizer.optimal("CBC").sequence == "AAAT"


@pytest.mark.parametrize("mapping", [{False: "R"}, {0.0: "R"}, [42], [" R "]])
def test_regulator_identity_is_not_coerced(mapping: object):
    optimizer = Optimizer(["AAA"], 3, "single")
    with pytest.raises(ValueError, match="regulator_by_index"):
        optimizer.add_regulator_constraints(mapping, min_required_regulators=1)


@pytest.mark.parametrize("library", [[""], ["N"], [42], []])
def test_direct_results_validate_motifs(library: list[object]):
    with pytest.raises(ValueError, match=r"motif|library"):
        DenseArray(library, 3, [0] * len(library), [None] * len(library))


@pytest.mark.parametrize("offset", [True, 0.5, "0"])
def test_direct_results_reject_nondiscrete_offsets(offset: object):
    with pytest.raises(ValueError, match="offset"):
        DenseArray(["AAA"], 3, [offset], [None])


def test_one_orientation_per_library_entry():
    with pytest.raises(ValueError, match="orientation"):
        DenseArray(["AT"], 2, [0], [0])


def test_optimizer_configuration_has_defensive_views():
    library = ["AAA", "AAT"]
    optimizer = Optimizer(library, 4, "single")
    library[1] = "CCC"
    optimizer.library[1] = "CCC"
    optimizer.adjacency_matrix[0][1] = 99
    assert optimizer.library == ["AAA", "AAT"]
    assert optimizer.optimal("CBC").sequence == "AAAT"
    for name, value in [
        ("library", ["CCC"]),
        ("sequence_length", 99),
        ("strands", "double"),
    ]:
        with pytest.raises(AttributeError):
            setattr(optimizer, name, value)


def test_bias_update_is_atomic_and_views_are_defensive():
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    optimizer.add_side_biases(left=["CCC"], right=["AAA"])
    with pytest.raises(ValueError, match="All motifs"):
        optimizer.add_side_biases(left=["AAA"], right=["XXX"])
    assert optimizer.ilefts == [1]
    assert optimizer.irights == [0]
    optimizer.ilefts.clear()
    assert optimizer.ilefts == [1]


def test_promoter_configuration_cannot_be_changed_through_views():
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    optimizer.add_promoter_constraints(upstream="AAA", downstream="CCC")
    constraint = optimizer.promoters[0]
    with pytest.raises(AttributeError):
        constraint.upstream_index = 1
    optimizer.promoters.clear()
    assert len(optimizer.promoters) == 1


def test_densearray_is_an_immutable_snapshot():
    library = ["AAA"]
    offsets = [0]
    solution = DenseArray(library, 3, offsets, [None])
    library[0] = "CCC"
    offsets[0] = 2
    solution.library[0] = "CCC"
    solution.offsets_fwd[0] = 2
    assert solution.library == ["AAA"]
    assert solution.offsets_fwd == [0]
    assert solution.sequence == "AAA"
    with pytest.raises(AttributeError):
        solution.sequence = "CCC"


@pytest.mark.parametrize("index", [-1, True, 1.2, 2])
def test_bad_weight_index_leaves_model_unchanged(index: object):
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    optimizer.build_model("CBC")
    before = optimizer.model.ExportModelAsLpFormat(obfuscate=False)
    with pytest.raises(ValueError, match="imotif"):
        optimizer.set_motif_weight(index, 2.0)
    assert optimizer.model.ExportModelAsLpFormat(obfuscate=False) == before


@pytest.mark.parametrize("weight", [True, float("nan"), float("inf"), "2"])
def test_bad_weight_leaves_model_unchanged(weight: object):
    optimizer = Optimizer(["AAA", "CCC"], 6, "single")
    optimizer.build_model("CBC")
    before = optimizer.model.ExportModelAsLpFormat(obfuscate=False)
    with pytest.raises(ValueError, match="weight"):
        optimizer.set_motif_weight(0, weight)
    assert optimizer.model.ExportModelAsLpFormat(obfuscate=False) == before


@pytest.mark.parametrize(
    "solution",
    [
        DenseArray(["CCC"], 3, [0], [None]),
        DenseArray(["AAA"], 4, [0], [None]),
        DenseArray(["AAA"], 3, [None], [0]),
    ],
)
def test_foreign_forbid_is_atomic(solution: DenseArray):
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.build_model("CBC")
    before = optimizer.model.NumConstraints()
    with pytest.raises(ValueError, match=r"problem|strand"):
        optimizer.forbid(solution)
    assert optimizer.model.NumConstraints() == before


def test_forbid_rejects_nonpath_placements():
    optimizer = Optimizer(["AAAA", "AA"], 4, "single")
    optimizer.build_model("CBC")
    before = optimizer.model.NumConstraints()
    direct = DenseArray(["AAAA", "AA"], 4, [0, 1], [None, None])
    with pytest.raises(ValueError, match="path"):
        optimizer.forbid(direct)
    assert optimizer.model.NumConstraints() == before


@pytest.mark.parametrize("required", ["R", ["R"]])
def test_required_regulator_collection_is_explicit(required: object):
    optimizer = Optimizer(["AAA"], 3, "single")
    with pytest.raises(ValueError, match="required"):
        optimizer.add_regulator_constraints(["R"], required=required)


def test_minimum_count_collection_is_explicit():
    optimizer = Optimizer(["AAA"], 3, "single")
    with pytest.raises(ValueError, match="min_count_by_regulator"):
        optimizer.add_regulator_constraints(["R"], min_count_by_regulator=[("R", 1)])


def test_bias_collection_does_not_interpret_string_as_entries():
    optimizer = Optimizer(["A", "T"], 2, "single")
    with pytest.raises(ValueError, match="list"):
        optimizer.add_side_biases(left="AT")
    assert optimizer.ilefts == []


def test_unrepresentable_weight_is_rejected_before_mutation():
    optimizer = Optimizer(["AAA"], 3, "single")
    optimizer.build_model("CBC")
    before = optimizer.model.ExportModelAsLpFormat(obfuscate=False)
    with pytest.raises(ValueError, match="weight"):
        optimizer.set_motif_weight(0, 10**400)
    assert optimizer.model.ExportModelAsLpFormat(obfuscate=False) == before
