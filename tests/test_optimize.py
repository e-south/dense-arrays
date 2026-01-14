"""Tests of the optimize.py module."""

import itertools as it
from collections import Counter

import pytest

import dense_arrays as da


def _assert_promoter_constraint(
    sol: da.DenseArray,
    library: list[str],
    constraint: dict[str, object],
) -> None:
    upstream = constraint["upstream"]
    downstream = constraint["downstream"]
    upstream_pos = constraint["upstream_pos"]
    spacer_length = constraint["spacer_length"]
    upstream_idx = library.index(upstream)
    downstream_idx = library.index(downstream)
    upstream_offset = sol.offsets_fwd[upstream_idx]
    downstream_offset = sol.offsets_fwd[downstream_idx]
    assert upstream_offset is not None
    assert downstream_offset is not None
    assert upstream_pos[0] <= upstream_offset <= upstream_pos[1]
    spacer = downstream_offset - upstream_offset - len(upstream)
    assert spacer_length[0] <= spacer <= spacer_length[1]


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


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("ATGC", "GCAT"),
        ("A-T", "A-T"),
    ],
)
def test_reverse_complement(sequence: str, expected: str):
    assert da.optimize.reverse_complement(sequence) == expected


def test_adjacency_matrix():
    motifs = ["AAA", "AAT"]
    adj = da.optimize.adjacency_matrix(motifs)
    assert len(adj) == 2
    assert all(len(row) == 2 for row in adj)
    assert adj[0][1] == da.optimize.shift_metric("AAA", "AAT")


def test_dispatch_labels_simple():
    library = ["ATG", "CCC"]
    offsets = [0, 3]
    lines = da.optimize.dispatch_labels(library, offsets, rev=False)
    assert lines == ["ATG", "   CCC"]

    lines_rev = da.optimize.dispatch_labels(library, offsets, rev=True)
    assert lines_rev == ["GTA", "   CCC"]


def test_densearray_sequence_and_ordering():
    library = ["AAA", "CC"]
    sol = da.DenseArray(
        library=library,
        sequence_length=4,
        offsets_fwd=[0, 0],
        offsets_rev=[None, None],
    )
    assert sol.offset_indices_in_order() == [(0, 1), (0, 0)]
    assert sol.sequence.startswith("AAA")


def test_densearray_str_and_compression_ratio():
    library = ["AAA", "CCC"]
    sol = da.DenseArray(
        library=library,
        sequence_length=6,
        offsets_fwd=[0, 3],
        offsets_rev=[None, None],
    )
    assert sol.compression_ratio == 1.0
    rendered = str(sol)
    assert "-->" in rendered
    assert "<--" in rendered
    assert "AAACCC" in rendered


def test_invalid_motif_inputs():
    with pytest.raises(ValueError, match="non-empty"):
        da.Optimizer([""], sequence_length=5, strands="single")
    with pytest.raises(ValueError, match="invalid bases"):
        da.Optimizer(["ATGN"], sequence_length=5, strands="single")


@pytest.mark.parametrize(
    "sequence_length, strands, library, match",
    [
        (0, "single", ["AAA"], "sequence_length"),
        (-1, "single", ["AAA"], "sequence_length"),
        (5, "invalid", ["AAA"], "strands"),
        (5, "single", [], "library"),
        (5, "single", ["aaa"], "invalid bases"),
    ],
)
def test_optimizer_validation(
    sequence_length: int, strands: str, library: list[str], match: str
):
    with pytest.raises(ValueError, match=match):
        da.Optimizer(library, sequence_length=sequence_length, strands=strands)


def test_build_model_invalid_solver():
    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(RuntimeError, match="Could not create model"):
        opt.build_model(solver="NOT_A_SOLVER")


def test_requires_model_for_advanced_api():
    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(RuntimeError, match="Model not built"):
        opt.solve()
    with pytest.raises(RuntimeError, match="Model not built"):
        opt.forbid(
            da.DenseArray(
                library=["AAA", "CCC"],
                sequence_length=6,
                offsets_fwd=[0, None],
                offsets_rev=[None, None],
            )
        )
    with pytest.raises(RuntimeError, match="Model not built"):
        opt.set_motif_weight(0, 1.0)


def test_add_side_biases_invalid_motif():
    opt = da.Optimizer(["AAA"], sequence_length=3, strands="single")
    with pytest.raises(ValueError, match="All motifs must belong"):
        opt.add_side_biases(left=["AAA", "CCC"])


def test_regulator_constraints_invalid_inputs():
    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="length"):
        opt.add_regulator_constraints(["R1"], required={"R1"})

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="cover all motif indices"):
        opt.add_regulator_constraints({0: "R1"}, required={"R1"})

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="non-empty"):
        opt.add_regulator_constraints(["R1", ""], required={"R1"})

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="At least one regulator constraint"):
        opt.add_regulator_constraints(["R1", "R2"])

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="missing"):
        opt.add_regulator_constraints(["R1", "R2"], required={"R3"})

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="min_count_by_regulator"):
        opt.add_regulator_constraints(["R1", "R2"], min_count_by_regulator={"R1": 0})

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="not in mapping"):
        opt.add_regulator_constraints(["R1", "R2"], min_count_by_regulator={"R9": 1})

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="cannot satisfy min_count"):
        opt.add_regulator_constraints(["R1", "R2"], min_count_by_regulator={"R1": 2})

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="min_required_regulators must be > 0"):
        opt.add_regulator_constraints(["R1", "R2"], min_required_regulators=0)

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    with pytest.raises(ValueError, match="exceeds available regulators"):
        opt.add_regulator_constraints(["R1", "R2"], min_required_regulators=3)

    opt = da.Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
    opt.add_regulator_constraints(["R1", "R2"], required={"R1"})
    with pytest.raises(RuntimeError, match="already set"):
        opt.add_regulator_constraints(["R1", "R2"], required={"R2"})


def test_regulator_constraints_double_strand():
    motifs = ["AAA", "CCC"]
    regulators = ["R1", "R2"]
    opt = da.Optimizer(motifs, sequence_length=6, strands="double")
    opt.add_regulator_constraints(regulators, min_required_regulators=1)
    sol = opt.optimal()
    used = _used_indices(sol)
    assert len(used) >= 1


def test_approximate_solution_basic():
    opt = da.Optimizer(["AAA", "CCC", "GGG"], sequence_length=6, strands="single")
    sol = opt.approximate()
    assert len(sol.sequence) == 6
    for motif, offset in zip(opt.library, sol.offsets_fwd, strict=True):
        if offset is not None:
            assert sol.sequence[offset : offset + len(motif)] == motif


def test_take_best_run_matches_maximum():
    library = ["AAA", "CCC"]
    sequence = "AAAXXCCC"
    seq_len = 5
    best = da.optimize.take_best_run(sequence, seq_len, library, strands="single")

    def count_motifs(subseq: str) -> int:
        return sum(motif in subseq for motif in library)

    all_subseq = [sequence[i : i + seq_len] for i in range(len(sequence) - seq_len + 1)]
    max_count = max(count_motifs(s) for s in all_subseq)
    assert count_motifs(best) == max_count


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
    "strands, sequence_length, noprom",
    [
        ("single", 10, 4),
        ("double", 8, 4),
    ],
)
def test_promoter_constraints(strands: str, sequence_length: int, noprom: int):
    opt = da.Optimizer(
        ["GCA", "CCC", "ATGC", "CATT"], sequence_length=sequence_length, strands=strands
    )
    sol_noprom = opt.optimal()
    assert sol_noprom.nb_motifs == noprom
    opt.add_promoter_constraints(
        upstream="ATGC", downstream="CCC", upstream_pos=(0, 2), spacer_length=(0, 3)
    )
    sol_prom = opt.optimal()
    _assert_promoter_constraint(
        sol_prom,
        ["GCA", "CCC", "ATGC", "CATT"],
        {
            "upstream": "ATGC",
            "downstream": "CCC",
            "upstream_pos": (0, 2),
            "spacer_length": (0, 3),
        },
    )


def test_multiple_promoter_constraints():
    opt = da.Optimizer(
        ["GCA", "CCC", "ACACT", "ATGC", "CATT", "ACACT", "GCA", "AAA"],
        sequence_length=10,
        strands="double",
    )
    opt.add_promoter_constraints(
        upstream="ATGC", downstream="CCC", upstream_pos=(0, 2), spacer_length=(0, 3)
    )
    with pytest.raises(ValueError, match="must appear several times"):
        opt.add_promoter_constraints(
            upstream="ATGC", downstream="GCA", upstream_pos=(0, 2), spacer_length=(0, 3)
        )
    with pytest.raises(ValueError, match="must appear several times"):
        opt.add_promoter_constraints(
            upstream="GCA", downstream="ATGC", upstream_pos=(0, 2), spacer_length=(0, 3)
        )
    with pytest.raises(ValueError, match="must appear several times"):
        opt.add_promoter_constraints(
            upstream="CATT", downstream="CCC", upstream_pos=(0, 2), spacer_length=(0, 3)
        )
    with pytest.raises(ValueError, match="must appear several times"):
        opt.add_promoter_constraints(
            upstream="CCC", downstream="CATT", upstream_pos=(0, 2), spacer_length=(0, 3)
        )
    opt.add_promoter_constraints(
        upstream="GCA", downstream="CATT", upstream_pos=(0, 2), spacer_length=(0, 3)
    )
    opt.add_promoter_constraints(
        upstream="ACACT", downstream="ACACT", upstream_pos=(0, 2), spacer_length=(0, 3)
    )
    opt.add_promoter_constraints(
        upstream="GCA", downstream="AAA", upstream_pos=(0, 2), spacer_length=(0, 3)
    )
    assert len(opt.promoters) == 4
    upstream = {p.upstream_index for p in opt.promoters}
    downstream = {p.downstream_index for p in opt.promoters}
    assert len(upstream | downstream) == 2 * len(opt.promoters)


def test_multiple_constraints_optimization():
    opt = da.Optimizer(
        ["AAA", "CCC", "AAA", "CCC"], sequence_length=12, strands="double"
    )
    opt.add_promoter_constraints(
        upstream="AAA", downstream="CCC", upstream_pos=(0, 0), spacer_length=(0, 0)
    )
    opt.add_promoter_constraints(
        upstream="AAA", downstream="CCC", upstream_pos=(6, 6), spacer_length=(0, 0)
    )
    sol = opt.optimal()
    assert sol.offset_indices_in_order() == [(0, 0), (3, 1), (6, 2), (9, 3)]


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

    opt.add_promoter_constraints(
        upstream="GGGT", downstream="CTTC", upstream_pos=(0, 3), spacer_length=(1, 4)
    )

    opt.add_side_biases(left=library[::2], right=library[1::2])
    sol_left = opt.optimal()
    assert sol_left.nb_motifs >= 3
    _assert_promoter_constraint(
        sol_left,
        library,
        {
            "upstream": "GGGT",
            "downstream": "CTTC",
            "upstream_pos": (0, 3),
            "spacer_length": (1, 4),
        },
    )

    opt.add_side_biases(left=library[1::2], right=library[::2])
    sol_right = opt.optimal()
    assert sol_right.nb_motifs >= 3
    _assert_promoter_constraint(
        sol_right,
        library,
        {
            "upstream": "GGGT",
            "downstream": "CTTC",
            "upstream_pos": (0, 3),
            "spacer_length": (1, 4),
        },
    )


def _used_indices(sol: da.DenseArray) -> set[int]:
    return {
        i
        for i, (fwd, rev) in enumerate(
            zip(sol.offsets_fwd, sol.offsets_rev, strict=True),
        )
        if fwd is not None or rev is not None
    }


def test_regulator_constraints_required_feasible():
    motifs = ["AAA", "CCC"]
    regulators = ["R1", "R2"]
    opt = da.Optimizer(motifs, sequence_length=6, strands="single")
    opt.add_regulator_constraints(regulators, required={"R1", "R2"})
    sol = opt.optimal()
    used = _used_indices(sol)
    assert used == {0, 1}


def test_regulator_constraints_required_infeasible():
    motifs = ["AAA", "CCC"]
    regulators = ["R1", "R2"]
    opt = da.Optimizer(motifs, sequence_length=3, strands="single")
    opt.add_regulator_constraints(regulators, required={"R1", "R2"})
    with pytest.raises(ValueError, match="feasible"):
        opt.optimal()


def test_regulator_constraints_k_of_n():
    motifs = ["AAA", "CCC", "GGG"]
    regulators = ["R1", "R2", "R3"]
    opt = da.Optimizer(motifs, sequence_length=6, strands="single")
    opt.add_regulator_constraints(regulators, min_required_regulators=2)
    sol = opt.optimal()
    used = _used_indices(sol)
    covered = {regulators[i] for i in used}
    assert len(covered) >= 2

    opt2 = da.Optimizer(motifs, sequence_length=3, strands="single")
    opt2.add_regulator_constraints(regulators, min_required_regulators=2)
    with pytest.raises(ValueError, match="feasible"):
        opt2.optimal()


def test_regulator_constraints_min_count():
    motifs = ["AAA", "AAT", "CCC"]
    regulators = ["R1", "R1", "R2"]
    opt = da.Optimizer(motifs, sequence_length=6, strands="single")
    opt.add_regulator_constraints(regulators, min_count_by_regulator={"R1": 2})
    sol = opt.optimal()
    used = _used_indices(sol)
    counts = Counter(regulators[i] for i in used)
    assert counts["R1"] >= 2
