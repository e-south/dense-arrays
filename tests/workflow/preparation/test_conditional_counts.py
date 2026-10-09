"""Independent sequence oracles for conditional background proposals.

Author: Eric J. South.
"""

from collections import Counter
from itertools import product
from math import prod

import pytest
from ortools.sat.python import cp_model

from dense_arrays.constraints import GC, Avoid
from dense_arrays.parts.background import ConditionalLimits, compile_background


def test_uniform_mass_and_each_rank_match_the_complete_valid_language():
    rules = (GC("gc", "sequence", 0.25, 0.75), Avoid("words", ("ACA", "CC")))
    expected = [
        seq
        for bases in product("ACGT", repeat=4)
        if 1 <= (seq := "".join(bases)).count("C") + seq.count("G") <= 3
        and not any(word in seq for word in ("ACA", "TGT", "CC", "GG"))
    ]
    result = compile_background(
        minimum=4, maximum=4, probabilities=(0.25,) * 4, screening=rules
    )
    assert result.report.status == "feasible"
    assert result.sampler.mass == len(expected)
    assert [result.sampler.unrank(i) for i in range(len(expected))] == expected


@pytest.mark.parametrize(
    "weights,length,gc,words",
    [
        ((1, 2, 3, 4), 3, (1, 2), ("AA",)),
        ((2, 1, 1, 0), 4, (1, 3), ("CG", "TA")),
        ((0, 1, 1, 0), 4, (0, 4), ("CC", "GG")),
    ],
)
def test_every_weighted_rank_has_exact_background_multiplicity(
    weights: tuple[int, ...], length: int, gc: tuple[int, int], words: tuple[str, ...]
):
    expected = Counter()
    for bases in product("ACGT", repeat=length):
        seq = "".join(bases)
        count = seq.count("C") + seq.count("G")
        mass = prod(weights["ACGT".index(b)] for b in seq)
        if gc[0] <= count <= gc[1] and not any(w in seq for w in words) and mass:
            expected[seq] = mass
    result = compile_background(
        minimum=length,
        maximum=length,
        probabilities=tuple(w / sum(weights) for w in weights),
        screening=(
            GC("gc", "sequence", gc[0] / length, gc[1] / length),
            Avoid("words", words, strands="forward"),
        ),
    )
    assert result.sampler.mass == sum(expected.values())
    assert (
        Counter(result.sampler.unrank(i) for i in range(result.sampler.mass))
        == expected
    )


def test_ranged_lengths_condition_the_original_uniform_length_prior():
    result = compile_background(
        minimum=1,
        maximum=2,
        probabilities=(0.25,) * 4,
        screening=(Avoid("word", ("AA",), strands="forward"),),
    )
    expected = Counter(dict.fromkeys("ACGT", 4))
    expected.update("".join(b) for b in product("ACGT", repeat=2) if b != ("A", "A"))
    assert result.sampler.mass == 31
    assert Counter(result.sampler.unrank(i) for i in range(31)) == expected


def test_gc_integer_boundary_matches_the_declared_fraction():
    result = compile_background(
        minimum=50,
        maximum=50,
        probabilities=(0.25,) * 4,
        screening=(GC("gc", "sequence", 0.14, 0.14),),
    )
    assert result.report.status == "feasible"
    for index in range(1, 11):
        seq = result.sampler.draw(seed=7, index=index)
        assert (seq.count("C") + seq.count("G")) / len(seq) == 0.14


@pytest.mark.parametrize(
    "probabilities,rules",
    [
        ((1, 0, 0, 0), (GC("gc", "sequence", 0.1, 1),)),
        ((0.25,) * 4, (Avoid("none", tuple("ACGT")),)),
        ((0.25,) * 4, (GC("a", "sequence", 0, 0.25), GC("b", "sequence", 0.75, 1))),
    ],
)
def test_zero_mass_proves_impossible_support(
    probabilities: tuple[float, ...], rules: tuple
):
    result = compile_background(
        minimum=2, maximum=4, probabilities=probabilities, screening=rules
    )
    assert result.sampler is None
    assert result.report.status == "infeasible"
    assert result.report.mass == "0x0"


@pytest.mark.parametrize(
    "limits,reason",
    [
        (ConditionalLimits(states=1), "states"),
        (ConditionalLimits(automaton_states=1), "automaton_states"),
        (ConditionalLimits(mass_bits=1), "mass_bits"),
    ],
)
def test_resource_limits_do_not_claim_infeasibility(
    limits: ConditionalLimits, reason: str
):
    result = compile_background(
        minimum=20,
        maximum=20,
        probabilities=(0.25,) * 4,
        screening=(Avoid("word", ("AAA",)),),
        limits=limits,
    )
    assert result.sampler is None
    assert result.report.status == "limited"
    assert result.report.reason == reason
    assert result.report.mass is None


def test_expired_cooperative_deadline_does_not_claim_infeasibility():
    result = compile_background(
        minimum=20,
        maximum=20,
        probabilities=(0.25,) * 4,
        screening=(),
        deadline=0,
    )
    assert result.report.status == "limited"
    assert result.report.reason == "time_budget"


def test_candidate_stream_is_independent_of_draw_order_and_effort_caps():
    args = {
        "minimum": 4,
        "maximum": 8,
        "probabilities": (0.1, 0.2, 0.3, 0.4),
        "screening": (),
    }
    first = compile_background(**args).sampler
    second = compile_background(
        **args, limits=ConditionalLimits(states=100_000)
    ).sampler
    sequences = [first.draw(seed=7, index=i) for i in range(1, 21)]
    assert [second.draw(seed=7, index=i) for i in range(20, 0, -1)] == sequences[::-1]
    assert [first.draw(seed=8, index=i) for i in range(1, 21)] != sequences


def test_long_sequences_do_not_depend_on_python_recursion_limit():
    result = compile_background(
        minimum=1500, maximum=1500, probabilities=(1, 0, 0, 0), screening=()
    )
    assert result.sampler.unrank(0) == "A" * 1500


@pytest.mark.parametrize("cap,expected", [(5, "limited"), (6, "feasible")])
def test_shared_suffixes_are_admitted_once(cap: int, expected: str):
    result = compile_background(
        minimum=4,
        maximum=4,
        probabilities=(0.25,) * 4,
        screening=(),
        limits=ConditionalLimits(states=cap),
    )
    # One length choice plus five suffix states; repeated bases share suffixes.
    assert result.report.states == cap
    assert result.report.status == expected
    if result.sampler:
        assert result.sampler.mass == 256
        assert result.sampler.unrank(0) == "AAAA"
        assert result.sampler.unrank(255) == "TTTT"


@pytest.mark.parametrize(
    "strands,expected", [("forward", "feasible"), ("both", "infeasible")]
)
def test_strand_scope_and_zero_probability_bases(strands: str, expected: str):
    result = compile_background(
        minimum=2,
        maximum=2,
        probabilities=(0, 0, 0, 1),
        screening=(Avoid("a", ("A",), strands=strands),),
    )
    assert result.report.status == expected
    if result.sampler:
        assert result.sampler.unrank(0) == "TT"


def test_cp_sat_clause_encoding_agrees_with_native_exact_support():
    model = cp_model.CpModel()
    bases = [
        [model.new_bool_var(f"{pos}_{base}") for base in "ACGT"] for pos in range(4)
    ]
    for row in bases:
        model.add_exactly_one(row)
    model.add_linear_constraint(sum(row[1] + row[2] for row in bases), 1, 3)
    for word in ("ACA", "TGT", "CC", "GG"):
        for start in range(5 - len(word)):
            model.add(
                sum(bases[start + i]["ACGT".index(b)] for i, b in enumerate(word))
                < len(word)
            )
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 5
    solver.parameters.num_search_workers = 1
    solver.parameters.enumerate_all_solutions = True

    class Collect(cp_model.CpSolverSolutionCallback):
        def __init__(self) -> None:
            super().__init__()
            self.sequences = set()

        def on_solution_callback(self) -> None:
            self.sequences.add(
                "".join(
                    "ACGT"[next(i for i, b in enumerate(row) if self.value(b))]
                    for row in bases
                )
            )

    found = Collect()
    assert solver.solve(model, found) == cp_model.OPTIMAL
    native = compile_background(
        minimum=4,
        maximum=4,
        probabilities=(0.25,) * 4,
        screening=(GC("gc", "sequence", 0.25, 0.75), Avoid("words", ("ACA", "CC"))),
    ).sampler
    assert {native.unrank(i) for i in range(native.mass)} == found.sequences


@pytest.mark.parametrize("value", [0, -1, True, 1.5])
def test_limits_reject_nonpositive_or_noninteger_work_caps(value: object):
    with pytest.raises((TypeError, ValueError)):
        ConditionalLimits(states=value)


def test_later_gc_rule_cannot_override_an_earlier_rule():
    first = GC("first", "sequence", 0.25, 0.75)
    second = GC("second", "sequence", 0.5, 1)
    result = compile_background(
        minimum=4, maximum=4, probabilities=(0.25,) * 4, screening=(first, second)
    )
    observed = {result.sampler.unrank(i) for i in range(result.sampler.mass)}
    expected = {
        "".join(seq)
        for seq in product("ACGT", repeat=4)
        if 2 <= sum(b in "GC" for b in seq) <= 3
    }
    assert observed == expected


def test_patterns_longer_than_every_candidate_need_no_automaton_states():
    result = compile_background(
        minimum=4,
        maximum=4,
        probabilities=(0.25,) * 4,
        screening=(Avoid("long", ("A" * 1000,)),),
        limits=ConditionalLimits(automaton_states=1),
    )
    assert result.report.status == "feasible"
    assert result.sampler.mass == 256
    assert result.report.automaton_states == 1
