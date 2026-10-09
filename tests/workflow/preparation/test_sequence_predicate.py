"""Sequence acceptance shares detailed-screen semantics without collecting hits.

Author: Eric J. South.
"""

from collections.abc import Iterator
from itertools import product

import pytest

from dense_arrays.constraints import GC, Avoid
from dense_arrays.parts import Part
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.eligibility import rejection_reasons
from dense_arrays.parts.sampling import Eligibility
from dense_arrays.parts.screening import sequence as screens


@pytest.mark.parametrize(
    "rule",
    [
        Avoid("forward", ("AA", "ACG"), strands="forward"),
        Avoid("both", ("AA", "ACG"), strands="both"),
        Avoid("palindrome", ("AT", "CG"), strands="both"),
        Avoid("long", ("AAAAAA",)),
        GC("gc", "sequence", 0.25, 0.75),
        GC("exact", "sequence", 0.5, 0.5),
    ],
)
def test_sequence_predicate_matches_full_observations(rule: Avoid | GC):
    for length in range(6):
        for bases in product("ACGT", repeat=length):
            sequence = "".join(bases)
            assert (
                screens.passes_sequence(rule, sequence)
                == screens.evaluate_sequence(rule, sequence)["passed"]
            )


def test_eligibility_does_not_materialize_discarded_match_observations(
    monkeypatch: pytest.MonkeyPatch,
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("eligibility collected all match observations")

    monkeypatch.setattr(screens, "_violations", forbidden)
    assert rejection_reasons(
        Candidate(1, Part("candidate_1", "AAAAA", "background")),
        Eligibility(),
        (Avoid("overlap", ("AA",)),),
        requires_hit=False,
    ) == ("overlap",)


def test_predicate_stops_at_first_violation(monkeypatch: pytest.MonkeyPatch):
    def first_only(*_args: object, **_kwargs: object) -> Iterator[dict[str, int]]:
        yield {"start": 0}
        pytest.fail("boolean acceptance requested another match")

    monkeypatch.setattr(screens, "_iter_violations", first_only)
    assert not screens.passes_sequence(Avoid("word", ("AA",)), "AAAAA")


def test_detailed_screen_retains_overlaps_sorted_and_palindromes_once():
    rule = Avoid("words", ("AAA", "AA", "AT"))
    observed = screens.evaluate_sequence(rule, "AAAATTTT")["observed"]
    expected = [
        ("AA", "forward", 0, 2),
        ("AAA", "forward", 0, 3),
        ("AA", "forward", 1, 3),
        ("AAA", "forward", 1, 4),
        ("AA", "forward", 2, 4),
        ("AT", "forward", 3, 5),
        ("AA", "reverse", 4, 6),
        ("AAA", "reverse", 4, 7),
        ("AA", "reverse", 5, 7),
        ("AAA", "reverse", 5, 8),
        ("AA", "reverse", 6, 8),
    ]
    assert [
        (match["pattern"], match["strand"], match["start"], match["end"])
        for match in observed
    ] == expected
    assert all(match["intersections"] == [] for match in observed)


@pytest.mark.parametrize(
    "rule,message",
    [
        (Avoid("except", ("AA",), except_placements=("fixed",)), "exceptions"),
        (GC("padding", "padding", 0, 1), "sequence scope"),
    ],
)
def test_predicate_keeps_sequence_scope_errors(rule: Avoid | GC, message: str):
    for operation in (screens.passes_sequence, screens.evaluate_sequence):
        with pytest.raises(ValueError, match=message):
            operation(rule, "AAAA")


def test_gc_exact_fraction_and_empty_sequence_boundaries():
    rule = GC("gc", "sequence", 0.14, 0.14)
    assert screens.passes_sequence(rule, "G" * 7 + "A" * 43)
    assert not screens.passes_sequence(rule, "G" * 8 + "A" * 42)
    assert screens.passes_sequence(rule, "")
    assert screens.evaluate_sequence(rule, "")["status"] == "not_applicable"
