"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/background/compiler.py

Compile native sequence constraints into a bounded conditional distribution.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from dataclasses import dataclass
from fractions import Fraction
from math import gcd, lcm

from dense_arrays._record_validation import integer, semantic_digest
from dense_arrays.constraints import GC, Avoid
from dense_arrays.parts.motifs.models import numeric_row
from dense_arrays.sequence import reverse_complement

from .automaton import transitions
from .contracts import (
    CONDITIONAL_POLICY,
    ConditionalLimits,
    ConstructionLimitError,
    ConstructionReport,
    ConstructionWork,
)
from .counting import CompletionTable, CountedBackground, state


def integer_weights(probabilities: tuple[float, ...]) -> tuple[int, ...]:
    """Preserve normalized decimal probability ratios with exact integer masses."""
    fractions = tuple(
        Fraction(str(p)) for p in numeric_row(probabilities, probability=True)
    )
    denominator = lcm(*(p.denominator for p in fractions))
    weights = tuple(p.numerator * (denominator // p.denominator) for p in fractions)
    divisor = gcd(*weights)
    return tuple(w // divisor for w in weights)


def model_identity(
    *,
    minimum: int,
    maximum: int,
    probabilities: tuple[float, ...],
    screening: tuple[Avoid | GC, ...],
) -> str:
    """Bind distribution and constraints without including labels or effort caps."""
    integer(minimum, field_name="conditional minimum", minimum=1)
    integer(maximum, field_name="conditional maximum", minimum=minimum)
    for rule in screening:
        if (
            not isinstance(rule, (Avoid, GC))
            or (isinstance(rule, Avoid) and rule.except_placements)
            or (isinstance(rule, GC) and rule.scope != "sequence")
        ):
            msg = (
                "conditional sampling compiles sequence GC and Avoid without exceptions"
            )
            raise ValueError(msg)
    return semantic_digest(
        {
            "policy": CONDITIONAL_POLICY,
            "minimum": minimum,
            "maximum": maximum,
            "weights": integer_weights(probabilities),
            "gc": [(r.min, r.max) for r in screening if isinstance(r, GC)],
            "avoid": [
                (r.patterns, r.strands) for r in screening if isinstance(r, Avoid)
            ],
        }
    )


def _gc_bounds(length: int, rules: tuple[Avoid | GC, ...]) -> tuple[int, int]:
    """Match the native count/length predicate without rounded multiplication."""
    minimum = max((r.min for r in rules if isinstance(r, GC)), default=0)
    maximum = min((r.max for r in rules if isinstance(r, GC)), default=1)
    low, high = 0, length + 1
    while low < high:
        middle = (low + high) // 2
        if middle / length < minimum:
            low = middle + 1
        else:
            high = middle
    first = low
    low, high = 0, length + 1
    while low < high:
        middle = (low + high) // 2
        if middle / length <= maximum:
            low = middle + 1
        else:
            high = middle
    return first, low - 1


@dataclass(frozen=True)
class BackgroundCompilation:
    """A qualified sampler or an explicit infeasible/limited construction result."""

    report: ConstructionReport
    sampler: CountedBackground | None


def compile_background(  # noqa: PLR0913 - explicit pure model and execution limits
    *,
    minimum: int,
    maximum: int,
    probabilities: tuple[float, ...],
    screening: tuple[Avoid | GC, ...],
    limits: ConditionalLimits | None = None,
    deadline: float | None = None,
) -> BackgroundCompilation:
    """Count valid suffixes once, conditioning the uniform prior over lengths."""
    model_id = model_identity(
        minimum=minimum,
        maximum=maximum,
        probabilities=probabilities,
        screening=screening,
    )
    limits = limits or ConditionalLimits()
    if not isinstance(limits, ConditionalLimits):
        msg = "conditional construction requires ConditionalLimits"
        raise TypeError(msg)
    expires = time.monotonic() + limits.seconds
    work = ConstructionWork(
        limits, min(expires, deadline) if deadline is not None else expires
    )
    sampler, mass, status, reason = None, None, "limited", "time_budget"
    try:
        patterns = set()
        for rule in screening:
            work.check_time()
            if isinstance(rule, Avoid):
                for word in rule.patterns:
                    work.check_time()
                    if len(word) > maximum:
                        continue
                    patterns.add(word)
                    if rule.strands == "both":
                        patterns.add(reverse_complement(word))
        weights = integer_weights(probabilities)
        table = CompletionTable(
            transitions(tuple(sorted(patterns)), work), weights, work
        )
        lengths, total = [], 0
        for length in range(minimum, maximum + 1):
            work.admit("states")
            root = state(length, 0, *_gc_bounds(length, screening))
            count = table.count(root)
            scale = work.power(sum(weights), maximum - length) if count else 1
            weighted = work.multiply(scale, count)
            total += weighted
            work.admit(
                "mass_bits", max(1, scale.bit_length()) + max(1, weighted.bit_length())
            )
            lengths.append((root, scale, weighted))
        work.admit("mass_bits", max(1, total.bit_length()))
        mass, status, reason = (
            hex(total),
            "feasible" if total else "infeasible",
            "completed",
        )
        if total:
            sampler = CountedBackground(model_id, table, tuple(lengths), total)
    except ConstructionLimitError as err:
        reason = err.reason
    return BackgroundCompilation(
        ConstructionReport(
            model_id,
            status,
            reason,
            work.states,
            work.automaton_states,
            work.mass_bits,
            mass,
        ),
        sampler,
    )
