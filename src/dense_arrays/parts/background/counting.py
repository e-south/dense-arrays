"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/background/counting.py

Exact completion masses and integer unranking for supported DNA strings.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
import time
from dataclasses import dataclass

from dense_arrays._record_validation import canonical_json, integer

from .contracts import CONDITIONAL_POLICY, ConstructionLimitError, ConstructionWork

State = tuple[int, int, int, int]


def state(remaining: int, node: int, low: int, high: int) -> State:
    """Coalesce equivalent suffix intervals while retaining impossible bounds."""
    return remaining, node, max(0, low), min(remaining, high)


class CompletionTable:
    """One execution's memoized masses, shared across length choices and draws."""

    def __init__(
        self,
        rows: tuple[tuple[int, ...], ...],
        weights: tuple[int, ...],
        work: ConstructionWork,
    ) -> None:
        """Keep construction state private to the owning preparation execution."""
        self.rows, self.weights, self.work = rows, weights, work
        self.masses: dict[State, int] = {}

    def children(self, key: State) -> tuple[tuple[int, int, State], ...]:
        """Return positive-support transitions in canonical ACGT order."""
        remaining, node, low, high = key
        return tuple(
            (
                base,
                weight,
                state(
                    remaining - 1,
                    child,
                    low - int(base in (1, 2)),
                    high - int(base in (1, 2)),
                ),
            )
            for base, (child, weight) in enumerate(
                zip(self.rows[node], self.weights, strict=True)
            )
            if child >= 0 and weight
        )

    def count(self, root: State) -> int:
        """Evaluate an acyclic suffix graph without recursion-depth dependence."""
        pending = [(root, False)]
        while pending:
            self.work.check_time()
            key, expanded = pending.pop()
            if key in self.masses:
                continue
            # Children reduce remaining length. DFS completes a state before
            # another queued copy reaches it, so only its first frame is new.
            if not expanded:
                self.work.admit("states")
            remaining, _, low, high = key
            if low > high:
                mass = 0
            elif remaining == 0:
                mass = 1
            elif not expanded:
                pending.append((key, True))
                pending.extend(
                    (child, False)
                    for _, _, child in self.children(key)
                    if child not in self.masses
                )
                continue
            else:
                mass = sum(
                    self.work.multiply(weight, self.masses[child])
                    for _, weight, child in self.children(key)
                )
                self.work.check_number(mass)
            self.work.admit("mass_bits", max(1, mass.bit_length()))
            self.masses[key] = mass
        return self.masses[root]


@dataclass(frozen=True)
class CountedBackground:
    """An exact conditional distribution ready for independent candidate draws."""

    model_id: str
    table: CompletionTable
    lengths: tuple[tuple[State, int, int], ...]
    mass: int

    def unrank(self, rank: int, *, deadline: float | None = None) -> str:
        """Map each unit of integer mass to its sequence without sampling bias."""
        integer(rank, field_name="conditional rank", minimum=0)
        if rank >= self.mass:
            msg = "conditional rank exceeds total completion mass"
            raise ValueError(msg)
        for root, scale, mass in self.lengths:
            if rank < mass:
                key = root
                rank //= scale
                break
            rank -= mass
        sequence = []
        while key[0]:
            if deadline is not None and time.monotonic() >= deadline:
                reason = "time_budget"
                raise ConstructionLimitError(reason)
            for base, weight, child in self.table.children(key):
                branch = weight * self.table.masses[child]
                if rank < branch:
                    rank //= weight
                    sequence.append("ACGT"[base])
                    key = child
                    break
                rank -= branch
        return "".join(sequence)

    def draw(self, *, seed: int, index: int, deadline: float | None = None) -> str:
        """Draw an unbiased mass rank with candidate-local SHAKE256 entropy."""
        integer(seed, field_name="conditional seed", minimum=0)
        integer(index, field_name="conditional candidate", minimum=1)
        bits = (self.mass - 1).bit_length()
        size = (bits + 7) // 8
        for attempt in range(128):
            if deadline is not None and time.monotonic() >= deadline:
                reason = "time_budget"
                raise ConstructionLimitError(reason)
            encoded = canonical_json(
                {
                    "policy": CONDITIONAL_POLICY,
                    "model_id": self.model_id,
                    "seed": seed,
                    "candidate": index,
                    "draw": attempt,
                }
            )
            rank = int.from_bytes(
                hashlib.shake_256(encoded.encode()).digest(size), "big"
            ) & ((1 << bits) - 1)
            if rank < self.mass:
                return self.unrank(rank, deadline=deadline)
        msg = "conditional sampling exhausted its bounded entropy draws"
        raise RuntimeError(msg)
