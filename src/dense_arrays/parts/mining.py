"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/mining.py

Versioned candidate-local random streams for DNA and PWM draws.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
import math
from bisect import bisect_right
from dataclasses import dataclass
from itertools import accumulate
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json, integer, semantic_digest
from dense_arrays.parts.background.contracts import CONDITIONAL_POLICY
from dense_arrays.parts.motifs.models import numeric_row

if TYPE_CHECKING:
    from dense_arrays.parts.motifs import Motif

POLICY = "candidate_shake256.v1"
POLICIES = {
    "stochastic": POLICY,
    "consensus": "consensus_embed_shake256.v1",
    "background": "background_draw_shake256.v1",
    "conditional": CONDITIONAL_POLICY,
}


@dataclass(frozen=True)
class Proposal:
    """One constructed sequence and its intended forward motif interval, if any."""

    sequence: str
    start: int | None
    end: int | None

    def evidence(self, strategy: str) -> dict[str, object]:
        """Record construction geometry independently of later scoring results."""
        return {
            "schema": "dense_arrays.proposal.v1",
            "strategy": strategy,
            "policy": POLICIES[strategy],
            "start": self.start,
            "end": self.end,
        }


def consensus(motif: Motif) -> str:
    """Choose the highest-probability base, breaking exact ties in ACGT order."""
    return "".join(
        "ACGT"[max(range(4), key=row.__getitem__)] for row in motif.probabilities
    )


def sample_sequence(
    *,
    length: int,
    seed: int,
    index: int,
    probabilities: tuple[float, ...],
    motif: Motif | None = None,
) -> str:
    """Preserve the original stochastic interface and candidate-local stream."""
    return propose_sequence(
        length=length, seed=seed, index=index, probabilities=probabilities, motif=motif
    ).sequence


def propose_sequence(  # noqa: PLR0913 - explicit pure proposal inputs
    *,
    length: int,
    seed: int,
    index: int,
    probabilities: tuple[float, ...],
    motif: Motif | None = None,
    strategy: str = "stochastic",
) -> Proposal:
    """Draw a bounded proposal independently of eligibility, retention or scoring."""
    integer(length, field_name="proposal.length", minimum=1)
    integer(seed, field_name="proposal.seed", minimum=0)
    integer(index, field_name="proposal.index", minimum=1)
    probabilities = numeric_row(probabilities, probability=True)
    if strategy == "conditional":
        msg = "conditional proposals require a compiled background distribution"
        raise ValueError(msg)
    if strategy not in POLICIES or (strategy != "stochastic" and motif is None):
        msg = "consensus and background proposal strategies require a motif"
        raise ValueError(msg)
    if motif is not None and length < motif.width:
        msg = "proposal length cannot be shorter than the motif"
        raise ValueError(msg)
    identity = {
        "policy": POLICIES[strategy],
        "seed": seed,
        "candidate": index,
        "motif": semantic_digest(
            {"probabilities": motif.probabilities, "background": motif.background}
        )
        if motif
        else None,
    }
    stream = hashlib.shake_256(canonical_json(identity).encode()).digest(
        8 * (length + 1)
    )
    values = [
        (int.from_bytes(stream[offset : offset + 8], "big") >> 11) / 2**53
        for offset in range(0, len(stream), 8)
    ]
    embedded = motif is not None and strategy != "background"
    offset = (
        min(int(values[-1] * (length - motif.width + 1)), length - motif.width)
        if embedded
        else 0
    )
    core = consensus(motif) if strategy == "consensus" else None
    background_distribution = _distribution(probabilities)
    sequence = []
    for position, value in enumerate(values[:-1]):
        in_motif = embedded and offset <= position < offset + motif.width
        if in_motif and core is not None:
            sequence.append(core[position - offset])
        else:
            bases, cumulative = (
                _distribution(motif.probabilities[position - offset])
                if in_motif
                else background_distribution
            )
            sequence.append(bases[bisect_right(cumulative, value)])
    return Proposal(
        "".join(sequence),
        offset if embedded else None,
        offset + motif.width if embedded else None,
    )


def _distribution(weights: tuple[float, ...]) -> tuple[tuple[str, ...], list[float]]:
    """Prepare one categorical distribution for reuse within a proposal."""
    positive = [
        (base, weight)
        for base, weight in zip("ACGT", weights, strict=True)
        if weight > 0
    ]
    total = math.fsum(weight for _, weight in positive)
    cumulative = list(accumulate(weight / total for _, weight in positive))
    cumulative[-1] = 1.0
    return tuple(base for base, _ in positive), cumulative


LENGTH_POLICY = "candidate_length_shake256.v1"
_MAX_LENGTH_DRAWS = 128


def sample_length(*, minimum: int, maximum: int, seed: int, index: int) -> int:
    """Draw uniformly from inclusive integer bounds with candidate-local entropy."""
    for name, value, floor in (
        ("minimum", minimum, 1),
        ("maximum", maximum, minimum),
        ("seed", seed, 0),
        ("index", index, 1),
    ):
        integer(value, field_name=f"sampled_length.{name}", minimum=floor)
    span = maximum - minimum + 1
    if span == 1:
        return minimum
    bits = (span - 1).bit_length()
    size = (bits + 7) // 8
    for draw in range(_MAX_LENGTH_DRAWS):
        identity = {
            "policy": LENGTH_POLICY,
            "seed": seed,
            "candidate": index,
            "minimum": minimum,
            "maximum": maximum,
            "draw": draw,
        }
        value = int.from_bytes(
            hashlib.shake_256(canonical_json(identity).encode()).digest(size), "big"
        ) & ((1 << bits) - 1)
        if value < span:
            return minimum + value
    msg = "sampled length exhausted its bounded entropy draws"
    raise RuntimeError(msg)
