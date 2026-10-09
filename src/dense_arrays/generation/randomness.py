"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/randomness.py

Versioned random DNA streams bound to logical work, never execution ordering.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

import hashlib

from dense_arrays._record_validation import canonical_json, integer, semantic_digest

PADDING_POLICY = "single_side_uniform_shake256.v1"


def padding_dna(length: int, *, seed: int, attempt: int, trial: int) -> tuple[str, str]:
    """Expand a domain-separated SHAKE-256 stream with four ordered bases per byte."""
    integer(length, field_name="padding length", minimum=0)
    identity = {
        "policy": PADDING_POLICY,
        "seed": seed,
        "cell": "default",
        "batch": 1,
        "attempt": attempt,
        "trial": trial,
    }
    stream = hashlib.shake_256(canonical_json(identity).encode()).digest(
        (length + 3) // 4
    )
    bases = "".join(
        "ACGT"[(byte >> shift) & 3] for byte in stream for shift in (6, 4, 2, 0)
    )
    return bases[:length], semantic_digest(identity)
