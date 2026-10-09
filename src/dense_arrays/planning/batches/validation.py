"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/batches/validation.py

Eligibility and membership checks shared by preparation and artifact readers.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from collections import Counter
from collections.abc import Sequence

from dense_arrays.parts import Part

from .models import BatchSampling


def validate_eligible(parts: Sequence[Part], policy: BatchSampling) -> None:
    """Require known annotations wherever a declared policy relies on them."""
    if (
        policy.unique_cores
        or policy.max_per_group is not None
        or policy.strategy == "group_balanced"
    ) and any(p.group is None for p in parts):
        msg = "this batch policy requires a group for every eligible part"
        raise ValueError(msg)
    if policy.unique_cores and any(p.core_sequence is None for p in parts):
        msg = "core uniqueness requires core annotations for every eligible part"
        raise ValueError(msg)


def validate_selected(parts: Sequence[Part], policy: BatchSampling) -> None:
    """Check hard membership limits without drawing or optimizing again."""
    if policy.unique_sequences and len({p.sequence for p in parts}) != len(parts):
        msg = "batch membership violates sequence uniqueness"
        raise ValueError(msg)
    if policy.unique_cores and len({(p.group, p.core_sequence) for p in parts}) != len(
        parts
    ):
        msg = "batch membership violates group-scoped core uniqueness"
        raise ValueError(msg)
    if policy.max_per_group is not None and any(
        count > policy.max_per_group
        for count in Counter(p.group for p in parts).values()
    ):
        msg = "batch membership exceeds max_per_group"
        raise ValueError(msg)
