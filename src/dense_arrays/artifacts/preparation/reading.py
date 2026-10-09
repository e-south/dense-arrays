"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/reading.py

Read verified sampled evidence once for bounded report construction.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays.artifacts.pools import POOL_DATABASE, stored_preparation
from dense_arrays.artifacts.preparation.verification import verify_sampled
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.artifacts.store import reader

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.pool_records import PoolSummary
    from dense_arrays.artifacts.reading import ReadLimits, Verification
    from dense_arrays.parts.candidates import Candidate
    from dense_arrays.planning.preparation import PreparationPlan


@dataclass(frozen=True)
class VerifiedPreparation:
    """Saved decisions and the shared allowance used to verify their native joins."""

    plan: PreparationPlan
    candidates: tuple[Candidate, ...]
    verification: Verification
    budget: ReadBudget


def read_verified_sampled(
    path: Path, summary: PoolSummary, limits: ReadLimits
) -> VerifiedPreparation:
    """Retain bounded candidates after verification so reports need no second scan."""
    budget = ReadBudget(limits)
    with reader(path, filename=POOL_DATABASE) as connection:
        row = connection.execute(
            "SELECT payload FROM preparation WHERE id=1"
        ).fetchone()
        if row is None:
            msg = "sampled pool has no preparation plan"
            raise ValueError(msg)
        budget.examine(row[0])
        plan = stored_preparation(connection, max_identities=limits.identities)
        verification, candidates = verify_sampled(
            connection, path, plan, summary, budget
        )
    return VerifiedPreparation(plan, candidates, verification, budget)
