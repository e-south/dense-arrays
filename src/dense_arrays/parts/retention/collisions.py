"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/retention/collisions.py

Compare retained sequence and observed core identities across named recipes.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dense_arrays.parts.candidates import Candidate


def validate_collisions(
    candidates: tuple[Candidate, ...], *, sequences: str, cores: str
) -> None:
    """Reject declared cross-recipe collisions without changing any selection.

    Cores use their recorded motif orientation. Missing cores are excluded;
    reverse complements are not collapsed and scores are not compared.
    """
    comparisons = []
    if sequences == "error":
        comparisons.append(("sequence", {}))
    if cores == "error":
        comparisons.append(("core", {}))
    if not comparisons:
        return
    for candidate in candidates:
        if not candidate.retained:
            continue
        for kind, seen in comparisons:
            key = (
                candidate.part.sequence
                if kind == "sequence"
                else candidate.part.core_sequence
            )
            if key is None:
                continue
            prior = seen.get(key)
            if prior is not None and prior.recipe_id != candidate.recipe_id:
                msg = (
                    f"retained {kind} occurs in multiple recipes: "
                    f"{prior.recipe_id!r} ({prior.part.part_id}) and "
                    f"{candidate.recipe_id!r} ({candidate.part.part_id}); "
                    f"revise the recipes or set {kind}_collisions=preserve "
                    "to keep distinct occurrences"
                )
                raise ValueError(msg)
            seen[key] = candidate
