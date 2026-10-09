"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/search.py

Recorded search-method evidence, distinct from backend proof reports.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import object_fields

if TYPE_CHECKING:
    from collections.abc import Mapping


@dataclass(frozen=True)
class HeuristicEvidence:
    """One deterministic greedy result; no optimality or infeasibility proof."""

    status: str

    def __post_init__(self) -> None:
        """Keep a finished proposal, finite exhaustion and a time stop distinct."""
        if self.status not in {"candidate", "exhausted", "time_limit"}:
            msg = "unsupported greedy search status"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, str]:
        """Bind evidence to its method version."""
        return {
            "schema": "dense_arrays.heuristic.v1",
            "method": "greedy_multistart.v1",
            "status": self.status,
        }

    @classmethod
    def from_dict(cls, value: object) -> HeuristicEvidence:
        """Require complete, recognized method evidence."""
        data = object_fields(value, {"schema", "method", "status"}, "heuristic")
        result = cls(data.get("status"))
        if data != result.to_dict():
            msg = "unsupported or incomplete heuristic evidence"
            raise ValueError(msg)
        return result


def validate_candidate_search(evidence: Mapping[str, object]) -> None:
    """Require either proven exact packing or explicitly unproven greedy evidence."""
    if "heuristic" in evidence:
        value = HeuristicEvidence.from_dict(evidence["heuristic"])
        if (
            value.status == "candidate"
            and evidence.get("proof_scope") is None
            and not {"solver_status", "backend_status"}.intersection(evidence)
        ):
            return
    elif (
        evidence.get("solver_status") == "optimal"
        and evidence.get("proof_scope") == "offered_packing_model"
    ):
        return
    msg = "candidate requires optimal packing evidence or explicit greedy evidence"
    raise ValueError(msg)


def search_exhausted(evidence: Mapping[str, object]) -> bool:
    """Whether this offered search ended, without broadening its proof scope."""
    return evidence.get("solver_status") == "infeasible" or (
        "heuristic" in evidence
        and HeuristicEvidence.from_dict(evidence["heuristic"]).status == "exhausted"
    )


def batch_search_finished(
    evidence: Mapping[str, object], *, on_unproven: str = "stop"
) -> bool:
    """Recognize exhaustion or a declared unproven-search batch boundary."""
    return search_exhausted(evidence) or (
        on_unproven == "next_batch"
        and evidence.get("solver_status") in {"unknown", "unproven"}
    )


def validate_heuristic_attempt(
    evidence: Mapping[str, object], *, outcome: str, candidate: bool
) -> None:
    """Reject a false proof, mixed methods or incomplete heuristic result."""
    if "heuristic" not in evidence:
        return
    value = HeuristicEvidence.from_dict(evidence["heuristic"])
    if (
        {"solver_status", "backend_status", "packing_objective"}.intersection(evidence)
        or evidence.get("proof_scope") is not None
        or evidence.get("termination_reason") != f"heuristic_{value.status}"
        or candidate != (value.status == "candidate")
        or (not candidate and outcome != "no_candidate")
    ):
        msg = "heuristic evidence disagrees with the recorded outcome or proof"
        raise ValueError(msg)
