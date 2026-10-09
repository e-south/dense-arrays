"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/candidates.py

Saved packing and final candidate evidence, separate from accepted designs.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import mutable_json, object_fields
from dense_arrays.artifacts.search import validate_candidate_search
from dense_arrays.playback.serialization import (
    realized_array_from_dict,
    realized_array_to_dict,
)
from dense_arrays.realized import RealizedArray

if TYPE_CHECKING:
    from collections.abc import Mapping


@dataclass(frozen=True, repr=False)
class CandidateEvidence:
    """One packing and the last evaluated final sequence, if evaluation began."""

    packed: RealizedArray
    final: RealizedArray | None

    def __repr__(self) -> str:
        """Describe saved evidence without dumping sequence or placement records."""
        final_length = None if self.final is None else len(self.final.sequence)
        return (
            f"CandidateEvidence(packed_length={len(self.packed.sequence)}, "
            f"final_length={final_length}, placements={len(self.packed.placements)})"
        )

    def __post_init__(self) -> None:
        """Preserve realized-array authority and require one candidate identity."""
        if not isinstance(self.packed, RealizedArray) or (
            self.final is not None and not isinstance(self.final, RealizedArray)
        ):
            msg = "candidate evidence requires realized arrays"
            raise TypeError(msg)
        if self.final is not None and self.final.source_id != self.packed.source_id:
            msg = "packed and final candidate identities differ"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Encode supplied placements; no generation or acceptance is implied."""
        return {
            "schema": "dense_arrays.candidate.v1",
            "packed": realized_array_to_dict(self.packed),
            "final": None if self.final is None else realized_array_to_dict(self.final),
        }

    def validate_outcome(self, outcome: str, evidence: Mapping[str, object]) -> None:
        """Check structural joins without making a plan-dependent acceptance claim."""
        if outcome not in {"accepted", "rejected", "duplicate", "no_candidate"}:
            msg = "candidate evidence requires a completed candidate evaluation"
            raise ValueError(msg)
        validate_candidate_search(evidence)
        if self.final is None and (
            outcome != "no_candidate" or evidence.get("assembly_trials") != 0
        ):
            msg = "only an unevaluated candidate may omit its final sequence"
            raise ValueError(msg)

    @classmethod
    def from_dict(cls, value: object) -> CandidateEvidence:
        """Reject unknown versions and malformed saved candidate geometry."""
        keys = {"schema", "packed", "final"}
        data = object_fields(value, keys, "candidate")
        if set(data) != keys or data["schema"] != "dense_arrays.candidate.v1":
            msg = "unsupported or incomplete candidate schema"
            raise ValueError(msg)
        return cls(
            realized_array_from_dict(mutable_json(data["packed"])),
            None
            if data["final"] is None
            else realized_array_from_dict(mutable_json(data["final"])),
        )
