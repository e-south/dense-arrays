"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/candidates.py

Pool-scoped candidate records preserve saved preparation decisions.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass

from dense_arrays._record_validation import digest, object_fields
from dense_arrays.parts.candidates import Candidate

CANDIDATE_RECORD_SCHEMA = "dense_arrays.pool_candidate.v1"


@dataclass(frozen=True, repr=False)
class PoolCandidate:
    """A candidate decision identified by its pool and original mining index."""

    pool_id: str
    candidate: Candidate

    def __post_init__(self) -> None:
        """Require a complete saved decision, including non-retained representatives."""
        digest(self.pool_id, field_name="pool_id")
        if not isinstance(self.candidate, Candidate):
            msg = "pool candidate requires Candidate"
            raise TypeError(msg)
        value = self.candidate
        if not value.error and not value.reasons and value.representative is None:
            msg = "pool candidate requires a recorded representative decision"
            raise ValueError(msg)

    @property
    def outcome(self) -> str:
        """Name one mutually exclusive terminal stage from the saved decision."""
        value = self.candidate
        if value.error is not None:
            return "execution_error"
        if value.reasons:
            return "eligibility_rejected"
        if value.representative != value.index:
            return "duplicate_discarded"
        return "retained" if value.retained else "not_selected"

    def to_dict(self) -> dict[str, object]:
        """Preserve the complete decision with its collection identity and outcome."""
        return {
            "schema": CANDIDATE_RECORD_SCHEMA,
            "pool_id": self.pool_id,
            "outcome": self.outcome,
            "candidate": self.candidate.to_dict(),
        }

    def __repr__(self) -> str:
        """Keep notebook output independent of sequence and observation size."""
        return (
            f"PoolCandidate({self.pool_id[:12]}, index={self.candidate.index}, "
            f"outcome={self.outcome})"
        )

    @classmethod
    def from_dict(cls, value: object) -> "PoolCandidate":
        """Reject contradictory outcome labels and incomplete native records."""
        fields = {"schema", "pool_id", "outcome", "candidate"}
        data = object_fields(value, fields, "pool candidate")
        if set(data) != fields or data["schema"] != CANDIDATE_RECORD_SCHEMA:
            msg = "unsupported or incomplete pool candidate schema"
            raise ValueError(msg)
        result = cls(data["pool_id"], Candidate.from_dict(data["candidate"]))
        if result.to_dict() != value:
            msg = "pool candidate outcome or contents disagree with saved decision"
            raise ValueError(msg)
        return result
