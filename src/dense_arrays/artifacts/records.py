"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/records.py

Logical native records; runtime paths are separate from content identity.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from dense_arrays._record_validation import (
    digest,
    immutable_json_mapping,
    integer,
    mutable_json,
    object_fields,
    records,
    required_text,
    semantic_digest,
)
from dense_arrays.artifacts.candidates import CandidateEvidence
from dense_arrays.artifacts.objectives import PackingObjective
from dense_arrays.artifacts.search import (
    validate_candidate_search,
    validate_heuristic_attempt,
)
from dense_arrays.playback.serialization import (
    realized_array_from_dict,
    realized_array_to_dict,
)
from dense_arrays.realized import RealizedArray

COMPOSITION_POLICY = "library_composition.v2"
DESIGN_RECORD_SCHEMA = "dense_arrays.design_record.v1"
RUN_SCHEMA = "dense_arrays.run.v2"
OUTCOMES = (
    "accepted",
    "rejected",
    "duplicate",
    "no_candidate",
    "error",
    "interrupted_unresolved",
    "in_progress",
)


def requirement_evidence(value: object) -> tuple[Mapping[str, object], ...]:
    """Freeze declared check results without coercing truth values or missing data."""
    required = {"id", "observed", "passed"}
    results = []
    seen = set()
    for record in records(value, Mapping, field_name="requirements"):
        data = object_fields(record, required | {"status"}, "requirement evidence")
        if not required <= data.keys():
            msg = "incomplete requirement evidence"
            raise ValueError(msg)
        required_text(data["id"], field_name="requirement.id")
        if not isinstance(data["passed"], bool):
            msg = "requirement.passed must be a boolean"
            raise TypeError(msg)
        if data["id"] in seen:
            msg = "duplicate requirement evidence"
            raise ValueError(msg)
        if "status" in data and data["status"] != "not_applicable":
            msg = "unsupported requirement evidence status"
            raise ValueError(msg)
        seen.add(data["id"])
        results.append(immutable_json_mapping(data))
    return tuple(results)


def _validate_attempt_effort(evidence: Mapping[str, object]) -> None:
    """Validate optional effort fields without inventing absent historical values."""
    for name, minimum in (("assembly_trials", 0), ("cell_attempt", 1)):
        if name in evidence:
            integer(evidence[name], field_name=name, minimum=minimum)


@dataclass(frozen=True, repr=False)
class Attempt:
    """One reserved attempt's status and typed solver or execution evidence."""

    attempt_id: int
    cell_id: str
    outcome: str
    evidence: Mapping[str, object]
    design_ref: str | None = None
    candidate: CandidateEvidence | None = field(init=False, repr=False, compare=False)

    def __repr__(self) -> str:
        """Show the outcome without expanding saved candidate or screening records."""
        return (
            f"Attempt(attempt_id={self.attempt_id}, cell_id={self.cell_id!r}, "
            f"outcome={self.outcome!r}, candidate={self.candidate!r})"
        )

    def __post_init__(self) -> None:
        """Validate a complete outcome without guessing backend termination causes."""
        integer(self.attempt_id, field_name="attempt_id", minimum=1)
        required_text(self.cell_id, field_name="cell_id")
        if self.outcome not in OUTCOMES:
            msg = f"unknown attempt outcome {self.outcome!r}"
            raise ValueError(msg)
        allowed = {
            "solver_status",
            "backend_status",
            "proof_scope",
            "termination_reason",
            "detail",
            "code",
            "assembly_trials",
            "requirements",
            "matched_design_ref",
            "candidate_sequence_id",
            "candidate",
            "cell_attempt",
            "batch_id",
            "batch_index",
            "batch_attempt",
            "packing_objective",
            "heuristic",
        }
        evidence = object_fields(self.evidence, allowed, "attempt.evidence")
        if "packing_objective" in evidence:
            PackingObjective.from_dict(evidence["packing_objective"])
        candidate = (
            CandidateEvidence.from_dict(evidence["candidate"])
            if "candidate" in evidence
            else None
        )
        object.__setattr__(self, "candidate", candidate)
        validate_heuristic_attempt(
            evidence, outcome=self.outcome, candidate=candidate is not None
        )
        if candidate is not None:
            candidate.validate_outcome(self.outcome, evidence)
        if evidence.get("code") == "parent_duplicate" or {
            "matched_design_ref",
            "candidate_sequence_id",
        }.intersection(evidence):
            if (
                self.outcome != "duplicate"
                or evidence.get("code") != "parent_duplicate"
            ):
                msg = "parent duplicate evidence requires a duplicate outcome"
                raise ValueError(msg)
            required_text(
                evidence.get("matched_design_ref"), field_name="matched_design_ref"
            )
            digest(
                evidence.get("candidate_sequence_id"),
                field_name="candidate_sequence_id",
            )
        if "requirements" in evidence:
            evidence["requirements"] = requirement_evidence(evidence["requirements"])
        _validate_batch_context(evidence)
        _validate_attempt_effort(evidence)
        object.__setattr__(self, "evidence", immutable_json_mapping(evidence))
        if self.outcome == "accepted":
            required_text(self.design_ref, field_name="design_ref")
            validate_candidate_search(self.evidence)
        elif self.design_ref is not None:
            msg = "only accepted attempts may reference a published design"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Serialize one attempt event, omitting absent design references."""
        value = {
            "schema": "dense_arrays.attempt.v1",
            "attempt_id": self.attempt_id,
            "cell_id": self.cell_id,
            "outcome": self.outcome,
            "evidence": mutable_json(self.evidence),
        }
        if self.design_ref is not None:
            value["design_ref"] = self.design_ref
        return value

    @classmethod
    def from_dict(cls, value: object) -> Attempt:
        """Reject unsupported attempt versions and unknown fields."""
        required = {"schema", "attempt_id", "cell_id", "outcome", "evidence"}
        data = object_fields(value, required | {"design_ref"}, "attempt")
        if not required <= set(data) or data.pop("schema") != "dense_arrays.attempt.v1":
            msg = "unsupported or incomplete attempt schema"
            raise ValueError(msg)
        return cls(**data)


def _validate_batch_context(evidence: dict[str, object]) -> None:
    if "batch_index" in evidence or "batch_attempt" in evidence:
        integer(evidence.get("batch_index"), field_name="batch_index", minimum=1)
        integer(evidence.get("batch_attempt"), field_name="batch_attempt", minimum=1)
        digest(evidence.get("batch_id"), field_name="batch_id")
    if "batch_id" in evidence:
        digest(evidence["batch_id"], field_name="attempt.batch_id")


@dataclass(frozen=True)
class RunHandle:
    """A lightweight run reference; constructing it opens no files."""

    path: Path
    run_id: str
    revision: int | None = None

    def __post_init__(self) -> None:
        """Freeze location and validate operational identity."""
        object.__setattr__(self, "path", Path(self.path).absolute())
        required_text(self.run_id, field_name="run_id")
        if self.revision is not None:
            integer(self.revision, field_name="revision", minimum=0)


@dataclass(frozen=True)
class Design:
    """One realized array and its plan, attempt and requirement evidence."""

    run_id: str
    cell_id: str
    design_id: str
    plan_id: str
    attempt_id: int
    realized: RealizedArray
    requirements: tuple[Mapping[str, object], ...]
    batch_id: str | None = None

    def __post_init__(self) -> None:
        """Freeze evidence while preserving the existing realized-array authority."""
        for name in ("run_id", "cell_id", "design_id"):
            required_text(getattr(self, name), field_name=name)
        digest(self.plan_id, field_name="plan_id")
        if self.batch_id is not None:
            digest(self.batch_id, field_name="batch_id")
        integer(self.attempt_id, field_name="attempt_id", minimum=1)
        if not isinstance(self.realized, RealizedArray):
            msg = "realized must be RealizedArray"
            raise TypeError(msg)
        object.__setattr__(
            self,
            "requirements",
            requirement_evidence(self.requirements),
        )

    @property
    def sequence_id(self) -> str:
        """Versioned equivalence of exact final DNA, independent of its history."""
        return semantic_digest(
            {"schema": "dense_arrays.sequence.v1", "sequence": self.realized.sequence}
        )

    @property
    def reference(self) -> str:
        """Full operational identity; local IDs cannot collapse different runs."""
        return f"{self.run_id}/{self.cell_id}/{self.design_id}"

    def to_dict(self) -> dict[str, object]:
        """Serialize the logical design without duplicating authoritative DNA."""
        return {
            "schema": DESIGN_RECORD_SCHEMA,
            "run_id": self.run_id,
            "cell_id": self.cell_id,
            "design_id": self.design_id,
            "plan_id": self.plan_id,
            "attempt_id": self.attempt_id,
            "sequence_id": self.sequence_id,
            "realized": realized_array_to_dict(self.realized),
            "requirements": [mutable_json(r) for r in self.requirements],
            **({"batch_id": self.batch_id} if self.batch_id is not None else {}),
        }

    @classmethod
    def from_dict(cls, value: object) -> Design:
        """Reject unsupported schemas, noncanonical geometry and changed digests."""
        keys = {
            "schema",
            "run_id",
            "cell_id",
            "design_id",
            "plan_id",
            "attempt_id",
            "sequence_id",
            "realized",
            "requirements",
        }
        data = object_fields(value, keys | {"batch_id"}, "design_record")
        if not keys <= set(data) or data.pop("schema") != DESIGN_RECORD_SCHEMA:
            msg = "unsupported or incomplete design record schema"
            raise ValueError(msg)
        sequence_id = data.pop("sequence_id")
        realized = data.pop("realized")
        result = cls(realized=realized_array_from_dict(realized), **data)
        if result.sequence_id != sequence_id or result.to_dict() != value:
            msg = "design digest or canonical representation mismatch"
            raise ValueError(msg)
        return result
