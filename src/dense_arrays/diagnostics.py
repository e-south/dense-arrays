"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/diagnostics.py

Immutable diagnostic records shared by planning and result inspection.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    mutable_json,
    object_fields,
    records,
    required_text,
)

if TYPE_CHECKING:
    from collections.abc import Mapping


@dataclass(frozen=True)
class Diagnostic:
    """A scoped observation with stable code, evidence and a suggested next action."""

    code: str
    severity: str
    stage: str
    observed: Mapping[str, object]
    expected: Mapping[str, object]
    evidence_refs: tuple[str, ...]
    next_action: str
    requirement_id: str | None = None
    proof_scope: str | None = None

    def __post_init__(self) -> None:
        """Freeze evidence without interpreting an unknown value as zero."""
        for name in ("code", "stage", "next_action"):
            required_text(getattr(self, name), field_name=name)
        if self.severity not in {"info", "warning", "error"}:
            msg = "unknown diagnostic severity"
            raise ValueError(msg)
        for name in ("observed", "expected"):
            object.__setattr__(self, name, immutable_json_mapping(getattr(self, name)))
        references = records(self.evidence_refs, str, field_name="evidence_refs")
        for reference in references:
            required_text(reference, field_name="evidence_refs")
        for name in ("requirement_id", "proof_scope"):
            if getattr(self, name) is not None:
                required_text(getattr(self, name), field_name=name)
        object.__setattr__(self, "evidence_refs", references)

    def to_dict(self) -> dict[str, object]:
        """Return evidence with its requirement and proof boundary explicit."""
        return {
            "schema": "dense_arrays.diagnostic.v1",
            "code": self.code,
            "severity": self.severity,
            "stage": self.stage,
            "requirement_id": self.requirement_id,
            "observed": mutable_json(self.observed),
            "expected": mutable_json(self.expected),
            "evidence_refs": list(self.evidence_refs),
            "proof_scope": self.proof_scope,
            "next_action": self.next_action,
        }

    @classmethod
    def from_dict(cls, value: object) -> Diagnostic:
        """Read one supported record without accepting future fields or versions."""
        keys = {
            "schema",
            "code",
            "severity",
            "stage",
            "requirement_id",
            "observed",
            "expected",
            "evidence_refs",
            "proof_scope",
            "next_action",
        }
        data = object_fields(value, keys, "diagnostic")
        if set(data) != keys or data.pop("schema") != "dense_arrays.diagnostic.v1":
            msg = "unsupported or incomplete diagnostic schema"
            raise ValueError(msg)
        return cls(**data)
