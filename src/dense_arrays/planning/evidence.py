"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/evidence.py

Resolved generation meaning, independent of executable input locations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING

from dense_arrays._record_validation import digest, object_fields, semantic_digest
from dense_arrays.parts import BoundParts, PartTable, PoolSource
from dense_arrays.parts.bound import collection_identity
from dense_arrays.parts.provenance import ImportReport
from dense_arrays.planning.batches.validation import (
    validate_eligible,
    validate_selected,
)
from dense_arrays.planning.libraries import ParentLibrary
from dense_arrays.planning.models import DesignSpec
from dense_arrays.planning.requirements import Fixed
from dense_arrays.planning.serialization import (
    PLAN_SCHEMA,
    policies_for,
    request_from_dict,
    request_to_dict,
)
from dense_arrays.planning.validation import validate_requirements

if TYPE_CHECKING:
    from dense_arrays.planning.libraries import ExcludedDesign

EVIDENCE_SCHEMA = "dense_arrays.plan_evidence.v1"


@dataclass(frozen=True, repr=False)
class PlanEvidence:
    """Frozen parts, rules and provenance sufficient to verify accepted designs.

    Input fingerprints describe original sources. This record has no source
    locators and cannot be submitted to run as an executable generation plan.
    """

    request: DesignSpec
    input_digests: tuple[str, ...] = ()
    import_report: ImportReport | None = None
    parent: ParentLibrary | None = None
    plan_id: str = field(init=False)

    def __post_init__(self) -> None:
        """Validate one normalized packing contract and retain its original identity."""
        if not isinstance(self.request, DesignSpec) or isinstance(
            self.request.parts, (PartTable, PoolSource, BoundParts)
        ):
            msg = "generation plans require a DesignSpec with resolved parts"
            raise TypeError(msg)
        if not isinstance(self.input_digests, (list, tuple)):
            msg = "input_digests must be an ordered array"
            raise TypeError(msg)
        for value in self.input_digests:
            digest(value, field_name="input.sha256")
        object.__setattr__(self, "input_digests", tuple(self.input_digests))
        if self.request.length.exact is not None and self.request.assembly is None:
            msg = "exact length requires an explicit assembly policy"
            raise ValueError(msg)
        if (
            self.request.length.maximum is not None
            and self.request.assembly is not None
        ):
            msg = "assembly requires length.exact; omit assembly for length.maximum"
            raise ValueError(msg)
        validate_requirements(self.request)
        self._validate_exclusions()
        if self.import_report is None:
            object.__setattr__(
                self, "import_report", ImportReport("inline", len(self.request.parts))
            )
        if not isinstance(
            self.import_report, ImportReport
        ) or self.import_report.rows != len(self.request.parts):
            msg = "import report must describe the resolved parts"
            raise ValueError(msg)
        self._validate_batch()
        object.__setattr__(self, "plan_id", semantic_digest(self.content()))

    def _validate_batch(self) -> None:
        """Bind prepared membership to this exact eligible collection."""
        if self.request.resampling is not None:
            sampling = self.request.resampling.sampling
            fixed = {
                r.part_id for r in self.request.requirements if isinstance(r, Fixed)
            }
            if not len(fixed) <= sampling.size <= len(self.request.parts):
                msg = (
                    "batch size must cover fixed occurrences "
                    "and not exceed eligible parts"
                )
                raise ValueError(msg)
            validate_eligible(self.request.parts, sampling)
        batches = (
            self.request.schedule.batches
            if self.request.schedule is not None
            else (() if self.request.batch is None else (self.request.batch,))
        )
        eligible = {p.part_id: p for p in self.request.parts} if batches else {}
        validated = set()
        for batch in batches:
            if batch.collection_id != self.collection_id:
                msg = "candidate batch belongs to a different eligible collection"
                raise ValueError(msg)
            if set(batch.part_ids) - eligible.keys():
                msg = "candidate batch references unknown part identities"
                raise ValueError(msg)
            if batch.sampling is not None:
                if batch.sampling not in validated:
                    validate_eligible(self.request.parts, batch.sampling)
                    validated.add(batch.sampling)
                validate_selected([eligible[p] for p in batch.part_ids], batch.sampling)
            if batch.feedback is not None:
                batch.feedback.weights(self.request.parts)

    def _validate_exclusions(self) -> None:
        """Require frozen source evidence and one exclusion owner per plan."""
        if self.parent is not None and not isinstance(self.parent, ParentLibrary):
            msg = "plan parent must be a frozen ParentLibrary"
            raise TypeError(msg)
        if self.request.exclude is not None:
            if not isinstance(self.request.exclude.source, ParentLibrary):
                msg = "generation plans require a frozen exclusion library"
                raise TypeError(msg)
            if self.parent is not None:
                msg = "extension plans cannot also declare request exclusions"
                raise ValueError(msg)

    def content(self) -> dict[str, object]:
        """Canonical semantic payload shared with executable generation plans."""
        return {
            "schema": PLAN_SCHEMA,
            "request": request_to_dict(self.request),
            "input_digests": list(self.input_digests),
            "policies": policies_for(self.request),
            "import_report": self.import_report.to_dict(),
            **({"parent": self.parent.to_dict()} if self.parent is not None else {}),
        }

    @property
    def exclusions(self) -> tuple[ExcludedDesign, ...]:
        """Resolved sequence exclusions shared by execution and verification."""
        if self.parent is not None:
            return self.parent.exclusions
        return (
            ()
            if self.request.exclude is None
            else self.request.exclude.source.exclusions
        )

    @cached_property
    def collection_id(self) -> str:
        """Keep part collection identity independent of target, lineage and paths."""
        return collection_identity(
            self.request.parts, self.import_report, self.input_digests
        )

    def to_dict(self) -> dict[str, object]:
        """Encode a verification record without inventing missing input locators."""
        return {
            "schema": EVIDENCE_SCHEMA,
            "plan_id": self.plan_id,
            "content": self.content(),
        }

    def __repr__(self) -> str:
        """Show bounded identity and scope without traversing the resolved records."""
        return (
            f"PlanEvidence(plan={self.plan_id[:12]}, parts={len(self.request.parts)}, "
            f"requirements={len(self.request.requirements)})"
        )

    @classmethod
    def from_dict(cls, value: object) -> PlanEvidence:
        """Reject changed policies or identities and omitted normalized fields."""
        data = object_fields(value, {"schema", "plan_id", "content"}, "plan evidence")
        if data.get("schema") != EVIDENCE_SCHEMA:
            msg = "unsupported plan evidence schema"
            raise ValueError(msg)
        content = object_fields(
            data.get("content"),
            {
                "schema",
                "request",
                "input_digests",
                "policies",
                "import_report",
                "parent",
            },
            "plan evidence content",
        )
        request = request_from_dict(content.get("request"))
        if content.get("schema") != PLAN_SCHEMA or content.get(
            "policies"
        ) != policies_for(request):
            msg = "unsupported generation evidence policies"
            raise ValueError(msg)
        result = cls(
            request,
            content["input_digests"],
            ImportReport.from_dict(content["import_report"]),
            ParentLibrary.from_dict(content["parent"]) if "parent" in content else None,
        )
        if result.plan_id != data.get("plan_id"):
            msg = "plan evidence digest mismatch"
            raise ValueError(msg)
        if result.content() != content:
            msg = "plan evidence requires complete normalized fields"
            raise ValueError(msg)
        return result
