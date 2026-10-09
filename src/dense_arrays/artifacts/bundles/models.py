"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/bundles/models.py

Bundle summaries retain collection scope and original source accounting.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import immutable_json_mapping, mutable_json
from dense_arrays.artifacts.reading import ReadCost, ReadLimits

if TYPE_CHECKING:
    from collections.abc import Mapping

RUNTIME_BUNDLE_SCHEMA = "dense_arrays.bundle.v2"
BUNDLE_SCHEMA = "dense_arrays.bundle.v1"
BUNDLE_DATABASE = "bundle.sqlite3"
BUNDLE_MANIFEST = "bundle.json"
EVIDENCE_BOUNDARY = {
    "verification": "selected_designs_and_resolved_plans",
    "original_inputs": "fingerprints_only",
    "attempt_history": "not_included",
    "ancestor_records": "identity_references_only",
}

RUNTIME_EVIDENCE_BOUNDARY = {
    **EVIDENCE_BOUNDARY,
    "runtime_batches": "membership_and_recorded_feedback",
    "feedback_history": "not_included",
}


@dataclass(frozen=True)
class BundleSummary:
    """Committed selected-collection metadata; source attainment retains its scope."""

    manifest: Mapping[str, object]
    read_limits: ReadLimits = field(default_factory=ReadLimits, repr=False)
    verification: Mapping[str, object] | None = None

    def __post_init__(self) -> None:
        """Detach metadata from mutable decoded input."""
        object.__setattr__(self, "manifest", immutable_json_mapping(self.manifest))
        if self.verification is not None:
            object.__setattr__(
                self, "verification", immutable_json_mapping(self.verification)
            )

    @property
    def bundle_id(self) -> str:
        """Identify the committed evidence independently of its directory."""
        return self.manifest["bundle_id"]

    @property
    def designs(self) -> int:
        """Number of included full design identities."""
        return self.manifest["designs"]

    @property
    def verified(self) -> bool:
        """Whether this read independently checked all included evidence."""
        return self.verification is not None

    @property
    def cost(self) -> ReadCost:
        """Describe a manifest read without scanning selected designs."""
        return ReadCost(self.bundle_id, 0, "manifest", "summary", 1, self.read_limits)

    @property
    def verification_cost(self) -> ReadCost:
        """Declare the complete selected evidence scan, plus physical file hashing."""
        return ReadCost(
            self.bundle_id,
            0,
            "scan",
            "verification",
            1
            + len(self.manifest["plans"])
            + self.designs
            + self.manifest.get("batches", 0),
            self.read_limits,
            bytes_estimate=self.manifest["file"]["bytes"],
        )

    def to_dict(self) -> dict[str, object]:
        """Present source summaries as original context, never subset completion."""
        from dense_arrays.reporting.summary import RunSummary  # noqa: PLC0415

        return {
            "schema": "dense_arrays.bundle_summary.v2"
            if "batches" in self.manifest
            else "dense_arrays.bundle_summary.v1",
            **(
                {"batches": self.manifest["batches"]}
                if "batches" in self.manifest
                else {}
            ),
            "bundle_id": self.bundle_id,
            "scope": "selected_collection",
            "designs": self.designs,
            "selection": mutable_json(self.manifest["selection"]),
            "sources": mutable_json(self.manifest["sources"]),
            "source_runs": [
                RunSummary.from_manifest(mutable_json(s)).to_dict()
                for s in self.manifest["source_runs"]
            ],
            "plans": list(self.manifest["plans"]),
            "evidence": mutable_json(self.manifest["evidence"]),
            "verified": self.verified,
            **(
                {"verification": mutable_json(self.verification)}
                if self.verified
                else {}
            ),
        }
