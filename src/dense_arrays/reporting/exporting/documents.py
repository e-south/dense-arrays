"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/exporting/documents.py

Native JSON documents preserve input schemas and declared report scope.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
from contextlib import nullcontext
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.bundles.models import BundleSummary
from dense_arrays.artifacts.pool_records import PoolSummary
from dense_arrays.artifacts.publication import new_text_file
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits
from dense_arrays.artifacts.receipts import ExportReceipt
from dense_arrays.parts import PreparationSet, PreparationSpec
from dense_arrays.planning import (
    DesignSpec,
    ExtensionSpec,
    GenerationPlan,
    MatrixPlan,
    MatrixSpec,
    PlanEvidence,
    PreparationPlan,
)
from dense_arrays.planning.preparation import preparation_to_dict
from dense_arrays.planning.serialization import request_to_dict
from dense_arrays.reporting.diagnostics import DiagnosticReport
from dense_arrays.reporting.plans import (
    PlanComparison,
    check_plan_limits,
    editable_request,
)
from dense_arrays.reporting.plans.requests import RequestReport
from dense_arrays.reporting.pools import PoolQualityReport, PoolQualitySnapshot
from dense_arrays.reporting.quality import (
    QualityComparison,
    QualityReport,
    QualitySnapshot,
)
from dense_arrays.reporting.selections import SelectionSnapshot
from dense_arrays.reporting.summary import RunSummary

if TYPE_CHECKING:
    from typing import TextIO

DOCUMENT_TYPES = {
    "selection": (SelectionSnapshot,),
    "plan": (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan),
    "request": (
        RequestReport,
        GenerationPlan,
        PlanEvidence,
        PreparationPlan,
        DesignSpec,
        PreparationSpec,
        PreparationSet,
        ExtensionSpec,
        MatrixPlan,
        MatrixSpec,
    ),
    "quality": (QualityReport, QualitySnapshot, PoolQualityReport, PoolQualitySnapshot),
    "diagnostics": (DiagnosticReport,),
    "summary": (RunSummary, PoolSummary, BundleSummary),
    "comparison": (PlanComparison, QualityComparison),
}


def document_view(value: object) -> str | None:
    """Infer a document projection from an explicitly typed Python value."""
    return next(
        (view for view, types in DOCUMENT_TYPES.items() if isinstance(value, types)),
        None,
    )


@dataclass(frozen=True)
class DocumentView:
    """A typed document and its export projection, without performing report scans."""

    value: (
        GenerationPlan
        | PlanEvidence
        | PreparationPlan
        | PreparationSpec
        | PreparationSet
        | DesignSpec
        | ExtensionSpec
        | QualityReport
        | PoolQualityReport
        | PoolQualitySnapshot
        | QualitySnapshot
        | QualityComparison
        | DiagnosticReport
        | RunSummary
        | PoolSummary
        | BundleSummary
        | PlanComparison
        | SelectionSnapshot
        | RequestReport
    )
    view: str
    read_limits: ReadLimits = field(default_factory=ReadLimits)

    def __post_init__(self) -> None:
        """Reject mismatched projections before opening any output."""
        if self.view not in DOCUMENT_TYPES or not isinstance(
            self.value, DOCUMENT_TYPES[self.view]
        ):
            msg = "document does not support the requested export view"
            raise TypeError(msg)
        if not isinstance(self.read_limits, ReadLimits):
            msg = "document exports require ReadLimits"
            raise TypeError(msg)
        if isinstance(self.value, SelectionSnapshot):
            ReadBudget(self.read_limits).retain(
                self.value.selected + len(self.value.sources) + len(self.value.counts)
            )
        if isinstance(
            self.value,
            (RequestReport, DesignSpec, PreparationSpec, PreparationSet, MatrixSpec),
        ):
            request = (
                self.value
                if isinstance(self.value, RequestReport)
                else RequestReport(self.value)
            )
            ReadBudget(self.read_limits).retain(request.identities)
        if isinstance(
            self.value, (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan)
        ):
            check_plan_limits(self.read_limits, self.value)
        elif isinstance(self.value, PlanComparison):
            check_plan_limits(self.read_limits, self.value.before, self.value.after)
        if self.view == "request" and isinstance(
            self.value, (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan)
        ):
            editable_request(self.value)

    @property
    def cost(self) -> ReadCost:
        """Reuse the report's exact cost; native input serialization is one document."""
        if isinstance(
            self.value,
            (
                QualityReport,
                PoolQualityReport,
                PoolQualitySnapshot,
                QualitySnapshot,
                QualityComparison,
                DiagnosticReport,
                RunSummary,
                PoolSummary,
                BundleSummary,
                SelectionSnapshot,
            ),
        ):
            return self.value.cost
        identity = (
            self.value.plan_id
            if isinstance(
                self.value, (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan)
            )
            else semantic_digest(
                request_to_dict(self.value)
                if isinstance(self.value, DesignSpec)
                else preparation_to_dict(self.value)
                if isinstance(self.value, (PreparationSpec, PreparationSet))
                else self.value.to_dict()
            )
        )
        return ReadCost(identity, 0, "manifest", self.view, 1, self.read_limits)

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Bind report sources or the native input identity outside output bytes."""
        value = self.value
        if isinstance(
            value,
            (
                QualityComparison,
                QualitySnapshot,
                SelectionSnapshot,
                QualityReport,
                PoolQualityReport,
                PoolQualitySnapshot,
            ),
        ):
            return (
                tuple(s.content() for s in value.sources)
                if isinstance(value, SelectionSnapshot)
                else value.sources
            )
        if isinstance(value, BundleSummary):
            return ({"bundle_id": value.bundle_id, "revision": 0},)
        if isinstance(value, (DiagnosticReport, RunSummary)):
            summary = value.summary if isinstance(value, DiagnosticReport) else value
            return (
                {
                    "run_id": summary.run_id,
                    "revision": summary.revision,
                    "plan_id": summary.plan_id,
                },
            )
        if isinstance(value, PoolSummary):
            return (
                {"pool_id": value.pool_id, "plan_id": value.plan_id, "revision": 0},
            )
        if isinstance(value, PlanComparison):
            return tuple(
                {"plan_id": plan.plan_id} for plan in (value.before, value.after)
            )
        return (
            {
                "plan_id"
                if isinstance(
                    value, (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan)
                )
                else "request_id": self.cost.source_id
            },
        )

    def to_dict(self, *, base: Path) -> dict[str, object]:
        """Keep reports unwrapped and encode input paths relative to their new file."""
        if self.view == "request":
            request = (
                editable_request(self.value)
                if isinstance(
                    self.value,
                    (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan),
                )
                else self.value
            )
            if isinstance(request, RequestReport):
                return request.to_dict(base=base)
            if isinstance(request, (PreparationSpec, PreparationSet)):
                return preparation_to_dict(request, base=base)
            if isinstance(request, (ExtensionSpec, MatrixSpec)):
                return request.to_dict(base=base)
            return request_to_dict(request, base=base)
        if isinstance(
            self.value, (GenerationPlan, PreparationPlan, MatrixPlan, SelectionSnapshot)
        ):
            return self.value.to_dict(base=base)
        return self.value.to_dict()


def export_document(query: DocumentView, *, out: str | Path | TextIO) -> ExportReceipt:
    """Compute a complete declared report before publishing one create-only file."""
    path = Path(out).absolute() if isinstance(out, (str, Path)) else None
    if path is not None and (path.exists() or path.is_symlink()):
        msg = f"output destination already exists: {path}"
        raise FileExistsError(msg)
    if path is None and not callable(getattr(out, "write", None)):
        msg = "out must be a file path or writable text stream"
        raise TypeError(msg)
    value = query.to_dict(base=path.parent if path is not None else Path.cwd())
    payload = canonical_json(value) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    sources = tuple(
        {
            **s,
            "document_schema": value["schema"],
            "document_sha256": digest,
            "scope": "declared_document",
        }
        for s in query.sources
    )
    with new_text_file(path) if path is not None else nullcontext(out) as stream:
        stream.write(payload)
    files = (
        ()
        if path is None
        else ({"name": path.name, "bytes": len(payload.encode()), "sha256": digest},)
    )
    return ExportReceipt(
        str(path) if path is not None else "<stream>",
        "json",
        query.view,
        1,
        sources,
        (),
        files,
    )
