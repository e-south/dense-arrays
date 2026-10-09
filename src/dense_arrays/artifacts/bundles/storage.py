"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/bundles/storage.py

Read bundle commits and checksummed plan/design records without external paths.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import json
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    canonical_json,
    digest,
    integer,
    object_fields,
    semantic_digest,
)
from dense_arrays.artifacts.bundles.models import (
    BUNDLE_DATABASE,
    BUNDLE_MANIFEST,
    BUNDLE_SCHEMA,
    EVIDENCE_BOUNDARY,
    RUNTIME_BUNDLE_SCHEMA,
    RUNTIME_EVIDENCE_BOUNDARY,
    BundleSummary,
)
from dense_arrays.artifacts.errors import integrity_boundary
from dense_arrays.artifacts.reading import ReadBudget, ReadLimitError, ReadLimits
from dense_arrays.artifacts.records import COMPOSITION_POLICY
from dense_arrays.artifacts.run_state import manifest_cell_ids
from dense_arrays.artifacts.store import checked_payload
from dense_arrays.planning import PlanEvidence
from dense_arrays.planning.batches.bindings import encoded_batch_size

if TYPE_CHECKING:
    import sqlite3
    from pathlib import Path


def is_bundle(path: Path) -> bool:
    """Recognize completed and interrupted bundle destinations."""
    return any(
        (path / name).exists()
        for name in (BUNDLE_MANIFEST, BUNDLE_DATABASE, ".bundle-pending")
    )


def read_summary(path: Path, limits: ReadLimits) -> BundleSummary:
    """Validate the commit and bounded metadata without scanning design evidence."""
    with integrity_boundary(path):
        if any((path / name).exists() for name in ("run.sqlite3", "pool.sqlite3")):
            msg = "ambiguous directory contains both bundle and native artifacts"
            raise ValueError(msg)
        manifest, database = path / BUNDLE_MANIFEST, path / BUNDLE_DATABASE
        if (
            (path / ".bundle-pending").exists()
            or not manifest.is_file()
            or not database.is_file()
        ):
            msg = "bundle is incomplete: missing committed manifest or database"
            raise ValueError(msg)
        if manifest.is_symlink() or database.is_symlink():
            msg = "bundle evidence must be contained files, not symbolic links"
            raise ValueError(msg)
        payload = manifest.read_text(encoding="utf-8")
        data = object_fields(
            json.loads(payload),
            {
                "schema",
                "bundle_id",
                "scope",
                "designs",
                "sources",
                "source_runs",
                "selection",
                "plans",
                "evidence",
                "metric_policy",
                "file",
                "batches",
            },
            "bundle",
        )
        identity = data.pop("bundle_id")
        if (
            data.get("schema") not in {BUNDLE_SCHEMA, RUNTIME_BUNDLE_SCHEMA}
            or data.get("scope") != "selected_collection"
        ):
            msg = "unsupported bundle schema or collection scope"
            raise ValueError(msg)
        if (
            identity != semantic_digest(data)
            or canonical_json({**data, "bundle_id": identity}) + "\n" != payload
        ):
            msg = "bundle manifest checksum or canonical encoding mismatch"
            raise ValueError(msg)
        _validate_metadata(data, limits)
        return BundleSummary({**data, "bundle_id": identity}, limits)


def _validate_metadata(data: dict[str, object], limits: ReadLimits) -> None:
    from dense_arrays.reporting.summary import RunSummary  # noqa: PLC0415

    _validate_runtime_schema(data)
    integer(data["designs"], field_name="designs", minimum=0)
    for name in ("sources", "source_runs", "plans"):
        if not isinstance(data[name], list):
            msg = "bundle source and plan inventories must be arrays"
            raise TypeError(msg)
    if (
        sum(len(data[name]) for name in ("sources", "source_runs", "plans"))
        > limits.identities
    ):
        msg = "read_limits.identities cannot hold bundle source metadata"
        raise ReadLimitError(msg)
    if not data["sources"] or not data["source_runs"] or not data["plans"]:
        msg = "bundle requires declared sources and bound plans"
        raise ValueError(msg)
    for source in data["source_runs"]:
        RunSummary.from_manifest(source)
    for plan_id in data["plans"]:
        digest(plan_id, field_name="plan_id")
    if len(set(data["plans"])) != len(data["plans"]):
        msg = "bundle plan inventory repeats an identity"
        raise ValueError(msg)
    _validate_selection(data)
    if (
        data["evidence"]
        != (
            RUNTIME_EVIDENCE_BOUNDARY
            if data["schema"] == RUNTIME_BUNDLE_SCHEMA
            else EVIDENCE_BOUNDARY
        )
        or data["metric_policy"] != COMPOSITION_POLICY
    ):
        msg = "unsupported bundle evidence or metric policy"
        raise ValueError(msg)
    _validate_file(data["file"])


def _validate_file(value: object) -> None:
    """Require one contained database with an explicit byte fingerprint."""
    file = object_fields(value, {"name", "bytes", "sha256"}, "bundle file")
    if file.get("name") != BUNDLE_DATABASE:
        msg = "bundle database must use the contained canonical filename"
        raise ValueError(msg)
    integer(file["bytes"], field_name="file.bytes", minimum=1)
    digest(file["sha256"], field_name="file.sha256")


def read_evidence(
    connection: sqlite3.Connection, plan_id: str, budget: ReadBudget
) -> PlanEvidence:
    """Check one bound plan before constructing its part/requirement state."""
    row = connection.execute(
        "SELECT payload,digest FROM plans WHERE plan_id=?", (plan_id,)
    ).fetchone()
    budget.examine(None if row is None else row[0])
    value = checked_payload(row)
    content = value.get("content", {})
    request = content.get("request", {})
    parts, requirements = request.get("parts"), request.get("requirements")
    if not isinstance(parts, list) or not isinstance(requirements, list):
        msg = "bundle plans require resolved part and requirement arrays"
        raise TypeError(msg)
    exclusions = content.get("parent", {}).get("exclusions", [])
    from dense_arrays.planning.libraries import exclusion_count  # noqa: PLC0415

    budget.retain(
        len(parts)
        + encoded_batch_size(request.get("batch"))
        + sum(
            encoded_batch_size(b)
            for b in request.get("schedule", {}).get("batches", [])
        )
        + len(requirements)
        + len(exclusions)
        + (0 if request.get("exclude") is None else exclusion_count(request["exclude"]))
    )
    plan = PlanEvidence.from_dict(value)
    if plan.plan_id != plan_id:
        msg = "bundle plan identity does not match its index"
        raise ValueError(msg)
    return plan


def _validate_selection(data: dict[str, object]) -> None:
    """Validate the declared predicate or saved allocation metadata."""
    from dense_arrays.reporting.design_filters import DesignFilter  # noqa: PLC0415
    from dense_arrays.reporting.selections.snapshots import (  # noqa: PLC0415
        validate_summary,
    )

    if data["selection"] is not None:
        if not isinstance(data["selection"], dict):
            msg = "bundle selection must be a declared object"
            raise TypeError(msg)
        if data["selection"].get("schema") == "dense_arrays.selection_summary.v1":
            selected = validate_summary(
                data["selection"],
                cells={
                    f"{s['run_id']}/{c}"
                    for s in data["source_runs"]
                    for c in manifest_cell_ids(s)
                },
            )
            if selected["selected"] != data["designs"]:
                msg = "bundle design count disagrees with its saved selection"
                raise ValueError(msg)
        else:
            DesignFilter.from_dict(data["selection"])


def _validate_runtime_schema(data: dict[str, object]) -> None:
    """Keep runtime evidence explicit in the declared bundle version."""
    if data["schema"] == RUNTIME_BUNDLE_SCHEMA:
        integer(data.get("batches"), field_name="batches", minimum=0)
    elif "batches" in data:
        msg = "runtime batch evidence requires bundle schema v2"
        raise ValueError(msg)
