"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/bundles/verification.py

Independently recount contained designs without claiming search-history replay.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
from typing import TYPE_CHECKING

from dense_arrays._record_validation import mutable_json
from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE, BUNDLE_MANIFEST
from dense_arrays.artifacts.reading import ReadBudget, Verification
from dense_arrays.artifacts.run_state import origin_bindings
from dense_arrays.artifacts.store import reader
from dense_arrays.generation.acceptance import evaluate
from dense_arrays.reporting.bundles.batches import load_batches
from dense_arrays.reporting.bundles.reading import decode_design, load_plans
from dense_arrays.reporting.bundles.views import BundleView

if TYPE_CHECKING:
    from collections.abc import Collection
    from pathlib import Path

    from dense_arrays.artifacts.bundles.models import BundleSummary
    from dense_arrays.artifacts.records import Design
    from dense_arrays.planning import PlanEvidence


def verify_design(
    design: Design,
    plan: PlanEvidence,
    excluded: Collection[str],
    *,
    batches: dict | None = None,
) -> None:
    """Recount final geometry and requirements against its resolved plan."""
    if design.plan_id != plan.plan_id or design.realized.source_id != design.reference:
        msg = "bundle design provenance does not match its plan or source identity"
        raise ValueError(msg)
    if plan.request.resampling is not None:
        decision = (batches or {}).get((design.run_id, design.cell_id, design.batch_id))
        if decision is None:
            msg = "bundle design has no runtime membership evidence"
            raise ValueError(msg)
        plan = decision.bind(plan)
    evidence = evaluate(design.realized, plan, batch_id=design.batch_id)
    if tuple(mutable_json(r) for r in design.requirements) != evidence or not all(
        r["passed"] for r in evidence
    ):
        msg = "bundle requirement evidence does not match independent recount"
        raise ValueError(msg)
    if design.sequence_id in excluded:
        msg = "bundle design duplicates a frozen parent or ancestor sequence"
        raise ValueError(msg)


def exclusion_index(plans: dict, budget: ReadBudget) -> dict[str, set[str]]:
    """Build bounded membership indexes once for all included plan records."""
    result = {}
    for identity, plan in plans.items():
        exclusions = plan.exclusions
        budget.retain(len(exclusions))
        result[identity] = {e.sequence_id for e in exclusions}
    return result


def verify_bundle(path: Path, summary: BundleSummary) -> dict[str, object]:
    """Check physical bytes, plan identities and recounted requirements."""
    file = summary.manifest["file"]
    database = path / BUNDLE_DATABASE
    with database.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if file["bytes"] != database.stat().st_size or file["sha256"] != digest:
        msg = "bundle database byte count or checksum mismatch"
        raise ValueError(msg)
    budget = ReadBudget(summary.read_limits)
    budget.examine((path / BUNDLE_MANIFEST).read_text())
    query = BundleView(
        path, summary, "designs", limit=None, read_limits=summary.read_limits
    )
    sources = set(origin_bindings(summary.manifest["source_runs"]))
    budget.retain(len(sources))
    seen = set()
    count = 0
    with reader(path, filename=BUNDLE_DATABASE) as connection:
        plans = load_plans(connection, query, budget)
        batches = load_batches(connection, summary.manifest, plans, budget)
        excluded = exclusion_index(plans, budget)
        for row in connection.execute(
            "SELECT ordinal,design_ref,run_id,cell_id,local_id,plan_id,payload,digest "
            "FROM designs ORDER BY ordinal"
        ):
            ordinal, design = decode_design(row, budget)
            count += 1
            if (
                ordinal != count
                or (design.run_id, design.cell_id, design.plan_id) not in sources
                or design.realized.source_id != design.reference
            ):
                msg = "bundle design ordering or source provenance mismatch"
                raise ValueError(msg)
            plan = plans[design.plan_id]
            verify_design(design, plan, excluded[design.plan_id], batches=batches)
            key = (design.run_id, design.cell_id, design.sequence_id)
            if key in seen:
                msg = (
                    "bundle violates declared sequence uniqueness or parent exclusions"
                )
                raise ValueError(msg)
            budget.retain()
            seen.add(key)
    if count != summary.designs:
        msg = "bundle design count does not match its manifest"
        raise ValueError(msg)
    return {
        **Verification(
            "selected_designs_and_resolved_plans", budget.examined, budget.bytes_checked
        ).to_dict(),
        "file_bytes_checked": file["bytes"],
        "search_history": "not_included",
    }
