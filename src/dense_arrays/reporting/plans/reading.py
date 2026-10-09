"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/plans/reading.py

Read resolved plans and compare semantic generation choices without solving.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.bundles.storage import is_bundle, read_summary
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.pools import POOL_DATABASE, pool_summary, stored_preparation
from dense_arrays.artifacts.reading import ReadLimitError
from dense_arrays.artifacts.store import checked_payload, latest, reader, stored_plan
from dense_arrays.parts import PoolHandle
from dense_arrays.planning import (
    GenerationPlan,
    LibraryExclusion,
    Lineage,
    MatrixPlan,
    PlanEvidence,
    PreparationPlan,
    RunReference,
)
from dense_arrays.planning.batches.bindings import membership_size
from dense_arrays.reporting.plans.requests import RequestReport

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.reading import ReadLimits


def read_plan(
    artifact: GenerationPlan
    | PlanEvidence
    | PreparationPlan
    | RunHandle
    | PoolHandle
    | Path,
    limits: ReadLimits,
    *,
    as_request: bool = False,
) -> GenerationPlan | PlanEvidence | PreparationPlan | MatrixPlan | RequestReport:
    """Read a native immutable plan without reopening any original source inputs."""
    if isinstance(
        artifact, (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan)
    ):
        check_plan_limits(limits, artifact)
        return editable_request(artifact) if as_request else artifact
    path = artifact.path if isinstance(artifact, (RunHandle, PoolHandle)) else artifact
    if is_bundle(path):
        from .bundles import read_bundle_plan  # noqa: PLC0415

        if isinstance(artifact, (RunHandle, PoolHandle)):
            msg = "run and pool handles cannot identify a portable bundle"
            raise TypeError(msg)
        plan = read_bundle_plan(path, read_summary(path, limits), limits)
        return editable_request(plan) if as_request else plan
    if (path / POOL_DATABASE).exists() and (path / "run.sqlite3").exists():
        msg = "artifact contains both run.sqlite3 and pool.sqlite3"
        raise ValueError(msg)
    if (path / POOL_DATABASE).exists():
        plan = _read_pool_plan(artifact, path, limits)
        return editable_request(plan) if as_request else plan
    with reader(path) as connection:
        manifest = (
            checked_payload(
                connection.execute(
                    "SELECT payload,digest FROM commits WHERE revision=?",
                    (artifact.revision,),
                ).fetchone()
            )
            if isinstance(artifact, RunHandle) and artifact.revision is not None
            else latest(connection)
        )
        if isinstance(artifact, RunHandle) and manifest["run_id"] != artifact.run_id:
            msg = "RunHandle identity does not match the stored plan"
            raise ValueError(msg)
        plan = stored_plan(connection, max_identities=limits.identities)
        if isinstance(artifact, PoolHandle):
            msg = "PoolHandle does not identify a prepared pool"
            raise TypeError(msg)
        if plan.plan_id != manifest["plan_id"]:
            msg = "stored generation does not match the run plan identity"
            raise ArtifactIntegrityError(msg, artifact=path)
        if as_request:
            return _request_at_revision(plan, manifest)
        return plan


def _request_at_revision(
    plan: GenerationPlan | MatrixPlan, manifest: dict
) -> RequestReport:
    """Preserve native lineage where the editable request has one design recipe."""
    report = editable_request(plan)
    if isinstance(plan, MatrixPlan):
        return report
    return RequestReport(
        replace(
            report.request,
            lineage=Lineage(
                RunReference(
                    manifest["run_id"],
                    manifest["plan_id"],
                    manifest["revision"],
                )
            ),
        )
    )


def _read_pool_plan(
    artifact: object, path: Path, limits: ReadLimits
) -> PreparationPlan:
    """Bind a preparation plan to its pool identity before exposing its request."""
    summary = pool_summary(path)
    if isinstance(artifact, RunHandle) or (
        isinstance(artifact, PoolHandle) and artifact.pool_id != summary.pool_id
    ):
        msg = "artifact handle does not match the stored pool"
        raise ValueError(msg)
    with reader(path, filename=POOL_DATABASE) as connection:
        plan = stored_preparation(connection, max_identities=limits.identities)
    if plan.plan_id != summary.plan_id:
        msg = "stored preparation does not match the pool plan identity"
        raise ArtifactIntegrityError(msg, artifact=path)
    return plan


def check_plan_limits(
    limits: ReadLimits,
    *plans: GenerationPlan | PlanEvidence | PreparationPlan | MatrixPlan,
) -> None:
    """Cap embedded identities before materializing a plan export or comparison."""
    count = sum(plan_identities(plan) for plan in plans)
    if count > limits.identities:
        msg = "read_limits.identities cannot hold the plan identities"
        raise ReadLimitError(msg)
    if len(plans) > limits.records:
        msg = "read_limits.records cannot hold the plan documents"
        raise ReadLimitError(msg)


def plan_identities(
    plan: GenerationPlan | PlanEvidence | PreparationPlan | MatrixPlan,
) -> int:
    """Count retained part, requirement and lineage identities under one policy."""
    if isinstance(plan, MatrixPlan):
        return (
            len(plan.cells)
            + plan_identities(plan.base)
            + sum(plan_identities(c.plan) for c in plan.cells)
            + sum(source.identities for source in plan.request.sources.values())
        )
    return (
        plan.identity_count
        if isinstance(plan, PreparationPlan)
        else len(plan.request.parts)
        + len(plan.request.requirements)
        + len(plan.exclusions)
        + membership_size(plan.request)
    )


def editable_request(
    plan: GenerationPlan | PlanEvidence | PreparationPlan | MatrixPlan,
) -> RequestReport:
    """Expose normalized rules without silently discarding frozen lineage."""
    from dense_arrays.parts import BoundParts  # noqa: PLC0415

    if isinstance(plan, MatrixPlan):
        return RequestReport(
            plan.request.with_changes(base=editable_request(plan.base).request)
        )
    request = plan.request
    if isinstance(plan, (GenerationPlan, PlanEvidence)):
        evidence = plan.evidence if isinstance(plan, GenerationPlan) else plan
        if evidence.input_digests:
            request = replace(
                request,
                parts=BoundParts(
                    request.parts,
                    plan.import_report,
                    evidence.input_digests,
                    tuple(i.path for i in plan.inputs)
                    if isinstance(plan, GenerationPlan)
                    and plan.inputs
                    and plan.parent is None
                    else None,
                ),
            )
    if isinstance(plan, (GenerationPlan, PlanEvidence)) and plan.parent is not None:
        return RequestReport(
            replace(
                request,
                lineage=Lineage(
                    RunReference(
                        plan.parent.run_id, plan.parent.plan_id, plan.parent.revision
                    )
                ),
                exclude=LibraryExclusion(
                    plan.parent,
                    {"default": "default"},
                    "exact_sequence_per_cell.v1",
                ),
            )
        )
    return RequestReport(request)
