"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/evidence.py

Bind contained designs to original run summaries and resolved plan evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE, BundleSummary
from dense_arrays.artifacts.bundles.storage import read_evidence
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.artifacts.run_plans import cell_plans
from dense_arrays.artifacts.store import reader, stored_plan
from dense_arrays.reporting.bundles.views import BundleView
from dense_arrays.reporting.summary import RunSummary

from .accumulation import QualityAccumulator

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadLimits
    from dense_arrays.planning import PlanEvidence
    from dense_arrays.reporting.collections.sources import SourceView

    from .models import QualityReport


def run_summaries(summary: RunSummary | BundleSummary) -> Iterator[RunSummary]:
    """Retain original scope when the supplied artifact contains only a subset."""
    if isinstance(summary, BundleSummary):
        yield from (
            RunSummary.from_manifest(s) for s in summary.manifest["source_runs"]
        )
    else:
        yield replace(summary, verified=False, verification=None)


def origins(
    inputs: tuple[SourceView, ...],
    summaries: tuple[RunSummary | BundleSummary, ...],
    limits: ReadLimits,
) -> dict[str, RunSummary]:
    """Require one consistent committed origin revision across the input union."""
    budget = ReadBudget(limits)
    result = {}
    for source, summary in zip(inputs, summaries, strict=True):
        if isinstance(source, BundleView):
            valid = (
                isinstance(summary, BundleSummary)
                and source.summary.bundle_id == summary.bundle_id
            )
        else:
            valid = isinstance(summary, RunSummary) and (
                source.run_id,
                source.revision,
            ) == (summary.run_id, summary.revision)
        if not valid:
            msg = "quality source does not match its summary"
            raise ValueError(msg)
        for origin in run_summaries(summary):
            if origin.run_id not in result:
                budget.retain()
                result[origin.run_id] = origin
            elif result[origin.run_id] != origin:
                msg = "quality requires one consistent revision per run"
                raise ValueError(msg)
    return result


def source_qualities(
    report: QualityReport, budget: ReadBudget
) -> dict[str, QualityAccumulator]:
    """Read every cell's bound evidence while retaining original run ownership."""
    qualities = {}
    for source, summary in zip(report.inputs, report.summaries, strict=True):
        for origin in run_summaries(summary):
            if isinstance(source, BundleView):
                plans = None
            else:
                budget.examine()
                with reader(source.path) as connection:
                    saved = stored_plan(
                        connection,
                        max_identities=budget.limits.identities - budget.identities,
                    )
                if saved.plan_id != origin.plan_id:
                    msg = "quality source plan does not match its run"
                    raise ArtifactIntegrityError(msg)
                plans = cell_plans(saved)
            for cell_id, cell in origin.cell_summaries.items():
                before = budget.identities
                if plans is None:
                    if cell.plan_id not in summary.manifest["plans"]:
                        msg = "bundle origin references an undeclared plan"
                        raise ArtifactIntegrityError(msg)
                    with reader(source.path, filename=BUNDLE_DATABASE) as connection:
                        plan = read_evidence(connection, cell.plan_id, budget)
                else:
                    plan = plans[cell_id].evidence
                if (
                    plan.plan_id != cell.plan_id
                    or plan.request.target.count != cell.target
                ):
                    msg = "quality plan identity or target does not match its cell"
                    raise ArtifactIntegrityError(msg)
                key = f"{origin.run_id}/{cell_id}"
                budget.identities = before
                _retain_quality(qualities, key, plan, budget)
    return qualities


def _retain_quality(
    qualities: dict, key: str, plan: PlanEvidence, budget: ReadBudget
) -> None:
    """Deduplicate cell origins while keeping conflicting evidence visible."""
    if key not in qualities:
        budget.retain(1 + len(plan.exclusions))
        qualities[key] = QualityAccumulator(plan, budget)
    elif qualities[key].plan != plan:
        msg = "conflicting source plans for the same cell"
        raise ArtifactIntegrityError(msg)
