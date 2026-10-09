"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/bundle_inspection.py

Compose bundle metadata, verification and lazy record queries.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays.artifacts.bundles.storage import read_summary
from dense_arrays.artifacts.errors import integrity_boundary
from dense_arrays.reporting.bundles import BundleView
from dense_arrays.reporting.bundles.verification import verify_bundle
from dense_arrays.reporting.plans.bundles import read_bundle_plan
from dense_arrays.reporting.quality import QualityReport
from dense_arrays.workflow.plans import inspect_plan

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.artifacts.bundles.models import BundleSummary
    from dense_arrays.artifacts.cursors import Cursor
    from dense_arrays.artifacts.reading import ReadLimits
    from dense_arrays.planning import PlanEvidence
    from dense_arrays.reporting.design_filters import DesignFilter
    from dense_arrays.reporting.plans import PlanComparison, PlanFilter
    from dense_arrays.reporting.plans.requests import RequestReport


def inspect_bundle(  # noqa: PLR0913 - one shared inspect operation
    path: Path,
    *,
    view: str,
    verify: bool,
    limit: int | None,
    all_rows: bool,
    read_limits: ReadLimits,
    after: Cursor | None,
    selected: DesignFilter | PlanFilter | None,
    compare: object = None,
) -> (
    BundleSummary
    | BundleView
    | QualityReport
    | PlanEvidence
    | PlanComparison
    | RequestReport
):
    """Retain collection scope and reject unavailable search-history projections."""
    if compare is not None and view != "plan":
        msg = "bundle inspection does not compare generation plans"
        raise ValueError(msg)
    if view not in {
        "summary",
        "designs",
        "sequences",
        "placements",
        "quality",
        "plan",
        "plans",
        "request",
        "batches",
    }:
        msg = (
            "bundle views support summary, designs, sequences, placements, "
            "quality, plan, plans, batches and request"
        )
        raise ValueError(msg)
    if view == "summary" and (
        limit is not None or all_rows or after is not None or selected is not None
    ):
        msg = "bundle summary does not accept record filters or pagination"
        raise ValueError(msg)
    summary = read_summary(path, read_limits)
    if view in {"plan", "request"}:
        if verify or limit is not None or all_rows or after is not None:
            msg = "single-plan inspection does not accept verification or pagination"
            raise ValueError(msg)
        evidence = read_bundle_plan(path, summary, read_limits, selected)
        return inspect_plan(evidence, compare, read_limits, view=view)
    if verify:
        with integrity_boundary(path):
            summary = replace(summary, verification=verify_bundle(path, summary))
    if view == "summary":
        return summary
    if view == "quality":
        return QualityReport(
            BundleView(path, summary, "designs", None, selected, read_limits),
            (summary,),
            limit=100 if limit is None else limit,
            after=after,
        )
    return BundleView(
        path,
        summary,
        view,
        None if all_rows else (100 if limit is None else limit),
        selected,
        read_limits,
        after,
    )
