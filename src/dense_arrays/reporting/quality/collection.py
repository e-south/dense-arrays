"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/collection.py

Collect selected composition while retaining each source run's native status.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter
from contextlib import closing
from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays.artifacts.bundles.models import BundleSummary
from dense_arrays.artifacts.errors import ArtifactIntegrityError, integrity_boundary
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.reporting.collections import LibraryView
from dense_arrays.reporting.collections.filtering import resolve_filter
from dense_arrays.reporting.collections.reading import read_library
from dense_arrays.reporting.collections.sources import read_designs
from dense_arrays.reporting.readers import RecordView
from dense_arrays.reporting.selections.membership import Membership
from dense_arrays.reporting.selections.views import SelectionView, check_sources

from .accumulation import QualityAccumulator
from .evidence import source_qualities
from .models import QUALITY_POLICY, QUALITY_SCHEMA
from .search import collect_search

if TYPE_CHECKING:
    from dense_arrays.reporting.summary import RunSummary

    from .models import QualityReport


def collect(
    report: QualityReport, *, budget: ReadBudget | None = None
) -> dict[str, object]:
    """Scan complete source evidence and apply design filters only to composition."""
    with integrity_boundary(report.inputs[0].path if len(report.inputs) == 1 else None):
        return _collect(report, budget or ReadBudget(report.read_limits))


def _collect(report: QualityReport, budget: ReadBudget) -> dict[str, object]:
    """Accumulate one population while retaining the caller's work counters."""
    initial_examined, initial_entries = budget.examined, budget.identities
    membership = _membership(report, budget)
    runs = tuple(replace(r, select=None) for r in report.inputs)
    qualities = source_qualities(report, budget)
    selected = resolve_filter(runs, report.query.select, budget)
    raw = replace(report.query, select=None)
    iterator = (
        read_library(raw, budget)
        if isinstance(raw, LibraryView)
        else read_designs(raw, budget)
    )
    seen = Counter()
    source_sequences = {run_id: set() for run_id in report.origins}
    with closing(iterator) as designs:
        for design in designs:
            cell_ref = f"{design.run_id}/{design.cell_id}"
            if cell_ref not in qualities:
                msg = "unsupported run or cell in quality population"
                raise ArtifactIntegrityError(msg)
            quality = qualities[cell_ref]
            if design.plan_id != quality.plan.plan_id:
                msg = "quality design provenance does not match its source"
                raise ArtifactIntegrityError(msg)
            seen[design.run_id] += 1
            if design.sequence_id not in source_sequences[design.run_id]:
                budget.retain()
                source_sequences[design.run_id].add(design.sequence_id)
            matches = (
                membership.matches(design)
                if membership
                else selected is None
                or selected.matches(design, quality.parts, quality.collection_id)
            )
            if matches:
                quality.observe(design)
    if membership is not None:
        membership.finish()
    summaries = {
        s.run_id: s for s in report.summaries if not isinstance(s, BundleSummary)
    }
    _reconcile_designs(seen, summaries)
    search, source_runs = collect_search(
        report,
        qualities,
        budget,
        {run_id: len(values) for run_id, values in source_sequences.items()},
        seen,
    )
    combined = not isinstance(report.query, RecordView)
    selected_population = (
        combined or report.query.select is not None or report.snapshot is not None
    )
    multi_cell = len(qualities) > 1
    if combined or multi_cell:
        quality = QualityAccumulator(None, budget)
        for cell_ref, cell in qualities.items():
            quality.merge(cell, cell_ref=cell_ref)
    else:
        quality = next(iter(qualities.values()))
    first = next(iter(report.origins.values()))
    return {
        "schema": QUALITY_SCHEMA,
        "policy": QUALITY_POLICY,
        "run_id": None if combined else first.run_id,
        "revision": report.cursor.revision,
        "status": "exact",
        "population": "selected_accepted_designs_at_revisions"
        if combined
        else "filtered_accepted_designs_at_revision"
        if report.query.select or report.snapshot
        else "all_accepted_designs_at_revision",
        "cell_ref": None if combined or multi_cell else next(iter(qualities)),
        "attainment": None if combined else source_runs[0]["attainment"],
        "selection": {
            "query_id": report.query.cursor.query_id,
            "sources": list(report.sources),
            "filter": None
            if report.query.select is None
            else report.query.select.to_dict(),
            "designs": quality.accepted,
            "distinct_sequences": len(quality.sequences),
            **({"snapshot": report.snapshot.summary()} if report.snapshot else {}),
        },
        "source_runs": source_runs,
        **quality.to_dict(selected=selected_population),
        "cells": [
            {
                "cell_ref": cell_ref,
                "selected_designs": cell.accepted,
                "distinct_sequences": len(cell.sequences),
                **cell.to_dict(selected=selected_population),
            }
            for cell_ref, cell in qualities.items()
        ],
        "search": search,
        "examined": budget.examined - initial_examined,
        "state_entries": budget.identities - initial_entries,
    }


def _membership(report: QualityReport, budget: ReadBudget) -> Membership | None:
    """Bind a saved panel before collecting source-wide denominators."""
    membership = None
    if report.snapshot is not None:
        check_sources(SelectionView(report.query, report.snapshot), budget)
        membership = Membership(report.snapshot, budget)
    return membership


def _reconcile_designs(seen: Counter, summaries: dict[str, RunSummary]) -> None:
    """Require each native design population to match its committed counters."""
    for run_id, summary in summaries.items():
        if seen[run_id] != summary.accepted:
            msg = "quality designs do not reconcile with committed counters"
            raise ArtifactIntegrityError(msg)
