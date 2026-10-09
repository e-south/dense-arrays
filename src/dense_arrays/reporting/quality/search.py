"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/search.py

Describe checked search histories without inventing records absent from a bundle.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
from contextlib import closing
from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.reporting.accounting import AttemptTotals
from dense_arrays.reporting.bundles.views import BundleView
from dense_arrays.reporting.readers import read_records

if TYPE_CHECKING:
    from collections import Counter

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.reporting.readers import RecordView
    from dense_arrays.reporting.summary import RunSummary

    from .accumulation import QualityAccumulator
    from .models import QualityReport


def collect_search(
    report: QualityReport,
    qualities: dict[str, QualityAccumulator],
    budget: ReadBudget,
    distinct_sequences: dict[str, int],
    included: Counter,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Validate native histories once per origin, checking every supplied copy."""
    combined = AttemptTotals(budget=budget)
    histories, fingerprints = {}, {}
    for source, summary in zip(report.inputs, report.summaries, strict=True):
        if isinstance(source, BundleView):
            continue
        before = budget.identities
        totals, fingerprint = _history(source, summary, budget)
        if summary.run_id in fingerprints:
            if fingerprints[summary.run_id] != fingerprint:
                msg = f"conflicting search evidence for run {summary.run_id}"
                raise ArtifactIntegrityError(msg, artifact=source.path)
            budget.identities = before
            continue
        budget.retain()
        fingerprints[summary.run_id] = fingerprint
        combined.merge(totals)
        histories[summary.run_id] = {
            **totals.to_dict(),
            "active_seconds": summary.active_seconds,
            "termination_reason": summary.termination_reason,
        }
    sources = [
        _origin_result(
            summary,
            sum(
                qualities[f"{run_id}/{cell_id}"].accepted
                for cell_id in summary.cell_ids
            ),
            histories.get(run_id),
            included[run_id],
            distinct_sequences[run_id],
        )
        for run_id, summary in report.origins.items()
    ]
    unavailable = [
        f"{s.run_id}/{s.revision}"
        for s in report.origins.values()
        if s.run_id not in histories
    ]
    availability = (
        "complete" if not unavailable else "partial" if histories else "not_included"
    )
    search = {
        **(combined.to_dict() if histories else dict.fromkeys(combined.to_dict())),
        "availability": availability,
        "population": "all_attempts_in_source_snapshots"
        if not unavailable
        else "attempts_in_available_native_sources",
        "unavailable_source_refs": unavailable,
        "active_seconds": sum(h["active_seconds"] for h in histories.values())
        if histories
        else None,
        "termination_reason": next(iter(histories.values()))["termination_reason"]
        if len(report.origins) == 1 and histories
        else None,
    }
    return search, sources


def _origin_result(
    summary: RunSummary,
    selected: int,
    search: dict | None,
    included: int,
    distinct_sequences: int,
) -> dict[str, object]:
    """Keep original attainment and contained population as separate observations."""
    return {
        "run_id": summary.run_id,
        "revision": summary.revision,
        "plan_id": summary.plan_id,
        "state": summary.state,
        **(
            {"producer": summary.producer.to_dict()}
            if summary.producer is not None
            else {}
        ),
        "attainment": {
            "target": summary.target,
            "accepted": summary.accepted,
            "shortfall": summary.target - summary.accepted,
            "distinct_sequences": distinct_sequences if search is not None else None,
        },
        "included_designs": included,
        "selected_designs": selected,
        "search": search,
        "evidence": "native_records"
        if search is not None
        else "included_designs_and_origin_manifest",
    }


def _history(
    source: RecordView, summary: RunSummary, budget: ReadBudget
) -> tuple[AttemptTotals, bytes]:
    """Reconcile one native attempt population and fingerprint its ordered records."""
    fingerprint = hashlib.sha256()
    totals = AttemptTotals(budget=budget)
    attempts = replace(source, view="attempts", select=None)
    with closing(read_records(attempts, budget)) as records:
        for attempt in records:
            if attempt.cell_id not in summary.cell_ids:
                msg = "unsupported cell in source attempt population"
                raise ArtifactIntegrityError(msg)
            totals.observe(attempt)
            fingerprint.update((canonical_json(attempt.to_dict()) + "\n").encode())
    totals.reconcile(summary.counts)
    return totals, fingerprint.digest()
