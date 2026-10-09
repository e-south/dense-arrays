"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/exporting/source_bindings.py

Bind record publications to the exact native manifests they read.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dense_arrays._record_validation import mutable_json, semantic_digest
from dense_arrays.artifacts.bundles.storage import read_summary
from dense_arrays.artifacts.pools import pool_summary
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.artifacts.store import checked_payload, reader
from dense_arrays.reporting.bundles import BundleView
from dense_arrays.reporting.collections import LibraryView
from dense_arrays.reporting.readers import RecordView
from dense_arrays.reporting.selections.views import SelectionView
from dense_arrays.reporting.summary import RunSummary

type ExportQuery = RecordView | LibraryView | BundleView | SelectionView


def bind_sources(query: ExportQuery) -> tuple[dict[str, object], ...]:
    """Read bounded manifest metadata without charging projected data-row scans."""
    inputs = (
        query.inputs if isinstance(query, (LibraryView, SelectionView)) else (query,)
    )
    ReadBudget(query.read_limits).retain(len(inputs))
    bound = []
    for source, descriptor in zip(inputs, query.sources, strict=True):
        digest = _manifest_digest(source)
        if descriptor.get("manifest_digest", digest) != digest:
            msg = "saved selection source manifest changed"
            raise ValueError(msg)
        bound.append({**descriptor, "manifest_digest": digest})
    return tuple(bound)


def check_sources(query: ExportQuery, expected: tuple[dict[str, object], ...]) -> None:
    """Recheck the same revisions before publishing a completed record stream."""
    if bind_sources(query) != expected:
        msg = "export source manifest changed during publication"
        raise ValueError(msg)


def _manifest_digest(source: RecordView | BundleView) -> str:
    """Validate identity and hash the canonical committed manifest, including schema."""
    if isinstance(source, BundleView):
        summary = read_summary(source.path, source.read_limits)
        if summary.bundle_id != source.summary.bundle_id:
            msg = "export source bundle identity changed"
            raise ValueError(msg)
        value = mutable_json(summary.manifest)
    elif source.pool_id is not None:
        summary = pool_summary(source.path, read_limits=source.read_limits)
        if summary.pool_id != source.pool_id:
            msg = "export source pool identity changed"
            raise ValueError(msg)
        value = summary.to_dict()
    else:
        with reader(source.path) as connection:
            value = checked_payload(
                connection.execute(
                    "SELECT payload,digest FROM commits WHERE revision=?",
                    (source.revision,),
                ).fetchone()
            )
        summary = RunSummary.from_manifest(value)
        if summary.run_id != source.run_id or summary.revision != source.revision:
            msg = "export source run identity changed"
            raise ValueError(msg)
    return semantic_digest(value)
