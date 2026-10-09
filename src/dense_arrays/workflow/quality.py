"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/quality.py

Compose quality comparisons through the shared artifact and selection readers.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from pathlib import Path

from dense_arrays.artifacts.reading import ReadLimits
from dense_arrays.reporting.pools import PoolQualityReport, PoolQualitySnapshot
from dense_arrays.reporting.quality import (
    QualityComparison,
    QualityReport,
    QualitySnapshot,
)
from dense_arrays.reporting.selections import LibrarySelection


def compare_quality(
    before: object,
    after: object,
    *,
    selected: object,
    read_limits: ReadLimits | None,
) -> QualityComparison:
    """Bind two complete populations without modifying a supplied report's scope."""
    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    if isinstance(selected, LibrarySelection):
        msg = "quality comparison requires a saved selection or already-bound reports"
        raise TypeError(msg)
    reports = []
    for value in (before, after):
        source = value
        if isinstance(value, (str, Path)) and Path(value).is_file():
            from dense_arrays.workflow.inputs import read_quality  # noqa: PLC0415

            source = read_quality(Path(source), read_limits)
        if isinstance(
            source,
            (QualityReport, QualitySnapshot, PoolQualityReport, PoolQualitySnapshot),
        ):
            if selected is not None:
                msg = "a quality report already binds its selection"
                raise ValueError(msg)
            report = source
        else:
            report = inspect(
                source, view="quality", select=selected, read_limits=read_limits
            )
        reports.append(report)
    return QualityComparison(*reports, read_limits=read_limits or ReadLimits())


def saved_quality(
    source: object, limits: ReadLimits | None
) -> QualityReport | QualitySnapshot | PoolQualityReport | PoolQualitySnapshot:
    """Return already-bound values or read a native quality document without sources."""
    from dense_arrays.workflow.inputs import read_quality  # noqa: PLC0415

    if isinstance(
        source, (QualityReport, QualitySnapshot, PoolQualityReport, PoolQualitySnapshot)
    ):
        if limits is not None:
            msg = "a quality report already binds its read limits"
            raise ValueError(msg)
        return source
    return read_quality(Path(source), limits)


def handles_quality(source: object, view: str, compare: object) -> bool:
    """Route bound reports, saved reports and comparisons before native dispatch."""
    return view == "quality" and (
        compare is not None
        or isinstance(
            source,
            (QualityReport, QualitySnapshot, PoolQualityReport, PoolQualitySnapshot),
        )
        or (isinstance(source, (str, Path)) and Path(source).is_file())
    )


def inspect_quality(  # noqa: PLR0913 - existing inspect options, one quality owner
    source: object,
    *,
    compare: object,
    selected: object,
    read_limits: ReadLimits | None,
    verify: bool,
    limit: int | None,
    all_rows: bool,
    after: str | None,
) -> (
    QualityComparison
    | QualityReport
    | QualitySnapshot
    | PoolQualityReport
    | PoolQualitySnapshot
):
    """Keep comparison aggregates independent of report usage-table pagination."""
    if verify or limit is not None or all_rows or after is not None:
        msg = (
            "saved quality reports and comparisons do not accept "
            "verify, limit, all or after"
        )
        raise ValueError(msg)
    if compare is not None:
        return compare_quality(
            source, compare, selected=selected, read_limits=read_limits
        )
    if selected is not None:
        msg = "a saved quality report binds its scope; filter source artifacts instead"
        raise ValueError(msg)
    return saved_quality(source, read_limits)
