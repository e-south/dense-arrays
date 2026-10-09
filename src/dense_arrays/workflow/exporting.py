"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/exporting.py

Resolve one shared Python/CLI export query before disclosing cost and writing.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays.artifacts.reading import ReadCost
from dense_arrays.reporting import LibraryView, ReadLimits, RecordView
from dense_arrays.reporting.bundles import BundleView
from dense_arrays.reporting.exporting import export_records, validate_format
from dense_arrays.reporting.exporting.documents import (
    DOCUMENT_TYPES,
    DocumentView,
    document_view,
    export_document,
)
from dense_arrays.reporting.quality import QualityReport
from dense_arrays.reporting.selections import LibrarySelection, SelectionSnapshot
from dense_arrays.reporting.selections.views import SelectionView

if TYPE_CHECKING:
    from typing import TextIO

    from dense_arrays.artifacts import ExportReceipt
    from dense_arrays.parts import PartFilter
    from dense_arrays.reporting import (
        AttemptFilter,
        CandidateFilter,
        DesignFilter,
        PlanFilter,
    )


def _validate_scope(
    view: str, *, all_rows: bool, read_limits: ReadLimits | None
) -> None:
    """Distinguish a complete record selection from one bounded document."""
    if not isinstance(all_rows, bool):
        msg = "all must be a boolean"
        raise TypeError(msg)
    if read_limits is not None and not isinstance(read_limits, ReadLimits):
        msg = "read_limits must be ReadLimits"
        raise TypeError(msg)
    if view in DOCUMENT_TYPES and all_rows:
        msg = (
            "all applies to complete record exports; "
            "document exports retain their declared report scope"
        )
        raise ValueError(msg)
    if view not in DOCUMENT_TYPES and not all_rows:
        msg = "record export requires all=True to declare its complete scope"
        raise ValueError(msg)


def resolve_export(  # noqa: PLR0913 - shared operation and CLI options
    artifact: object,
    *,
    view: str | None,
    format_name: str,
    all_rows: bool,
    selected: PartFilter
    | CandidateFilter
    | AttemptFilter
    | DesignFilter
    | PlanFilter
    | LibrarySelection
    | SelectionSnapshot
    | None,
    read_limits: ReadLimits | None,
    limit: int | None = None,
    compare: object = None,
) -> RecordView | LibraryView | BundleView | SelectionView | DocumentView:
    """Bind the declared scope once; reports keep their explicit display bounds."""
    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    if not isinstance(all_rows, bool):
        msg = "all must be a boolean"
        raise TypeError(msg)
    bound = isinstance(
        artifact, (RecordView, LibraryView, BundleView, SelectionView, DocumentView)
    )
    view = view or (
        artifact.view
        if bound
        else "designs"
        if isinstance(artifact, SelectionSnapshot) and format_name != "selection"
        else document_view(artifact) or "designs"
    )
    validate_format(view, format_name)
    bounded = (
        isinstance(artifact, (SelectionSnapshot, SelectionView))
        or isinstance(selected, SelectionSnapshot)
        or (isinstance(selected, LibrarySelection) and selected.take is not None)
    )
    if bounded and all_rows:
        msg = "all contradicts an explicitly bounded or saved selection"
        raise ValueError(msg)
    if format_name == "selection" or view == "selection":
        _validate_scope(
            "designs", all_rows=all_rows or bounded, read_limits=read_limits
        )
        if limit is not None or compare is not None:
            msg = "selection export does not accept limit or compare"
            raise ValueError(msg)
        value = inspect(
            artifact, view="selection", select=selected, read_limits=read_limits
        )
        return DocumentView(value, "selection", read_limits or ReadLimits())
    _validate_scope(
        view,
        all_rows=all_rows or (bounded and view not in DOCUMENT_TYPES),
        read_limits=read_limits,
    )
    if compare is not None:
        artifact, view, read_limits = _comparison_source(
            artifact, view, compare, read_limits, selected
        )
        selected = None
    if (
        view == "request"
        and isinstance(artifact, (str, Path))
        and Path(artifact).is_file()
    ):
        from dense_arrays.workflow.inputs import read_source  # noqa: PLC0415

        artifact = read_source(Path(artifact))
    if bound:
        if artifact.view != view or any(
            v is not None for v in (selected, read_limits, limit)
        ):
            msg_0 = "an export view already binds its view, filter and read limits"
            raise ValueError(msg_0)
        return artifact
    if view in DOCUMENT_TYPES:
        return _resolve_document(
            artifact, view=view, selected=selected, limit=limit, read_limits=read_limits
        )
    return inspect(
        artifact,
        view=view,
        select=selected,
        all=True,
        read_limits=read_limits,
        limit=limit,
    )


def _comparison_source(
    artifact: object,
    view: str,
    compare: object,
    limits: ReadLimits | None,
    selected: object,
) -> tuple[object, str, ReadLimits | None]:
    """Reuse shared comparison reports and their already-bound work limits."""
    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    if view not in {"plan", "quality"}:
        msg = "comparison export requires view='plan' or 'quality'"
        raise ValueError(msg)
    value = inspect(
        artifact, view=view, compare=compare, read_limits=limits, select=selected
    )
    return value, "comparison", None if hasattr(value, "read_limits") else limits


def export_cost(
    query: RecordView | LibraryView | BundleView | DocumentView, format_name: str
) -> ReadCost:
    """Include evidence reads in the bundle publication estimate."""
    cost = query.cost
    if format_name != "bundle":
        return cost
    extra = (
        1 + 2 * len(query.summary.manifest["plans"])
        if isinstance(query, BundleView)
        else sum(
            1 + 2 * len(source.summary.manifest["plans"])
            if isinstance(source, BundleView)
            else 3
            for source in (
                query.inputs
                if isinstance(query, (LibraryView, SelectionView))
                else (query,)
            )
        )
    )
    return ReadCost(
        cost.source_id,
        cost.revision,
        "scan",
        "bundle",
        None if cost.records_estimate is None else cost.records_estimate + extra,
        cost.limits,
    )


def publish_export(
    query: RecordView | LibraryView | BundleView | DocumentView,
    *,
    format_name: str,
    out: str | Path | TextIO,
) -> ExportReceipt:
    """Write the same bound query described to either frontend."""
    if format_name == "bundle":
        from dense_arrays.reporting.exporting.bundles import (  # noqa: PLC0415
            export_bundle,
        )

        receipt = export_bundle(query, out=out)
    elif isinstance(query, DocumentView):
        receipt = export_document(query, out=out)
    else:
        receipt = export_records(query, format_name=format_name, out=out)
    snapshot = (
        query.snapshot
        if isinstance(query, SelectionView)
        else query.value
        if isinstance(query, DocumentView)
        and isinstance(query.value, SelectionSnapshot)
        else query.value.snapshot
        if isinstance(query, DocumentView) and isinstance(query.value, QualityReport)
        else None
    )
    if snapshot is not None:
        receipt = replace(
            receipt,
            format=format_name,
            design_refs=tuple(snapshot.references()),
            selection=snapshot.summary(),
        )
    return receipt


def _resolve_document(
    artifact: object,
    *,
    view: str,
    selected: object,
    limit: int | None,
    read_limits: ReadLimits | None,
) -> DocumentView:
    """Preserve an already bound report or compose one scoped document query."""
    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    if document_view(artifact) is not None and not isinstance(
        artifact, SelectionSnapshot
    ):
        if selected is not None or limit is not None:
            msg_0 = "a typed document already binds its filter and display limit"
            raise ValueError(msg_0)
        if hasattr(artifact, "read_limits") and read_limits is not None:
            msg_0 = "a typed report already binds its read limits"
            raise ValueError(msg_0)
        value = artifact
    else:
        value = inspect(
            artifact,
            view=view,
            select=selected,
            limit=limit,
            read_limits=read_limits,
        )
    return DocumentView(
        value, view, read_limits or getattr(value, "read_limits", ReadLimits())
    )
