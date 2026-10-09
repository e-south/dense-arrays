"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/selections.py

Compose selection materialization with the shared inspection operations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays.artifacts import RunHandle
from dense_arrays.artifacts.reading import ReadLimits

if TYPE_CHECKING:
    from dense_arrays.artifacts.cursors import Cursor
    from dense_arrays.artifacts.reading import ReadCost
    from dense_arrays.reporting.quality import QualityReport
from dense_arrays.reporting.design_filters import DesignFilter
from dense_arrays.reporting.selections import LibrarySelection, SelectionSnapshot
from dense_arrays.reporting.selections.materialization import materialize
from dense_arrays.reporting.selections.views import SelectionView


def inspect_selection(
    artifact: object, selected: object, limits: ReadLimits
) -> SelectionSnapshot:
    """Resolve eligibility once, then materialize allocation over the full query."""
    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    if isinstance(artifact, SelectionSnapshot):
        if selected is not None:
            msg = "a saved selection already declares its membership"
            raise ValueError(msg)
        return artifact
    if isinstance(selected, SelectionSnapshot):
        # Reuse requires availability of the exact saved revisions, never a new draw.
        query = selection_view(
            artifact, selected, view="designs", limit=None, limits=limits, after=None
        )
        with query.records() as rows:
            for _ in rows:
                pass
        return selected
    if selected is None or isinstance(selected, DesignFilter):
        selected = LibrarySelection(filter=selected or DesignFilter())
    if not isinstance(selected, LibrarySelection):
        msg = "selection inspection requires LibrarySelection or DesignFilter"
        raise TypeError(msg)
    query = inspect(
        artifact, view="designs", all=True, select=selected.filter, read_limits=limits
    )
    return materialize(query, selected)


def selection_view(  # noqa: PLR0913 - resolved inspection options
    artifact: object,
    selected: object,
    *,
    view: str,
    limit: int | None,
    limits: ReadLimits,
    after: Cursor | None,
) -> SelectionView | QualityReport:
    """Locate the saved source revisions before constructing lazy record reads."""
    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    if isinstance(artifact, SelectionSnapshot):
        if selected is not None:
            msg = "a saved selection already declares its membership"
            raise ValueError(msg)
        snapshot = artifact
        supplied = [s.path for s in snapshot.sources]
    else:
        snapshot = (
            selected
            if isinstance(selected, SelectionSnapshot)
            else inspect_selection(artifact, selected, limits)
        )
        supplied = list(artifact) if isinstance(artifact, (list, tuple)) else [artifact]
    if len(supplied) != len(snapshot.sources):
        msg = "saved selection requires its original source count and order"
        raise ValueError(msg)
    bound = []
    for source, binding in zip(supplied, snapshot.sources, strict=True):
        path = source.path if isinstance(source, RunHandle) else Path(source)
        if isinstance(source, RunHandle) and (
            source.run_id != binding.source_id
            or source.revision not in {None, binding.revision}
        ):
            msg = "RunHandle conflicts with the saved selection revision"
            raise ValueError(msg)
        bound.append(
            RunHandle(path, binding.source_id, binding.revision)
            if binding.kind == "run"
            else path
        )
    query = inspect(
        bound[0] if len(bound) == 1 else bound,
        view="designs",
        all=True,
        read_limits=limits,
    )
    if view == "quality":
        from dense_arrays.reporting.quality import QualityReport  # noqa: PLC0415

        summaries = tuple(inspect(source, read_limits=limits) for source in bound)
        return QualityReport(
            query, summaries, limit=limit, after=after, snapshot=snapshot
        )
    return SelectionView(query, snapshot, view, limit, after)


def selection_cost(
    artifact: object, selected: object, limits: ReadLimits
) -> ReadCost | None:
    """Estimate materialization before a frontend initiates its population scan."""
    from dataclasses import replace  # noqa: PLC0415

    from dense_arrays.workflow.operations import inspect  # noqa: PLC0415

    if isinstance(selected, SelectionSnapshot):
        return selection_view(
            artifact, selected, view="designs", limit=None, limits=limits, after=None
        ).cost
    predicate = selected.filter if isinstance(selected, LibrarySelection) else selected
    query = inspect(
        artifact, view="designs", all=True, select=predicate, read_limits=limits
    )
    count = len(query.inputs) if hasattr(query, "inputs") else 1
    estimate = query.cost.records_estimate
    return replace(
        query.cost,
        projection="selection",
        mode="scan",
        records_estimate=None if estimate is None else estimate + count,
    )


def inspect_selected(  # noqa: PLR0913 - shared public query options
    artifact: object,
    *,
    selected: object,
    view: str,
    verify: bool,
    limit: int | None,
    all_rows: bool,
    read_limits: ReadLimits | None,
    after: str | None,
    compare: object,
) -> SelectionSnapshot | SelectionView | QualityReport:
    """Validate query controls once before materializing or reusing membership."""
    from dense_arrays.workflow.operations import _read_options  # noqa: PLC0415

    cursor = _read_options(
        view=view,
        verify=verify,
        limit=limit,
        all_rows=all_rows,
        read_limits=read_limits,
        after=after,
    )
    if verify or compare is not None:
        msg = "selection views do not accept verify or compare"
        raise ValueError(msg)
    limits = read_limits or ReadLimits()
    if view == "selection":
        if limit is not None or all_rows or after is not None:
            msg = "selection materialization does not accept limit, all or after"
            raise ValueError(msg)
        return inspect_selection(artifact, selected, limits)
    return selection_view(
        artifact,
        selected,
        view=view,
        limit=None if all_rows else (100 if limit is None else limit),
        limits=limits,
        after=cursor,
    )
