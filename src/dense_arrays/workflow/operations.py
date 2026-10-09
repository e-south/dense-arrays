"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/operations.py

Thin application composition over semantic workflow owners.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import replace
from pathlib import Path
from typing import TextIO

from dense_arrays._record_validation import integer
from dense_arrays.artifacts import ExportReceipt, RunHandle
from dense_arrays.artifacts.bundles.storage import is_bundle
from dense_arrays.artifacts.cursors import Cursor
from dense_arrays.artifacts.errors import integrity_boundary
from dense_arrays.artifacts.pool_records import PoolSummary
from dense_arrays.artifacts.pools import (
    POOL_DATABASE,
    pool_summary,
    publish_pool,
    validate_filter,
    verify_pool,
)
from dense_arrays.artifacts.store import checked_payload, latest, reader
from dense_arrays.parts import PartFilter, PoolHandle, PreparationSet, PreparationSpec
from dense_arrays.planning import (
    BatchSampling,
    DesignSpec,
    ExtensionSpec,
    GenerationPlan,
    MatrixPlan,
    MatrixSpec,
    ParentRun,
    PlanEvidence,
    PreparationPlan,
)
from dense_arrays.planning.matrices.resolution import resolve_matrix
from dense_arrays.planning.preparation import resolve_preparation
from dense_arrays.planning.resolution import resolve_design
from dense_arrays.reporting import (
    AttemptFilter,
    BundleSummary,
    BundleView,
    CandidateFilter,
    DesignFilter,
    DiagnosticReport,
    LibrarySelection,
    LibraryView,
    PlanComparison,
    PlanFilter,
    PoolQualityReport,
    PoolQualitySnapshot,
    QualityComparison,
    QualityReport,
    QualitySnapshot,
    ReadLimitError,
    ReadLimits,
    RecordView,
    RequestReport,
    RunSummary,
    SelectionSnapshot,
)
from dense_arrays.reporting.verification import verify_run
from dense_arrays.workflow.plans import inspect_plan
from dense_arrays.workflow.quality import handles_quality, inspect_quality


def export(  # noqa: PLR0913 - paired Python/CLI query and output options
    artifact: RunHandle
    | SelectionSnapshot
    | PoolHandle
    | RecordView
    | LibraryView
    | BundleView
    | GenerationPlan
    | MatrixPlan
    | PlanEvidence
    | PreparationPlan
    | PreparationSpec
    | PreparationSet
    | DesignSpec
    | MatrixSpec
    | RequestReport
    | ExtensionSpec
    | QualityReport
    | PoolQualityReport
    | PoolQualitySnapshot
    | QualityComparison
    | QualitySnapshot
    | DiagnosticReport
    | PlanComparison
    | RunSummary
    | PoolSummary
    | BundleSummary
    | str
    | Path
    | list[RunHandle | str | Path]
    | tuple[RunHandle | str | Path, ...],
    *,
    out: str | Path | TextIO,
    format: str = "json",  # noqa: A002 - canonical public output option
    view: str | None = None,
    select: PartFilter
    | CandidateFilter
    | AttemptFilter
    | DesignFilter
    | PlanFilter
    | LibrarySelection
    | SelectionSnapshot
    | None = None,
    all: bool = False,  # noqa: A002 - explicit complete output scope
    read_limits: ReadLimits | None = None,
    compare: GenerationPlan
    | MatrixPlan
    | PlanEvidence
    | RunHandle
    | QualityReport
    | QualitySnapshot
    | SelectionSnapshot
    | str
    | Path
    | list[RunHandle | str | Path]
    | tuple[RunHandle | str | Path, ...]
    | None = None,
) -> ExportReceipt:
    """Export accepted evidence through the same readers used by inspect."""
    from dense_arrays.workflow.exporting import (  # noqa: PLC0415
        publish_export,
        resolve_export,
    )

    query = resolve_export(
        artifact,
        view=view,
        format_name=format,
        selected=select,
        all_rows=all,
        read_limits=read_limits,
        compare=compare,
    )
    return publish_export(query, format_name=format, out=out)


def render(
    artifact: RunHandle
    | PoolHandle
    | SelectionSnapshot
    | QualityReport
    | QualitySnapshot
    | PoolQualityReport
    | PoolQualitySnapshot
    | str
    | Path
    | list[RunHandle | str | Path]
    | tuple[RunHandle | str | Path, ...],
    *,
    out: str | Path,
    view: str = "design",
    read_limits: ReadLimits | None = None,
    select: DesignFilter | LibrarySelection | SelectionSnapshot | None = None,
) -> ExportReceipt:
    """Render persisted design or quality evidence without running generation."""
    from dense_arrays.reporting.rendering import (  # noqa: PLC0415
        render_design,
        render_quality,
        validate_design_destination,
    )

    report_types = (
        QualityReport,
        QualitySnapshot,
        PoolQualityReport,
        PoolQualitySnapshot,
    )
    if view in {"library-quality", "preparation-quality"}:
        if view == "preparation-quality" and select is not None:
            msg = (
                "preparation-quality reports cover complete recipes; "
                "filters do not apply"
            )
            raise ValueError(msg)
        report = (
            artifact
            if isinstance(artifact, report_types)
            else inspect(
                artifact, view="quality", read_limits=read_limits, select=select
            )
        )
        if isinstance(artifact, report_types) and (
            read_limits is not None or select is not None
        ):
            msg = "a quality report already binds its read limits and filter"
            raise ValueError(msg)
        expected = (
            "preparation-quality"
            if isinstance(report, (PoolQualityReport, PoolQualitySnapshot))
            else "library-quality"
        )
        if view != expected:
            msg = f"this report requires view={expected!r}"
            raise TypeError(msg)
        return render_quality(report, Path(out).absolute())
    if isinstance(artifact, report_types):
        expected = (
            "preparation-quality"
            if isinstance(artifact, (PoolQualityReport, PoolQualitySnapshot))
            else "library-quality"
        )
        msg = f"quality reports require view={expected!r}"
        raise TypeError(msg)
    destination = Path(out).absolute()
    validate_design_destination(destination, view=view)
    query = inspect(
        artifact, view="designs", select=select, limit=2, read_limits=read_limits
    )
    return render_design(query, destination)


def plan(
    request: DesignSpec | PreparationSpec | PreparationSet | ExtensionSpec | MatrixSpec,
    *,
    read_limits: ReadLimits | None = None,
) -> GenerationPlan | PreparationPlan | MatrixPlan:
    """Preview a design or preparation recipe without sampling or solving."""
    if read_limits is not None and not isinstance(read_limits, ReadLimits):
        msg = "read_limits must be ReadLimits"
        raise TypeError(msg)
    if isinstance(request, ExtensionSpec):
        from dense_arrays.workflow.extensions import resolve_extension  # noqa: PLC0415

        return resolve_extension(request, read_limits or ReadLimits())
    if (
        isinstance(request, DesignSpec)
        and request.exclude is not None
        and isinstance(request.exclude.source, ParentRun)
    ):
        from dense_arrays.workflow.extensions import (  # noqa: PLC0415
            resolve_exclusion_design,
        )

        return resolve_exclusion_design(request, read_limits or ReadLimits())
    if read_limits is not None:
        msg = "planning read_limits apply to parent-library resolution"
        raise ValueError(msg)
    if isinstance(request, MatrixSpec):
        return resolve_matrix(request)
    return (
        resolve_preparation(request)
        if isinstance(request, (PreparationSpec, PreparationSet))
        else resolve_design(request)
    )


def run(
    request: DesignSpec
    | ExtensionSpec
    | GenerationPlan
    | MatrixSpec
    | MatrixPlan
    | None = None,
    *,
    out: str | Path | None = None,
    resume: str | Path | None = None,
) -> RunHandle:
    """Generate a new run or exclusively resume its unchanged saved contract."""
    from dense_arrays.workflow.execution import execute  # noqa: PLC0415

    if resume is not None:
        from dense_arrays.workflow.recovery import resume_run  # noqa: PLC0415

        if request is not None or out is not None:
            msg = "resume is exclusive with request and out"
            raise ValueError(msg)
        return resume_run(Path(resume).absolute())
    if out is None:
        msg = "run requires an explicit out destination"
        raise ValueError(msg)
    if not isinstance(
        request, (DesignSpec, ExtensionSpec, GenerationPlan, MatrixSpec, MatrixPlan)
    ):
        msg = (
            "run requires a design request or GenerationPlan; "
            "use prepare for preparation"
        )
        raise TypeError(msg)
    resolved = (
        request if isinstance(request, (GenerationPlan, MatrixPlan)) else plan(request)
    )
    return execute(resolved, Path(out).absolute())


def prepare(  # noqa: PLR0913 - matching Python/CLI preparation options
    request: PreparationSpec
    | PreparationSet
    | PreparationPlan
    | GenerationPlan
    | MatrixPlan,
    *,
    out: str | Path,
    sampling: BatchSampling | None = None,
    batch_count: int = 1,
    attempts_per_batch: int | None = None,
    accepted_per_batch: int | None = None,
) -> PoolHandle | GenerationPlan | MatrixPlan:
    """Publish a prepared pool or save sampled batches in an executable plan."""
    if isinstance(request, (GenerationPlan, MatrixPlan)):
        from dense_arrays.workflow.batches import prepare_batches  # noqa: PLC0415

        target = Path(out).absolute()
        if target.exists() or target.is_symlink():
            msg = f"output destination already exists: {target}"
            raise FileExistsError(msg)
        return prepare_batches(
            request,
            sampling,
            target,
            batch_count=batch_count,
            attempts_per_batch=attempts_per_batch,
            accepted_per_batch=accepted_per_batch,
        )
    if (
        sampling is not None
        or batch_count != 1
        or attempts_per_batch is not None
        or accepted_per_batch is not None
    ):
        msg = "batch sampling requires a resolved generation or matrix plan"
        raise ValueError(msg)
    if not isinstance(request, (PreparationSpec, PreparationSet, PreparationPlan)):
        msg = (
            "prepare requires a preparation request or PreparationPlan; "
            "use run for designs"
        )
        raise TypeError(msg)
    target = Path(out).absolute()
    if target.exists() or target.is_symlink():
        msg = f"output destination already exists: {target}"
        raise FileExistsError(msg)
    resolved = (
        request
        if isinstance(request, PreparationPlan)
        else resolve_preparation(request)
    )
    if resolved.sampled:
        from dense_arrays.workflow.preparation import (  # noqa: PLC0415 - resolve owner after module initialization
            execute_preparation,
        )

        return execute_preparation(resolved, target)
    return publish_pool(resolved, target)


def inspect(  # noqa: PLR0913 - public query options share one operation
    artifact: RunHandle
    | DesignSpec
    | MatrixSpec
    | PreparationSpec
    | PreparationSet
    | SelectionSnapshot
    | PoolHandle
    | GenerationPlan
    | MatrixPlan
    | PlanEvidence
    | PreparationPlan
    | QualityReport
    | PoolQualityReport
    | PoolQualitySnapshot
    | QualitySnapshot
    | str
    | Path
    | list[RunHandle | str | Path]
    | tuple[RunHandle | str | Path, ...],
    *,
    view: str = "summary",
    verify: bool = False,
    limit: int | None = None,
    select: PartFilter
    | CandidateFilter
    | AttemptFilter
    | DesignFilter
    | PlanFilter
    | LibrarySelection
    | SelectionSnapshot
    | None = None,
    all: bool = False,  # noqa: A002 - public Python/CLI full-scope option
    read_limits: ReadLimits | None = None,
    after: str | None = None,
    compare: GenerationPlan
    | MatrixPlan
    | PlanEvidence
    | RunHandle
    | QualityReport
    | QualitySnapshot
    | SelectionSnapshot
    | str
    | Path
    | list[RunHandle | str | Path]
    | tuple[RunHandle | str | Path, ...]
    | None = None,
) -> (
    RunSummary
    | PoolSummary
    | BundleSummary
    | BundleView
    | SelectionSnapshot
    | RecordView
    | DiagnosticReport
    | QualityReport
    | PoolQualityReport
    | PoolQualitySnapshot
    | QualityComparison
    | QualitySnapshot
    | GenerationPlan
    | MatrixPlan
    | PlanEvidence
    | PlanComparison
    | LibraryView
    | RequestReport
    | PreparationPlan
):
    """Read requests or committed evidence without invoking generation or repairs."""
    cursor = _read_options(
        view=view,
        verify=verify,
        limit=limit,
        all_rows=all,
        read_limits=read_limits,
        after=after,
    )
    if handles_quality(artifact, view, compare):
        return inspect_quality(
            artifact,
            compare=compare,
            selected=select,
            read_limits=read_limits,
            verify=verify,
            limit=limit,
            all_rows=all,
            after=after,
        )
    if (
        view == "selection"
        or isinstance(artifact, SelectionSnapshot)
        or isinstance(select, (LibrarySelection, SelectionSnapshot))
    ):
        from dense_arrays.workflow.selections import inspect_selected  # noqa: PLC0415

        return inspect_selected(
            artifact,
            selected=select,
            view=view,
            verify=verify,
            limit=limit,
            all_rows=all,
            read_limits=read_limits,
            after=after,
            compare=compare,
        )
    if isinstance(artifact, (str, Path)) and is_bundle(Path(artifact)):
        from dense_arrays.workflow.bundle_inspection import (  # noqa: PLC0415
            inspect_bundle,
        )

        return inspect_bundle(
            Path(artifact),
            view=view,
            verify=verify,
            limit=limit,
            all_rows=all,
            read_limits=read_limits or ReadLimits(),
            after=cursor,
            selected=select,
            compare=compare,
        )
    if view in {"plan", "request"} or isinstance(
        artifact, (GenerationPlan, PlanEvidence, PreparationPlan, MatrixPlan)
    ):
        if (
            view not in {"plan", "request", "summary"}
            or select is not None
            or limit is not None
            or all
            or after is not None
            or verify
        ):
            msg = (
                "plan inspection supports comparison, not record filters, "
                "pagination or run verification"
            )
            raise ValueError(msg)
        return inspect_plan(artifact, compare, read_limits or ReadLimits(), view=view)
    if compare is not None:
        msg = "comparison currently requires view='plan'"
        raise ValueError(msg)
    if isinstance(artifact, (list, tuple)):
        return _inspect_collection(
            artifact,
            view=view,
            verify=verify,
            limit=None if all else (100 if limit is None else limit),
            selected=select,
            read_limits=read_limits or ReadLimits(),
            after=cursor,
        )
    return _inspect_native(
        artifact,
        view=view,
        verify=verify,
        limit=limit,
        selected=select,
        all_rows=all,
        read_limits=read_limits,
        after=cursor,
    )


def _inspect_native(  # noqa: PLR0913 - resolved public inspection options
    artifact: RunHandle | PoolHandle | str | Path,
    *,
    view: str,
    verify: bool,
    limit: int | None,
    selected: object,
    all_rows: bool,
    read_limits: ReadLimits | None,
    after: Cursor | None,
) -> RunSummary | PoolSummary | RecordView | QualityReport | DiagnosticReport:
    """Route one native directory to its pool or committed-run reader."""
    path = (
        artifact.path
        if isinstance(artifact, (RunHandle, PoolHandle))
        else Path(artifact)
    )

    if (path / POOL_DATABASE).exists():
        if (path / "run.sqlite3").exists():
            msg = "ambiguous directory contains both run and pool artifacts"
            raise ValueError(msg)
        return _inspect_pool(
            artifact,
            path,
            view=view,
            verify=verify,
            limit=limit,
            selected=selected,
            all_rows=all_rows,
            read_limits=read_limits or ReadLimits(),
            after=after,
        )
    if isinstance(artifact, PoolHandle):
        msg = "pool artifact is unavailable at the declared location"
        raise FileNotFoundError(msg)
    _validate_run_filter(view, selected)
    return _inspect_run(
        artifact,
        path,
        view=view,
        verify=verify,
        limit=limit,
        all_rows=all_rows,
        read_limits=read_limits or ReadLimits(),
        after=after,
        selected=selected,
    )


def _inspect_collection(  # noqa: PLR0913 - resolved public query options
    sources: list | tuple,
    *,
    view: str,
    verify: bool,
    limit: int | None,
    selected: DesignFilter | None,
    read_limits: ReadLimits,
    after: Cursor | None,
) -> LibraryView | QualityReport:
    """Bind every artifact before streaming the combined query."""
    if not sources:
        msg = "a combined library requires at least one source"
        raise ValueError(msg)
    if view not in {"designs", "sequences", "placements", "quality"} or verify:
        msg = (
            "combined sources support designs, sequences, placements and quality; "
            "verify artifacts individually"
        )
        raise ValueError(msg)
    if len(sources) > read_limits.identities:
        msg = "read_limits.identities cannot hold the source bindings"
        raise ReadLimitError(msg)
    if after is not None and len(after.bindings) != len(sources):
        msg = "cursor source count does not match the combined query"
        raise ValueError(msg)
    inputs = [
        _collection_input(
            source, None if after is None else after.bindings[index], read_limits
        )
        for index, source in enumerate(sources)
    ]
    queries, summaries = zip(*inputs, strict=True)
    if view == "quality":
        query = LibraryView(tuple(queries), "designs", None, selected, read_limits)
        return QualityReport(query, tuple(summaries), limit=limit, after=after)
    return LibraryView(tuple(queries), view, limit, selected, read_limits, after)


def _collection_input(
    source: object,
    binding: tuple[str, int] | None,
    limits: ReadLimits,
) -> tuple[RecordView | BundleView, RunSummary | BundleSummary]:
    """Bind a native prefix or immutable bundle to one union input position."""
    from dense_arrays.artifacts.bundles.storage import read_summary  # noqa: PLC0415

    if not isinstance(source, (str, Path, RunHandle)):
        msg = "combined libraries require run/bundle paths or RunHandle values"
        raise TypeError(msg)
    path = source.path if isinstance(source, RunHandle) else Path(source)
    if is_bundle(path):
        if isinstance(source, RunHandle):
            msg = "RunHandle cannot identify a portable bundle"
            raise ValueError(msg)
        summary = read_summary(path, limits)
        if binding is not None and binding != (summary.bundle_id, 0):
            msg = "bundle identity does not match the cursor source binding"
            raise ValueError(msg)
        return BundleView(path, summary, "designs", None, read_limits=limits), summary
    bound = source
    if binding is not None:
        run_id, revision = binding
        if isinstance(source, RunHandle) and (
            source.run_id != run_id or source.revision not in {None, revision}
        ):
            msg = "RunHandle conflicts with the cursor snapshot"
            raise ValueError(msg)
        bound = RunHandle(path, run_id, revision=revision)
    summary = _run_summary(bound, path, None, limits)
    return RecordView(
        path,
        summary.revision,
        "designs",
        None,
        run_id=summary.run_id,
        source_records=summary.accepted,
        read_limits=limits,
    ), summary


def _read_options(  # noqa: PLR0913 - shared query syntax
    *,
    view: str,
    verify: bool,
    limit: int | None,
    all_rows: bool,
    read_limits: ReadLimits | None,
    after: str | None,
) -> Cursor | None:
    """Validate shared query flags before native reads."""
    if not isinstance(all_rows, bool) or not isinstance(verify, bool):
        msg = "all_rows and verify must be booleans"
        raise TypeError(msg)
    if limit is not None:
        integer(limit, field_name="limit", minimum=1)
    if view == "summary" and (limit is not None or all_rows):
        msg = "limit and all apply to record views, not summary"
        raise ValueError(msg)
    if all_rows and limit is not None:
        msg = "all and limit are exclusive"
        raise ValueError(msg)
    if read_limits is not None and not isinstance(read_limits, ReadLimits):
        msg = "read_limits must be ReadLimits"
        raise TypeError(msg)
    cursor = None if after is None else Cursor.from_token(after)
    if cursor is not None and view in {"summary", "diagnostics"}:
        msg = "cursor pagination applies to record views, not aggregate reports"
        raise ValueError(msg)
    if view in {"diagnostics", "quality"} and all_rows:
        msg = (
            "aggregate reports require a bounded display limit; "
            "--all applies to records"
        )
        raise ValueError(msg)
    return cursor


def _inspect_run(  # noqa: PLR0913 - resolved read options
    artifact: object,
    path: Path,
    *,
    view: str,
    verify: bool,
    limit: int | None,
    all_rows: bool,
    read_limits: ReadLimits,
    after: Cursor | None,
    selected: AttemptFilter | DesignFilter | None,
) -> RunSummary | RecordView | DiagnosticReport | QualityReport:
    summary = _run_summary(artifact, path, after, read_limits)
    if selected is not None:
        if selected.identities > read_limits.identities:
            msg = "read_limits.identities cannot hold the attempt filter"
            raise ReadLimitError(msg)
        if isinstance(selected, AttemptFilter):
            selected.validate(summary)
    if verify:
        with integrity_boundary(path):
            evidence = verify_run(path, summary, read_limits)
        summary = replace(summary, verified=True, verification=evidence)
    if view == "summary":
        return summary
    if view == "diagnostics":
        return DiagnosticReport(
            path,
            summary,
            limit=20 if limit is None else limit,
            read_limits=read_limits,
            select=selected,
        )
    if view == "quality":
        return QualityReport(
            RecordView(
                path,
                summary.revision,
                "designs",
                None,
                run_id=summary.run_id,
                source_records=summary.accepted,
                read_limits=read_limits,
                select=selected,
            ),
            (summary,),
            limit=100 if limit is None else limit,
            after=after,
        )
    return RecordView(
        path,
        summary.revision,
        view,
        None if all_rows else (100 if limit is None else limit),
        run_id=summary.run_id,
        source_records=(
            (summary.batch_count or 0)
            if view == "batches"
            else summary.accepted
            if view in {"designs", "sequences", "placements"}
            else summary.counts["started"]
        ),
        read_limits=read_limits,
        after=after,
        select=selected,
    )


def _inspect_pool(  # noqa: PLR0913 - resolved read options
    artifact: object,
    path: Path,
    *,
    view: str,
    verify: bool,
    limit: int | None,
    selected: PartFilter | CandidateFilter | None,
    all_rows: bool,
    read_limits: ReadLimits,
    after: Cursor | None,
) -> PoolSummary | RecordView | PoolQualityReport:
    summary = replace(
        pool_summary(path, read_limits=read_limits), read_limits=read_limits
    )
    if isinstance(artifact, RunHandle) or (
        isinstance(artifact, PoolHandle) and artifact.pool_id != summary.pool_id
    ):
        msg = "artifact handle does not match the stored pool"
        raise ValueError(msg)
    if view == "quality" and summary.preparation is not None:
        if selected is not None or after is not None or limit is not None or all_rows:
            msg = "pool quality binds the complete candidate population"
            raise ValueError(msg)
        return PoolQualityReport(path, summary, read_limits)
    if view == "candidates" and summary.preparation is None:
        msg = "this pool has no saved candidate evidence; inspect parts instead"
        raise ValueError(msg)
    if view not in {"summary", "parts", "candidates"}:
        msg = (
            f"unsupported pool view {view!r}; "
            "supported: summary, parts, sampled candidates and quality"
        )
        raise ValueError(msg)
    if view == "summary" and (limit is not None or selected is not None or all_rows):
        msg = "record filters and pagination apply to parts, not pool summary"
        raise ValueError(msg)
    _validate_pool_filter(path, summary, view, selected, read_limits)
    if verify:
        with integrity_boundary(path):
            evidence = verify_pool(path, summary, read_limits)
        summary = replace(summary, verified=True, verification=evidence)
    if view == "summary":
        return replace(summary, verified=verify)
    return RecordView(
        path,
        0,
        view,
        None if all_rows else (100 if limit is None else limit),
        pool_id=summary.pool_id,
        select=selected,
        source_records=summary.source_parts
        if view == "candidates"
        else summary.retained_parts,
        read_limits=read_limits,
        after=after,
    )


def _validate_pool_filter(
    path: Path,
    summary: PoolSummary,
    view: str,
    selected: object,
    read_limits: ReadLimits,
) -> None:
    """Resolve the predicate against its declared pool population and state cap."""
    expected = CandidateFilter if view == "candidates" else PartFilter
    if selected is not None and not isinstance(selected, expected):
        msg = f"pool {view} require {expected.__name__}"
        raise TypeError(msg)
    if selected is not None:
        count = (
            selected.identities
            if isinstance(selected, CandidateFilter)
            else len(selected.part_ids) + len(selected.groups)
        )
        if count > read_limits.identities:
            msg = "read_limits.identities cannot hold the requested filter"
            raise ReadLimitError(msg)
        if isinstance(selected, CandidateFilter):
            selected.validate(summary)
        else:
            validate_filter(path, selected)


def _run_summary(
    artifact: object, path: Path, after: Cursor | None, read_limits: ReadLimits
) -> RunSummary:
    """Bind a path, handle or cursor to one validated committed revision."""
    revision = artifact.revision if isinstance(artifact, RunHandle) else None
    if after is not None:
        if revision is not None and revision != after.revision:
            msg = "cursor revision does not match the bound RunHandle"
            raise ValueError(msg)
        revision = after.revision
    with reader(path) as connection:
        value = (
            latest(connection)
            if revision is None
            else checked_payload(
                connection.execute(
                    "SELECT payload,digest FROM commits WHERE revision=?",
                    (revision,),
                ).fetchone()
            )
        )
        summary = replace(RunSummary.from_manifest(value), read_limits=read_limits)
    if isinstance(artifact, RunHandle) and summary.run_id != artifact.run_id:
        msg = "RunHandle identity does not match the artifact"
        raise ValueError(msg)
    return summary


def _validate_run_filter(view: str, select: object) -> None:
    """Require a filter whose population matches the requested native view."""
    expected = (
        AttemptFilter
        if view in {"attempts", "diagnostics"}
        else DesignFilter
        if view in {"designs", "sequences", "placements", "quality"}
        else None
    )
    if select is not None and (expected is None or not isinstance(select, expected)):
        msg = (
            "run queries require the matching AttemptFilter or DesignFilter; "
            "this view does not support the supplied filter"
        )
        raise TypeError(msg)
