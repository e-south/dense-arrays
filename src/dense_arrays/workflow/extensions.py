"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/extensions.py

Resolve an additional target from terminal native parent evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from contextlib import closing
from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays._record_validation import semantic_digest
from dense_arrays.artifacts.reading import ReadBudget, ReadLimits
from dense_arrays.artifacts.run_plans import RunPlan, cell_plans
from dense_arrays.artifacts.store import latest, reader, stored_plan
from dense_arrays.planning import (
    Allocation,
    DesignSpec,
    GenerationPlan,
    LibraryExclusion,
    Lineage,
    MatrixPlan,
    ParentRun,
    RunReference,
    Target,
)
from dense_arrays.planning.libraries import (
    ExcludedDesign,
    ParentLibrary,
    accepted_digest,
)
from dense_arrays.planning.resolution import resolve_design
from dense_arrays.reporting.readers import RecordView, read_records
from dense_arrays.reporting.summary import RunSummary
from dense_arrays.reporting.verification import verify_run

if TYPE_CHECKING:
    from dense_arrays.planning import ExtensionSpec


def resolve_extension(request: ExtensionSpec, limits: ReadLimits) -> RunPlan:
    """Freeze verified per-cell parent and ancestor exclusions under new targets."""
    parent, libraries = read_parent_libraries(request.parent, limits)
    if isinstance(parent, MatrixPlan):
        if not isinstance(request.additional, Mapping):
            msg = "matrix extension requires explicit per-cell additional targets"
            raise TypeError(msg)
        allocation = Allocation(counts=request.additional)
        allocation.resolve(tuple(cell.cell_id for cell in parent.cells))
        origin = next(iter(libraries.values()))
        base = GenerationPlan(
            replace(
                parent.base.request,
                limits=request.limits,
                seed=request.seed,
                lineage=Lineage(
                    RunReference(origin.run_id, origin.plan_id, origin.revision)
                ),
            ),
            import_report=parent.base.import_report,
            embedded_input_digests=parent.base.evidence.input_digests,
        )
        matrix = parent.request.with_changes(
            base=base.request,
            sources={
                cell: replace(source, locations=None)
                for cell, source in parent.request.sources.items()
            },
            allocation=allocation,
            exclude={
                cell: LibraryExclusion(
                    library, {cell: cell}, "exact_sequence_per_cell.v1"
                )
                for cell, library in libraries.items()
            },
        )
        return MatrixPlan(matrix, base)
    if isinstance(request.additional, Mapping):
        msg = "single-cell extension requires an integer additional target"
        raise TypeError(msg)
    return GenerationPlan(
        replace(
            parent.request,
            target=Target(request.additional),
            limits=request.limits,
            seed=request.seed,
            exclude=None,
        ),
        parent.inputs,
        parent.import_report,
        libraries["default"],
        parent.embedded_input_digests,
    )


def resolve_exclusion_design(request: DesignSpec, limits: ReadLimits) -> GenerationPlan:
    """Resolve an explicit exclusion without inheriting any of its design rules."""
    _, library = read_parent_library(request.exclude.source, limits)
    return resolve_design(
        replace(request, exclude=replace(request.exclude, source=library))
    )


def read_parent_libraries(
    source: ParentRun, limits: ReadLimits
) -> tuple[RunPlan, dict[str, ParentLibrary]]:
    """Freeze a terminal accepted library and its already-bound ancestor exclusions."""
    import fcntl  # noqa: PLC0415 - same native locking capability as run ownership

    path = source.run.absolute()
    budget = ReadBudget(limits)
    with (path / ".writer.lock").open("rb") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as err:
            msg = (
                "accepted-library resolution requires a terminal source "
                "with no active writer"
            )
            raise ValueError(msg) from err
        with reader(path) as connection:
            summary = RunSummary.from_manifest(latest(connection))
            if summary.state not in {"completed", "stopped", "failed"}:
                msg = (
                    "accepted-library resolution requires a terminal source; "
                    "stop generation first"
                )
                raise ValueError(msg)
            budget.examine()
            parent = stored_plan(connection, max_identities=limits.identities)
            plans = cell_plans(parent)
        verify_run(path, summary, budget=budget)
        query = RecordView(
            path, summary.revision, "designs", None, run_id=summary.run_id
        )
        accepted = {cell: [] for cell in plans}
        with closing(read_records(query, budget)) as designs:
            for design in designs:
                budget.retain()
                accepted[design.cell_id].append(
                    ExcludedDesign(
                        design.reference,
                        design.sequence_id,
                        semantic_digest(design.to_dict()),
                        cell_id=design.cell_id,
                    )
                )
        libraries = {}
        for name, plan in plans.items():
            ancestors = plan.exclusions
            budget.retain(len(ancestors))
            direct = tuple(accepted[name])
            cell = summary.cell_summaries[name]
            libraries[name] = ParentLibrary(
                summary.run_id,
                parent.plan_id,
                summary.revision,
                cell.state,
                cell.target,
                cell.accepted,
                accepted_digest(direct),
                (*ancestors, *direct),
                cell_id=name,
            )
    return parent, libraries


def read_parent_library(
    source: ParentRun, limits: ReadLimits
) -> tuple[GenerationPlan, ParentLibrary]:
    """Read the default cell for a single-cell accepted-library request."""
    parent, libraries = read_parent_libraries(source, limits)
    if isinstance(parent, MatrixPlan):
        msg = "matrix exclusion requires an explicit frozen per-cell source"
        raise TypeError(msg)
    return parent, libraries["default"]
