"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/collections/sources.py

Read native and bundled design evidence through explicit source operations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from contextlib import closing, contextmanager
from typing import TYPE_CHECKING

from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE
from dense_arrays.artifacts.bundles.storage import read_evidence
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.run_plans import cell_plans
from dense_arrays.artifacts.run_state import manifest_cell_ids
from dense_arrays.artifacts.store import reader, stored_plan
from dense_arrays.reporting.bundles.views import BundleView
from dense_arrays.reporting.design_queries import bundle_design_aliases, design_aliases
from dense_arrays.reporting.readers import RecordView, read_records
from dense_arrays.reporting.summary import RunSummary

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.artifacts.records import Design
    from dense_arrays.parts.models import Part

type SourceView = RecordView | BundleView
type Annotations = dict[str, tuple[dict[str, Part], str]]


def read_designs(source: SourceView, budget: ReadBudget) -> Iterator[Design]:
    """Open the contained designs at the bound source snapshot."""
    if isinstance(source, BundleView):
        from dense_arrays.reporting.bundles.reading import read_bundle  # noqa: PLC0415

        iterator = read_bundle(source, budget)
    else:
        iterator = read_records(source, budget)
    count = 0
    with closing(iterator) as records:
        for design in records:
            count += 1
            yield design
    expected = design_count(source)
    if expected is not None and count != expected:
        msg = "contained designs do not reconcile with the source snapshot"
        raise ArtifactIntegrityError(msg, artifact=source.path)


def cells(source: SourceView) -> Iterator[str]:
    """Expose declared cell namespaces, including cells with no included designs."""
    if isinstance(source, BundleView):
        yield from (
            f"{s['run_id']}/{cell}"
            for s in source.summary.manifest["source_runs"]
            for cell in manifest_cell_ids(s)
        )
    else:
        from dense_arrays.artifacts.store import checked_payload  # noqa: PLC0415

        with reader(source.path) as connection:
            value = checked_payload(
                connection.execute(
                    "SELECT payload,digest FROM commits WHERE revision=?",
                    (source.revision,),
                ).fetchone()
            )
        summary = RunSummary.from_manifest(value)
        yield from (f"{source.run_id}/{cell}" for cell in summary.cell_ids)


def aliases(
    source: SourceView, requested: tuple[str, ...]
) -> Iterator[tuple[str, str]]:
    """Resolve only requested design labels using the artifact's identity index."""
    if isinstance(source, BundleView):
        with reader(source.path, filename=BUNDLE_DATABASE) as connection:
            yield from bundle_design_aliases(connection, requested)
    else:
        with reader(source.path) as connection:
            yield from design_aliases(connection, source.revision, requested)


@contextmanager
def annotations(
    source: SourceView, budget: ReadBudget, *, needed: bool
) -> Iterator[Annotations]:
    """Hold one source's plan annotations only while its records need them."""
    lookup = {}
    before = budget.identities
    if needed:
        if isinstance(source, BundleView):
            with reader(source.path, filename=BUNDLE_DATABASE) as connection:
                for identity in source.summary.manifest["plans"]:
                    plan = read_evidence(connection, identity, budget)
                    budget.retain()
                    lookup[identity] = (
                        {p.part_id: p for p in plan.request.parts},
                        plan.collection_id,
                    )
        else:
            budget.examine()
            with reader(source.path) as connection:
                plan = stored_plan(
                    connection,
                    max_identities=budget.limits.identities - budget.identities,
                )
            for child in cell_plans(plan).values():
                budget.retain(1 + len(child.request.parts))
                lookup[child.plan_id] = (
                    {p.part_id: p for p in child.request.parts},
                    child.collection_id,
                )
    retained = budget.identities - before
    try:
        yield lookup
    finally:
        budget.identities -= retained


def annotation_cost(source: SourceView) -> int:
    """Count the plan documents needed for this source's annotations."""
    return (
        len(source.summary.manifest["plans"]) if isinstance(source, BundleView) else 1
    )


def design_count(source: SourceView) -> int | None:
    """Describe contained designs independently of original source attainment."""
    return (
        source.summary.designs
        if isinstance(source, BundleView)
        else source.source_records
    )
