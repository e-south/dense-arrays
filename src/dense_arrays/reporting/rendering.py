"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/rendering.py

Render native design evidence through the existing playback presentation owner.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
from contextlib import closing
from dataclasses import replace
from importlib import import_module
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE
from dense_arrays.artifacts.bundles.storage import read_evidence
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.artifacts.receipts import ExportReceipt
from dense_arrays.artifacts.run_plans import cell_plans, validate_run_binding
from dense_arrays.artifacts.store import checked_payload, reader, stored_plan
from dense_arrays.planning import (
    GC,
    Avoid,
    Fixed,
    GenerationPlan,
    GroupCoverage,
    Occurrences,
    PlanEvidence,
    Spacing,
)
from dense_arrays.playback.models import PlaybackNotice
from dense_arrays.playback.output import publish_exports
from dense_arrays.playback.presentation import PlaybackDocument
from dense_arrays.playback.reconstruction import reconstruct_playback
from dense_arrays.playback.theme import PlaybackPresentation
from dense_arrays.realized import DeclaredConstraint
from dense_arrays.reporting.bundles.reading import read_bundle
from dense_arrays.reporting.bundles.views import BundleView
from dense_arrays.reporting.collections import LibraryView
from dense_arrays.reporting.collections.reading import read_library
from dense_arrays.reporting.pools import PoolQualityReport, PoolQualitySnapshot
from dense_arrays.reporting.quality import QualitySnapshot
from dense_arrays.reporting.readers import RecordView, read_records
from dense_arrays.reporting.selections.views import SelectionView, read_selection

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path

    from dense_arrays.artifacts.records import Design
    from dense_arrays.reporting.quality import QualityReport


def render_quality(
    report: QualityReport | QualitySnapshot | PoolQualityReport | PoolQualitySnapshot,
    out: Path,
) -> ExportReceipt:
    """Publish an aggregate view of the same snapshot report available to inspection."""
    if out.exists() or out.is_symlink():
        msg = f"output destination already exists: {out}"
        raise FileExistsError(msg)
    pool = isinstance(report, (PoolQualityReport, PoolQualitySnapshot))
    view = "preparation-quality" if pool else "library-quality"
    if out.suffix.lower() != ".png":
        msg = f"{view} rendering requires a .png output"
        raise ValueError(msg)
    try:
        import_module("matplotlib.figure")
        if pool:
            from dense_arrays.playback.quality.preparation import (  # noqa: PLC0415
                preparation_figure as quality_figure,
            )
        else:
            from dense_arrays.playback.quality import quality_figure  # noqa: PLC0415
    except ImportError as err:
        msg = (
            "render requires optional playback dependencies; "
            "install dense-arrays[playback]"
        )
        raise ValueError(msg) from err
    value = report.to_dict()
    records = value["counts"]["retained"] if pool else value["selection"]["designs"]
    description = (
        f"Preparation quality for {records} retained parts ({value['state']})"
        if pool
        else f"Library quality for {records} selected designs"
    )
    report_digest = semantic_digest(value)
    rendered = []

    def publish(paths: Mapping[str, Path]) -> None:
        target = paths["quality.png"]
        figure = quality_figure(value)
        try:
            figure.savefig(
                target,
                dpi=160,
                metadata={
                    "Description": description,
                    "DenseArraysReport": canonical_json(value),
                    "ReportDigest": report_digest,
                },
            )
        finally:
            figure.clear()
        payload = target.read_bytes()
        rendered.append(
            {
                "name": out.name,
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        )

    publish_exports(
        None
        if isinstance(report, (QualitySnapshot, PoolQualitySnapshot))
        else report.path / "pool.sqlite3"
        if pool
        else report.inputs[0].path
        / (
            "bundle.sqlite3"
            if isinstance(report.inputs[0], BundleView)
            else "run.sqlite3"
        ),
        {"quality.png": out},
        publish,
        replace=False,
    )
    return ExportReceipt(
        str(out),
        "png",
        view,
        records,
        tuple({**source, "report_sha256": report_digest} for source in report.sources),
        (),
        tuple(rendered),
        selection=None if pool else value["selection"].get("snapshot"),
    )


def validate_design_destination(out: Path, *, view: str) -> None:
    """Reject unsupported outputs before reading source evidence."""
    if out.exists() or out.is_symlink():
        msg = f"output destination already exists: {out}"
        raise FileExistsError(msg)
    if view != "design" or out.suffix.lower() != ".png":
        msg = "design rendering requires view='design' with a .png output"
        raise ValueError(msg)


def render_design(
    query: RecordView | LibraryView | BundleView | SelectionView, out: Path
) -> ExportReceipt:
    """Render exactly one selected design at the query's bound source revisions."""
    validate_design_destination(out, view="design")
    budget = ReadBudget(query.read_limits)
    design = _selected_design(query, budget)
    plan, origin_plan = _selected_plan(query, design, budget)
    document = _document(design, plan)
    try:
        for module in ("matplotlib.pyplot", "networkx", "PIL.Image"):
            import_module(module)
        from dense_arrays.playback.matplotlib_renderer import (  # noqa: PLC0415 - optional renderer
            render_collection_poster_png,
        )
    except ImportError as err:
        msg = (
            "render requires optional playback dependencies; "
            "install dense-arrays[playback]"
        )
        raise ValueError(msg) from err

    rendered = []

    def publish(paths: Mapping[str, Path]) -> None:
        target = paths["design.png"]
        render_collection_poster_png((document,), target)
        payload = target.read_bytes()
        rendered.append(
            {
                "name": out.name,
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        )

    publish_exports(None, {"design.png": out}, publish, replace=False)
    return ExportReceipt(
        str(out),
        "png",
        "design",
        1,
        tuple(
            {
                **source,
                **(
                    {"plan_id": origin_plan}
                    if source.get("run_id") == design.run_id and origin_plan is not None
                    else {}
                ),
            }
            for source in query.sources
        ),
        (design.reference,),
        tuple(rendered),
        selection=query.snapshot.summary()
        if isinstance(query, SelectionView)
        else None,
    )


def _selected_design(
    query: RecordView | LibraryView | BundleView | SelectionView, budget: ReadBudget
) -> Design:
    """Read at most two matches; no implicit truncation establishes uniqueness."""
    iterator = (
        read_selection(query, budget)
        if isinstance(query, SelectionView)
        else read_library(query, budget)
        if isinstance(query, LibraryView)
        else read_bundle(query, budget)
        if isinstance(query, BundleView)
        else read_records(query, budget)
    )
    with closing(iterator):
        design = next(iterator, None)
        second = next(iterator, None)
    if design is None or second is not None:
        msg = (
            "design rendering requires exactly one selected design; "
            "use a design ID, filter or saved selection"
        )
        raise ValueError(msg)
    return design


def _selected_plan(
    query: RecordView | LibraryView | BundleView | SelectionView,
    design: Design,
    budget: ReadBudget,
) -> tuple[GenerationPlan | PlanEvidence, str | None]:
    """Bind geometry to its source plan using the same aggregate read allowance."""
    sources = (
        query.inputs if isinstance(query, (LibraryView, SelectionView)) else (query,)
    )
    for source in sources:
        if isinstance(source, BundleView):
            if design.plan_id not in source.summary.manifest["plans"]:
                continue
            with reader(source.path, filename=BUNDLE_DATABASE) as connection:
                return read_evidence(connection, design.plan_id, budget), None
        if source.run_id != design.run_id:
            continue
        with reader(source.path) as connection:
            row = connection.execute(
                "SELECT payload,digest FROM commits WHERE revision=?",
                (source.revision,),
            ).fetchone()
            budget.examine(None if row is None else row[0])
            manifest = checked_payload(row)
            budget.examine()
            plan = stored_plan(
                connection, max_identities=budget.limits.identities - budget.identities
            )
            validate_run_binding(plan, manifest)
        child = cell_plans(plan).get(design.cell_id)
        if (
            manifest["run_id"] != design.run_id
            or child is None
            or child.plan_id != design.plan_id
        ):
            msg = "stored design does not match its committed run and plan"
            raise ValueError(msg)
        return child, plan.plan_id
    msg = "selected design has no matching source plan"
    raise ValueError(msg)


def _document(design: Design, plan: GenerationPlan | PlanEvidence) -> PlaybackDocument:
    constraints = []
    by_part = {p.feature_id: p for p in design.realized.placements}
    for rule in plan.request.requirements:
        if isinstance(rule, Spacing):
            if rule.min < 0:
                msg = (
                    f"{rule.id}: playback v1 cannot render negative spacing; "
                    "inspect stored requirements"
                )
                raise ValueError(msg)
            constraints.append(
                DeclaredConstraint(
                    rule.id,
                    by_part[rule.upstream].placement_id,
                    by_part[rule.downstream].placement_id,
                    rule.min,
                    rule.max,
                )
            )
    evidence = {r["id"]: r for r in design.requirements}
    if set(evidence) != {r.id for r in plan.request.requirements}:
        msg = "stored requirement evidence does not match its plan"
        raise ValueError(msg)
    notices = tuple(
        PlaybackNotice(
            code="stored_requirement",
            message=_requirement_text(rule, evidence[rule.id]),
        )
        for rule in plan.request.requirements
    )
    realized = replace(design.realized, constraints=tuple(constraints))
    return PlaybackDocument(
        reconstruct_playback(realized, notices=notices),
        title=f"Dense array {design.design_id}",
        subtitle=f"Source: {design.reference}",
        presentation=PlaybackPresentation(show_authority_notice=True),
        label_overrides={p.placement_id: p.feature_id for p in realized.placements},
    )


def _requirement_text(rule: object, evidence: Mapping[str, object]) -> str:
    observed = evidence["observed"]
    status = "PASSED" if evidence["passed"] else "FAILED"
    if isinstance(rule, Fixed):
        detail = f"{rule.part_id} {rule.orientation}, start {observed['start']}"
        if rule.start is not None:
            upper = rule.start.max if rule.start.max is not None else "unbounded"
            detail += f" (allowed {rule.start.min or 0}..{upper})"
    elif isinstance(rule, Spacing):
        detail = f"spacing {observed} bp (required {rule.min}..{rule.max} bp)"
    elif isinstance(rule, GC):
        detail = (
            "padding GC not applicable: no added bases"
            if observed is None
            else (
                f"{rule.scope} GC {observed['fraction']:.1%} "
                f"({observed['gc_bases']}/{observed['bases']} bases; "
                f"required {rule.min:.1%}..{rule.max:.1%})"
            )
        )
    elif isinstance(rule, Avoid):
        detail = (
            f"{len(observed)} forbidden matches; "
            f"{', '.join(rule.patterns)}; {rule.strands} strands"
        )
        if rule.except_placements:
            detail += (
                f"; fixed interval exceptions: {', '.join(rule.except_placements)}"
            )
    elif isinstance(rule, GroupCoverage):
        detail = f"{observed} distinct groups (required at least {rule.min})"
    elif isinstance(rule, Occurrences):
        upper = rule.max if rule.max is not None else "unbounded"
        detail = f"{observed} supplied occurrences (required {rule.min or 0}..{upper})"
    else:
        msg = "unsupported requirement for presentation"
        raise TypeError(msg)
    return f"{status} {rule.id}: {detail}"
