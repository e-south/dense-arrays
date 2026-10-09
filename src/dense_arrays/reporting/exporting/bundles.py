"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/exporting/bundles.py

Export selected designs with path-free, bound verification evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
from contextlib import closing
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    canonical_json,
    mutable_json,
    semantic_digest,
)
from dense_arrays.artifacts.bundles.models import (
    BUNDLE_DATABASE,
    BUNDLE_MANIFEST,
    BUNDLE_SCHEMA,
    EVIDENCE_BOUNDARY,
    RUNTIME_BUNDLE_SCHEMA,
    RUNTIME_EVIDENCE_BOUNDARY,
)
from dense_arrays.artifacts.bundles.publication import destination
from dense_arrays.artifacts.bundles.storage import read_evidence
from dense_arrays.artifacts.publication import write_new
from dense_arrays.artifacts.reading import ReadBudget
from dense_arrays.artifacts.receipts import ExportReceipt
from dense_arrays.artifacts.records import COMPOSITION_POLICY
from dense_arrays.artifacts.run_plans import cell_plans, validate_run_binding
from dense_arrays.artifacts.run_state import origin_bindings
from dense_arrays.artifacts.store import checked_payload, reader, stored_plan
from dense_arrays.reporting.bundles import BundleView
from dense_arrays.reporting.bundles.batches import load_batches
from dense_arrays.reporting.bundles.reading import read_bundle
from dense_arrays.reporting.bundles.verification import exclusion_index, verify_design
from dense_arrays.reporting.collections import LibraryView
from dense_arrays.reporting.collections.reading import read_library
from dense_arrays.reporting.exporting.batch_evidence import copy_batches
from dense_arrays.reporting.readers import RecordView, read_records
from dense_arrays.reporting.selections.views import SelectionView, read_selection
from dense_arrays.reporting.summary import RunSummary

if TYPE_CHECKING:
    import sqlite3

    from dense_arrays.reporting.collections.sources import SourceView


def export_bundle(
    query: RecordView | LibraryView | BundleView | SelectionView, *, out: str | Path
) -> ExportReceipt:
    """Publish exactly the declared collection; original source locations stay local."""
    if not isinstance(out, (str, Path)):
        msg = "bundle export requires a new directory path"
        raise TypeError(msg)
    if query.view != "designs" or query.limit is not None or query.after is not None:
        msg = "bundle export requires a complete designs selection with all=True"
        raise ValueError(msg)
    path = Path(out).absolute()
    if path.exists() or path.is_symlink():
        msg = f"output destination already exists: {path}"
        raise FileExistsError(msg)
    sources = (
        query.inputs if isinstance(query, (LibraryView, SelectionView)) else (query,)
    )
    budget = ReadBudget(query.read_limits)
    references = []
    with destination(path) as connection:
        source_runs, plan_ids = _source_evidence(connection, sources, budget)
        plans = {
            identity: read_evidence(connection, identity, budget)
            for identity in plan_ids
        }
        batch_count = copy_batches(connection, sources, budget)
        batches = load_batches(
            connection,
            {"batches": batch_count, "source_runs": source_runs},
            plans,
            budget,
        )
        excluded = exclusion_index(plans, budget)
        source_bindings = set(origin_bindings(source_runs))
        seen_sequences = set()
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
            for ordinal, design in enumerate(iterator, 1):
                verify_design(
                    design,
                    plans[design.plan_id],
                    excluded[design.plan_id],
                    batches=batches,
                )
                key = (design.run_id, design.cell_id, design.sequence_id)
                if (
                    design.run_id,
                    design.cell_id,
                    design.plan_id,
                ) not in source_bindings or key in seen_sequences:
                    msg = "bundle design source or sequence uniqueness mismatch"
                    raise ValueError(msg)
                budget.retain(2)
                seen_sequences.add(key)
                references.append(design.reference)
                value = design.to_dict()
                connection.execute(
                    "INSERT INTO designs VALUES (?,?,?,?,?,?,?,?)",
                    (
                        ordinal,
                        design.reference,
                        design.run_id,
                        design.cell_id,
                        design.design_id,
                        design.plan_id,
                        canonical_json(value),
                        semantic_digest(value),
                    ),
                )
        connection.commit()
        database = path / BUNDLE_DATABASE
        with database.open("rb") as stream:
            sha256 = hashlib.file_digest(stream, "sha256").hexdigest()
        file = {
            "name": BUNDLE_DATABASE,
            "bytes": database.stat().st_size,
            "sha256": sha256,
        }
        manifest = {
            "schema": RUNTIME_BUNDLE_SCHEMA if batch_count else BUNDLE_SCHEMA,
            **({"batches": batch_count} if batch_count else {}),
            "scope": "selected_collection",
            "designs": len(references),
            "selection": query.snapshot.summary()
            if isinstance(query, SelectionView)
            else None
            if query.select is None
            else query.select.to_dict(),
            "sources": list(query.sources),
            "source_runs": source_runs,
            "plans": plan_ids,
            "evidence": RUNTIME_EVIDENCE_BOUNDARY if batch_count else EVIDENCE_BOUNDARY,
            "metric_policy": COMPOSITION_POLICY,
            "file": file,
        }
        manifest["bundle_id"] = semantic_digest(manifest)
        write_new(path / BUNDLE_MANIFEST, canonical_json(manifest) + "\n")
    manifest_path = path / BUNDLE_MANIFEST
    manifest_file = {
        "name": BUNDLE_MANIFEST,
        "bytes": manifest_path.stat().st_size,
        "sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
    }
    return ExportReceipt(
        str(path),
        "bundle",
        "designs",
        len(references),
        query.sources,
        tuple(references),
        (file, manifest_file),
    )


def _source_evidence(
    target: sqlite3.Connection, sources: tuple[SourceView, ...], budget: ReadBudget
) -> tuple[list[dict[str, object]], list[str]]:
    """Copy resolved evidence and original snapshot summaries, once per identity."""
    summaries = {}
    plans = {}
    for source in sources:
        if isinstance(source, BundleView):
            _copy_bundle_evidence(target, source, budget, summaries, plans)
            continue
        with reader(source.path) as connection:
            row = connection.execute(
                "SELECT payload,digest FROM commits WHERE revision=?",
                (source.revision,),
            ).fetchone()
            budget.examine(None if row is None else row[0])
            manifest = checked_payload(row)
            summary = RunSummary.from_manifest(manifest)
            if summary.run_id != source.run_id:
                msg = "source run identity changed during bundle export"
                raise ValueError(msg)
            key = (summary.run_id, summary.revision)
            if key in summaries and summaries[key] != manifest:
                msg = "conflicting source snapshot summaries"
                raise ValueError(msg)
            if key not in summaries:
                budget.retain()
                summaries[key] = manifest
            row = connection.execute("SELECT payload FROM plan WHERE id=1").fetchone()
            budget.examine(None if row is None else row[0])
            plan = stored_plan(
                connection, max_identities=budget.limits.identities - budget.identities
            )
            validate_run_binding(plan, manifest)
            for child in cell_plans(plan).values():
                if child.plan_id not in plans:
                    budget.retain()
                    value = child.evidence.to_dict()
                    plans[child.plan_id] = semantic_digest(value)
                    target.execute(
                        "INSERT INTO plans VALUES (?,?,?)",
                        (child.plan_id, canonical_json(value), plans[child.plan_id]),
                    )
    return list(summaries.values()), list(plans)


def _copy_bundle_evidence(
    target: sqlite3.Connection,
    source: BundleView,
    budget: ReadBudget,
    summaries: dict,
    plans: dict,
) -> None:
    """Carry included origin evidence forward when publishing a smaller collection."""
    budget.examine()
    for item in source.summary.manifest["source_runs"]:
        manifest = mutable_json(item)
        key = (manifest["run_id"], manifest["revision"])
        if key not in summaries:
            budget.retain()
            summaries[key] = manifest
        elif summaries[key] != manifest:
            msg = "conflicting source snapshot summaries"
            raise ValueError(msg)
    with reader(source.path, filename=BUNDLE_DATABASE) as connection:
        for identity in source.summary.manifest["plans"]:
            state = budget.identities
            plan = read_evidence(connection, identity, budget)
            value = plan.to_dict()
            budget.identities = state
            if identity not in plans:
                budget.retain()
                plans[identity] = semantic_digest(value)
                target.execute(
                    "INSERT INTO plans VALUES (?,?,?)",
                    (identity, canonical_json(value), plans[identity]),
                )
