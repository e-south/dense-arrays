"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/pools.py

Publish immutable prepared pools in one coherent local transaction.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
import sqlite3
from contextlib import closing
from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.pool_records import (
    POOL_IDENTITY_SCHEMA,
    PoolPart,
    PoolSummary,
)
from dense_arrays.artifacts.preparation.bands import accounting_bands_size
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.reading import (
    ReadBudget,
    ReadLimitError,
    ReadLimits,
    Verification,
)
from dense_arrays.artifacts.store import checked_payload, reader
from dense_arrays.parts.ingestion import ImportedParts, validate_parts
from dense_arrays.parts.pools import PoolHandle
from dense_arrays.parts.provenance import ImportReport
from dense_arrays.planning.preparation import PreparationPlan

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from dense_arrays.parts.filters import PartFilter
    from dense_arrays.parts.pools import PoolSource

POOL_DATABASE = "pool.sqlite3"


def pool_identity(plan: PreparationPlan) -> str:
    """Name the prepared collection independently of destination paths."""
    return semantic_digest(
        {"schema": POOL_IDENTITY_SCHEMA, "preparation_plan_id": plan.plan_id}
    )


def publish_pool(plan: PreparationPlan, out: Path) -> PoolHandle:
    """Recheck source bytes and exclusively publish retained parts and summary."""
    import fcntl  # noqa: PLC0415 - local lock capability before publication

    if out.exists() or out.is_symlink():
        msg = f"output destination already exists: {out}"
        raise FileExistsError(msg)
    plan.verify_inputs()
    producer = Producer.capture()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=False)
    handle = PoolHandle(out, pool_identity(plan))
    with (out / ".writer.lock").open("xb") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with closing(sqlite3.connect(out / POOL_DATABASE)) as connection:
            connection.execute("PRAGMA synchronous=FULL")
            connection.executescript("""
                CREATE TABLE manifest (id INTEGER PRIMARY KEY CHECK(id=1),
                    payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE TABLE preparation (id INTEGER PRIMARY KEY CHECK(id=1),
                    payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE TABLE parts (ordinal INTEGER PRIMARY KEY,
                    part_id TEXT UNIQUE NOT NULL, group_name TEXT,
                    payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE INDEX parts_group ON parts(group_name);
            """)
            summary = PoolSummary(
                handle.pool_id,
                plan.plan_id,
                len(plan.parts),
                len(plan.retained_indices),
                producer=producer,
            )
            with connection:
                for ordinal, index in enumerate(plan.retained_indices, 1):
                    part = plan.parts[index]
                    payload = PoolPart(handle.pool_id, ordinal, part).to_dict()
                    connection.execute(
                        "INSERT INTO parts VALUES (?,?,?,?,?)",
                        (
                            ordinal,
                            part.part_id,
                            part.group,
                            canonical_json(payload),
                            semantic_digest(payload),
                        ),
                    )
                payload = plan.to_dict()
                connection.execute(
                    "INSERT INTO preparation VALUES (1,?,?)",
                    (canonical_json(payload), semantic_digest(payload)),
                )
                payload = summary.to_dict()
                connection.execute(
                    "INSERT INTO manifest VALUES (1,?,?)",
                    (canonical_json(payload), semantic_digest(payload)),
                )
    return handle


def pool_summary(path: Path, *, read_limits: ReadLimits | None = None) -> PoolSummary:
    """Read exactly one committed manifest and validate its declared family."""
    with reader(path, filename=POOL_DATABASE) as connection:
        value = checked_payload(
            connection.execute(
                "SELECT payload,digest FROM manifest WHERE id=1"
            ).fetchone()
        )
        ReadBudget(read_limits or ReadLimits()).retain(
            accounting_bands_size(value.get("preparation"))
        )
        return PoolSummary.from_dict(value)


def stored_preparation(
    connection: sqlite3.Connection, *, max_identities: int | None = None
) -> PreparationPlan:
    """Read the bound snapshot without reopening the original source table."""
    value = checked_payload(
        connection.execute(
            "SELECT payload,digest FROM preparation WHERE id=1"
        ).fetchone()
    )
    return decode_preparation(value, max_identities=max_identities)


def decode_preparation(
    value: dict, *, max_identities: int | None = None, base: Path | None = None
) -> PreparationPlan:
    """Bound embedded model state for both native stores and saved plan files."""
    count = _preparation_identities(value)
    if max_identities is not None and count > max_identities:
        msg = "read_limits.identities cannot hold the stored preparation"
        raise ReadLimitError(msg)
    return PreparationPlan.from_dict(value, base=base)


def _preparation_identities(value: dict) -> int:
    """Admit embedded recipe models before constructing resolved objects."""
    if value.get("schema") == "dense_arrays.preparation_plan.v3":
        recipes = value.get("recipes")
        if not isinstance(recipes, list) or any(
            not isinstance(item, dict)
            or not isinstance(item.get("plan"), dict)
            or item["plan"].get("schema") != "dense_arrays.preparation_plan.v2"
            for item in recipes
        ):
            msg = "stored set requires complete sampled recipe plans"
            raise ValueError(msg)
        return sum(1 + _preparation_identities(item["plan"]) for item in recipes)
    if value.get("schema") == "dense_arrays.preparation_plan.v1":
        if not isinstance(value.get("parts"), list):
            msg = "stored preparation requires a parts array"
            raise TypeError(msg)
        return len(value["parts"])

    def motif_size(source: dict) -> int:
        return (
            1
            + len(source.get("motif", {}).get("probabilities", []))
            + (
                1
                + len(
                    source.get("scoring", {}).get("motif", {}).get("probabilities", [])
                )
                if "window" in source
                else 0
            )
        )

    return (
        motif_size(value.get("source", {}))
        + len(value.get("request", {}).get("screening", []))
        + len(
            value.get("request", {}).get("score_bands", {}).get("upper_fractions", [])
        )
        + sum(
            motif_size(motif)
            for screen in value.get("screens", [])
            for motif in screen["motifs"]
        )
    )


def validate_filter(path: Path, selected: PartFilter | None) -> None:
    """Resolve requested identities through indexes without loading part sequences."""
    if selected is None:
        return
    found = []
    with reader(path, filename=POOL_DATABASE) as connection:
        for query, labels in (
            ("SELECT 1 FROM parts WHERE part_id=? LIMIT 1", selected.part_ids),
            ("SELECT 1 FROM parts WHERE group_name=? LIMIT 1", selected.groups),
        ):
            values = {
                label
                for label in labels
                if connection.execute(query, (label,)).fetchone() is not None
            }
            found.append(values)
    selected.validate_available(*found)


def iter_parts(  # noqa: PLR0913 - explicit immutable query and work budget
    path: Path,
    *,
    pool_id: str,
    selected: PartFilter | None = None,
    limit: int | None = None,
    budget: ReadBudget | None = None,
    after: int = 0,
) -> Iterator[PoolPart]:
    """Stream a stable pool in retained order and validate every returned record."""
    with reader(path, filename=POOL_DATABASE) as connection:
        summary = PoolSummary.from_dict(
            checked_payload(
                connection.execute(
                    "SELECT payload,digest FROM manifest WHERE id=1"
                ).fetchone()
            )
        )
        if summary.pool_id != pool_id:
            msg = "pool identity changed since selecting the snapshot"
            raise ValueError(msg)
        emitted = 0
        for ordinal, part_id, group, payload, checksum in connection.execute(
            "SELECT ordinal,part_id,group_name,payload,digest FROM parts "
            "WHERE ordinal>? ORDER BY ordinal",
            (after,),
        ):
            if budget is not None:
                budget.examine(payload)
            record = PoolPart.from_dict(checked_payload((payload, checksum)))
            if (record.pool_id, record.ordinal, record.part_id, record.part.group) != (
                pool_id,
                ordinal,
                part_id,
                group,
            ):
                msg = "pool record does not match its identity or indexed fields"
                raise ValueError(msg)
            if selected is None or selected.matches(record.part):
                if budget is not None:
                    budget.position = ordinal
                yield record
                emitted += 1
                if limit is not None and emitted >= limit:
                    return


def verify_pool(
    path: Path, summary: PoolSummary, limits: ReadLimits | None = None
) -> Verification:
    """Verify preparation and retained records without reopening source inputs."""
    budget = ReadBudget(limits or ReadLimits())
    with reader(path, filename=POOL_DATABASE) as connection:
        payload = connection.execute(
            "SELECT payload FROM preparation WHERE id=1"
        ).fetchone()[0]
        budget.examine(payload)
        plan = stored_preparation(connection, max_identities=budget.limits.identities)
        if plan.sampled:
            from dense_arrays.artifacts.preparation.verification import (  # noqa: PLC0415 - resolve owner after module initialization
                verify_sampled,
            )

            verification, _ = verify_sampled(connection, path, plan, summary, budget)
            return verification
        budget.retain(len(plan.parts))
    expected = PoolSummary(
        pool_identity(plan),
        plan.plan_id,
        len(plan.parts),
        len(plan.retained_indices),
        producer=summary.producer,
    )
    if replace(summary, verified=False) != expected:
        msg = "pool summary does not match its preparation plan"
        raise ValueError(msg)
    records = iter_parts(path, pool_id=summary.pool_id, budget=budget)
    try:
        for ordinal, (index, record) in enumerate(
            zip(plan.retained_indices, records, strict=True), 1
        ):
            if record.ordinal != ordinal or record.part != plan.parts[index]:
                msg = "retained part does not match the bound preparation source"
                raise ValueError(msg)
    finally:
        records.close()
    return Verification(
        "preparation_and_retained_parts", budget.examined, budget.bytes_checked
    )


def read_pool_source(source: PoolSource) -> ImportedParts:
    """Resolve one verified immutable pool selection for generation planning."""
    path = source.path
    summary = pool_summary(path)
    if isinstance(source.pool, PoolHandle) and source.pool.pool_id != summary.pool_id:
        msg = "PoolHandle identity does not match the stored pool"
        raise ValueError(msg)
    if summary.state != "completed":
        msg = (
            f"pool is incomplete: {path}; inspect it with --view quality, "
            "then revise preparation and write a completed pool to a new destination"
        )
        raise ValueError(msg)
    before = hashlib.sha256((path / POOL_DATABASE).read_bytes()).hexdigest()
    verify_pool(path, summary)
    validate_filter(path, source.select)
    values = tuple(
        r.part
        for r in iter_parts(path, pool_id=summary.pool_id, selected=source.select)
    )
    after = hashlib.sha256((path / POOL_DATABASE).read_bytes()).hexdigest()
    if before != after:
        msg = "pool changed while resolving inputs; create a new plan"
        raise ValueError(msg)
    values = validate_parts(values)
    return ImportedParts(
        values,
        after,
        report=ImportReport(
            "pool", len(values), collection_id=summary.pool_id, selection=source.select
        ),
    )
