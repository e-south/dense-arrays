"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/reading.py

Bounded, closeable array reads with snapshot-bound continuation cursors.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import semantic_digest
from dense_arrays.artifacts.cursors import Cursor
from dense_arrays.artifacts.errors import integrity_boundary
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits
from dense_arrays.artifacts.store import checked_payload
from dense_arrays.playback.serialization import realized_array_from_dict
from dense_arrays.reporting.readers import Records

from .geometry import validate_array
from .models import ArrayFilter, ArrayRecord, CollectionSummary, sequence_identity
from .projections import ArrayPlacement, ArraySequence, CollectionPart, project
from .storage import read_parts, read_summary, reader

if TYPE_CHECKING:
    import sqlite3
    from collections.abc import Iterator, Mapping
    from pathlib import Path

    from dense_arrays.parts import Part
    from dense_arrays.realized import RealizedArray


@dataclass(frozen=True)
class CollectionView:
    """An immutable descriptor; each iterator owns its own read-only connection."""

    path: Path
    summary: CollectionSummary
    view: str
    limit: int | None = 100
    select: ArrayFilter = field(default_factory=ArrayFilter)
    read_limits: ReadLimits = field(default_factory=ReadLimits)
    after: Cursor | None = None
    revision: int = field(default=0, init=False)

    @property
    def cost(self) -> ReadCost:
        """Account for the catalog and the maximum number of examined arrays."""
        selected = (
            len(self.select.array_ids) if self.select.array_ids else self.summary.arrays
        )
        count = self.summary.parts + (0 if self.view == "parts" else selected)
        return ReadCost(
            self.summary.collection_id,
            0,
            "indexed" if self.select.array_ids else "scan",
            self.view,
            count,
            self.read_limits,
        )

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Describe source identity without exposing its filesystem location."""
        return ({"collection_id": self.summary.collection_id},)

    @property
    def cursor(self) -> Cursor:
        """Bind continuation to content, view and predicate."""
        return Cursor(
            self.summary.collection_id,
            0,
            semantic_digest({"view": self.view, "filter": self.select.to_dict()}),
        )

    def records(
        self,
    ) -> Records[ArrayRecord | ArraySequence | ArrayPlacement | CollectionPart]:
        """Return an independently closeable iterator without starting a scan."""
        budget = ReadBudget(self.read_limits)
        return Records(_read(self, budget), budget, self.cursor, self.limit)


def _validate_filter(
    connection: sqlite3.Connection,
    query: CollectionView,
    parts: Mapping[str, Part],
    budget: ReadBudget,
) -> None:
    selected = query.select
    budget.retain(
        sum(
            len(getattr(selected, name)) for name in ("array_ids", "part_ids", "groups")
        )
    )
    unknown = set(selected.part_ids) - parts.keys()
    if unknown or set(selected.groups) - {part.group for part in parts.values()}:
        msg = "array filter references an unknown part or group"
        raise ValueError(msg)
    for identity in selected.array_ids:
        if (
            connection.execute(
                "SELECT 1 FROM arrays WHERE array_id=?", (identity,)
            ).fetchone()
            is None
        ):
            msg = f"unknown array ID {identity!r}"
            raise ValueError(msg)


def _matches(
    array: RealizedArray, selected: ArrayFilter, parts: Mapping[str, Part]
) -> bool:
    identifiers = {placement.feature_id for placement in array.placements}
    return (
        not selected.part_ids or bool(identifiers.intersection(selected.part_ids))
    ) and (
        not selected.groups
        or any(parts[identity].group in selected.groups for identity in identifiers)
    )


def _query_arrays(
    connection: sqlite3.Connection, query: CollectionView, after: Cursor
) -> sqlite3.Cursor:
    sql = (
        "SELECT ordinal,array_id,sequence_id,payload,digest "
        "FROM arrays WHERE ordinal>=?"
    )
    params = [after.ordinal]
    if query.select.array_ids:
        sql += (
            " AND array_id IN (" + ",".join("?" for _ in query.select.array_ids) + ")"
        )
        params.extend(query.select.array_ids)
    return connection.execute(sql + " ORDER BY ordinal", params)


def _read(
    query: CollectionView, budget: ReadBudget
) -> Iterator[ArrayRecord | ArraySequence | ArrayPlacement | CollectionPart]:
    if (
        read_summary(query.path, query.read_limits).collection_id
        != query.summary.collection_id
    ):
        msg = "array collection changed after query binding"
        raise ValueError(msg)
    cursor = query.after or query.cursor
    if (cursor.source_id, cursor.revision, cursor.query_id, cursor.bindings) != (
        query.cursor.source_id,
        0,
        query.cursor.query_id,
        (),
    ):
        msg = "cursor does not match this array collection query"
        raise ValueError(msg)
    with reader(query.path) as connection, integrity_boundary(query.path):
        parts = read_parts(connection, query.summary, budget)
        _validate_filter(connection, query, parts, budget)
        if query.view == "parts":
            yield from _part_rows(query, budget, cursor, parts)
            return
        for ordinal, identity, sequence_id, payload, checksum in _query_arrays(
            connection, query, cursor
        ):
            budget.examine(payload)
            array = realized_array_from_dict(checked_payload((payload, checksum)))
            validate_array(array, parts)
            _check_identity(array, identity, sequence_id)
            if not _matches(array, query.select, parts):
                continue
            record = ArrayRecord(query.summary.collection_id, array)
            for offset, result in enumerate(project(record, query.view, parts), 1):
                if ordinal == cursor.ordinal and offset <= cursor.offset:
                    continue
                budget.position, budget.offset = ordinal, offset
                yield result
                if query.limit is not None and budget.returned >= query.limit:
                    return
    if (
        read_summary(query.path, query.read_limits).collection_id
        != query.summary.collection_id
    ):
        msg = "array collection changed during reading"
        raise ValueError(msg)


def _part_rows(
    query: CollectionView, budget: ReadBudget, cursor: Cursor, parts: Mapping[str, Part]
) -> Iterator[CollectionPart]:
    """Page the source catalog independently of realized-array membership."""
    for ordinal, part in enumerate(parts.values(), 1):
        if (
            ordinal <= cursor.ordinal
            or (query.select.part_ids and part.part_id not in query.select.part_ids)
            or (query.select.groups and part.group not in query.select.groups)
        ):
            continue
        budget.position = ordinal
        yield CollectionPart(query.summary.collection_id, part)
        if query.limit is not None and budget.returned >= query.limit:
            return


def _check_identity(array: RealizedArray, identity: str, sequence_id: str) -> None:
    if array.source_id != identity or sequence_identity(array.sequence) != sequence_id:
        msg = "collection array identity differs"
        raise ValueError(msg)
