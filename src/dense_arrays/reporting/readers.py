"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/readers.py

Snapshot-bound record iterators with explicit early-close ownership.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    integer,
    required_text,
    semantic_digest,
)
from dense_arrays.artifacts.batches import BatchDecision
from dense_arrays.artifacts.batches.records import DECISION_SCHEMA
from dense_arrays.artifacts.cursors import Cursor
from dense_arrays.artifacts.objectives import encoded_objective_size
from dense_arrays.artifacts.preparation.candidates import (
    CANDIDATE_RECORD_SCHEMA,
    PoolCandidate,
)
from dense_arrays.artifacts.reading import (
    ReadBudget,
    ReadCost,
    ReadLimitError,
    ReadLimits,
)
from dense_arrays.artifacts.records import DESIGN_RECORD_SCHEMA, Attempt, Design
from dense_arrays.artifacts.store import checked_payload, reader
from dense_arrays.parts.filters import PartFilter
from dense_arrays.planning.batches.bindings import encoded_batch_size
from dense_arrays.reporting.design_filters import DesignFilter
from dense_arrays.reporting.design_queries import design_context
from dense_arrays.reporting.filters import AttemptFilter
from dense_arrays.reporting.pools.filters import CandidateFilter
from dense_arrays.reporting.projections import PlacementRecord, SequenceRecord, project
from dense_arrays.reporting.summary import RunSummary

if TYPE_CHECKING:
    import sqlite3
    from pathlib import Path
    from types import TracebackType

    from dense_arrays.artifacts.pool_records import PoolPart


class Records[T](Iterator[T]):
    """An independent lazy iterator, closable explicitly or as a context manager."""

    def __init__(
        self,
        iterator: Iterator[T],
        budget: ReadBudget,
        cursor: Cursor,
        limit: int | None,
    ) -> None:
        """Own one lazy generator, opening no files until iteration begins."""
        self._iterator = iterator
        self._budget = budget
        self._cursor = cursor
        self._limit = limit
        self._exhausted = False
        self._closed = False

    def __next__(self) -> T:
        """Read one row; exhaustion closes the underlying generator's resources."""
        if self._closed:
            raise StopIteration
        try:
            value = next(self._iterator)
        except StopIteration:
            self._exhausted = True
            self.close()
            raise
        except BaseException:
            self.close()
            raise
        self._budget.returned += 1
        return value

    @property
    def examined(self) -> int:
        """Data records decoded, including records rejected by the predicate."""
        return self._budget.examined

    @property
    def returned(self) -> int:
        """Rows emitted by this independent iterator."""
        return self._budget.returned

    @property
    def next_cursor(self) -> str | None:
        """Continue after the last emitted row; a full page may be the final page."""
        if not self.returned or (
            self._exhausted and (self._limit is None or self.returned < self._limit)
        ):
            return None
        return self._cursor.advance(self._budget.position, self._budget.offset).token()

    def close(self) -> None:
        """Release the iterator's open reader after early termination."""
        self._iterator.close()
        self._closed = True

    def retain_identities(self, count: int = 1) -> None:
        """Charge consumer lookup state against this iterator's shared state cap."""
        integer(count, field_name="retained identities", minimum=1)
        self._budget.retain(count)

    def __enter__(self) -> Records[T]:
        """Return this iterator without starting a scan."""
        return self

    def __exit__(
        self,
        _type: type[BaseException] | None,
        _value: BaseException | None,
        _traceback: TracebackType | None,
    ) -> None:
        """Release the owned reader even when the caller breaks or raises."""
        self.close()


@dataclass(frozen=True, repr=False)
class RecordView:
    """A bounded descriptor; each records() call opens an independent snapshot."""

    path: Path
    revision: int
    view: str
    limit: int | None = 100
    pool_id: str | None = None
    select: PartFilter | CandidateFilter | AttemptFilter | DesignFilter | None = None
    run_id: str | None = None
    source_records: int | None = None
    read_limits: ReadLimits = field(default_factory=ReadLimits)
    after: Cursor | None = None

    def __post_init__(self) -> None:
        """Validate the bounded view without opening a database."""
        integer(self.revision, field_name="revision", minimum=0)
        if self.limit is not None:
            integer(self.limit, field_name="limit", minimum=1)
        if not isinstance(self.read_limits, ReadLimits):
            msg = "record views require ReadLimits"
            raise TypeError(msg)

        if self.source_records is not None:
            integer(self.source_records, field_name="source_records", minimum=0)
        self._validate_source()
        self._validate_filter()
        if self.after is not None and (
            not isinstance(self.after, Cursor) or self.after.advance(0) != self.cursor
        ):
            msg = "cursor does not match this query, source identity or revision"
            raise ValueError(msg)
        if self.after is not None and self.after.offset and self.view != "placements":
            msg = "cursor offset is only valid for placement rows"
            raise ValueError(msg)

    def _validate_source(self) -> None:
        if self.view not in {
            "designs",
            "sequences",
            "placements",
            "attempts",
            "parts",
            "candidates",
            "batches",
        }:
            msg = f"unsupported record view {self.view!r}"
            raise ValueError(msg)
        if (self.view in {"parts", "candidates"}) != (self.pool_id is not None):
            msg = "part views require a pool snapshot identity"
            raise ValueError(msg)
        if self.view in {"parts", "candidates"}:
            digest(self.pool_id, field_name="pool_id")
            if self.run_id is not None or self.revision != 0:
                msg = "immutable pool views require revision zero and no run ID"
                raise ValueError(msg)
        else:
            required_text(self.run_id, field_name="run_id")

    def _validate_filter(self) -> None:
        if self.select is not None:
            expected = {
                "parts": PartFilter,
                "candidates": CandidateFilter,
                "attempts": AttemptFilter,
                "designs": DesignFilter,
                "sequences": DesignFilter,
                "placements": DesignFilter,
            }.get(self.view)
            if expected is None or not isinstance(self.select, expected):
                msg = (
                    "record view requires its matching PartFilter, "
                    "DesignFilter or AttemptFilter"
                )
                raise TypeError(msg)
            if (
                isinstance(self.select, (AttemptFilter, DesignFilter, CandidateFilter))
                and self.select.identities > self.read_limits.identities
            ):
                msg = "read_limits.identities cannot hold the attempt filter"
                raise ReadLimitError(msg)

    @property
    def record_schema(self) -> str:
        """Declare the row wire schema even when a selection is empty."""
        return {
            "parts": "dense_arrays.pool_part.v1",
            "candidates": CANDIDATE_RECORD_SCHEMA,
            "designs": DESIGN_RECORD_SCHEMA,
            "sequences": "dense_arrays.sequence_record.v1",
            "placements": "dense_arrays.placement_record.v1",
            "attempts": "dense_arrays.attempt.v1",
            "batches": DECISION_SCHEMA,
        }[self.view]

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Describe the exact source and query used by portable exports."""
        return (
            {
                "pool_id" if self.pool_id else "run_id": self.pool_id or self.run_id,
                "revision": self.revision,
                "filter": None if self.select is None else self.select.to_dict(),
                "scope": "all_matching_records",
                "record_schema": self.record_schema,
            },
        )

    @property
    def cursor(self) -> Cursor:
        """Bind source, view, predicate, record version and retained ordering."""
        return Cursor(
            self.pool_id or self.run_id,
            self.revision,
            semantic_digest(
                {
                    "schema": "dense_arrays.record_query.v1",
                    "view": self.view,
                    "record_schema": self.record_schema,
                    "selection": None if self.select is None else self.select.to_dict(),
                    "order": "native_ordinal.v1",
                }
            ),
        )

    @property
    def start_ordinal(self) -> int:
        """Use a validated token's position, or begin at the first record."""
        return 0 if self.after is None else self.after.ordinal

    @property
    def cost(self) -> ReadCost:
        """Disclose the worst-case data scan before opening an iterator."""
        count = self.source_records
        if count is not None:
            count = max(
                0,
                count
                - self.start_ordinal
                + (1 if self.after and self.after.offset else 0),
            )
        if (
            self.select is None
            and count is not None
            and self.limit is not None
            and self.view != "placements"
        ):
            count = min(count, self.limit)
        if count is not None and (
            self.view == "placements"
            or (
                isinstance(self.select, DesignFilter)
                and (self.select.part_ids or self.select.groups)
            )
        ):
            count += 1
        return ReadCost(
            source_id=self.pool_id or self.run_id,
            revision=self.revision,
            mode="indexed" if self.select is None else "scan",
            projection=self.view,
            records_estimate=count,
            limits=self.read_limits,
        )

    def records(
        self,
    ) -> Records[
        Design
        | Attempt
        | PoolPart
        | PoolCandidate
        | SequenceRecord
        | PlacementRecord
        | BatchDecision
    ]:
        """Create a fresh iterator over the same revision and page."""
        budget = ReadBudget(self.read_limits)
        return Records(read_records(self, budget), budget, self.cursor, self.limit)

    def __repr__(self) -> str:
        """Describe the snapshot without expanding predicates or reading records."""
        cost = self.cost
        return (
            f"RecordView({self.view}, source={cost.source_id[:12]}, "
            f"revision={self.revision}, limit={self.limit}, mode={cost.mode})"
        )


def read_records(
    view: RecordView, budget: ReadBudget
) -> Iterator[
    Design
    | Attempt
    | PoolPart
    | PoolCandidate
    | SequenceRecord
    | PlacementRecord
    | BatchDecision
]:
    """Read records with a shared work budget across report projections."""
    if view.view == "candidates":
        from dense_arrays.reporting.pools.reading import (  # noqa: PLC0415 - avoid reader initialization cycle
            iter_candidates,
        )

        yield from iter_candidates(view, budget)
        return
    if view.view == "parts":
        from dense_arrays.artifacts.pools import iter_parts  # noqa: PLC0415

        yield from iter_parts(
            view.path,
            pool_id=view.pool_id,
            selected=view.select,
            limit=view.limit,
            budget=budget,
            after=view.start_ordinal,
        )
        return
    with reader(view.path) as connection:
        snapshot = _check_snapshot(connection, view)
        design_view = view.view in {"designs", "sequences", "placements"}
        annotations = (
            design_context(connection, view, budget, snapshot.cell_ids)
            if design_view
            else {}
        )
        start = view.start_ordinal - (1 if view.after and view.after.offset else 0)
        if view.view == "batches" and snapshot.batch_count is None:
            return
        query = _record_query(view.view)
        returned = 0
        for row in connection.execute(
            query,
            (
                view.revision,
                start,
                -1
                if view.limit is None
                or view.select is not None
                or view.view == "placements"
                else view.limit,
            ),
        ):
            budget.examine(row[1])
            value = checked_payload(row[1:])
            record = _decode_record(value, view.view, budget)
            parts, collection_id = (
                annotations.get(record.plan_id, ({}, "")) if design_view else ({}, "")
            )
            if view.select is not None:
                matches = (
                    view.select.matches(record, parts, collection_id)
                    if design_view
                    else view.select.matches(record, view.run_id)
                )
                if not matches:
                    continue
            projected = (
                project(record, view.view, parts, collection_id)
                if design_view
                else (record,)
            )
            for offset, item in enumerate(projected, 1):
                if (
                    view.after
                    and row[0] == view.after.ordinal
                    and offset <= view.after.offset
                ):
                    continue
                budget.position = row[0]
                budget.offset = offset if view.view == "placements" else 0
                yield item
                returned += 1
                if view.limit is not None and returned >= view.limit:
                    return


def _check_snapshot(connection: sqlite3.Connection, view: RecordView) -> RunSummary:
    """Validate the declared run snapshot before opening record iteration."""
    exists = connection.execute(
        "SELECT payload,digest FROM commits WHERE revision=?", (view.revision,)
    ).fetchone()
    if exists is None:
        msg = "requested snapshot revision is unavailable"
        raise ValueError(msg)
    snapshot = RunSummary.from_manifest(checked_payload(exists))
    if snapshot.run_id != view.run_id or snapshot.revision != view.revision:
        msg = "run identity changed since selecting the snapshot"
        raise ValueError(msg)
    if isinstance(view.select, AttemptFilter):
        view.select.validate(snapshot)
    return snapshot


def _decode_record(
    value: dict, view: str, budget: ReadBudget
) -> Design | Attempt | BatchDecision:
    """Enforce per-record state before constructing runtime membership."""
    if view == "batches":
        size = 2 + encoded_batch_size(value.get("batch"))
        budget.retain(size)
        budget.identities -= size
        return BatchDecision.from_dict(value)
    if view == "attempts":
        evidence = value.get("evidence")
        size = encoded_objective_size(
            evidence.get("packing_objective") if isinstance(evidence, dict) else None
        )
        budget.retain(size)
        budget.identities -= size
        return Attempt.from_dict(value)
    return Design.from_dict(value)


def _record_query(view: str) -> str:
    """Choose an indexed projection without mixing record parsing into iteration."""
    if view in {"batches", "designs", "sequences", "placements"}:
        return (
            "SELECT ordinal,payload,digest FROM batches WHERE revision<=? "
            "AND ordinal>? ORDER BY ordinal LIMIT ?"
            if view == "batches"
            else "SELECT ordinal,payload,digest FROM designs WHERE revision<=? "
            "AND ordinal>? ORDER BY ordinal LIMIT ?"
        )
    return """SELECT a.attempt,a.payload,a.digest FROM attempts a
        WHERE a.revision =
        (SELECT MAX(b.revision) FROM attempts b
         WHERE b.attempt=a.attempt AND b.revision<=?)
        AND a.attempt>? ORDER BY a.attempt LIMIT ?"""
