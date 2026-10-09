"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/storage.py

Exclusive sampled pool publication and candidate-evidence reads.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import sqlite3
from contextlib import closing, contextmanager
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.pool_records import PoolPart, PoolSummary
from dense_arrays.artifacts.preparation.sets import SetAccounting
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.store import checked_payload
from dense_arrays.parts.candidates import Candidate
from dense_arrays.parts.pools import PoolHandle

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from dense_arrays.artifacts.reading import ReadBudget
    from dense_arrays.planning.preparation import PreparationPlan

    from .records import PoolAccounting


@contextmanager
def pool_destination(out: Path) -> Iterator[sqlite3.Connection]:
    """Reserve an output directory before mining and commit one coherent pool."""
    import fcntl  # noqa: PLC0415 - local lock capability before output ownership

    if out.exists() or out.is_symlink():
        msg = f"output destination already exists: {out}"
        raise FileExistsError(msg)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=False)
    with (out / ".writer.lock").open("xb") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with closing(sqlite3.connect(out / "pool.sqlite3")) as connection:
            connection.execute("PRAGMA synchronous=FULL")
            connection.executescript("""
                CREATE TABLE manifest (id INTEGER PRIMARY KEY CHECK(id=1),
                    payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE TABLE preparation (id INTEGER PRIMARY KEY CHECK(id=1),
                    payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE TABLE candidates (ordinal INTEGER PRIMARY KEY,
                    payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE TABLE parts (ordinal INTEGER PRIMARY KEY,
                    part_id TEXT UNIQUE NOT NULL, group_name TEXT,
                    payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE INDEX parts_group ON parts(group_name);
            """)
            with connection:
                try:
                    yield connection
                except KeyboardInterrupt as err:
                    err.add_note(
                        f"Preparation interrupted; no pool was committed at {out}. "
                        "Restart preparation with a new destination."
                    )
                    raise


def result_identity(
    plan_id: str,
    candidates: tuple[Candidate, ...],
    accounting: PoolAccounting | SetAccounting,
) -> str:
    """Bind actual decisions as well as the recipe, independently of destination."""
    return semantic_digest(
        {
            "schema": "dense_arrays.prepared_pool.v1",
            "plan_id": plan_id,
            "candidates": [semantic_digest(c.to_dict()) for c in candidates],
            "accounting": accounting.to_dict(),
        }
    )


def publish(
    connection: sqlite3.Connection,
    out: Path,
    plan: PreparationPlan,
    candidates: tuple[Candidate, ...],
    accounting: PoolAccounting | SetAccounting,
) -> PoolHandle:
    """Write plan, decisions, retained parts and manifest in the caller transaction."""
    identity = result_identity(plan.plan_id, candidates, accounting)
    summary = PoolSummary(
        identity,
        plan.plan_id,
        accounting.counts["processed"],
        accounting.counts["retained"],
        state=accounting.state,
        producer=Producer.capture(),
        preparation=accounting,
    )
    for candidate in candidates:
        payload = candidate.to_dict()
        connection.execute(
            "INSERT INTO candidates VALUES (?,?,?)",
            (candidate.index, canonical_json(payload), semantic_digest(payload)),
        )
    order = (
        {name: index for index, name in enumerate(accounting.recipes)}
        if isinstance(accounting, SetAccounting)
        else {}
    )
    retained = sorted(
        (c for c in candidates if c.retained),
        key=lambda c: (order.get(c.recipe_id, 0), c.rank),
    )
    for ordinal, candidate in enumerate(retained, 1):
        part = candidate.part
        payload = PoolPart(identity, ordinal, part).to_dict()
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
    for statement, payload in (
        ("INSERT INTO preparation VALUES (1,?,?)", plan.to_dict()),
        ("INSERT INTO manifest VALUES (1,?,?)", summary.to_dict()),
    ):
        connection.execute(
            statement, (canonical_json(payload), semantic_digest(payload))
        )
    return PoolHandle(out, identity)


def read_candidates(
    connection: sqlite3.Connection, budget: ReadBudget
) -> tuple[Candidate, ...]:
    """Read indexed candidate evidence under explicit record and state caps."""
    result = []
    for ordinal, payload, checksum in connection.execute(
        "SELECT ordinal,payload,digest FROM candidates ORDER BY ordinal"
    ):
        budget.examine(payload)
        budget.retain()
        candidate = Candidate.from_dict(checked_payload((payload, checksum)))
        if ordinal != len(result) + 1 or ordinal != candidate.index:
            msg = "sampled candidate ordinal sequence is inconsistent"
            raise ValueError(msg)
        result.append(candidate)
    return tuple(result)
