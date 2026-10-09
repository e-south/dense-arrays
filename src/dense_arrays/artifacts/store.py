"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/store.py

Coherent local run commits with append-only revisions and exclusive ownership.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import json
import sqlite3
from contextlib import closing, contextmanager
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING
from uuid import uuid4

from dense_arrays._record_validation import canonical_json, semantic_digest
from dense_arrays.artifacts.errors import integrity_boundary
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.records import (
    OUTCOMES,
    RUN_SCHEMA,
    Attempt,
    Design,
    RunHandle,
)
from dense_arrays.artifacts.run_plans import (
    RunPlan,
    cell_plans,
    decode_plan,
    excluded_sequences,
    run_limits,
)
from dense_arrays.artifacts.run_state import MATRIX_RUN_SCHEMA, RESAMPLING_RUN_SCHEMA
from dense_arrays.planning import Limits, MatrixPlan
from dense_arrays.solver import SolverIdentity

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from dense_arrays.artifacts.batches import BatchDecision

DATABASE = "run.sqlite3"
_DDL = """
CREATE TABLE plan (id INTEGER PRIMARY KEY CHECK(id=1),
                   payload TEXT NOT NULL, digest TEXT NOT NULL);
CREATE TABLE commits (revision INTEGER PRIMARY KEY,
                      payload TEXT NOT NULL, digest TEXT NOT NULL);
CREATE TABLE attempts (attempt INTEGER NOT NULL, revision INTEGER NOT NULL,
                      payload TEXT NOT NULL,
                      digest TEXT NOT NULL, PRIMARY KEY(attempt,revision));
CREATE TABLE batches (ordinal INTEGER PRIMARY KEY, revision INTEGER NOT NULL,
                      cell_id TEXT NOT NULL, batch_index INTEGER NOT NULL,
                      payload TEXT NOT NULL, digest TEXT NOT NULL,
                      UNIQUE(cell_id, batch_index));
CREATE TABLE designs (ordinal INTEGER PRIMARY KEY, revision INTEGER NOT NULL,
                      sequence_id TEXT NOT NULL,
                      payload TEXT NOT NULL, digest TEXT NOT NULL,
                      cell_id TEXT NOT NULL DEFAULT 'default',
                      UNIQUE(cell_id, sequence_id));
"""


def checked_payload(row: tuple[str, str] | None) -> dict[str, object]:
    """Decode stored JSON only after checking its content fingerprint."""
    if row is None:
        msg = "native artifact is incomplete: no committed record"
        raise ValueError(msg)
    value = json.loads(row[0])
    if (
        not isinstance(value, dict)
        or canonical_json(value) != row[0]
        or semantic_digest(value) != row[1]
    ):
        msg = "native artifact checksum or canonical encoding mismatch"
        raise ValueError(msg)
    return value


@contextmanager
def reader(path: Path, *, filename: str = DATABASE) -> Iterator[sqlite3.Connection]:
    """Open a database strictly read-only and close it on every exit path."""
    database = path.absolute() / filename
    with closing(
        sqlite3.connect(database.as_uri() + "?mode=ro", uri=True)
    ) as connection:
        connection.execute("PRAGMA query_only=ON")
        with integrity_boundary(path):
            yield connection


def latest(connection: sqlite3.Connection) -> dict[str, object]:
    """Read the committed summary in constant record count using the revision index."""
    return checked_payload(
        connection.execute(
            "SELECT payload,digest FROM commits ORDER BY revision DESC LIMIT 1"
        ).fetchone()
    )


def stored_plan(
    connection: sqlite3.Connection, *, max_identities: int | None = None
) -> RunPlan:
    """Load the immutable bound plan without reopening original inputs."""
    value = checked_payload(
        connection.execute("SELECT payload,digest FROM plan WHERE id=1").fetchone()
    )
    return decode_plan(value, max_identities)


@dataclass
class RunWriter:
    """One owned connection; transaction completion precedes in-memory advancement."""

    connection: sqlite3.Connection
    handle: RunHandle
    state: dict[str, object]
    parent_exclusions: dict[tuple[str, str], str]
    limits: Limits

    def refresh(self) -> None:
        """Reload committed state after an exception at a transaction boundary."""
        self.state = latest(self.connection)

    def begin_resume(self, *, active_seconds: float) -> None:
        """Mark active execution before admitting work or rebuilding a model."""
        state = json.loads(canonical_json(self.state))
        state.update(
            state="running",
            termination_reason=None,
            resumable=False,
            active_seconds=active_seconds,
        )
        for cell in state.get("cells", {}).values():
            if (
                cell["state"] == "stopped"
                and cell["termination_reason"] == "interrupted"
            ):
                cell.update(state="running", termination_reason=None)
        self._commit(state)

    def _commit(
        self,
        next_state: dict[str, object],
        *,
        attempt: dict[str, object] | None = None,
        design: Design | None = None,
        batch_decision: BatchDecision | None = None,
    ) -> None:
        revision = self.state["revision"] + 1
        next_state["revision"] = revision
        with self.connection:
            if batch_decision is not None:
                payload = batch_decision.to_dict()
                self.connection.execute(
                    "INSERT INTO batches VALUES (?,?,?,?,?,?)",
                    (
                        next_state["batch_count"],
                        revision,
                        batch_decision.cell_id,
                        batch_decision.index,
                        canonical_json(payload),
                        semantic_digest(payload),
                    ),
                )
            if attempt is not None:
                self.connection.execute(
                    "INSERT INTO attempts VALUES (?,?,?,?)",
                    (
                        attempt["attempt_id"],
                        revision,
                        canonical_json(attempt),
                        semantic_digest(attempt),
                    ),
                )
            if design is not None:
                payload = design.to_dict()
                matrix = bool(self.state.get("cells"))
                self.connection.execute(
                    (
                        "INSERT INTO designs "
                        "(ordinal,revision,sequence_id,payload,digest,cell_id) "
                        "VALUES (?,?,?,?,?,?)"
                    )
                    if matrix
                    else (
                        "INSERT INTO designs "
                        "(ordinal,revision,sequence_id,payload,digest) "
                        "VALUES (?,?,?,?,?)"
                    ),
                    (
                        next_state["counts"]["accepted"],
                        revision,
                        design.sequence_id,
                        canonical_json(payload),
                        semantic_digest(payload),
                    )
                    + ((design.cell_id,) if matrix else ()),
                )
            self.connection.execute(
                "INSERT INTO commits VALUES (?,?,?)",
                (revision, canonical_json(next_state), semantic_digest(next_state)),
            )
        self.state = next_state

    def bind_solver(self, identity: SolverIdentity) -> None:
        """Record the actual model's backend before reserving any search work."""
        if not isinstance(identity, SolverIdentity):
            msg = "solver binding requires SolverIdentity"
            raise TypeError(msg)
        producer = Producer.from_dict(self.state["producer"])
        if producer.solver == identity:
            return
        if producer.solver is not None or self.state["counts"]["started"]:
            msg = "cannot change the solver identity after binding or starting search"
            raise ValueError(msg)
        self._commit(
            dict(self.state, producer=replace(producer, solver=identity).to_dict())
        )

    def reserve(  # noqa: PLR0913 - reservation and optional first batch commit
        self,
        *,
        active_seconds: float,
        cell_id: str = "default",
        batch_id: str | None = None,
        batch_index: int | None = None,
        batch_attempt: int | None = None,
        batch_decision: BatchDecision | None = None,
    ) -> int:
        """Durably reserve an attempt before invoking any candidate search."""
        state = json.loads(canonical_json(self.state))
        if self.state["counts"]["in_progress"]:
            msg = "cannot reserve while another attempt is in progress"
            raise ValueError(msg)
        if "cells" in state:
            cell = state["cells"][cell_id]
            if cell["counts"]["accepted"] >= cell["target"] or cell["state"] not in {
                "created",
                "running",
            }:
                msg = "cannot reserve work for a completed or stopped cell"
                raise ValueError(msg)
            cell["counts"]["started"] += 1
            cell["counts"]["in_progress"] += 1
            cell["state"] = "running"
        elif cell_id != "default":
            msg = "unknown cell identity"
            raise ValueError(msg)
        state["counts"]["started"] += 1
        state["counts"]["in_progress"] += 1
        state.update(state="running", active_seconds=active_seconds)
        ordinal = state["counts"]["started"]
        context = {
            **({"cell_attempt": cell["counts"]["started"]} if "cells" in state else {}),
            **({"batch_id": batch_id} if batch_id is not None else {}),
            **(
                {"batch_index": batch_index, "batch_attempt": batch_attempt}
                if batch_index is not None or batch_attempt is not None
                else {}
            ),
        }
        record = Attempt(
            attempt_id=ordinal, cell_id=cell_id, outcome="in_progress", evidence=context
        )
        if batch_decision is not None:
            progress = self.state.get("cells", {}).get(cell_id, self.state)
            if (
                "batch_count" not in state
                or batch_attempt != 1
                or batch_decision.run_id != self.handle.run_id
                or batch_decision.cell_id != cell_id
                or batch_decision.plan_id != progress["plan_id"]
                or batch_decision.index != batch_index
                or batch_decision.batch.batch_id != batch_id
                or batch_decision.after_attempt != ordinal - 1
            ):
                msg = "batch decision does not match its first reservation"
                raise ValueError(msg)
            state["batch_count"] += 1
        self._commit(state, attempt=record.to_dict(), batch_decision=batch_decision)
        return ordinal

    def _reserved_evidence(
        self, attempt: int, evidence: dict[str, object]
    ) -> dict[str, object]:
        """Preserve the reserved batch through interruptions and outcomes."""
        reserved = checked_payload(
            self.connection.execute(
                "SELECT payload,digest FROM attempts WHERE attempt=? "
                "ORDER BY revision DESC LIMIT 1",
                (attempt,),
            ).fetchone()
        )
        for key in ("batch_id", "batch_index", "batch_attempt"):
            if key in reserved["evidence"]:
                value = reserved["evidence"][key]
                if evidence.get(key, value) != value:
                    msg = "attempt cannot change its reserved batch identity"
                    raise ValueError(msg)
                evidence = {**evidence, key: value}
        return evidence

    def publish(  # noqa: PLR0913 - one attempt publication and terminal cell state
        self,
        attempt: int,
        outcome: str,
        evidence: dict[str, object],
        *,
        active_seconds: float,
        design: Design | None = None,
        terminal: tuple[str, str] | None = None,
    ) -> str:
        """Commit an outcome, design, sequence uniqueness and counters together."""
        if outcome not in OUTCOMES or outcome == "in_progress":
            msg = "publish requires a terminal attempt outcome"
            raise ValueError(msg)
        if (
            attempt != self.state["counts"]["started"]
            or self.state["counts"]["in_progress"] != 1
        ):
            msg = "publish requires the currently reserved attempt"
            raise ValueError(msg)
        if (design is not None) != (outcome == "accepted"):
            msg = "an accepted outcome requires exactly one design"
            raise ValueError(msg)
        cell_id = next(
            (
                name
                for name, cell in self.state.get("cells", {}).items()
                if cell["counts"]["in_progress"]
            ),
            "default",
        )
        evidence = self._reserved_evidence(attempt, evidence)
        if design is not None and design.cell_id != cell_id:
            msg = "design does not belong to the reserved cell"
            raise ValueError(msg)
        if (
            design is not None
            and (cell_id, design.sequence_id) in self.parent_exclusions
        ):
            evidence = {
                **evidence,
                "code": "parent_duplicate",
                "candidate_sequence_id": design.sequence_id,
                "matched_design_ref": self.parent_exclusions[
                    (cell_id, design.sequence_id)
                ],
            }
            outcome, design = "duplicate", None
        if (
            design is not None
            and self.connection.execute(
                "SELECT 1 FROM designs WHERE sequence_id=? AND cell_id=?"
                if "cells" in self.state
                else "SELECT 1 FROM designs WHERE sequence_id=?",
                (design.sequence_id,) + ((cell_id,) if "cells" in self.state else ()),
            ).fetchone()
        ):
            outcome, design = "duplicate", None
        state = json.loads(canonical_json(self.state))
        state["counts"]["in_progress"] -= 1
        state["counts"][outcome] += 1
        state["active_seconds"] = active_seconds
        if "cells" in state:
            cell = state["cells"][cell_id]
            evidence = {**evidence, "cell_attempt": cell["counts"]["started"]}
            cell["counts"]["in_progress"] -= 1
            cell["counts"][outcome] += 1
            _finish_cell_outcome(cell, terminal)
        record = {
            "schema": "dense_arrays.attempt.v1",
            "attempt_id": attempt,
            "cell_id": cell_id,
            "outcome": outcome,
            "evidence": evidence,
        }
        if design is not None:
            record["design_ref"] = design.reference
        Attempt.from_dict(record)
        self._commit(state, attempt=record, design=design)
        return outcome

    def stop_cell(
        self, cell_id: str, state: str, reason: str, *, active_seconds: float
    ) -> None:
        """Record a schedule stop after its committed prefix reaches a boundary."""
        updated = json.loads(canonical_json(self.state))
        cell = updated["cells"][cell_id]
        if cell["state"] not in {"created", "running"} or cell["counts"]["in_progress"]:
            msg = "only an open cell without an active attempt can stop"
            raise ValueError(msg)
        _finish_cell_outcome(cell, (state, reason))
        updated["active_seconds"] = active_seconds
        self._commit(updated)

    def finish(self, state: str, reason: str, *, active_seconds: float) -> None:
        """Publish terminal run status without changing its original target."""
        if state not in {"completed", "stopped", "failed"}:
            msg = "finish requires a terminal run state"
            raise ValueError(msg)
        attained = self.state["counts"]["accepted"] == self.state["target"]
        if (state == "completed") != attained:
            msg = "completed means exactly the original target was attained"
            raise ValueError(msg)
        updated = dict(
            self.state,
            state=state,
            termination_reason=reason,
            active_seconds=active_seconds,
            resumable=(
                state == "stopped"
                and reason == "interrupted"
                and self.state["counts"]["started"] < self.limits.attempts
                and active_seconds < self.limits.active_seconds
                and (
                    "cells" not in self.state
                    or any(
                        c["state"] in {"created", "running"}
                        for c in self.state["cells"].values()
                    )
                )
            ),
        )
        if "cells" in updated:
            updated = json.loads(canonical_json(updated))
            for cell in updated["cells"].values():
                if cell["state"] in {"created", "running"}:
                    cell.update(
                        state="failed" if state == "failed" else "stopped",
                        termination_reason=reason,
                    )
        self._commit(updated)


def _finish_cell_outcome(
    cell: dict[str, object], terminal: tuple[str, str] | None
) -> None:
    """Close a cell in the same transaction as its final attempt evidence."""
    if terminal is not None:
        state, reason = terminal
        if (
            state not in {"stopped", "failed"}
            or not isinstance(reason, str)
            or not reason
        ):
            msg = "cell termination requires a stopped/failed state and reason"
            raise ValueError(msg)
        cell.update(state=state, termination_reason=reason)
    elif cell["counts"]["accepted"] == cell["target"]:
        cell.update(state="completed", termination_reason="target_attained")


@contextmanager
def create_run(plan: RunPlan, out: Path) -> Iterator[RunWriter]:
    """Own a new destination and stable lock inode for the entire writer lifetime."""
    import fcntl  # noqa: PLC0415 - Unix capability checked before destination creation

    producer = Producer.capture()
    cells = cell_plans(plan)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=False)
    with (out / ".writer.lock").open("xb") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with closing(sqlite3.connect(out / DATABASE)) as connection:
            connection.execute("PRAGMA synchronous=FULL")
            connection.executescript(_DDL)
            run_id = str(uuid4())
            state = {
                "schema": MATRIX_RUN_SCHEMA
                if isinstance(plan, MatrixPlan)
                else RUN_SCHEMA,
                "run_id": run_id,
                "plan_id": plan.plan_id,
                "revision": 0,
                "state": "created",
                "target": sum(p.request.target.count for p in cells.values()),
                "counts": {"started": 0, **dict.fromkeys(OUTCOMES, 0)},
                "termination_reason": None,
                "active_seconds": 0.0,
                "resumable": False,
                "producer": producer.to_dict(),
            }
            if any(p.request.resampling is not None for p in cells.values()):
                state["schema"] = RESAMPLING_RUN_SCHEMA
                state["batch_count"] = 0
            if isinstance(plan, MatrixPlan):
                state["cells"] = {
                    name: {
                        "cell_id": name,
                        "plan_id": p.plan_id,
                        "target": p.request.target.count,
                        "counts": {"started": 0, **dict.fromkeys(OUTCOMES, 0)},
                        "state": "created" if p.request.target.count else "inactive",
                        "termination_reason": None
                        if p.request.target.count
                        else "zero_target",
                    }
                    for name, p in cells.items()
                }
            encoded_plan = plan.to_dict()
            with connection:
                connection.execute(
                    "INSERT INTO plan VALUES (1,?,?)",
                    (
                        canonical_json(encoded_plan),
                        semantic_digest(encoded_plan),
                    ),
                )
                connection.execute(
                    "INSERT INTO commits VALUES (0,?,?)",
                    (
                        canonical_json(state),
                        semantic_digest(state),
                    ),
                )
            yield RunWriter(
                connection,
                RunHandle(out, run_id),
                state,
                excluded_sequences(plan),
                run_limits(plan),
            )
