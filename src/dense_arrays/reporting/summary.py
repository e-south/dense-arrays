"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/summary.py

Run summaries distinguish attainment, effort and search termination.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

from dense_arrays._record_validation import (
    digest,
    immutable_json_mapping,
    integer,
    object_fields,
)
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.reading import ReadCost, ReadLimits, Verification
from dense_arrays.artifacts.records import RUN_SCHEMA
from dense_arrays.artifacts.run_state import (
    MATRIX_RUN_SCHEMA,
    RESAMPLING_RUN_SCHEMA,
    CellSummary,
    validate_counts,
)


@dataclass(frozen=True)
class RunSummary:
    """Bounded metadata for one committed run revision, without scanning designs."""

    run_id: str
    plan_id: str
    revision: int
    state: str
    target: int
    counts: Mapping[str, int]
    termination_reason: str | None
    active_seconds: float
    resumable: bool
    verified: bool = False
    read_limits: ReadLimits = field(
        default_factory=ReadLimits, compare=False, repr=False
    )
    verification: Verification | None = field(default=None, compare=False)
    producer: Producer | None = None
    cells: Mapping[str, CellSummary] = field(default_factory=dict, repr=False)

    batch_count: int | None = None

    @property
    def cost(self) -> ReadCost:
        """Describe the bounded manifest read independently of verification."""
        return ReadCost(
            self.run_id, self.revision, "manifest", "summary", 1, self.read_limits
        )

    @property
    def verification_cost(self) -> ReadCost:
        """Preview the complete committed plan, attempt and accepted design scope."""
        return ReadCost(
            self.run_id,
            self.revision,
            "scan",
            "verification",
            1 + self.counts["started"] + self.accepted + (self.batch_count or 0),
            self.read_limits,
        )

    def __post_init__(self) -> None:
        """Validate stored counters independently of presentation."""
        if self.producer is not None and not isinstance(self.producer, Producer):
            msg = "producer must be a recorded Producer or null"
            raise TypeError(msg)
        digest(self.plan_id, field_name="plan_id")
        integer(self.revision, field_name="revision", minimum=0)
        integer(self.target, field_name="target", minimum=0)
        validate_counts(self.counts)
        if self.state not in {"created", "running", "completed", "stopped", "failed"}:
            msg = "unknown run state"
            raise ValueError(msg)
        if self.counts["accepted"] > self.target or (
            self.state == "completed" and self.counts["accepted"] != self.target
        ):
            msg = "run completion does not match its original target"
            raise ValueError(msg)
        if (
            isinstance(self.active_seconds, bool)
            or not isinstance(self.active_seconds, (int, float))
            or not math.isfinite(self.active_seconds)
            or self.active_seconds < 0
        ):
            msg = "active_seconds must be a finite nonnegative number"
            raise ValueError(msg)
        if not isinstance(self.resumable, bool):
            msg = "resumable must be a boolean"
            raise TypeError(msg)
        object.__setattr__(self, "counts", immutable_json_mapping(self.counts))
        if self.batch_count is not None:
            integer(self.batch_count, field_name="batch_count", minimum=0)
            if self.batch_count > self.counts["started"]:
                msg = "batch count exceeds reserved attempts"
                raise ValueError(msg)
        self._validate_cells()

    def _validate_cells(self) -> None:
        if not isinstance(self.cells, Mapping) or any(
            not isinstance(c, CellSummary) or name != c.cell_id
            for name, c in self.cells.items()
        ):
            msg = "run cells must map cell identities to CellSummary records"
            raise TypeError(msg)
        if self.cells and (
            sum(c.target for c in self.cells.values()) != self.target
            or any(
                sum(c.counts[name] for c in self.cells.values()) != count
                for name, count in self.counts.items()
            )
        ):
            msg = "cell targets or counters do not reconcile with the run"
            raise ValueError(msg)
        object.__setattr__(self, "cells", MappingProxyType(dict(self.cells)))

    @property
    def cell_ids(self) -> tuple[str, ...]:
        """Return declared cell identities, including inactive cells."""
        return tuple(self.cells) if self.cells else ("default",)

    @property
    def cell_summaries(self) -> Mapping[str, CellSummary]:
        """Expose cell attainment for both original and multi-cell run schemas."""
        return self.cells or MappingProxyType(
            {
                "default": CellSummary(
                    "default",
                    self.plan_id,
                    self.target,
                    self.counts,
                    self.state,
                    self.termination_reason,
                )
            }
        )

    @property
    def accepted(self) -> int:
        """Accepted design count, separate from the original target."""
        return self.counts["accepted"]

    def to_dict(self) -> dict[str, object]:
        """Return a versioned report, suitable for pure machine stdout."""
        return {
            "schema": "dense_arrays.run_summary.v4"
            if self.batch_count is not None
            else "dense_arrays.run_summary.v3"
            if self.cells
            else "dense_arrays.run_summary.v1"
            if self.producer is None
            else "dense_arrays.run_summary.v2",
            **(
                {"batch_count": self.batch_count}
                if self.batch_count is not None
                else {}
            ),
            "run_id": self.run_id,
            "plan_id": self.plan_id,
            "revision": self.revision,
            "state": self.state,
            "target": self.target,
            "accepted": self.accepted,
            "counts": dict(self.counts),
            "termination_reason": self.termination_reason,
            "active_seconds": self.active_seconds,
            "resumable": self.resumable,
            "verified": self.verified,
            **(
                {"cells": {name: c.to_dict() for name, c in self.cells.items()}}
                if self.cells
                else {}
            ),
            **(
                {"producer": self.producer.to_dict()}
                if self.producer is not None
                else {}
            ),
            **(
                {"verification": self.verification.to_dict()}
                if self.verification is not None
                else {}
            ),
        }

    @classmethod
    def from_manifest(cls, value: object) -> RunSummary:
        """Read a stored manifest without accepting future fields or versions."""
        expected = {
            "schema",
            "run_id",
            "plan_id",
            "revision",
            "state",
            "target",
            "counts",
            "termination_reason",
            "active_seconds",
            "resumable",
        }
        data = object_fields(
            value, expected | {"producer", "cells", "batch_count"}, "run"
        )
        schema = data.get("schema")
        if schema not in {
            "dense_arrays.run.v1",
            RUN_SCHEMA,
            MATRIX_RUN_SCHEMA,
            RESAMPLING_RUN_SCHEMA,
        }:
            msg = (
                f"unsupported run schema {schema!r}; "
                f"supported: dense_arrays.run.v1, {RUN_SCHEMA}"
            )
            raise ValueError(msg)
        if schema in {RUN_SCHEMA, MATRIX_RUN_SCHEMA, RESAMPLING_RUN_SCHEMA}:
            expected.add("producer")
        if schema == RESAMPLING_RUN_SCHEMA:
            expected.add("batch_count")
            if "cells" in data:
                expected.add("cells")
        if schema == MATRIX_RUN_SCHEMA:
            expected.add("cells")
        if set(data) != expected:
            msg = "unsupported or incomplete run schema"
            raise ValueError(msg)
        data.pop("schema")
        if "producer" in data:
            data["producer"] = Producer.from_dict(data["producer"])
        if "cells" in data:
            if not isinstance(data["cells"], Mapping) or not data["cells"]:
                msg = "matrix run requires nonempty cell summaries"
                raise ValueError(msg)
            data["cells"] = {
                name: CellSummary.from_dict(c) for name, c in data["cells"].items()
            }
        return cls(**data)
