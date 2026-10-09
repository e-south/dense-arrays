"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/run_state.py

Validated native cell attainment and reconciled attempt counters.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    immutable_json_mapping,
    integer,
    object_fields,
    required_text,
)
from dense_arrays.artifacts.records import OUTCOMES

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping

RESAMPLING_RUN_SCHEMA = "dense_arrays.run.v4"
MATRIX_RUN_SCHEMA = "dense_arrays.run.v3"


def manifest_cell_ids(manifest: Mapping[str, object]) -> tuple[str, ...]:
    """Expose cell labels from an already validated native origin manifest."""
    return tuple(manifest["cells"]) if "cells" in manifest else ("default",)


def origin_bindings(
    manifests: Iterable[Mapping[str, object]],
) -> Iterator[tuple[str, str, str]]:
    """Join declared run, cell and plan identities without inspecting designs."""
    for manifest in manifests:
        for cell in manifest_cell_ids(manifest):
            plan_id = (
                manifest["cells"][cell]["plan_id"]
                if "cells" in manifest
                else manifest["plan_id"]
            )
            yield manifest["run_id"], cell, plan_id


def validate_counts(counts: Mapping[str, int]) -> None:
    """Require nonnegative exclusive outcomes to reconcile with started work."""
    if set(counts) != {"started", *OUTCOMES}:
        msg = "run accounting contains missing or unknown outcome categories"
        raise ValueError(msg)
    for name, value in counts.items():
        integer(value, field_name=f"counts.{name}", minimum=0)
    if counts["started"] != sum(counts[name] for name in OUTCOMES):
        msg = "run accounting does not reconcile"
        raise ValueError(msg)


@dataclass(frozen=True)
class CellSummary:
    """One cell's immutable target, outcomes and termination at a run revision."""

    cell_id: str
    plan_id: str
    target: int
    counts: Mapping[str, int]
    state: str
    termination_reason: str | None

    def __post_init__(self) -> None:
        """Validate attainment independently from the aggregate run counters."""
        required_text(self.cell_id, field_name="cell_id")
        digest(self.plan_id, field_name="cell.plan_id")
        integer(self.target, field_name="cell.target", minimum=0)
        validate_counts(self.counts)
        if self.state not in {
            "created",
            "running",
            "completed",
            "stopped",
            "failed",
            "inactive",
        }:
            msg = "unknown cell state"
            raise ValueError(msg)
        if self.accepted > self.target or (
            self.state == "completed" and self.accepted != self.target
        ):
            msg = "cell completion does not match its original target"
            raise ValueError(msg)
        if self.state == "inactive" and (self.target or self.counts["started"]):
            msg = "inactive cells cannot request or start work"
            raise ValueError(msg)
        object.__setattr__(self, "counts", immutable_json_mapping(self.counts))

    @property
    def accepted(self) -> int:
        """Return accepted designs in this cell only."""
        return self.counts["accepted"]

    def to_dict(self) -> dict[str, object]:
        """Encode cell state inside the versioned run manifest."""
        return {
            "cell_id": self.cell_id,
            "plan_id": self.plan_id,
            "target": self.target,
            "counts": dict(self.counts),
            "state": self.state,
            "termination_reason": self.termination_reason,
        }

    @classmethod
    def from_dict(cls, value: object) -> CellSummary:
        """Reject missing or unknown native cell fields."""
        fields = {
            "cell_id",
            "plan_id",
            "target",
            "counts",
            "state",
            "termination_reason",
        }
        data = object_fields(value, fields, "cell summary")
        if set(data) != fields:
            msg = "incomplete cell summary"
            raise ValueError(msg)
        return cls(**data)
