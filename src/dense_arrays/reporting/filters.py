"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/filters.py

Typed predicates over persisted attempt identity and outcome evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    integer,
    object_fields,
    records,
    required_text,
)
from dense_arrays.artifacts.records import OUTCOMES

if TYPE_CHECKING:
    from dense_arrays.artifacts.records import Attempt
    from dense_arrays.reporting.summary import RunSummary


@dataclass(frozen=True)
class AttemptFilter:
    """OR within each field and AND across attempt IDs, cell references and outcomes."""

    attempt_ids: tuple[int, ...] = ()
    cells: tuple[str, ...] = ()
    outcomes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Freeze selectors without coercing boolean IDs or unknown outcome codes."""
        for name in ("attempt_ids", "cells", "outcomes"):
            values = records(
                getattr(self, name),
                int if name == "attempt_ids" else str,
                field_name=name,
            )
            for value in values:
                if name == "attempt_ids":
                    integer(value, field_name=name, minimum=1)
                else:
                    required_text(value, field_name=name)
            if len(set(values)) != len(values):
                msg = f"{name} must not repeat labels"
                raise ValueError(msg)
            object.__setattr__(self, name, values)
        if unknown := set(self.outcomes) - set(OUTCOMES):
            msg = f"unknown attempt outcomes: {sorted(unknown)}"
            raise ValueError(msg)

    @property
    def identities(self) -> int:
        """Number of explicit identity selectors held for this query."""
        return len(self.attempt_ids) + len(self.cells)

    def validate(self, summary: RunSummary) -> None:
        """Resolve IDs against the selected snapshot's contiguous attempt ordinals."""
        if any(value > summary.counts["started"] for value in self.attempt_ids):
            msg = "unknown attempt IDs at the selected revision"
            raise ValueError(msg)
        if set(self.cells) - (
            set(summary.cell_ids) | {f"{summary.run_id}/{c}" for c in summary.cell_ids}
        ):
            msg = "unknown cell references at the selected revision"
            raise ValueError(msg)

    def matches(self, attempt: Attempt, run_id: str) -> bool:
        """Match recorded fields; no evidence is regenerated."""
        return (
            (not self.attempt_ids or attempt.attempt_id in self.attempt_ids)
            and (not self.outcomes or attempt.outcome in self.outcomes)
            and (
                not self.cells
                or attempt.cell_id in self.cells
                or f"{run_id}/{attempt.cell_id}" in self.cells
            )
        )

    def to_dict(self) -> dict[str, object]:
        """Serialize a declared predicate for reuse across Python and CLI."""
        return {
            "schema": "dense_arrays.attempt-filter.v1",
            "attempt_ids": list(self.attempt_ids),
            "cells": list(self.cells),
            "outcomes": list(self.outcomes),
        }

    @classmethod
    def from_dict(cls, value: object) -> AttemptFilter:
        """Reject unknown schemas and fields before selecting any evidence."""
        data = object_fields(
            value, {"schema", "attempt_ids", "cells", "outcomes"}, "attempt filter"
        )
        if data.pop("schema", None) != "dense_arrays.attempt-filter.v1":
            msg = "unsupported attempt-filter schema"
            raise ValueError(msg)
        return cls(**data)
