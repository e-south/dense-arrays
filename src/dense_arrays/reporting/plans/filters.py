"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/plans/filters.py

Select included plans by their full semantic identity.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass

from dense_arrays._record_validation import digest, object_fields, records

PLAN_FILTER_SCHEMA = "dense_arrays.plan-filter.v1"


@dataclass(frozen=True)
class PlanFilter:
    """Choose explicit plan IDs; an empty predicate includes the whole inventory."""

    plan_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Require complete, unique semantic digests without prefix matching."""
        values = records(self.plan_ids, str, field_name="plan_ids")
        for value in values:
            digest(value, field_name="plan_id")
        if len(set(values)) != len(values):
            msg = "plan_ids must not repeat identities"
            raise ValueError(msg)
        object.__setattr__(self, "plan_ids", values)

    @property
    def identities(self) -> int:
        """Count explicit selectors retained by this predicate."""
        return len(self.plan_ids)

    def to_dict(self) -> dict[str, object]:
        """Serialize the same predicate accepted by the CLI selection file."""
        return {"schema": PLAN_FILTER_SCHEMA, "plan_ids": list(self.plan_ids)}

    @classmethod
    def from_dict(cls, value: object) -> PlanFilter:
        """Reject unknown schemas and fields before reading plan evidence."""
        data = object_fields(value, {"schema", "plan_ids"}, "plan filter")
        if data.get("schema") != PLAN_FILTER_SCHEMA:
            msg = "unsupported plan filter schema"
            raise ValueError(msg)
        return cls(data["plan_ids"])
