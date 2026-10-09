"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/pools/snapshots.py

Portable pool metrics with explicit recorded-evidence scope.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

from dense_arrays._record_validation import (
    digest,
    immutable_json_mapping,
    mutable_json,
    object_fields,
    semantic_digest,
)
from dense_arrays.artifacts.preparation.bands import band_report_size
from dense_arrays.artifacts.preparation.sets import (
    SET_ACCOUNTING_SCHEMA,
    read_accounting,
)
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits
from dense_arrays.reporting.pools.diversity import validate_diversity

POOL_QUALITY_SCHEMA = "dense_arrays.pool_quality.v1"


@dataclass(frozen=True, repr=False)
class PoolQualitySnapshot:
    """Recorded pool totals; checks consistency without claiming source verification."""

    data: Mapping[str, object]
    read_limits: ReadLimits = field(
        default_factory=ReadLimits, repr=False, compare=False
    )
    snapshot_id: str = field(init=False)

    def __post_init__(self) -> None:
        """Validate schema, identities, stage equations and declared completion."""
        if not isinstance(self.read_limits, ReadLimits):
            msg = "pool quality snapshots require ReadLimits"
            raise TypeError(msg)
        fields = {
            "schema",
            "pool_id",
            "plan_id",
            "state",
            "counts",
            "requested_retention",
            "candidate_budget",
            "stop_reason",
            "rejections",
        }
        data = object_fields(
            self.data,
            fields
            | {
                "retention",
                "recipes",
                "mining_target",
                "score_bands",
                "construction",
                "diversity",
            },
            "pool quality",
        )
        if not fields <= data.keys() or data["schema"] != POOL_QUALITY_SCHEMA:
            msg = "unsupported or incomplete pool quality schema"
            raise ValueError(msg)
        for key in ("counts", "rejections"):
            if not isinstance(data[key], Mapping):
                msg = f"pool quality {key} requires a mapping"
                raise TypeError(msg)
        retention = data.get("retention")
        budget = ReadBudget(self.read_limits)
        budget.retain(
            len(data)
            + len(data["counts"])
            + len(data["rejections"])
            + (len(retention) if isinstance(retention, Mapping) else 0)
            + len(data.get("mining_target", {}))
            + len(data.get("construction", {}))
            + (
                len(retention.get("sizing", {}))
                if isinstance(retention, Mapping)
                else 0
            )
            + band_report_size(data.get("score_bands"))
        )
        recipes = data.get("recipes", [])
        if not isinstance(recipes, list):
            msg = "pool quality recipes must be an ordered array"
            raise TypeError(msg)
        for raw in recipes:
            item = object_fields(raw, {"id", "accounting"}, "recipe accounting")
            record = item.get("accounting")
            if not isinstance(record, Mapping):
                msg = "recipe accounting requires a complete record"
                raise TypeError(msg)
            budget.retain(
                len(item)
                + len(record)
                + len(record.get("retention", {}).get("sizing", {}))
                + band_report_size(record.get("score_bands"))
                + sum(
                    len(record.get(key, {}))
                    for key in (
                        "counts",
                        "rejections",
                        "retention",
                        "mining_target",
                        "construction",
                    )
                )
            )
        for key in ("pool_id", "plan_id"):
            digest(data[key], field_name=key)
        accounting = read_accounting(
            {
                **{
                    k: v
                    for k, v in data.items()
                    if k not in {"pool_id", "plan_id", "state", "diversity"}
                },
                "schema": SET_ACCOUNTING_SCHEMA
                if "recipes" in data
                else "dense_arrays.pool_accounting.v1",
            }
        )
        if data["state"] != accounting.state:
            msg = "pool quality state disagrees with retained count or execution errors"
            raise ValueError(msg)
        validate_diversity(data, accounting, budget)
        object.__setattr__(self, "data", immutable_json_mapping(data))
        object.__setattr__(self, "snapshot_id", semantic_digest(data))

    @classmethod
    def from_dict(
        cls, value: object, *, read_limits: ReadLimits | None = None
    ) -> PoolQualitySnapshot:
        """Open recorded counts without reading source pools or invoking preparation."""
        return cls(value, read_limits or ReadLimits())

    def to_dict(self) -> dict[str, object]:
        """Preserve the native report schema and recorded population."""
        return mutable_json(self.data)

    @property
    def sources(self) -> tuple[dict[str, object], ...]:
        """Bind recorded origin identities independently of their locations."""
        return (
            {
                "pool_id": self.data["pool_id"],
                "plan_id": self.data["plan_id"],
                "revision": 0,
            },
        )

    @property
    def cost(self) -> ReadCost:
        """Saved aggregates need one document read, with no candidate scan."""
        return ReadCost(self.snapshot_id, 0, "manifest", "quality", 1, self.read_limits)

    def __repr__(self) -> str:
        """Display a bounded description without expanding saved reason counts."""
        return (
            f"PoolQualitySnapshot({self.snapshot_id[:12]}, "
            f"retained={self.data['counts']['retained']})"
        )
