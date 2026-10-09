"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/objectives.py

Recorded packing objectives, independently interpretable from proposed usage.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from dense_arrays._record_validation import (
    immutable_json_mapping,
    integer,
    mutable_json,
    object_fields,
    required_text,
)
from dense_arrays.model import part_usage_weights


def encoded_objective_size(value: object) -> int:
    """Bound both per-part maps before constructing typed objective evidence."""
    if value is None:
        return 0
    if not isinstance(value, Mapping) or any(
        not isinstance(value.get(key), Mapping) for key in ("usage", "weights")
    ):
        msg = "packing objective requires usage and weight maps"
        raise TypeError(msg)
    return len(value["usage"]) + len(value["weights"])


@dataclass(frozen=True)
class PackingObjective:
    """Occurrence maximization with a strictly subordinate part-usage bonus."""

    usage: Mapping[str, int]
    proposed_packings: int

    def __post_init__(self) -> None:
        """Require nonempty occurrence identities and achievable prior counts."""
        integer(self.proposed_packings, field_name="proposed_packings", minimum=0)
        counts = immutable_json_mapping(self.usage)
        for name, count in counts.items():
            required_text(name, field_name="part_id")
            integer(count, field_name=f"usage.{name}", minimum=0)
            if count > self.proposed_packings:
                msg = "part usage exceeds the number of proposed packings"
                raise ValueError(msg)
        part_usage_weights(tuple(counts.values()))
        object.__setattr__(self, "usage", counts)

    def to_dict(self) -> dict[str, object]:
        """Record exact weights and the population from which they were derived."""
        return {
            "schema": "dense_arrays.packing_objective.v1",
            "primary": "selected_occurrences.v1",
            "preference": "underused_parts.v1",
            "population": "proposed_packings_in_batch",
            "proposed_packings": self.proposed_packings,
            "usage": mutable_json(self.usage),
            "weights": dict(
                zip(
                    self.usage,
                    part_usage_weights(tuple(self.usage.values())),
                    strict=True,
                )
            ),
        }

    @classmethod
    def from_dict(cls, value: object) -> PackingObjective:
        """Reject changed algorithms, incomplete evidence and inconsistent weights."""
        keys = {
            "schema",
            "primary",
            "preference",
            "population",
            "proposed_packings",
            "usage",
            "weights",
        }
        data = object_fields(value, keys, "packing objective")
        if set(data) != keys:
            msg = "incomplete packing objective"
            raise ValueError(msg)
        result = cls(data["usage"], data["proposed_packings"])
        weights = object_fields(data["weights"], set(result.usage), "packing weights")
        if (
            any(
                isinstance(w, bool) or not isinstance(w, (int, float))
                for w in weights.values()
            )
            or data != result.to_dict()
        ):
            msg = "packing objective weights or policies disagree with proposed usage"
            raise ValueError(msg)
        return result
