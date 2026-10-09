"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/batches/feedback.py

Versioned sampling weights and immutable observation snapshots.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from numbers import Real
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    integer,
    object_fields,
    required_text,
)

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from dense_arrays.parts import Part

FEEDBACK_POLICY = "coverage_failure_priority.v1"
FEEDBACK_SCHEMA = "dense_arrays.feedback_snapshot.v1"
_PARAMETERS = {"coverage_alpha", "coverage_power", "failure_alpha", "failure_power"}


@dataclass(frozen=True)
class FeedbackPolicy:
    """Favor underused group/sequence pairs and optionally penalize failures."""

    coverage_alpha: float = 1
    coverage_power: float = 1
    failure_alpha: float = 0
    failure_power: float = 1

    def __post_init__(self) -> None:
        """Require finite nonnegative parameters without implicit boolean values."""
        for name in _PARAMETERS:
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(value)
                or value < 0
            ):
                msg = f"feedback.{name} must be a finite nonnegative number"
                raise ValueError(msg)
            object.__setattr__(self, name, float(value))

    def weight(self, used: int, failed: int) -> float:
        """Compute a positive weight; reject overflow or underflow explicitly."""
        integer(used, field_name="feedback.used", minimum=0)
        integer(failed, field_name="feedback.failed", minimum=0)
        try:
            coverage = (
                1 + self.coverage_alpha / (1 + used) ** self.coverage_power
                if self.coverage_alpha
                else 1
            )
            penalty = (
                (1 + self.failure_alpha * failed) ** self.failure_power
                if self.failure_alpha and failed
                else 1
            )
            weight = coverage / penalty
        except OverflowError as err:
            msg = "feedback weights exceed the supported numeric range"
            raise ValueError(msg) from err
        if not math.isfinite(weight) or weight <= 0:
            msg = "feedback weights exceed the supported numeric range"
            raise ValueError(msg)
        return weight

    def to_dict(self) -> dict[str, object]:
        """Bind formula, ranking version and parameter values."""
        return {
            "policy": FEEDBACK_POLICY,
            **{name: getattr(self, name) for name in sorted(_PARAMETERS)},
        }

    @classmethod
    def from_dict(cls, value: object) -> FeedbackPolicy:
        """Reject unsupported weight/ranking policies and missing parameters."""
        keys = _PARAMETERS | {"policy"}
        data = object_fields(value, keys, "feedback policy")
        if set(data) != keys or data.pop("policy") != FEEDBACK_POLICY:
            msg = "unsupported or incomplete feedback policy"
            raise ValueError(msg)
        return cls(**data)


@dataclass(frozen=True)
class FeedbackSnapshot:
    """Per-part observations, aggregated by group and sequence for weighting."""

    policy: FeedbackPolicy
    used: Mapping[str, int] = field(default_factory=dict)
    failed: Mapping[str, int] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Detach caller counts and require nonnegative integral observations."""
        if not isinstance(self.policy, FeedbackPolicy):
            msg = "feedback snapshot requires FeedbackPolicy"
            raise TypeError(msg)
        for name in ("used", "failed"):
            values = immutable_json_mapping(getattr(self, name))
            for key, count in values.items():
                required_text(key, field_name=f"feedback.{name}.part_id")
                integer(count, field_name=f"feedback.{name}.{key}", minimum=0)
            object.__setattr__(self, name, values)

    def weights(self, parts: Sequence[Part]) -> dict[str, float]:
        """Resolve aliases explicitly; unknown observation identities are errors."""
        known = {p.part_id for p in parts}
        if (self.used.keys() | self.failed.keys()) - known:
            msg = "unknown part identities in sampling feedback"
            raise ValueError(msg)
        used, failed = Counter(), Counter()
        for part in parts:
            key = (part.group, part.sequence)
            used[key] += self.used.get(part.part_id, 0)
            failed[key] += self.failed.get(part.part_id, 0)
        weights = {
            key: self.policy.weight(count, failed[key]) for key, count in used.items()
        }
        return {p.part_id: weights[(p.group, p.sequence)] for p in parts}

    def to_dict(self) -> dict[str, object]:
        """Persist observations without claiming biological interpretation."""
        return {
            "schema": FEEDBACK_SCHEMA,
            "policy": self.policy.to_dict(),
            "used": dict(self.used),
            "failed": dict(self.failed),
        }

    @classmethod
    def from_dict(cls, value: object) -> FeedbackSnapshot:
        """Decode the complete immutable context of a weighted selection."""
        keys = {"schema", "policy", "used", "failed"}
        data = object_fields(value, keys, "feedback snapshot")
        if set(data) != keys or data.pop("schema") != FEEDBACK_SCHEMA:
            msg = "unsupported or incomplete feedback snapshot"
            raise ValueError(msg)
        return cls(
            FeedbackPolicy.from_dict(data["policy"]), data["used"], data["failed"]
        )
