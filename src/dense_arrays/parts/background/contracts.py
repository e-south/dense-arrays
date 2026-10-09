"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/background/contracts.py

Resource and outcome contracts for exact conditional background sampling.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import time
from dataclasses import asdict, dataclass

from dense_arrays._record_validation import digest, integer, object_fields
from dense_arrays.parts.scoring.configuration import finite

CONDITIONAL_POLICY = "conditional_background_shake256.v1"


@dataclass(frozen=True)
class ConditionalLimits:
    """Bound table construction; seconds is a cooperative allowance, not a deadline."""

    states: int = 250_000
    automaton_states: int = 4096
    mass_bits: int = 32_000_000
    seconds: float = 30.0

    def __post_init__(self) -> None:
        """Require positive admission caps for tables and exact integer storage."""
        for name in ("states", "automaton_states", "mass_bits"):
            integer(getattr(self, name), field_name=f"conditional.{name}", minimum=1)
        if finite(self.seconds, "conditional.seconds") <= 0:
            msg = "conditional.seconds must be positive"
            raise ValueError(msg)
        object.__setattr__(self, "seconds", float(self.seconds))


class ConstructionLimitError(Exception):
    """A counting resource limit carries no feasibility conclusion."""

    def __init__(self, reason: str) -> None:
        """Retain the exact exhausted allowance."""
        self.reason = reason
        super().__init__(reason)


@dataclass
class ConstructionWork:
    """Execution-local work counters shared by automaton and completion tables."""

    limits: ConditionalLimits
    deadline: float
    states: int = 0
    automaton_states: int = 0
    mass_bits: int = 0

    def check_time(self) -> None:
        """Observe the shared cooperative deadline before the next unit of work."""
        if time.monotonic() >= self.deadline:
            reason = "time_budget"
            raise ConstructionLimitError(reason)

    def admit(self, name: str, amount: int = 1) -> None:
        """Stop before allocating a unit beyond the declared work allowance."""
        self.check_time()
        value = getattr(self, name) + amount
        if value > getattr(self.limits, name):
            raise ConstructionLimitError(name)
        setattr(self, name, value)

    def multiply(self, first: int, second: int) -> int:
        """Bound temporary integers as well as retained completion masses."""
        if (
            first
            and second
            and first.bit_length() + second.bit_length() - 1 > self.limits.mass_bits
        ):
            reason = "mass_bits"
            raise ConstructionLimitError(reason)
        result = first * second
        self.check_number(result)
        return result

    def check_number(self, value: int) -> None:
        """Reject a single oversized temporary mass before retaining it."""
        if value.bit_length() > self.limits.mass_bits:
            reason = "mass_bits"
            raise ConstructionLimitError(reason)

    def power(self, base: int, exponent: int) -> int:
        """Compute a bounded power without first allocating an enormous integer."""
        result = 1
        while exponent:
            self.check_time()
            if exponent & 1:
                result = self.multiply(result, base)
            exponent >>= 1
            if exponent:
                base = self.multiply(base, base)
        return result


@dataclass(frozen=True)
class ConstructionReport:
    """Recorded exact-support result, separate from candidate and retention counts."""

    model_id: str
    status: str
    reason: str
    states: int
    automaton_states: int
    mass_bits: int
    mass: str | None

    def __post_init__(self) -> None:
        """Reject contradictory completion, zero-support and limited outcomes."""
        digest(self.model_id, field_name="conditional model_id")
        for name in ("states", "automaton_states", "mass_bits"):
            integer(getattr(self, name), field_name=f"conditional.{name}", minimum=0)
        if self.status not in {"feasible", "infeasible", "limited"}:
            msg = "unsupported conditional construction status"
            raise ValueError(msg)
        if self.status == "limited":
            if (
                self.reason
                not in {"states", "automaton_states", "mass_bits", "time_budget"}
                or self.mass is not None
            ):
                msg = (
                    "limited construction requires its resource cause and unknown mass"
                )
                raise ValueError(msg)
        else:
            if self.reason != "completed" or not isinstance(self.mass, str):
                msg = "completed construction requires exact mass"
                raise ValueError(msg)
            try:
                value = int(self.mass, 16)
            except ValueError as err:
                msg = "conditional mass requires canonical hexadecimal"
                raise ValueError(msg) from err
            if (
                value < 0
                or hex(value) != self.mass
                or (value > 0) != (self.status == "feasible")
                or value.bit_length() > self.mass_bits
            ):
                msg = "conditional mass disagrees with its construction outcome"
                raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Encode exact masses without decimal-integer serialization limits."""
        return {
            "schema": "dense_arrays.conditional_construction.v1",
            "policy": CONDITIONAL_POLICY,
            **asdict(self),
        }

    @classmethod
    def from_dict(cls, value: object) -> ConstructionReport:
        """Restore complete construction evidence without rebuilding the table."""
        keys = {
            "model_id",
            "status",
            "reason",
            "states",
            "automaton_states",
            "mass_bits",
            "mass",
        }
        data = object_fields(
            value, keys | {"schema", "policy"}, "conditional construction"
        )
        if (
            data.pop("schema", None) != "dense_arrays.conditional_construction.v1"
            or data.pop("policy", None) != CONDITIONAL_POLICY
            or set(data) != keys
        ):
            msg = "unsupported or incomplete conditional construction record"
            raise ValueError(msg)
        return cls(**data)
