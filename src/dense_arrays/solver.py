"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/solver.py

Supported solver controls and evidence from one packing solve.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum
from numbers import Real
from typing import TYPE_CHECKING

from ._record_validation import object_fields, required_text

if TYPE_CHECKING:
    from .solution import DenseArray


@dataclass(frozen=True)
class SolverIdentity:
    """Requested backend name and version reported by its actual built model."""

    name: str
    version: str

    def __post_init__(self) -> None:
        """Reject missing backend observations instead of inventing version values."""
        required_text(self.name, field_name="solver.name")
        required_text(self.version, field_name="solver.version")

    def to_dict(self) -> dict[str, str]:
        """Return the observed backend identity without creating another model."""
        return {"name": self.name, "version": self.version}

    @classmethod
    def from_dict(cls, value: object) -> SolverIdentity:
        """Read the exact backend fields within a versioned producer record."""
        return cls(**object_fields(value, {"name", "version"}, "solver identity"))


@dataclass(frozen=True)
class SolverControls:
    """Cooperative backend limits, without a hard wall-clock guarantee.

    Time limits round up to the backend's millisecond resolution. Explicit
    thread controls are supported for SCIP only; CBC's OR-Tools interface
    does not honor them. A new backend needs its own qualification.
    """

    time_limit_seconds: float | None = None
    threads: int | None = None

    def __post_init__(self) -> None:
        """Reject controls that cannot be represented by the backend."""
        value = self.time_limit_seconds
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, Real)
            or not math.isfinite(value)
            or not 0 < value * 1000 < 2**63
        ):
            msg = "time_limit_seconds must be positive, finite and fit milliseconds"
            raise ValueError(msg)
        if self.threads is not None and (
            isinstance(self.threads, bool)
            or not isinstance(self.threads, int)
            or not 0 < self.threads < 2**31
        ):
            msg = "threads must be a positive integer"
            raise ValueError(msg)

    @property
    def time_limit_ms(self) -> int | None:
        """The positive millisecond limit passed to the backend."""
        if self.time_limit_seconds is None:
            return None
        return math.ceil(self.time_limit_seconds * 1000)


class SolveStatus(StrEnum):
    """Evidence categories without an inferred reason for backend termination."""

    OPTIMAL = "optimal"
    INFEASIBLE = "infeasible"
    UNPROVEN = "unproven"
    UNKNOWN = "unknown"
    BACKEND_ERROR = "backend_error"
    INVALID_RESULT = "invalid_result"


@dataclass(frozen=True)
class SolveReport:
    """One backend outcome; only optimal, validated results expose a solution.

    Proof applies to the offered packing model, including its current path
    exclusions. It does not establish feasibility of all possible DNA designs.
    Unknown causes remain unknown even when a time limit was configured.
    """

    status: SolveStatus
    backend_status: int | None
    solution: DenseArray | None = None
    proof_scope: str | None = None
    termination_reason: str = "unknown"
    detail: str = ""
