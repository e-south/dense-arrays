"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/scoring/configuration.py

Explicit optional scorer settings and finite work limits.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import dataclass, field
from numbers import Real
from pathlib import Path

from dense_arrays._record_validation import integer


def finite(value: object, name: str) -> float:
    """Accept finite real numbers without coercing booleans or strings."""
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
    ):
        msg = f"{name} must be a finite number"
        raise ValueError(msg)
    return float(value)


@dataclass(frozen=True)
class ScoringLimits:
    """Per-call wall time, oriented window count and combined tool output cap."""

    seconds: float = 60.0
    windows: int = 1_000_000
    output_bytes: int = 64 * 1024 * 1024

    def __post_init__(self) -> None:
        """Reject disabled bounds; the window count includes calibration work."""
        if finite(self.seconds, "scoring seconds") <= 0:
            msg = "scoring seconds must be positive"
            raise ValueError(msg)
        integer(self.windows, field_name="scoring windows", minimum=1)
        integer(self.output_bytes, field_name="scoring output_bytes", minimum=1)
        object.__setattr__(self, "seconds", float(self.seconds))


@dataclass(frozen=True)
class FimoScoring:
    """FIMO log-odds scoring with an explicit background and p-value threshold."""

    hit_pvalue_max: float = 1e-4
    background: Path | str | None = None
    strands: str = "double"
    pseudocount: float = 0.1
    executable: Path | str | None = None
    limits: ScoringLimits = field(default_factory=ScoringLimits)

    def __post_init__(self) -> None:
        """Validate settings without discovering or invoking an external tool."""
        if not 0 < finite(self.hit_pvalue_max, "hit_pvalue_max") <= 1:
            msg = "hit_pvalue_max must be in (0, 1]"
            raise ValueError(msg)
        if finite(self.pseudocount, "pseudocount") < 0:
            msg = "pseudocount must be nonnegative"
            raise ValueError(msg)
        if self.strands not in {"single", "double"}:
            msg = "scoring strands must be single or double"
            raise ValueError(msg)
        if not isinstance(self.limits, ScoringLimits):
            msg = "scoring limits must be ScoringLimits"
            raise TypeError(msg)
        object.__setattr__(self, "hit_pvalue_max", float(self.hit_pvalue_max))
        object.__setattr__(self, "pseudocount", float(self.pseudocount))
        for name in ("background", "executable"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, Path(value))
