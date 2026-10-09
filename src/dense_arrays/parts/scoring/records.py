"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/scoring/records.py

FIMO score evidence keeps p-values, score units and geometry explicit.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass

from dense_arrays._record_validation import digest, integer, object_fields
from dense_arrays.problem import motif_library

from .configuration import finite


@dataclass(frozen=True)
class FimoHit:
    """One motif-oriented core at zero-based half-open candidate coordinates."""

    start: int
    end: int
    strand: str
    core: str
    raw: float
    pvalue: float
    theoretical_max: float

    def __post_init__(self) -> None:
        """Reject inconsistent geometry, score bounds and probability values."""
        integer(self.start, field_name="hit.start", minimum=0)
        integer(self.end, field_name="hit.end", minimum=1)
        motif_library((self.core,))
        if self.end - self.start != len(self.core):
            msg = "FIMO hit interval must match core length"
            raise ValueError(msg)
        if self.strand not in {"forward", "reverse"}:
            msg = "FIMO hit strand must be forward or reverse"
            raise ValueError(msg)
        for name in ("raw", "pvalue", "theoretical_max"):
            object.__setattr__(self, name, finite(getattr(self, name), f"FIMO {name}"))
        if not 0 <= self.pvalue <= 1:
            msg = "FIMO p-value must be in [0, 1]"
            raise ValueError(msg)
        if self.raw > self.theoretical_max:
            msg = "FIMO score exceeds calibrated theoretical maximum"
            raise ValueError(msg)
        if self.fraction_of_max is not None:
            finite(self.fraction_of_max, "FIMO score ratio")

    @property
    def per_base(self) -> float:
        """Normalize by core width, never by flanking candidate length."""
        return self.raw / len(self.core)

    @property
    def fraction_of_max(self) -> float | None:
        """Expose a ratio only for a positive calibrated denominator."""
        return self.raw / self.theoretical_max if self.theoretical_max > 0 else None

    @property
    def units(self) -> str:
        """Name the backend's reported log-odds units, including its rounding."""
        return "fimo_log2_odds"

    def to_dict(self) -> dict[str, object]:
        """Serialize geometric and statistical evidence without conflating scales."""
        return {
            "schema": "dense_arrays.fimo_hit.v1",
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "core": self.core,
            "raw": self.raw,
            "pvalue": self.pvalue,
            "theoretical_max": self.theoretical_max,
            "per_base": self.per_base,
            "fraction_of_max": self.fraction_of_max,
            "units": self.units,
        }

    @classmethod
    def from_dict(cls, value: object) -> FimoHit:
        """Recompute normalized quantities when loading saved score evidence."""
        fields = {"start", "end", "strand", "core", "raw", "pvalue", "theoretical_max"}
        keys = fields | {"schema", "per_base", "fraction_of_max", "units"}
        data = object_fields(value, keys, "FIMO hit")
        if set(data) != keys or data["schema"] != "dense_arrays.fimo_hit.v1":
            msg = "unsupported or incomplete FIMO hit"
            raise ValueError(msg)
        result = cls(**{k: data[k] for k in fields})
        if data != result.to_dict():
            msg = "FIMO normalized score fields disagree"
            raise ValueError(msg)
        return result


@dataclass(frozen=True)
class FimoResult:
    """One best qualifying hit per candidate, with explicit scoring effort."""

    binding_id: str
    hits: tuple[FimoHit | None, ...]
    candidate_windows: int
    calibration_windows: int
    reported_hits: int
    theoretical_max: float

    def __post_init__(self) -> None:
        """Keep work accounting and the score denominator internally consistent."""
        digest(self.binding_id, field_name="FIMO binding_id")
        if not isinstance(self.hits, (tuple, list)) or not self.hits:
            msg = "FIMO results require a nonempty ordered candidate batch"
            raise ValueError(msg)
        object.__setattr__(self, "hits", tuple(self.hits))
        for name in ("candidate_windows", "calibration_windows", "reported_hits"):
            integer(getattr(self, name), field_name=f"FIMO {name}", minimum=0)
        if (
            self.candidate_windows < self.processed
            or self.calibration_windows not in {1, 2}
            or self.reported_hits > self.candidate_windows
        ):
            msg = "FIMO work counts disagree"
            raise ValueError(msg)
        object.__setattr__(
            self,
            "theoretical_max",
            finite(self.theoretical_max, "FIMO theoretical_max"),
        )
        for hit in self.hits:
            if hit is not None and (
                not isinstance(hit, FimoHit)
                or hit.theoretical_max != self.theoretical_max
            ):
                msg = "FIMO hit and result denominators disagree"
                raise ValueError(msg)
        if sum(hit is not None for hit in self.hits) > self.reported_hits:
            msg = "FIMO hit counts disagree"
            raise ValueError(msg)

    @classmethod
    def from_dict(cls, value: object) -> FimoResult:
        """Read persisted evidence with independent accounting validation."""
        keys = {
            "schema",
            "binding_id",
            "hits",
            "processed",
            "candidate_windows",
            "calibration_windows",
            "reported_hits",
            "theoretical_max",
            "maximum_policy",
        }
        data = object_fields(value, keys, "FIMO result")
        if (
            set(data) != keys
            or data["schema"] != "dense_arrays.fimo_result.v1"
            or not isinstance(data["hits"], list)
        ):
            msg = "unsupported or incomplete FIMO result"
            raise ValueError(msg)
        result = cls(
            data["binding_id"],
            tuple(
                FimoHit.from_dict(hit) if hit is not None else None
                for hit in data["hits"]
            ),
            data["candidate_windows"],
            data["calibration_windows"],
            data["reported_hits"],
            data["theoretical_max"],
        )
        if data != result.to_dict():
            msg = "FIMO result fields disagree"
            raise ValueError(msg)
        return result

    @property
    def processed(self) -> int:
        """Count all scored candidates, including candidates without a hit."""
        return len(self.hits)

    def to_dict(self) -> dict[str, object]:
        """Record absence as null, never as a synthetic zero score."""
        return {
            "schema": "dense_arrays.fimo_result.v1",
            "binding_id": self.binding_id,
            "hits": [hit.to_dict() if hit is not None else None for hit in self.hits],
            "processed": self.processed,
            "candidate_windows": self.candidate_windows,
            "calibration_windows": self.calibration_windows,
            "reported_hits": self.reported_hits,
            "theoretical_max": self.theoretical_max,
            "maximum_policy": "maximizing_core.v1",
        }
