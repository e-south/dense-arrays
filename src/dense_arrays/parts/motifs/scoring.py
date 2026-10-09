"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/scoring.py

Exact supplied-matrix scores and oriented best-hit geometry, without p-values.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import dataclass
from numbers import Real
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    integer,
    object_fields,
)
from dense_arrays.problem import motif_library
from dense_arrays.sequence import reverse_complement

from .models import BASES, Motif

if TYPE_CHECKING:
    from collections.abc import Iterable


@dataclass(frozen=True)
class MotifScore:
    """Explicit score denominators; values describe a matrix, not binding affinity."""

    model_id: str
    raw: float
    width: int
    theoretical_max: float
    units: str = "declared_log_odds"

    def __post_init__(self) -> None:
        """Reject invalid denominators, nonfinite values and mislabeled units."""
        digest(self.model_id, field_name="score.model_id")
        integer(self.width, field_name="score.width", minimum=1)
        if self.units != "declared_log_odds":
            msg = "unsupported motif score units"
            raise ValueError(msg)
        for name in ("raw", "theoretical_max"):
            value = getattr(self, name)
            if (
                isinstance(value, bool)
                or not isinstance(value, Real)
                or not math.isfinite(value)
            ):
                msg = "motif scores must be finite numbers"
                raise ValueError(msg)
            object.__setattr__(self, name, float(value))
        if self.raw > self.theoretical_max:
            msg = "motif score exceeds its theoretical maximum"
            raise ValueError(msg)
        if self.fraction_of_max is not None and not math.isfinite(self.fraction_of_max):
            msg = "motif score ratio exceeds the supported numeric range"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Publish labeled score values and their exact denominators."""
        return {
            "schema": "dense_arrays.motif_score.v1",
            "model_id": self.model_id,
            "raw": self.raw,
            "width": self.width,
            "per_base": self.per_base,
            "theoretical_max": self.theoretical_max,
            "fraction_of_max": self.fraction_of_max,
            "units": self.units,
        }

    @classmethod
    def from_dict(cls, value: object) -> MotifScore:
        """Independently recompute normalized fields when loading evidence."""
        keys = {
            "schema",
            "model_id",
            "raw",
            "width",
            "per_base",
            "theoretical_max",
            "fraction_of_max",
            "units",
        }
        data = object_fields(value, keys, "motif score")
        if set(data) != keys or data["schema"] != "dense_arrays.motif_score.v1":
            msg = "unsupported or incomplete motif score"
            raise ValueError(msg)
        result = cls(
            data["model_id"],
            data["raw"],
            data["width"],
            data["theoretical_max"],
            data["units"],
        )
        if result.to_dict() != data:
            msg = "normalized motif score fields disagree"
            raise ValueError(msg)
        return result

    @property
    def per_base(self) -> float:
        """Normalize by the scored motif width, never by flanking sequence length."""
        return self.raw / self.width

    @property
    def fraction_of_max(self) -> float | None:
        """Return a ratio only when a strictly positive denominator exists."""
        return self.raw / self.theoretical_max if self.theoretical_max > 0 else None


@dataclass(frozen=True)
class MotifHit:
    """Zero-based half-open source coordinates and the motif-oriented core."""

    start: int
    end: int
    strand: str
    core: str
    score: MotifScore

    def __post_init__(self) -> None:
        """Require an oriented core matching its scored source interval."""
        integer(self.start, field_name="hit.start", minimum=0)
        integer(self.end, field_name="hit.end", minimum=1)
        motif_library((self.core,))
        if not isinstance(self.score, MotifScore):
            msg = "hit.score must be MotifScore"
            raise TypeError(msg)
        if (
            self.end - self.start != self.score.width
            or len(self.core) != self.score.width
        ):
            msg = "hit interval and core must match the scored width"
            raise ValueError(msg)
        if self.strand not in {"forward", "reverse"}:
            msg = "hit strand must be forward or reverse"
            raise ValueError(msg)

    def to_dict(self) -> dict[str, object]:
        """Serialize source coordinates separately from motif-oriented sequence."""
        return {
            "schema": "dense_arrays.motif_hit.v1",
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "core": self.core,
            "score": self.score.to_dict(),
        }

    @classmethod
    def from_dict(cls, value: object) -> MotifHit:
        """Decode complete geometry and validated score evidence."""
        keys = {"schema", "start", "end", "strand", "core", "score"}
        data = object_fields(value, keys, "motif hit")
        if set(data) != keys or data.pop("schema") != "dense_arrays.motif_hit.v1":
            msg = "unsupported or incomplete motif hit"
            raise ValueError(msg)
        data["score"] = MotifScore.from_dict(data["score"])
        return cls(**data)


def _sum(values: Iterable[float]) -> float:
    try:
        value = math.fsum(values)
    except OverflowError as err:
        msg = "motif score exceeds the supported numeric range"
        raise ValueError(msg) from err
    if not math.isfinite(value):
        msg = "motif score exceeds the supported numeric range"
        raise ValueError(msg)
    return value


def _raw(motif: Motif, sequence: str) -> float:
    return _sum(
        row[BASES.index(base)]
        for row, base in zip(motif.log_odds, sequence, strict=True)
    )


def score_core(motif: Motif, sequence: str) -> MotifScore:
    """Score exactly one oriented core; never truncate mismatched inputs."""
    _require_scores(motif)
    motif_library((sequence,))
    if len(sequence) != motif.width:
        msg = "scored core length must equal motif width"
        raise ValueError(msg)
    maximum = _sum(max(row) for row in motif.log_odds)
    return MotifScore(motif.model_id, _raw(motif, sequence), motif.width, maximum)


def best_hit(motif: Motif, sequence: str, *, strands: str = "double") -> MotifHit:
    """Choose highest score, then earliest interval, then forward orientation."""
    _require_scores(motif)
    motif_library((sequence,))
    if strands not in {"single", "double"}:
        msg = "scan strands must be single or double"
        raise ValueError(msg)
    if len(sequence) < motif.width:
        msg = "scan sequence is shorter than the motif"
        raise ValueError(msg)
    maximum = _sum(max(row) for row in motif.log_odds)
    best = None
    for start in range(len(sequence) - motif.width + 1):
        core = sequence[start : start + motif.width]
        candidates = (
            (("forward", core),)
            if strands == "single"
            else (("forward", core), ("reverse", reverse_complement(core)))
        )
        for strand, oriented in candidates:
            raw = _raw(motif, oriented)
            if best is None or raw > best.score.raw:
                best = MotifHit(
                    start,
                    start + motif.width,
                    strand,
                    oriented,
                    MotifScore(motif.model_id, raw, motif.width, maximum),
                )
    return best


def _require_scores(motif: Motif) -> None:
    if motif.log_odds is None:
        msg = "motif has no declared score matrix; use an explicit scoring backend"
        raise ValueError(msg)
