"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/models.py

Pure request contracts shared by Python, planning and file parsers.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from numbers import Real

from dense_arrays._record_validation import integer
from dense_arrays.constraints import Length
from dense_arrays.parts import BoundParts, Part, PartTable, PoolSource
from dense_arrays.parts.ingestion import validate_parts
from dense_arrays.planning.batches import BatchSchedule, CandidateBatch, Resampling
from dense_arrays.planning.libraries import LibraryExclusion
from dense_arrays.planning.lineage import Lineage
from dense_arrays.planning.requirements import REQUIREMENT_TYPES, Requirement

_DEFAULT_MODEL_PAIRS = 250_000


@dataclass(frozen=True)
class Target:
    """Number of accepted unique designs; zero declares an inactive cell."""

    count: int = 1

    def __post_init__(self) -> None:
        """Require a nonnegative integral target."""
        integer(self.count, field_name="target.count", minimum=0)


@dataclass(frozen=True)
class Limits:
    """Finite effort limits; time bounds are cooperative, not hard deadlines."""

    attempts: int = 1000
    active_seconds: float = 300
    solver_seconds: float = 30
    model_pairs: int = field(default=_DEFAULT_MODEL_PAIRS, kw_only=True)

    def __post_init__(self) -> None:
        """Reject nonfinite, fractional-count or boolean effort limits."""
        integer(self.attempts, field_name="limits.attempts", minimum=1)
        integer(self.model_pairs, field_name="limits.model_pairs", minimum=1)
        for name in ("active_seconds", "solver_seconds"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Real):
                msg = f"limits.{name} must be a positive finite number"
                raise TypeError(msg)
            if not math.isfinite(value) or value <= 0:
                msg = f"limits.{name} must be a positive finite number"
                raise ValueError(msg)

    def to_dict(self) -> dict[str, int | float]:
        """Preserve default request identities while binding explicit size limits."""
        return {
            "attempts": self.attempts,
            "active_seconds": self.active_seconds,
            "solver_seconds": self.solver_seconds,
            **(
                {"model_pairs": self.model_pairs}
                if self.model_pairs != _DEFAULT_MODEL_PAIRS
                else {}
            ),
        }

    def admit_model(self, oriented_nodes: int) -> None:
        """Check quadratic work before adjacency or solver allocation."""
        integer(oriented_nodes, field_name="oriented_nodes", minimum=0)
        pairs = oriented_nodes**2
        if pairs > self.model_pairs:
            msg = (
                f"packing model requires {pairs} oriented pairs, exceeding "
                f"limits.model_pairs={self.model_pairs}; reduce the offered batch "
                "or explicitly raise limits.model_pairs"
            )
            raise ValueError(msg)


@dataclass(frozen=True)
class Padding:
    """Bounded uniform DNA proposals appended to one declared side."""

    side: str
    max_trials: int

    def __post_init__(self) -> None:
        """Reject ambiguous placement or unbounded inner search."""
        if self.side not in {"left", "right"}:
            msg = "padding.side must be left or right"
            raise ValueError(msg)
        integer(self.max_trials, field_name="padding.max_trials", minimum=1)


@dataclass(frozen=True)
class Assembly:
    """Explicit final-length policy; no padding unless requested."""

    padding: Padding | None = None

    def __post_init__(self) -> None:
        """Keep assembly choices typed and separate from packing requirements."""
        if self.padding is not None and not isinstance(self.padding, Padding):
            msg = "assembly.padding must be Padding"
            raise TypeError(msg)


@dataclass(frozen=True)
class DesignSpec:
    """A bounded design request; unspecified policies resolve once in planning."""

    parts: PartTable | PoolSource | BoundParts | tuple[Part, ...]
    length: Length
    requirements: tuple[Requirement, ...] = ()
    target: Target = field(default_factory=Target)
    seed: int = 0
    limits: Limits = field(default_factory=Limits)
    strands: str = "double"
    assembly: Assembly | None = None
    lineage: Lineage | None = None
    exclude: LibraryExclusion | None = None
    batch: CandidateBatch | None = None
    schedule: BatchSchedule | None = None
    resampling: Resampling | None = None
    packing_preference: str | None = None
    search: str = "exact"

    def __post_init__(self) -> None:
        """Freeze caller collections and validate the shared request syntax."""
        if not isinstance(self.parts, (PartTable, PoolSource, BoundParts)):
            object.__setattr__(self, "parts", validate_parts(self.parts))
        for name, expected in (
            ("length", Length),
            ("target", Target),
            ("limits", Limits),
        ):
            if not isinstance(getattr(self, name), expected):
                msg = f"{name} must be {expected.__name__}"
                raise TypeError(msg)
        integer(self.seed, field_name="seed", minimum=0)
        if self.search not in {"exact", "greedy"}:
            msg = "search must be exact or greedy"
            raise ValueError(msg)
        if self.packing_preference not in {None, "underused_parts"}:
            msg = "packing_preference must be underused_parts or None"
            raise ValueError(msg)
        for name, expected in (
            ("batch", CandidateBatch),
            ("schedule", BatchSchedule),
            ("resampling", Resampling),
            ("assembly", Assembly),
            ("lineage", Lineage),
            ("exclude", LibraryExclusion),
        ):
            value = getattr(self, name)
            if value is not None and not isinstance(value, expected):
                msg = f"{name} must be {expected.__name__}"
                raise TypeError(msg)
        if sum(v is not None for v in (self.batch, self.schedule, self.resampling)) > 1:
            msg = "batch, schedule and resampling are mutually exclusive"
            raise ValueError(msg)
        if self.strands not in {"single", "double"}:
            msg = "strands must be single or double"
            raise ValueError(msg)
        self._validate_requirements()

    def _validate_requirements(self) -> None:
        """Freeze typed rules with unique identities independently of search policy."""
        if not isinstance(self.requirements, (list, tuple)) or any(
            not isinstance(r, REQUIREMENT_TYPES) for r in self.requirements
        ):
            msg = "requirements must contain supported typed requirements"
            raise TypeError(msg)
        object.__setattr__(self, "requirements", tuple(self.requirements))
        ids = [r.id for r in self.requirements]
        if len(ids) != len(set(ids)):
            msg = "requirement IDs must be unique"
            raise ValueError(msg)

    def with_changes(self, **changes: object) -> DesignSpec:
        """Return a validated copy without reading inputs, planning or executing."""
        return replace(self, **changes)
