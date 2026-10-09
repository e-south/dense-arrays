"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/requirements.py

Typed hard requirements for supplied occurrences and final geometry.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass

from dense_arrays._record_validation import integer, required_text
from dense_arrays.constraints import GC, Avoid
from dense_arrays.parts import PartSelector


@dataclass(frozen=True)
class Occurrences:
    """Independent inclusive bounds on selected supplied identities."""

    id: str
    select: PartSelector
    min: int | None = None
    max: int | None = None

    def __post_init__(self) -> None:
        """Validate count syntax before resolving eligible identities."""
        required_text(self.id, field_name="requirement.id")
        if not isinstance(self.select, PartSelector):
            msg = "select must be PartSelector"
            raise TypeError(msg)
        if self.min is None and self.max is None:
            msg = "occurrences requires a minimum or maximum bound"
            raise ValueError(msg)
        for name in ("min", "max"):
            if (value := getattr(self, name)) is not None:
                integer(value, field_name=f"{self.id}.{name}", minimum=0)
        if self.min is not None and self.max is not None and self.min > self.max:
            msg = f"{self.id}: minimum exceeds maximum"
            raise ValueError(msg)


@dataclass(frozen=True)
class GroupCoverage:
    """Minimum distinct represented caller-defined groups."""

    id: str
    groups: tuple[str, ...]
    min: int

    def __post_init__(self) -> None:
        """Freeze group labels and reject static coverage contradictions."""
        required_text(self.id, field_name="requirement.id")
        selected = PartSelector(groups=self.groups)
        object.__setattr__(self, "groups", selected.groups)
        integer(self.min, field_name=f"{self.id}.min", minimum=1)
        if self.min > len(self.groups):
            msg = f"{self.id}: minimum exceeds available groups"
            raise ValueError(msg)


@dataclass(frozen=True)
class StartWindow:
    """Inclusive bounds on a placement start in final-sequence coordinates."""

    min: int | None = None
    max: int | None = None

    def __post_init__(self) -> None:
        """Require at least one nonnegative bound in increasing order."""
        if self.min is None and self.max is None:
            msg = "start window requires min or max"
            raise ValueError(msg)
        for name in ("min", "max"):
            if (value := getattr(self, name)) is not None:
                integer(value, field_name=f"start.{name}", minimum=0)
        if self.min is not None and self.max is not None and self.min > self.max:
            msg = "start minimum exceeds maximum"
            raise ValueError(msg)


@dataclass(frozen=True)
class Fixed:
    """Require one named occurrence with an explicit orientation and start window."""

    id: str
    part_id: str
    orientation: str
    start: StartWindow | None = None

    def __post_init__(self) -> None:
        """Check syntax without looking up a sequence or guessing an orientation."""
        required_text(self.id, field_name="requirement.id")
        required_text(self.part_id, field_name="fixed.part_id")
        if self.orientation not in {"forward", "reverse"}:
            msg = "fixed orientation must be forward or reverse"
            raise ValueError(msg)
        if self.start is not None and not isinstance(self.start, StartWindow):
            msg = "fixed start must be StartWindow"
            raise TypeError(msg)


@dataclass(frozen=True)
class Spacing:
    """Signed downstream-start minus upstream-end distance in bases."""

    id: str
    upstream: str
    downstream: str
    min: int
    max: int

    def __post_init__(self) -> None:
        """Require distinct named parts and inclusive signed integer bounds."""
        for name in ("id", "upstream", "downstream"):
            required_text(getattr(self, name), field_name=f"spacing.{name}")
        for name in ("min", "max"):
            integer(getattr(self, name), field_name=f"spacing.{name}")
        if self.upstream == self.downstream or self.min > self.max:
            msg = "spacing requires distinct parts and min <= max"
            raise ValueError(msg)


type Requirement = Occurrences | GroupCoverage | Fixed | Spacing | Avoid | GC
REQUIREMENT_TYPES = (Occurrences, GroupCoverage, Fixed, Spacing, Avoid, GC)
