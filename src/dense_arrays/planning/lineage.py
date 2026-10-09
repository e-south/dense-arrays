"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/lineage.py

Descriptive origin identities for editable design requests.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import asdict, dataclass
from pathlib import Path

from dense_arrays._record_validation import (
    digest,
    integer,
    object_fields,
    required_text,
)


@dataclass(frozen=True)
class RunReference:
    """One committed origin revision; no input locator or exclusion policy."""

    run_id: str
    plan_id: str
    revision: int

    def __post_init__(self) -> None:
        """Require explicit origin identity without reading a runtime workspace."""
        required_text(self.run_id, field_name="lineage.parent.run_id")
        digest(self.plan_id, field_name="lineage.parent.plan_id")
        integer(self.revision, field_name="lineage.parent.revision", minimum=0)


@dataclass(frozen=True)
class Lineage:
    """Record where a request originated without excluding accepted sequences."""

    parent: RunReference

    def __post_init__(self) -> None:
        """Keep origin evidence typed and immutable."""
        if not isinstance(self.parent, RunReference):
            msg = "lineage.parent must be RunReference"
            raise TypeError(msg)

    def to_dict(self) -> dict[str, object]:
        """Serialize a location-independent origin reference."""
        return asdict(self)

    @classmethod
    def from_dict(cls, value: object) -> "Lineage":
        """Reject omitted or unknown origin fields."""
        data = object_fields(value, {"parent"}, "lineage")
        parent = object_fields(
            data.get("parent"), {"run_id", "plan_id", "revision"}, "lineage.parent"
        )
        return cls(RunReference(**parent))


@dataclass(frozen=True)
class ParentRun:
    """An explicit native parent location, resolved once into a terminal snapshot."""

    run: Path

    def __post_init__(self) -> None:
        """Normalize location without opening or changing the parent."""
        object.__setattr__(self, "run", Path(self.run))
