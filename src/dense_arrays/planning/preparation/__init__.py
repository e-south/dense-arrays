"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/preparation/__init__.py

One preparation-plan interface over curated and sampled source evidence.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays._record_validation import canonical_json
from dense_arrays.artifacts.publication import write_new
from dense_arrays.parts import Part, PartTable, PreparationSet, PreparationSpec
from dense_arrays.parts.ingestion import read_parts
from dense_arrays.planning.resolution import InputBinding

from .curated import PREPARATION_PLAN_SCHEMA, CuratedPreparation
from .requests import (
    PREPARE_SCHEMA,
    PREPARE_SET_SCHEMA,
    PREPARE_WINDOWS_SCHEMA,
    preparation_from_dict,
    preparation_to_dict,
)
from .sampled import SAMPLED_PLAN_SCHEMA, SampledPreparation, resolve_sampled
from .sets import SET_PLAN_SCHEMA, SetPreparation, resolve_set

if TYPE_CHECKING:
    from dense_arrays.parts.provenance import ImportReport

__all__ = [
    "PREPARATION_PLAN_SCHEMA",
    "PREPARE_SCHEMA",
    "PREPARE_SET_SCHEMA",
    "PREPARE_WINDOWS_SCHEMA",
    "SAMPLED_PLAN_SCHEMA",
    "SET_PLAN_SCHEMA",
    "PreparationPlan",
    "preparation_from_dict",
    "preparation_to_dict",
    "resolve_preparation",
]


@dataclass(frozen=True, repr=False)
class PreparationPlan:
    """Immutable source-specific evidence with shared preview and file operations."""

    resolved: CuratedPreparation | SampledPreparation | SetPreparation

    def __post_init__(self) -> None:
        """Require a complete resolved source family."""
        if not isinstance(
            self.resolved, (CuratedPreparation, SampledPreparation, SetPreparation)
        ):
            msg = "PreparationPlan requires resolved curated or sampled evidence"
            raise TypeError(msg)

    @property
    def sampled(self) -> bool:
        """Distinguish candidate mining from deterministic curated filtering."""
        return isinstance(self.resolved, (SampledPreparation, SetPreparation))

    @property
    def request(self) -> PreparationSpec | PreparationSet:
        """Expose the immutable effective request."""
        return self.resolved.request

    @property
    def plan_id(self) -> str:
        """Preserve each schema's location-independent semantic identity."""
        return self.resolved.plan_id

    @property
    def preview(self) -> object:
        """Expose exact curated counts or bounded sampled effort and unknown yield."""
        return self.resolved.preview

    @property
    def identity_count(self) -> int:
        """Count retained source state without creating any sampled candidates."""
        return (
            self.resolved.identity_count if self.sampled else len(self.resolved.parts)
        )

    def _curated(self) -> CuratedPreparation:
        if self.sampled:
            msg = "sampled preparation has no curated parts before execution"
            raise TypeError(msg)
        return self.resolved

    @property
    def parts(self) -> tuple[Part, ...]:
        """Expose supplied curated parts; sampled candidates do not exist in plans."""
        return self._curated().parts

    @property
    def retained_indices(self) -> tuple[int, ...]:
        """Expose exact curated retention only."""
        return self._curated().retained_indices

    @property
    def input(self) -> InputBinding:
        """Expose the curated table's input binding."""
        return self._curated().input

    @property
    def import_report(self) -> ImportReport:
        """Expose the curated table's import transformations."""
        return self._curated().import_report

    def __repr__(self) -> str:
        """Show a bounded plan summary without embedded parts or matrices."""
        return (
            f"PreparationPlan({self.plan_id[:12]}, sampled={self.sampled}, "
            f"retained={self.preview['retained_parts']}, "
            f"budget={self.preview['candidate_budget']})"
        )

    def verify_inputs(self) -> None:
        """Revalidate resolved source bytes before executing the recipe."""
        self.resolved.verify_inputs()

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Serialize the complete source family's versioned plan."""
        return self.resolved.to_dict(base=base)

    def write(self, path: str | Path) -> None:
        """Publish create-only JSON using destination-relative locators."""
        path = Path(path).absolute()
        write_new(path, canonical_json(self.to_dict(base=path.parent)) + "\n")

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> PreparationPlan:
        """Read supported schemas without reapplying defaults or running tools."""
        owner = (
            SetPreparation
            if isinstance(value, dict) and value.get("schema") == SET_PLAN_SCHEMA
            else SampledPreparation
            if isinstance(value, dict) and value.get("schema") == SAMPLED_PLAN_SCHEMA
            else CuratedPreparation
        )
        return cls(owner.from_dict(value, base=base))


def resolve_preparation(request: PreparationSpec | PreparationSet) -> PreparationPlan:
    """Resolve one typed source through its owner without executing preparation."""
    if isinstance(request, PreparationSet):
        return PreparationPlan(resolve_set(request))
    if not isinstance(request.source, PartTable):
        return PreparationPlan(resolve_sampled(request))
    imported = read_parts(request.source)
    return PreparationPlan(
        CuratedPreparation(
            replace(
                request,
                source=replace(request.source, table=request.source.table.absolute()),
            ),
            imported.parts,
            InputBinding(request.source.table, imported.source_digest),
            imported.report,
        )
    )
