"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/plans/comparison.py

Compare semantic fields by domain identity, retaining ordering and locations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from functools import cached_property

from dense_arrays._record_validation import (
    canonical_json,
    immutable_json_mapping,
    mutable_json,
)
from dense_arrays.planning import GenerationPlan, MatrixPlan, PlanEvidence

from .fields import plan_locations, semantic_fields


@dataclass(frozen=True)
class PlanChange:
    """One added, removed or changed value at an identity-based semantic path."""

    path: tuple[str, ...]
    kind: str
    before: object
    after: object

    def __post_init__(self) -> None:
        """Freeze exact JSON values and reject ambiguous change descriptions."""
        if self.kind not in {"added", "removed", "changed"}:
            msg = "unsupported semantic change kind"
            raise ValueError(msg)
        if (
            not isinstance(self.path, (list, tuple))
            or not self.path
            or not self.path[0]
            or any(not isinstance(part, str) for part in self.path)
        ):
            msg = "semantic change paths require a named root and string components"
            raise ValueError(msg)
        values = immutable_json_mapping({"before": self.before, "after": self.after})
        if self.kind == "added" and values["before"] is not None:
            msg = "an added value requires an absent before value"
            raise ValueError(msg)
        if self.kind == "removed" and values["after"] is not None:
            msg = "a removed value requires an absent after value"
            raise ValueError(msg)
        if self.kind == "changed" and _same_value(values["before"], values["after"]):
            msg = "a changed value requires unequal before and after values"
            raise ValueError(msg)
        object.__setattr__(self, "path", tuple(self.path))
        object.__setattr__(self, "before", values["before"])
        object.__setattr__(self, "after", values["after"])

    @property
    def pointer(self) -> str:
        """Escape each identity as a JSON Pointer for unambiguous display."""
        return "/" + "/".join(
            p.replace("~", "~0").replace("/", "~1") for p in self.path
        )

    def to_dict(self) -> dict[str, object]:
        """Keep missing-side nulls distinct from changed nullable values by kind."""
        return {
            "path": list(self.path),
            "kind": self.kind,
            "before": mutable_json(self.before),
            "after": mutable_json(self.after),
        }


@dataclass(frozen=True, repr=False)
class PlanComparison:
    """Actionable generation changes, separate from original input locations."""

    before: GenerationPlan | PlanEvidence | MatrixPlan
    after: GenerationPlan | PlanEvidence | MatrixPlan

    def __post_init__(self) -> None:
        """Accept executable plans and their location-free evidence records."""
        generation = all(
            isinstance(p, (GenerationPlan, PlanEvidence))
            for p in (self.before, self.after)
        )
        matrix = all(isinstance(p, MatrixPlan) for p in (self.before, self.after))
        if not (generation or matrix):
            msg = (
                "comparison requires matching plan kinds: two matrix plans "
                "or two generation plans/evidence records"
            )
            raise TypeError(msg)

    @cached_property
    def changes(self) -> tuple[PlanChange, ...]:
        """Expose immutable changes, using IDs for parts, rules and exclusions."""
        return tuple(
            _differences(semantic_fields(self.before), semantic_fields(self.after))
        )

    @property
    def changed_fields(self) -> tuple[str, ...]:
        """Name changed categories in the normalized plan's field order."""
        return tuple(dict.fromkeys(change.path[0] for change in self.changes))

    @property
    def unchanged_fields(self) -> tuple[str, ...]:
        """Name inherited categories without confusing a reorder with value edits."""
        changed = set(self.changed_fields)
        return tuple(
            name for name in semantic_fields(self.before) if name not in changed
        )

    def to_dict(self) -> dict[str, object]:
        """Serialize exact changes without implying causality or significance."""
        return {
            "schema": "dense_arrays.matrix_plan_comparison.v1"
            if isinstance(self.before, MatrixPlan)
            else "dense_arrays.plan_comparison.v2",
            "before_plan_id": self.before.plan_id,
            "after_plan_id": self.after.plan_id,
            "changed_fields": list(self.changed_fields),
            "unchanged_fields": list(self.unchanged_fields),
            "changes": [change.to_dict() for change in self.changes],
            "locations": {
                name: plan_locations(plan)
                for name, plan in (("before", self.before), ("after", self.after))
            },
        }

    def __repr__(self) -> str:
        """Describe identities without traversing parts or materializing differences."""
        return (
            f"PlanComparison(before={self.before.plan_id[:12]}, "
            f"after={self.after.plan_id[:12]})"
        )


def _differences(
    before: object, after: object, path: tuple[str, ...] = ()
) -> Iterator[PlanChange]:
    """Compare named values recursively; retain ordered arrays as explicit values."""
    if isinstance(before, Mapping) and isinstance(after, Mapping):
        for name in dict.fromkeys((*before, *after)):
            child = (*path, name)
            if name not in before:
                yield PlanChange(child, "added", None, after[name])
            elif name not in after:
                yield PlanChange(child, "removed", before[name], None)
            else:
                yield from _differences(before[name], after[name], child)
    elif not _same_value(before, after):
        yield PlanChange(path, "changed", before, after)


def _same_value(before: object, after: object) -> bool:
    """Match the canonical wire identity, including booleans and numeric encoding."""
    return canonical_json(mutable_json(before)) == canonical_json(mutable_json(after))
