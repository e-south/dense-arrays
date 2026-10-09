"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/matrices/resolution.py

Resolved combinations bind part sources and independent named random streams.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import InitVar, dataclass, field, replace
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    canonical_json,
    object_fields,
    semantic_digest,
)
from dense_arrays.artifacts.publication import write_new
from dense_arrays.parts import BoundParts
from dense_arrays.planning.batches import BatchSchedule, CandidateBatch, Resampling
from dense_arrays.planning.models import Target
from dense_arrays.planning.resolution import (
    GenerationPlan,
    bind_design,
    resolve_design,
    resolve_source,
)

from .expansion import cell_id, expand, expansion_size
from .requests import MatrixSpec

if TYPE_CHECKING:
    from collections.abc import Mapping

    from dense_arrays.planning.requirements import Requirement

MATRIX_PLAN_SCHEMA = "dense_arrays.matrix_plan.v1"
SEED_POLICY = "matrix_cell_sha256.v1"


@dataclass(frozen=True)
class MatrixCell:
    """One design combination, identified by its named axis choices within a matrix."""

    cell_id: str
    choices: Mapping[str, str]
    plan: GenerationPlan

    def __post_init__(self) -> None:
        """Detach choice metadata from mutable caller mappings."""
        object.__setattr__(self, "choices", MappingProxyType(dict(self.choices)))

    @property
    def target(self) -> int:
        """Read the sole target owner in this resolved cell request."""
        return self.plan.request.target.count

    @property
    def active(self) -> bool:
        """Zero targets remain listed without admitting generation."""
        return self.target > 0

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Encode the complete cell request and its explicit allocation."""
        return {
            "cell_id": self.cell_id,
            "choices": dict(self.choices),
            "target": self.target,
            "active": self.active,
            "plan": self.plan.to_dict(base=base),
        }


@dataclass(frozen=True, repr=False)
class MatrixPlan:
    """Resolved cell recipes with explicit targets and one shared run budget."""

    request: MatrixSpec
    base: GenerationPlan
    cells: tuple[MatrixCell, ...] = field(init=False)
    plan_id: str = field(init=False)
    _execution_supported: bool = field(
        default=True, repr=False, compare=False, kw_only=True
    )

    _encoded_cells: InitVar[list[dict[str, object]] | None] = field(
        default=None, kw_only=True
    )

    def admit_work(self) -> None:
        """Check active cell models without charging the base or inactive cells."""
        for cell in self.cells:
            cell.plan.admit_work()

    def __post_init__(self, _encoded_cells: list[dict[str, object]] | None) -> None:
        """Compile immutable cells and bind their ordered semantic identities."""
        if not isinstance(self._execution_supported, bool):
            msg = "matrix preview.execution_supported must be a boolean"
            raise TypeError(msg)
        if not isinstance(self.request, MatrixSpec) or not isinstance(
            self.base, GenerationPlan
        ):
            msg = "matrix plan requires a MatrixSpec and resolved GenerationPlan"
            raise TypeError(msg)
        if self.base.parent is not None or self.base.request.exclude is not None:
            msg = "matrix base exclusions require explicit per-cell declarations"
            raise ValueError(msg)
        if self.request.base != self.base.request:
            msg = "matrix base must match its resolved source recipe"
            raise ValueError(msg)
        if any(not isinstance(s, BoundParts) for s in self.request.sources.values()):
            msg = (
                "matrix plan sources require frozen BoundParts; "
                "use plan to resolve sources"
            )
            raise TypeError(msg)
        combinations = expand(self.request)
        targets = _targets(self.request, combinations)
        if _encoded_cells is not None:
            _validate_encoded_cells(self.request, combinations, _encoded_cells)
        cells = tuple(
            _cell(self.request, self.base, choices, target)
            for choices, target in zip(combinations, targets, strict=True)
        )
        object.__setattr__(self, "cells", cells)
        # Input locations do not participate in matrix identity.
        shape = self.request.to_dict()
        shape.pop("base")
        if self.request.sources:
            shape["sources"] = {
                cell: source.snapshot_id
                for cell, source in self.request.sources.items()
            }
        object.__setattr__(
            self,
            "plan_id",
            semantic_digest(
                {
                    "schema": MATRIX_PLAN_SCHEMA,
                    "base_plan_id": self.base.plan_id,
                    "matrix": shape,
                    "seed_policy": SEED_POLICY,
                    "cells": [
                        {"cell_id": c.cell_id, "plan_id": c.plan.plan_id} for c in cells
                    ],
                }
            ),
        )

    def __repr__(self) -> str:
        """Show shape and identity without expanding embedded requests."""
        return (
            f"MatrixPlan({self.plan_id[:12]}, "
            f"cells={len(self.cells)}, target={self.total})"
        )

    @property
    def total(self) -> int:
        """Return the exact sum of all explicit cell targets."""
        return sum(cell.target for cell in self.cells)

    @property
    def preview(self) -> Mapping[str, object]:
        """Expose allocation, stream policy and current execution capability."""
        return MappingProxyType(
            {
                "cells": len(self.cells),
                "active_cells": sum(c.active for c in self.cells),
                "total": self.total,
                "allocation_policy": self.request.allocation.policy_version,
                "seed_policy": SEED_POLICY,
                "execution_supported": self._execution_supported,
                "limits_scope": "matrix_run",
                "uniqueness": "exact_sequence_per_cell.v1",
                **(
                    {
                        "excluded_sequences": sum(
                            len(c.plan.exclusions) for c in self.cells
                        )
                    }
                    if self.request.exclude
                    else {}
                ),
            }
        )

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Serialize a complete immutable preview with document-relative locators."""
        return {
            "schema": MATRIX_PLAN_SCHEMA,
            "plan_id": self.plan_id,
            "request": self.request.to_dict(base=base),
            "base": self.base.to_dict(base=base),
            "cells": [cell.to_dict(base=base) for cell in self.cells],
            "preview": dict(self.preview),
        }

    def write(self, path: str | Path) -> None:
        """Save a new matrix plan without replacing an existing file."""
        path = Path(path).absolute()
        write_new(path, canonical_json(self.to_dict(base=path.parent)) + "\n")

    def verify_inputs(self) -> None:
        """Check each distinct source binding once before any combination executes."""
        bindings = dict.fromkeys(
            item
            for plan in (self.base, *(cell.plan for cell in self.cells))
            for item in plan.inputs
        )
        for item in bindings:
            item.verify()

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> MatrixPlan:
        """Reconstruct from frozen source records, never reopening inputs."""
        keys = {"schema", "plan_id", "request", "base", "cells", "preview"}
        data = object_fields(value, keys, "matrix plan")
        if set(data) != keys or data["schema"] != MATRIX_PLAN_SCHEMA:
            msg = "unsupported or incomplete matrix plan schema"
            raise ValueError(msg)
        preview = object_fields(
            data["preview"],
            {
                "cells",
                "active_cells",
                "total",
                "allocation_policy",
                "seed_policy",
                "execution_supported",
                "limits_scope",
                "uniqueness",
                "excluded_sequences",
            },
            "matrix preview",
        )
        if "execution_supported" not in preview:
            msg = "matrix preview requires execution_supported"
            raise ValueError(msg)
        request = MatrixSpec.from_dict(data["request"], base=base)
        encoded_cells = data["cells"]
        if not isinstance(encoded_cells, list) or len(encoded_cells) != (
            expansion_size(request)
        ):
            msg = "saved matrix cell count does not match request cardinality"
            raise ValueError(msg)
        result = cls(
            request,
            GenerationPlan.from_dict(data["base"], base=base),
            _execution_supported=preview["execution_supported"],
            _encoded_cells=encoded_cells,
        )
        if result.to_dict(base=base) != data:
            msg = "matrix plan does not match its saved identities, targets or policies"
            raise ValueError(msg)
        return result


def resolve_matrix(request: MatrixSpec) -> MatrixPlan:
    """Validate combination references, then bind default and overridden sources."""
    choices = expand(request)
    _targets(request, choices)
    if request.base.exclude is not None:
        msg = "matrix exclusions require explicit per-cell declarations"
        raise ValueError(msg)
    base = resolve_design(request.base)
    sources = {cell: resolve_source(source) for cell, source in request.sources.items()}
    result = MatrixPlan(replace(request, base=base.request, sources=sources), base)
    result.verify_inputs()
    return result


def _targets(
    request: MatrixSpec, combinations: tuple[Mapping[str, str], ...]
) -> tuple[int, ...]:
    """Check allocation and exclusion identities before reading source parts."""
    identities = tuple(cell_id(c) for c in combinations)
    if set(request.exclude) - set(identities):
        msg = "matrix exclusions reference unknown cells"
        raise ValueError(msg)
    if set(request.batches) - set(identities):
        msg = "matrix batches reference unknown cells"
        raise ValueError(msg)
    if set(request.sources) - set(identities):
        msg = "matrix sources reference unknown cells"
        raise ValueError(msg)
    return request.allocation.resolve(identities)


def _cell(
    request: MatrixSpec, base: GenerationPlan, choices: Mapping[str, str], target: int
) -> MatrixCell:
    identity = cell_id(choices)
    source = request.sources.get(identity)
    available = base.request.parts if source is None else source.parts
    parts, rules = {}, {}
    for axis, choice in choices.items():
        variant = request.axes[axis][choice]
        for item in variant.parts:
            if item.part_id in parts:
                msg = f"axes both replace part {item.part_id!r}"
                raise ValueError(msg)
            parts[item.part_id] = item
        for rule in variant.requirements:
            if rule.id in rules:
                msg = f"axes both replace requirement {rule.id!r}"
                raise ValueError(msg)
            rules[rule.id] = rule
    if set(parts) - {p.part_id for p in available} or set(rules) - {
        r.id for r in base.request.requirements
    }:
        msg = "matrix variants must replace known part and requirement identities"
        raise ValueError(msg)
    additions = _added_rules(request, choices, rules)
    seed = int(
        semantic_digest(
            {"schema": SEED_POLICY, "seed": base.request.seed, "cell_id": identity}
        ),
        16,
    )
    resolved = base.request.with_changes(
        parts=tuple(parts.get(p.part_id, p) for p in available),
        requirements=tuple(rules.get(r.id, r) for r in base.request.requirements)
        + additions,
        target=Target(count=target),
        seed=seed,
        **_cell_policies(request, identity),
    )
    compiled = (
        replace(base, request=resolved)
        if source is None
        else bind_design(resolved, replace(source, parts=resolved.parts))
    )
    return MatrixCell(identity, choices, compiled)


def _added_rules(
    request: MatrixSpec,
    choices: Mapping[str, str],
    replacements: Mapping[str, Requirement],
) -> tuple[Requirement, ...]:
    """Append declared rules in axis order, rejecting every identity collision."""
    seen = {r.id for r in request.base.requirements} | set(replacements)
    additions = []
    for axis, choice in choices.items():
        for rule in request.axes[axis][choice].add_requirements:
            if rule.id in seen:
                msg = (
                    f"added requirement {rule.id!r} conflicts with another requirement"
                )
                raise ValueError(msg)
            seen.add(rule.id)
            additions.append(rule)
    return tuple(additions)


def _cell_policies(request: MatrixSpec, identity: str) -> dict[str, object]:
    """Resolve the one declared owner of each cell's batch and exclusion policy."""
    batch = request.batches.get(identity)
    return {
        "exclude": request.exclude.get(identity),
        "batch": batch if isinstance(batch, CandidateBatch) else None,
        "schedule": batch if isinstance(batch, BatchSchedule) else None,
        "resampling": (batch if isinstance(batch, Resampling) else None)
        if identity in request.batches
        else request.base.resampling,
    }


def _validate_encoded_cells(
    request: MatrixSpec,
    combinations: tuple[Mapping[str, str], ...],
    encoded: list[dict[str, object]],
) -> None:
    """Bind admitted child dimensions before reconstructing any child plan.

    Use the already expanded canonical choices. No part/rule copies, child
    DesignSpecs or GenerationPlans are created during this structural pass.
    Full content and identity equality remain checked after reconstruction.
    """
    for choices, cell in zip(combinations, encoded, strict=True):
        identity = cell_id(choices)
        if (
            not isinstance(cell, dict)
            or cell.get("cell_id") != identity
            or cell.get("choices") != choices
        ):
            msg = (
                "matrix plan does not match its saved identities, targets or policies: "
                "cell order or choices disagree with the request"
            )
            raise ValueError(msg)
        plan = cell.get("plan")
        child = plan.get("request") if isinstance(plan, dict) else None
        if not isinstance(child, dict):
            msg = "saved matrix cell requires a plan with a request"
            raise TypeError(msg)
        source = request.sources.get(identity)
        sizes = {
            "parts": len(request.base.parts if source is None else source.parts),
            "requirements": len(request.base.requirements)
            + sum(
                len(request.axes[axis][choice].add_requirements)
                for axis, choice in choices.items()
            ),
        }
        for name, size in sizes.items():
            if not isinstance(child.get(name), list) or len(child[name]) != size:
                msg = f"saved matrix cell {identity!r} has inconsistent {name} count"
                raise ValueError(msg)
        for name, policy in _cell_policies(request, identity).items():
            expected = None if policy is None else policy.to_dict()
            if child.get(name) != expected:
                msg = f"saved matrix cell {identity!r} has inconsistent {name} policy"
                raise ValueError(msg)
