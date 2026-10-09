"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/run_plans.py

Decode native execution plans and expose their declared cell recipes.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.reading import ReadLimitError
from dense_arrays.parts.bound import encoded_bound_size
from dense_arrays.planning import GenerationPlan, Limits, MatrixPlan
from dense_arrays.planning.batches.bindings import encoded_batch_size
from dense_arrays.planning.libraries import exclusion_count
from dense_arrays.planning.matrices.resolution import MATRIX_PLAN_SCHEMA

if TYPE_CHECKING:
    from pathlib import Path

type RunPlan = GenerationPlan | MatrixPlan


def validate_run_binding(plan: RunPlan, manifest: dict[str, object]) -> None:
    """Require the origin's full cell inventory to agree with its saved recipes."""
    actual = {
        name: (p.plan_id, p.request.target.count)
        for name, p in cell_plans(plan).items()
    }
    if (
        plan.plan_id != manifest["plan_id"]
        or sum(target for _, target in actual.values()) != manifest["target"]
    ):
        msg = "run plan identity or target mismatch"
        raise ValueError(msg)
    declared = (
        {name: (c["plan_id"], c["target"]) for name, c in manifest["cells"].items()}
        if "cells" in manifest
        else {"default": (manifest["plan_id"], manifest["target"])}
    )
    if actual != declared:
        msg = "cell plan identities or targets do not match the saved run plan"
        raise ValueError(msg)


def cell_plans(plan: RunPlan) -> dict[str, GenerationPlan]:
    """Expose cell recipes without replacing matrix identity with a child plan."""
    plans = (
        {cell.cell_id: cell.plan for cell in plan.cells}
        if isinstance(plan, MatrixPlan)
        else {"default": plan}
    )
    for cell, recipe in plans.items():
        if recipe.request.exclude is not None and set(
            recipe.request.exclude.cell_mapping
        ) != {cell}:
            msg = "exclusion target does not match the execution cell"
            raise ValueError(msg)
    return plans


def excluded_sequences(plan: RunPlan) -> dict[tuple[str, str], str]:
    """Index parent/ancestor sequences by destination cell and full source reference."""
    return {
        (cell, e.sequence_id): e.design_ref
        for cell, recipe in cell_plans(plan).items()
        for e in recipe.exclusions
    }


def run_limits(plan: RunPlan) -> Limits:
    """Return one effort allowance shared by the entire execution."""
    return (
        plan.base.request.limits
        if isinstance(plan, MatrixPlan)
        else plan.request.limits
    )


def decode_plan(
    value: dict[str, object], max_identities: int | None, *, base: Path | None = None
) -> RunPlan:
    """Apply input state bounds before reconstructing saved execution recipes."""
    check_plan_size(value, max_identities)
    cls = MatrixPlan if value.get("schema") == MATRIX_PLAN_SCHEMA else GenerationPlan
    return cls.from_dict(value, base=base)


def check_plan_size(value: dict[str, object], max_identities: int | None) -> None:
    """Admit the encoded state before any typed plan or evidence construction."""
    if not isinstance(value, dict):
        msg = "stored plan must be an object"
        raise TypeError(msg)
    matrix = value.get("schema") == MATRIX_PLAN_SCHEMA
    if max_identities is not None:
        if matrix:
            cells = value.get("cells")
            if not isinstance(cells, list):
                msg = "matrix plan cells must be an array"
                raise ValueError(msg)
            documents = [
                value.get("base"),
                *(c.get("plan") for c in cells if isinstance(c, dict)),
            ]
            size = len(cells) + sum(_plan_size(document) for document in documents)
            request = value.get("request")
            if not isinstance(request, dict) or not isinstance(
                request.get("sources", {}), dict
            ):
                msg = "matrix plan requires a request with a source mapping"
                raise TypeError(msg)
            size += sum(
                encoded_bound_size(s) for s in request.get("sources", {}).values()
            )
        else:
            size = _plan_size(value)
        if size > max_identities:
            msg = (
                f"read_limits.identities={max_identities} "
                "cannot hold the plan identities"
            )
            raise ReadLimitError(msg)


def _plan_size(value: object) -> int:
    if not isinstance(value, dict):
        msg = "stored plan must be an object"
        raise TypeError(msg)
    request = value.get("request")
    if not isinstance(request, dict) or not all(
        isinstance(request.get(name), list) for name in ("parts", "requirements")
    ):
        msg = "stored plan requires complete part and requirement arrays"
        raise TypeError(msg)
    parent = value.get("parent", {})
    if not isinstance(parent, dict) or not isinstance(
        parent.get("exclusions", []), list
    ):
        msg = "plan parent must contain an array of exclusions"
        raise TypeError(msg)
    return (
        encoded_batch_size(request.get("batch"))
        + _schedule_size(request.get("schedule"))
        + len(request["parts"])
        + len(request["requirements"])
        + len(parent.get("exclusions", []))
        + (0 if request.get("exclude") is None else exclusion_count(request["exclude"]))
    )


def _schedule_size(value: object) -> int:
    if value is None:
        return 0
    if not isinstance(value, dict) or not isinstance(value.get("batches"), list):
        msg = "prepared schedule requires an ordered batch array"
        raise TypeError(msg)
    return sum(encoded_batch_size(batch) for batch in value["batches"])
