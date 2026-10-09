"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/plans/fields.py

Semantic projections for saved recipes and named matrix combinations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from collections import Counter
from dataclasses import asdict

from dense_arrays.planning import GenerationPlan, MatrixPlan, PlanEvidence
from dense_arrays.planning.serialization import policies_for, request_to_dict


def semantic_fields(
    plan: GenerationPlan | PlanEvidence | MatrixPlan,
) -> dict[str, object]:
    """Expose resolved effects and their declarations without source locators."""
    if isinstance(plan, MatrixPlan):
        return {
            "cells": {
                cell.cell_id: generation_fields(cell.plan) for cell in plan.cells
            },
            "cell_order": [cell.cell_id for cell in plan.cells],
            "base": generation_fields(plan.base),
            "matrix": _matrix_fields(plan),
            "policies": {
                key: plan.preview[key]
                for key in (
                    "allocation_policy",
                    "seed_policy",
                    "limits_scope",
                    "uniqueness",
                )
            },
        }
    return generation_fields(plan)


def plan_locations(plan: GenerationPlan | PlanEvidence | MatrixPlan) -> object:
    """Report file bindings separately from design meaning, scoped to each cell."""
    if isinstance(plan, MatrixPlan):
        return {
            "base": plan_locations(plan.base),
            "cells": {cell.cell_id: plan_locations(cell.plan) for cell in plan.cells},
        }
    return (
        [str(item.path) for item in plan.inputs]
        if isinstance(plan, GenerationPlan)
        else None
    )


def _matrix_fields(plan: MatrixPlan) -> dict[str, object]:
    request = plan.request
    sources = {}
    for cell, source in request.sources.items():
        value = source.to_dict()
        value.pop("snapshot_id")
        value.pop("schema")
        value["part_order"] = [p["part_id"] for p in value["parts"]]
        value["parts"] = {p["part_id"]: p for p in value["parts"]}
        value["inputs"] = dict(Counter(source.input_digests))
        value["input_order"] = list(source.input_digests)
        sources[cell] = value
    axes = {}
    for axis, options in request.axes.items():
        axes[axis] = {}
        for choice, variant in options.items():
            value = variant.to_dict()
            value["part_order"] = [p["part_id"] for p in value["parts"]]
            value["requirement_order"] = [r["id"] for r in value["requirements"]]
            value["parts"] = {p["part_id"]: p for p in value["parts"]}
            value["requirements"] = {r["id"]: r for r in value["requirements"]}
            if "add_requirements" in value:
                value["added_requirement_order"] = [
                    r["id"] for r in value["add_requirements"]
                ]
                value["add_requirements"] = {
                    r["id"]: r for r in value["add_requirements"]
                }
            axes[axis][choice] = value
    return {
        "axes": axes,
        "axis_order": list(request.axes),
        "choice_order": {axis: list(options) for axis, options in request.axes.items()},
        "pairing": request.pairing,
        "pairs": [dict(pair) for pair in request.pairs],
        "allocation": request.allocation.to_dict(),
        "max_cells": request.max_cells,
        "sources": sources,
        "batches": {cell: value.to_dict() for cell, value in request.batches.items()},
        "exclude": {cell: value.to_dict() for cell, value in request.exclude.items()},
    }


def generation_fields(plan: GenerationPlan | PlanEvidence) -> dict[str, object]:
    """Normalize identities while preserving order, provenance and policy versions."""
    evidence = plan.evidence if isinstance(plan, GenerationPlan) else plan
    request = request_to_dict(evidence.request)
    request.pop("schema")
    request["part_order"] = [p["part_id"] for p in request["parts"]]
    request["requirement_order"] = [r["id"] for r in request["requirements"]]
    request["parts"] = {p["part_id"]: p for p in request["parts"]}
    request["requirements"] = {r["id"]: r for r in request["requirements"]}
    request["policies"] = policies_for(evidence.request)
    request["inputs"] = dict(Counter(evidence.input_digests))
    request["input_order"] = list(evidence.input_digests)
    request["import_report"] = evidence.import_report.to_dict()
    request["target_scope"] = "initial" if evidence.parent is None else "additional"
    request["exclusions"] = {e.sequence_id: asdict(e) for e in evidence.exclusions}
    if "exclude" in request:
        request["exclude"]["source"].pop("exclusions")
    request["lineage"] = (
        request.get("lineage")
        if evidence.parent is None
        else {
            "run_id": evidence.parent.run_id,
            "revision": evidence.parent.revision,
            "plan_id": evidence.parent.plan_id,
            "accepted_digest": evidence.parent.accepted_digest,
            "request": request.get("lineage"),
        }
    )
    return request
