"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/acceptance.py

Recompute requirement evidence from final placements independently of the solver.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.generation.geometry import validate_assembly
from dense_arrays.generation.screening import evaluate_screen
from dense_arrays.planning import (
    GC,
    Avoid,
    Fixed,
    GenerationPlan,
    GroupCoverage,
    Occurrences,
    PlanEvidence,
    Spacing,
)
from dense_arrays.planning.batches.bindings import offered_batch
from dense_arrays.realized import Orientation, Placement, PlacementKind, RealizedArray
from dense_arrays.sequence import reverse_complement

if TYPE_CHECKING:
    from dense_arrays.solution import DenseArray


def realize(
    solution: DenseArray, plan: GenerationPlan, *, source_id: str
) -> RealizedArray:
    """Use original occurrence IDs and final, half-open coordinates."""
    placements = []
    for index, part in enumerate(plan.request.parts):
        for orientation, offset in (
            (Orientation.FORWARD, solution.offsets_fwd[index]),
            (Orientation.REVERSE, solution.offsets_rev[index]),
        ):
            if offset is None:
                continue
            sequence = (
                part.sequence
                if orientation is Orientation.FORWARD
                else reverse_complement(part.sequence)
            )
            placements.append(
                Placement(
                    placement_id=f"p{index + 1}",
                    feature_id=part.part_id,
                    kind=PlacementKind.OTHER,
                    sequence=sequence,
                    start=offset,
                    orientation=orientation,
                    label=part.group,
                )
            )
    return RealizedArray(source_id, solution.sequence, tuple(placements))


def _validate_identity(
    realized: RealizedArray, plan: GenerationPlan | PlanEvidence, batch_id: str | None
) -> None:
    """Recount supplied occurrences and reject inconsistent identity/strand evidence."""
    request = plan.request
    by_id = {part.part_id: part for part in request.parts}
    selected = [placement.feature_id for placement in realized.placements]
    batch = offered_batch(request, batch_id)
    if batch is not None and set(selected) - set(batch.part_ids):
        msg = "placements reference a part outside the offered batch"
        raise ValueError(msg)
    if len(selected) != len(set(selected)) or set(selected) - set(by_id):
        msg = "placements repeat a supplied identity or reference an unknown part"
        raise ValueError(msg)
    if (
        request.length.maximum is not None
        and len(realized.sequence) > request.length.maximum
    ):
        msg = "realized sequence exceeds length.maximum"
        raise ValueError(msg)
    if (
        request.length.exact is not None
        and len(realized.sequence) != request.length.exact
    ):
        msg = "realized sequence does not match length.exact"
        raise ValueError(msg)
    for placement in realized.placements:
        expected = by_id[placement.feature_id].sequence
        if placement.orientation is Orientation.REVERSE:
            if request.strands == "single":
                msg = "reverse placement violates single-strand eligibility"
                raise ValueError(msg)
            expected = reverse_complement(expected)
        elif placement.orientation is not Orientation.FORWARD:
            msg = "generated placements require a declared orientation"
            raise ValueError(msg)
        if placement.sequence != expected:
            msg = "placement sequence does not match its supplied part and orientation"
            raise ValueError(msg)


def evaluate(
    realized: RealizedArray,
    plan: GenerationPlan | PlanEvidence,
    *,
    batch_id: str | None = None,
) -> tuple[dict[str, object], ...]:
    """Recount each requirement from final placements, independently of packing."""
    _validate_identity(realized, plan, batch_id)
    validate_assembly(realized, plan)
    return tuple(
        _evaluate_rule(rule, realized, plan) for rule in plan.request.requirements
    )


def _evaluate_rule(
    rule: object, realized: RealizedArray, plan: GenerationPlan | PlanEvidence
) -> dict[str, object]:
    if isinstance(rule, (Avoid, GC)):
        return evaluate_screen(rule, realized)
    request = plan.request
    by_id = {p.part_id: p for p in request.parts}
    selected = {p.feature_id for p in realized.placements}
    placements = {p.feature_id: p for p in realized.placements}
    if isinstance(rule, Occurrences):
        eligible = {
            request.parts[i].part_id for i in rule.select.indices(request.parts)
        }
        observed = len(eligible & set(selected))
        passed = (rule.min is None or observed >= rule.min) and (
            rule.max is None or observed <= rule.max
        )
    elif isinstance(rule, GroupCoverage):
        observed = len(
            {by_id[part_id].group for part_id in selected} & set(rule.groups)
        )
        passed = observed >= rule.min
    elif isinstance(rule, Fixed):
        placement = placements.get(rule.part_id)
        observed = (
            None
            if placement is None
            else {
                "start": placement.start,
                "orientation": "forward"
                if placement.orientation is Orientation.FORWARD
                else "reverse",
            }
        )
        passed = placement is not None and observed["orientation"] == rule.orientation
        if passed and rule.start is not None:
            passed = (rule.start.min is None or placement.start >= rule.start.min) and (
                rule.start.max is None or placement.start <= rule.start.max
            )
    elif isinstance(rule, Spacing):
        upstream, downstream = (
            placements.get(rule.upstream),
            placements.get(rule.downstream),
        )
        observed = (
            None
            if upstream is None or downstream is None
            else downstream.start - upstream.end
        )
        passed = observed is not None and rule.min <= observed <= rule.max
    else:
        msg = f"unsupported requirement type: {type(rule).__name__}"
        raise TypeError(msg)
    return {"id": rule.id, "observed": observed, "passed": passed}
