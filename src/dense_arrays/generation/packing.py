"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/generation/packing.py

Translate supported resolved requirements into the existing packing owner.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dense_arrays.generation.acceptance import realize
from dense_arrays.generation.batches.packing import BatchOptimizer
from dense_arrays.generation.heuristic import GreedySearch
from dense_arrays.generation.objectives import UsageBalancedOptimizer
from dense_arrays.optimizer import Optimizer
from dense_arrays.parts import PartSelector
from dense_arrays.planning import (
    Fixed,
    GenerationPlan,
    GroupCoverage,
    Occurrences,
    Spacing,
)
from dense_arrays.planning.batches.bindings import offered_batch
from dense_arrays.realized import Orientation, RealizedArray
from dense_arrays.sequence import reverse_complement, shift_metric
from dense_arrays.solution import DenseArray
from dense_arrays.solver import SolverControls

PackingEngine = Optimizer | BatchOptimizer | UsageBalancedOptimizer | GreedySearch


def restore_packing(
    packed: RealizedArray, plan: GenerationPlan, *, batch_id: str | None = None
) -> DenseArray:
    """Reconstruct identity-specific offsets without constructing a solver model."""
    request = plan.request
    batch = offered_batch(request, batch_id)
    if batch is not None and {p.feature_id for p in packed.placements} - set(
        batch.part_ids
    ):
        msg = "candidate packing references a part outside its offered batch"
        raise ValueError(msg)
    indices = {part.part_id: i for i, part in enumerate(request.parts)}
    forward, reverse = [None] * len(indices), [None] * len(indices)
    seen = set()
    for placement in packed.placements:
        if placement.feature_id not in indices or placement.feature_id in seen:
            msg = "candidate packing repeats or names an unknown part"
            raise ValueError(msg)
        seen.add(placement.feature_id)
        offsets = forward if placement.orientation is Orientation.FORWARD else reverse
        offsets[indices[placement.feature_id]] = placement.start
    result = DenseArray(
        [part.sequence for part in request.parts],
        request.length.maximum or request.length.exact,
        forward,
        reverse,
    )
    if packed != realize(result, plan, source_id=packed.source_id) or (
        request.strands == "single" and any(offset is not None for offset in reverse)
    ):
        msg = "candidate packing does not match its declared parts and orientations"
        raise ValueError(msg)
    _validate_packing_path(result)
    return result


def _validate_packing_path(solution: DenseArray) -> None:
    """Check adjacent shifts in linear path space, without an all-pairs model."""
    library = solution.library
    oriented = library + [reverse_complement(motif) for motif in library]
    expected_offset = 0
    previous = None
    for offset, index in solution.offset_indices_in_order():
        if previous is not None:
            expected_offset += shift_metric(oriented[previous], oriented[index])
        if offset != expected_offset:
            msg = "candidate offsets do not describe an exact packing path"
            raise ValueError(msg)
        previous = index


def build_optimizer(plan: GenerationPlan, *, seconds: float) -> PackingEngine:
    """Construct the declared search with supported controls and requirements."""
    plan.admit_work()
    if plan.request.batch is not None:
        return BatchOptimizer(plan, seconds=seconds)
    if plan.request.search == "greedy":
        return GreedySearch(plan)
    request = plan.request
    padding = None if request.assembly is None else request.assembly.padding
    optimizer = Optimizer(
        [part.sequence for part in request.parts],
        request.length.maximum or request.length.exact,
        request.strands,
        length_mode="exact"
        if request.length.exact is not None and padding is None
        else "maximum",
    )
    indices = {p.part_id: i for i, p in enumerate(request.parts)}
    for requirement in request.requirements:
        if isinstance(requirement, Fixed):
            window = requirement.start
            bounds = None if window is None else (window.min, window.max)
            origin = "start"
            if bounds is not None and padding is not None and padding.side == "left":
                origin = "end"
                bounds = tuple(
                    None if v is None else v - request.length.exact for v in bounds
                )
            optimizer.add_fixed_occurrence(
                indices[requirement.part_id],
                orientation=requirement.orientation,
                start=bounds,
                origin=origin,
            )
    for requirement in request.requirements:
        if isinstance(requirement, Occurrences):
            optimizer.add_count_constraint(
                list(requirement.select.indices(request.parts)),
                minimum=requirement.min,
                maximum=requirement.max,
            )
        elif isinstance(requirement, GroupCoverage):
            groups = [
                list(PartSelector(groups=(group,)).indices(request.parts))
                for group in requirement.groups
            ]
            optimizer.add_group_coverage(groups, minimum=requirement.min)
        elif isinstance(requirement, Spacing):
            optimizer.add_spacing_constraint(
                indices[requirement.upstream],
                indices[requirement.downstream],
                minimum=requirement.min,
                maximum=requirement.max,
            )
    optimizer.build_model(controls=SolverControls(time_limit_seconds=seconds))
    return (
        UsageBalancedOptimizer(optimizer, tuple(p.part_id for p in request.parts))
        if request.packing_preference is not None
        else optimizer
    )
