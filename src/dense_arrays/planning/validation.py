"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/validation.py

Resolve cross-requirement references and static contradictions before search.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dense_arrays.constraints import count_bounds
from dense_arrays.parts import PartSelector
from dense_arrays.planning.diagnostics import insufficient_parts
from dense_arrays.planning.models import DesignSpec
from dense_arrays.planning.requirements import (
    GC,
    Avoid,
    Fixed,
    GroupCoverage,
    Occurrences,
    Spacing,
)


def validate_requirements(request: DesignSpec) -> None:
    """Reject unsupported geometry and impossible identity/count declarations."""
    _validate_search(request)
    fixed_ids = _validate_fixed(request)
    _validate_spacing(request, fixed_ids)
    _validate_counts(request, fixed_ids)
    _validate_screens(request, fixed_ids)


def _validate_search(request: DesignSpec) -> None:
    """Admit only requirements enforced by the selected packing method."""
    if request.search != "greedy":
        return
    if request.packing_preference is not None:
        msg = "greedy search does not support packing_preference"
        raise ValueError(msg)
    for rule in request.requirements:
        if not isinstance(rule, (Avoid, GC)):
            msg = f"{rule.id}: greedy search does not support {type(rule).__name__}"
            raise ValueError(msg)  # noqa: TRY004 - valid rule, unsupported method
    if request.length.exact is not None and (
        request.assembly is None or request.assembly.padding is None
    ):
        msg = "greedy search requires padding for length.exact"
        raise ValueError(msg)


def _validate_screens(request: DesignSpec, fixed_ids: set[str]) -> None:
    for rule in request.requirements:
        if isinstance(rule, Avoid) and set(rule.except_placements) - fixed_ids:
            msg = f"{rule.id}: exceptions must name fixed part IDs"
            raise ValueError(msg)
        if (
            isinstance(rule, GC)
            and rule.scope == "padding"
            and (request.assembly is None or request.assembly.padding is None)
        ):
            msg = f"{rule.id}: padding GC requires an explicit padding policy"
            raise ValueError(msg)


def _validate_fixed(request: DesignSpec) -> set[str]:
    fixed = [r for r in request.requirements if isinstance(r, Fixed)]
    fixed_ids = {r.part_id for r in fixed}
    by_id = {p.part_id: p for p in request.parts}
    if len(fixed_ids) != len(fixed):
        msg = "only one fixed declaration is allowed per part"
        raise ValueError(msg)
    if fixed_ids - set(by_id):
        missing = sorted(fixed_ids - set(by_id))
        msg = f"fixed requirements reference unknown parts: {missing}"
        raise ValueError(msg)
    for rule in fixed:
        if rule.orientation == "reverse" and request.strands == "single":
            msg = f"{rule.id}: reverse fixed occurrence requires double strands"
            raise ValueError(msg)
        bound = request.length.maximum or request.length.exact
        if (
            rule.start is not None
            and (rule.start.min or 0) + len(by_id[rule.part_id].sequence) > bound
        ):
            msg = f"{rule.id}: fixed start cannot fit within the declared length"
            raise ValueError(msg)
    return fixed_ids


def _validate_spacing(request: DesignSpec, fixed_ids: set[str]) -> None:
    spacing = [r for r in request.requirements if isinstance(r, Spacing)]
    if len(spacing) > 1:
        msg = "only one fixed pair spacing requirement is supported"
        raise ValueError(msg)
    for rule in spacing:
        if {rule.upstream, rule.downstream} - fixed_ids:
            msg = f"{rule.id}: spacing requires both parts to be fixed"
            raise ValueError(msg)


def _validate_counts(request: DesignSpec, fixed_ids: set[str]) -> None:
    for rule in request.requirements:
        if isinstance(rule, Occurrences):
            try:
                selected = rule.select.indices(request.parts)
            except ValueError as error:
                msg = f"{rule.id}: {error}"
                raise ValueError(msg) from error
            if rule.min is not None and rule.min > len(selected):
                raise insufficient_parts(rule, selected)
            count_bounds(rule.min, rule.max, len(selected))
            mandatory = sum(request.parts[i].part_id in fixed_ids for i in selected)
            if rule.max is not None and mandatory > rule.max:
                msg = f"{rule.id}: fixed occurrences exceed the maximum count"
                raise ValueError(msg)
        elif isinstance(rule, GroupCoverage):
            try:
                PartSelector(groups=rule.groups).indices(request.parts)
            except ValueError as error:
                msg = f"{rule.id}: {error}"
                raise ValueError(msg) from error
