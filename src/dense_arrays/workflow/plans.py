"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/plans.py

Compose saved-plan inspection, editable requests and semantic comparisons.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import replace
from pathlib import Path

from dense_arrays.parts import PreparationSet, PreparationSpec
from dense_arrays.planning import (
    DesignSpec,
    GenerationPlan,
    MatrixPlan,
    MatrixSpec,
    PlanEvidence,
    PreparationPlan,
)
from dense_arrays.reporting import (
    PlanComparison,
    ReadLimitError,
    ReadLimits,
    RequestReport,
)
from dense_arrays.reporting.plans.reading import plan_identities


def _resolve_saved_plan(
    source: object, limits: ReadLimits, *, as_request: bool = False
) -> GenerationPlan | PlanEvidence | PreparationPlan | MatrixPlan | RequestReport:
    """Read saved plan files or native stores without compiling new requests."""
    from dense_arrays.reporting.plans import read_plan  # noqa: PLC0415
    from dense_arrays.workflow.inputs import read_source  # noqa: PLC0415

    if isinstance(source, (str, Path)):
        source = Path(source)
        if source.is_file():
            source = read_source(source, max_identities=limits.identities)
            if not isinstance(
                source,
                (
                    GenerationPlan,
                    PlanEvidence,
                    PreparationPlan,
                    MatrixPlan,
                    DesignSpec,
                    PreparationSpec,
                    PreparationSet,
                    MatrixSpec,
                ),
            ):
                msg = (
                    "plan inspection requires a resolved generation or preparation plan"
                )
                raise TypeError(msg)
    if isinstance(source, (DesignSpec, PreparationSpec, PreparationSet, MatrixSpec)):
        if not as_request:
            msg = "plan inspection requires a resolved generation or preparation plan"
            raise TypeError(msg)
        report = RequestReport(source)
        if report.identities > limits.identities:
            msg = "read_limits.identities cannot hold the request identities"
            raise ReadLimitError(msg)
        return report
    return read_plan(source, limits, as_request=as_request)


def inspect_plan(
    artifact: object, compare: object, limits: ReadLimits, *, view: str = "plan"
) -> (
    GenerationPlan
    | PlanEvidence
    | PreparationPlan
    | MatrixPlan
    | RequestReport
    | PlanComparison
):
    """Dispatch saved plan files or native runs without resolving new requests."""
    from dense_arrays.reporting.plans import (  # noqa: PLC0415
        check_plan_limits,
    )

    if not isinstance(limits, ReadLimits):
        msg = "plan inspection requires ReadLimits"
        raise TypeError(msg)
    comparison_records = 2
    if compare is not None and limits.records < comparison_records:
        msg = "read_limits.records cannot hold two plan records for comparison"
        raise ReadLimitError(msg)
    if view == "request":
        if compare is not None:
            msg = "request inspection does not support comparison; use view='plan'"
            raise ValueError(msg)
        return _resolve_saved_plan(artifact, limits, as_request=True)
    before = _resolve_saved_plan(artifact, limits)
    if compare is None:
        return before
    remaining = limits.identities - plan_identities(before)
    if remaining <= 0:
        msg = "read_limits.identities cannot retain the second comparison plan"
        raise ReadLimitError(msg)
    after = _resolve_saved_plan(
        compare, replace(limits, identities=remaining, records=limits.records - 1)
    )
    check_plan_limits(limits, before, after)
    return PlanComparison(before, after)
