"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/plans/requests.py

Typed editable requests with native file serialization.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from dense_arrays.parts import BoundParts, PreparationSet, PreparationSpec
from dense_arrays.planning.batches import BatchSchedule, CandidateBatch
from dense_arrays.planning.batches.bindings import batch_state_size, membership_size
from dense_arrays.planning.libraries import ParentLibrary
from dense_arrays.planning.matrices import MatrixSpec
from dense_arrays.planning.models import DesignSpec
from dense_arrays.planning.preparation import preparation_to_dict
from dense_arrays.planning.serialization import request_to_dict

if TYPE_CHECKING:
    from pathlib import Path


@dataclass(frozen=True, repr=False)
class RequestReport:
    """An editable request; serialization is a native input, not a report envelope."""

    request: DesignSpec | PreparationSpec | PreparationSet | MatrixSpec

    def __post_init__(self) -> None:
        """Reject reports that cannot be submitted to the matching planner."""
        if not isinstance(
            self.request, (DesignSpec, PreparationSpec, PreparationSet, MatrixSpec)
        ):
            msg = "RequestReport requires a design, matrix or preparation request"
            raise TypeError(msg)

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Encode source locations relative to the destination when supplied."""
        if isinstance(self.request, MatrixSpec):
            return self.request.to_dict(base=base)
        return (
            request_to_dict(self.request, base=base)
            if isinstance(self.request, DesignSpec)
            else preparation_to_dict(self.request, base=base)
        )

    @property
    def identities(self) -> int:
        """Count embedded parts, requirements and exclusions without opening sources."""
        request = self.request
        if isinstance(request, MatrixSpec):
            return (
                sum(
                    sum(batch_state_size(item) for item in b.batches)
                    if isinstance(b, BatchSchedule)
                    else batch_state_size(b)
                    if isinstance(b, CandidateBatch)
                    else 0
                    for b in request.batches.values()
                )
                + RequestReport(request.base).identities
                + sum(
                    source.identities
                    if isinstance(source, BoundParts)
                    else len(source)
                    if isinstance(source, tuple)
                    else 0
                    for source in request.sources.values()
                )
                + sum(
                    1 + len(v.parts) + len(v.requirements) + len(v.add_requirements)
                    for options in request.axes.values()
                    for v in options.values()
                )
                + sum(
                    len(scope.source.exclusions) for scope in request.exclude.values()
                )
            )
        if isinstance(request, PreparationSet):
            return len(request.recipes) + sum(
                RequestReport(item).identities for item in request.recipes.values()
            )
        if isinstance(request, PreparationSpec):
            return 0
        excluded = (
            len(request.exclude.source.exclusions)
            if request.exclude is not None
            and isinstance(request.exclude.source, ParentLibrary)
            else 0
        )
        return (
            (
                request.parts.identities
                if isinstance(request.parts, BoundParts)
                else len(request.parts)
                if isinstance(request.parts, tuple)
                else 0
            )
            + membership_size(request)
            + len(request.requirements)
            + excluded
        )

    def __repr__(self) -> str:
        """Summarize the request kind without traversing embedded parts."""
        return f"RequestReport({type(self.request).__name__})"
