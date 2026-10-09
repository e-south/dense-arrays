"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/diagnostics.py

Explain static requirement contradictions with counts and source references.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

from dense_arrays.diagnostics import Diagnostic

if TYPE_CHECKING:
    from dense_arrays.parts import BoundParts
    from dense_arrays.planning.requirements import Occurrences


class PlanningError(ValueError):
    """An invalid request with structured evidence and a suggested correction."""

    def __init__(self, message: str, diagnostic: Diagnostic) -> None:
        """Keep the human explanation and machine-readable observation together."""
        super().__init__(message)
        if not isinstance(diagnostic, Diagnostic):
            msg = "planning errors require a Diagnostic record"
            raise TypeError(msg)
        self.diagnostic = diagnostic

    def with_source(self, source: BoundParts) -> PlanningError:
        """Attach one-based data rows only for a directly bound part table."""
        if source.import_report.kind != "table" or not source.locations:
            return self
        if len(source.locations) != 1:
            return self
        references = self.diagnostic.evidence_refs
        referenced = set(references)
        rows = tuple(
            f"{source.locations[0].as_uri()}#row={index + 1}"
            for index in range(len(source.parts))
            if f"parts/{index}" in referenced
        )
        return PlanningError(
            str(self),
            replace(self.diagnostic, evidence_refs=(*references, *rows)),
        )


def insufficient_parts(rule: Occurrences, selected: tuple[int, ...]) -> PlanningError:
    """Explain a minimum exceeding the number of eligible supplied identities."""
    available = len(selected)
    message = (
        f"{rule.id}: requires at least {rule.min} occurrences, "
        f"but only {available} are available"
    )
    return PlanningError(
        message,
        Diagnostic(
            code="insufficient_parts",
            severity="error",
            stage="planning",
            requirement_id=rule.id,
            observed={"available": available},
            expected={"minimum": rule.min},
            evidence_refs=tuple(f"parts/{index}" for index in selected),
            proof_scope="declared_part_counts",
            next_action="Add matching parts or lower the occurrence minimum.",
        ),
    )
