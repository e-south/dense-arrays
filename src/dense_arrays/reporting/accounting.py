"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/accounting.py

Reconcile mutually exclusive attempt outcomes and overlapping evidence reasons.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.records import OUTCOMES, Attempt

if TYPE_CHECKING:
    from dense_arrays.artifacts.reading import ReadBudget


@dataclass
class AttemptTotals:
    """One exact accumulation over the declared attempt population."""

    counts: Counter = field(
        default_factory=lambda: Counter({"started": 0, **dict.fromkeys(OUTCOMES, 0)})
    )
    reasons: Counter = field(default_factory=Counter)
    proof_scopes: Counter = field(default_factory=Counter)
    budget: ReadBudget | None = field(default=None, repr=False)

    def observe(self, attempt: Attempt) -> None:
        """Count each outcome once; requirement failures may overlap other reasons."""
        self.counts["started"] += 1
        self.counts[attempt.outcome] += 1
        self._add(self.proof_scopes, attempt.evidence.get("proof_scope") or "unknown")
        code = attempt.evidence.get("code")
        if code:
            self._add(self.reasons, code)
        failed = sum(
            not result["passed"] for result in attempt.evidence.get("requirements", ())
        )
        if failed:
            self._add(self.reasons, "requirement_failed", failed)
        elif attempt.outcome != "accepted" and not code:
            self._add(
                self.reasons,
                "duplicate_sequence"
                if attempt.outcome == "duplicate"
                else attempt.outcome,
            )

    def _add(self, counter: Counter, key: str, count: int = 1) -> None:
        if key not in counter and self.budget is not None:
            self.budget.retain()
        counter[key] += count

    def merge(self, other: AttemptTotals) -> None:
        """Combine distinct source-run effort without treating filters as attempts."""
        self.counts.update(other.counts)
        for destination, source in (
            (self.reasons, other.reasons),
            (self.proof_scopes, other.proof_scopes),
        ):
            for key, count in source.items():
                self._add(destination, key, count)

    def reconcile(self, expected: object) -> None:
        """Reject an exact report whose source counters disagree."""
        if dict(self.counts) != dict(expected):
            msg = "attempts do not reconcile with committed counters"
            raise ArtifactIntegrityError(msg)

    def to_dict(self) -> dict[str, object]:
        """Keep exclusive outcomes distinct from multi-label reason totals."""
        return {
            "attempt_counts": dict(self.counts),
            "reason_counts": dict(self.reasons),
            "proof_scopes": dict(self.proof_scopes),
        }
