"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/diagnostics.py

Explain persisted outcomes without solving, screening again or changing requests.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from contextlib import closing
from dataclasses import dataclass, field
from functools import cached_property
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    immutable_json_mapping,
    integer,
)
from dense_arrays.artifacts.errors import ArtifactIntegrityError, integrity_boundary
from dense_arrays.artifacts.reading import ReadBudget, ReadCost, ReadLimits
from dense_arrays.artifacts.run_plans import cell_plans, run_limits
from dense_arrays.artifacts.store import reader, stored_plan
from dense_arrays.diagnostics import Diagnostic
from dense_arrays.planning.serialization import requirement_to_dict
from dense_arrays.reporting.accounting import AttemptTotals
from dense_arrays.reporting.filters import AttemptFilter
from dense_arrays.reporting.readers import RecordView, read_records

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping
    from pathlib import Path

    from dense_arrays.artifacts.records import Attempt
    from dense_arrays.reporting.summary import RunSummary


@dataclass(frozen=True)
class _DiagnosticData:
    records: tuple[Diagnostic, ...]
    omitted: int
    counts: Mapping[str, int]
    reasons: Mapping[str, int]
    proof_scopes: Mapping[str, int]
    examined: int


@dataclass(frozen=True, repr=False)
class DiagnosticReport:
    """A lazy snapshot report whose cost is available before evidence is scanned."""

    path: Path
    summary: RunSummary
    limit: int = 20
    read_limits: ReadLimits = field(default_factory=ReadLimits)
    select: AttemptFilter | None = None

    def __post_init__(self) -> None:
        """Validate display/work bounds without opening source records."""
        integer(self.limit, field_name="limit", minimum=1)
        if not isinstance(self.read_limits, ReadLimits):
            msg = "diagnostics require ReadLimits"
            raise TypeError(msg)
        if self.select is not None:
            if not isinstance(self.select, AttemptFilter):
                msg = "diagnostics requires AttemptFilter"
                raise TypeError(msg)
            self.select.validate(self.summary)

    def __repr__(self) -> str:
        """Expose source and scope without triggering a scan."""
        return (
            f"DiagnosticReport({self.summary.run_id}, "
            f"revision={self.summary.revision}, "
            f"shortfall={self.shortfall}, limit={self.limit})"
        )

    @property
    def cost(self) -> ReadCost:
        """Count the stored plan row and all attempt records at this revision."""
        return ReadCost(
            self.summary.run_id,
            self.summary.revision,
            "scan",
            "diagnostics",
            self.summary.counts["started"] + 1,
            self.read_limits,
        )

    @property
    def shortfall(self) -> int:
        """Keep target attainment separate from the diagnostic display limit."""
        return self.summary.target - self.summary.accepted

    @cached_property
    def _data(self) -> _DiagnosticData:
        with integrity_boundary(self.path):
            return _collect(self)

    @property
    def diagnostics(self) -> tuple[Diagnostic, ...]:
        """Return the bounded display; totals cover the complete snapshot."""
        return self._data.records

    @property
    def omitted(self) -> int:
        """Number of additional diagnostics omitted from this display."""
        return self._data.omitted

    @property
    def attempt_counts(self) -> Mapping[str, int]:
        """Mutually exclusive outcome totals independently read from attempts."""
        return self._data.counts

    @property
    def reason_counts(self) -> Mapping[str, int]:
        """Multi-label reasons; their sum is not the number of failed attempts."""
        return self._data.reasons

    @property
    def proof_scopes(self) -> Mapping[str, int]:
        """Observed proof scopes; missing scope is explicitly unknown."""
        return self._data.proof_scopes

    def to_dict(self) -> dict[str, object]:
        """Compute one exact report or raise; a truncated scan is never exact."""
        data = self._data
        return {
            "schema": "dense_arrays.diagnostics.v1",
            "policy": "stored_outcomes.v1",
            "run_id": self.summary.run_id,
            "revision": self.summary.revision,
            "status": "exact",
            "population": "all_attempts_at_revision"
            if self.select is None
            else "matching_attempts_at_revision",
            "filter": None if self.select is None else self.select.to_dict(),
            "target": self.summary.target,
            "accepted": self.summary.accepted,
            "shortfall": self.shortfall,
            "attempt_counts": dict(data.counts),
            "reason_counts": dict(data.reasons),
            "proof_scopes": dict(data.proof_scopes),
            "diagnostics": [d.to_dict() for d in data.records],
            "omitted": data.omitted,
            "examined": data.examined,
            "cost": self.cost.to_dict(),
        }


def _run_diagnostic(
    report: DiagnosticReport, limits: Mapping[str, object]
) -> Diagnostic | None:
    summary = report.summary
    if not report.shortfall:
        return None
    reason = summary.termination_reason
    limited = reason in {"attempt_limit", "active_time_limit"}
    return Diagnostic(
        "limit_reached" if limited else reason or "run_in_progress",
        "error" if summary.state == "failed" else "warning",
        "generation",
        {
            "accepted": summary.accepted,
            "attempts": summary.counts["started"],
            "active_seconds": summary.active_seconds,
            "termination_reason": reason,
        },
        {"target": summary.target, "limits": dict(limits)},
        (f"{summary.run_id}/revision/{summary.revision}",),
        "Review persisted outcomes before making a new request with explicit limits.",
        proof_scope="observed_search" if limited else None,
    )


def _attempt_diagnostics(
    attempt: Attempt, run_id: str, rules: Mapping[str, dict[str, object]]
) -> Iterator[Diagnostic]:
    reference = f"{run_id}/{attempt.cell_id}/attempt/{attempt.attempt_id}"
    scope = attempt.evidence.get("proof_scope")
    failed = [r for r in attempt.evidence.get("requirements", ()) if not r["passed"]]
    for result in failed:
        rule = rules.get(result["id"])
        if rule is None:
            msg = "attempt evidence names a requirement absent from its plan"
            raise ArtifactIntegrityError(msg)
        observed = result["observed"]
        yield Diagnostic(
            "requirement_failed",
            "warning",
            "screening",
            {"matches": observed} if rule["kind"] == "avoid" else {"value": observed},
            rule,
            (reference,),
            "Review this requirement and the recorded final intervals.",
            requirement_id=result["id"],
            proof_scope="observed_candidate",
        )
    if failed or attempt.outcome == "accepted":
        return
    code = attempt.evidence.get("code") or (
        "duplicate_sequence" if attempt.outcome == "duplicate" else attempt.outcome
    )
    yield Diagnostic(
        code,
        "error" if attempt.outcome == "error" else "info",
        "generation",
        {
            "outcome": attempt.outcome,
            "solver_status": attempt.evidence.get("solver_status"),
            **(
                {"heuristic": dict(attempt.evidence["heuristic"])}
                if "heuristic" in attempt.evidence
                else {}
            ),
            "termination_reason": attempt.evidence.get("termination_reason"),
            "detail": attempt.evidence.get("detail"),
        },
        {"outcome": "accepted"},
        (reference,),
        "Inspect the offered batch and recorded outcome; "
        "no global feasibility claim follows.",
        proof_scope=scope,
    )


def _collect(report: DiagnosticReport) -> _DiagnosticData:
    budget = ReadBudget(report.read_limits)
    budget.examine()
    with reader(report.path) as connection:
        plan = stored_plan(connection, max_identities=report.read_limits.identities)
    if plan.plan_id != report.summary.plan_id:
        msg = "diagnostic plan identity does not match the snapshot"
        raise ArtifactIntegrityError(msg, artifact=report.path)
    rules = {
        cell: {r.id: requirement_to_dict(r) for r in child.request.requirements}
        for cell, child in cell_plans(plan).items()
    }
    limits = run_limits(plan)
    first = _run_diagnostic(
        report,
        {
            "attempts": limits.attempts,
            "active_seconds": limits.active_seconds,
            "solver_seconds": limits.solver_seconds,
        },
    )
    selected = [] if first is None else [first]
    total = len(selected)
    totals = AttemptTotals()
    view = RecordView(
        report.path,
        report.summary.revision,
        "attempts",
        None,
        run_id=report.summary.run_id,
        read_limits=report.read_limits,
        select=report.select,
    )
    with closing(read_records(view, budget)) as attempts:
        for attempt in attempts:
            totals.observe(attempt)
            for diagnostic in _attempt_diagnostics(
                attempt, report.summary.run_id, rules[attempt.cell_id]
            ):
                if len(selected) < report.limit:
                    selected.append(diagnostic)
                total += 1
    if report.select is None:
        totals.reconcile(report.summary.counts)
    return _DiagnosticData(
        tuple(selected),
        total - len(selected),
        immutable_json_mapping(totals.counts),
        immutable_json_mapping(totals.reasons),
        immutable_json_mapping(totals.proof_scopes),
        budget.examined,
    )
