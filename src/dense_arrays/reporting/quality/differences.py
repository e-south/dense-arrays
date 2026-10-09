"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/differences.py

Descriptive numeric differences with explicit denominators and availability.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer
from dense_arrays.reporting.metrics import COMPOSITION_METRICS

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping


@dataclass(frozen=True)
class MetricDifference:
    """An after-minus-before difference; unavailable values never become zero."""

    path: tuple[str, ...]
    before: int | float | None
    after: int | float | None
    before_denominator: int | None
    after_denominator: int | None
    status: str = "comparable"
    reason: str | None = None
    delta: int | float | None = field(init=False)

    def __post_init__(self) -> None:
        """Validate finite values and prevent arithmetic across incompatible metrics."""
        if (
            not isinstance(self.path, (list, tuple))
            or not self.path
            or any(not isinstance(part, str) or not part for part in self.path)
        ):
            msg = "metric path must contain nonempty strings"
            raise TypeError(msg)
        object.__setattr__(self, "path", tuple(self.path))
        for name in ("before", "after"):
            value = getattr(self, name)
            if value is not None and (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                msg = f"{name} metric must be a finite number or null"
                raise ValueError(msg)
        for name in ("before_denominator", "after_denominator"):
            if (value := getattr(self, name)) is not None:
                integer(value, field_name=name, minimum=0)
        if self.status not in {"comparable", "incomparable", "unavailable"} or (
            (self.status == "comparable") != (self.reason is None)
        ):
            msg = "metric status requires an incompatibility or availability reason"
            raise ValueError(msg)
        if self.status == "comparable" and (self.before is None or self.after is None):
            msg = "comparable metrics require two observed values"
            raise ValueError(msg)
        object.__setattr__(
            self,
            "delta",
            self.after - self.before if self.status == "comparable" else None,
        )

    def to_dict(self) -> dict[str, object]:
        """Serialize numeric observations separately from the derived difference."""
        return {**asdict(self), "path": list(self.path)}


def metric_differences(
    before: Mapping,
    after: Mapping,
    *,
    policy: str,
) -> tuple[MetricDifference, ...]:
    """Compare aggregate metrics; usage pages never constrain their populations."""
    reason = (
        "metric_policy_mismatch"
        if before["policy"] != after["policy"]
        else "unsupported_metric_policy"
        if before["policy"] != policy
        else None
    )
    return tuple(_differences(before, after, reason))


def _differences(
    before: Mapping, after: Mapping, reason: str | None
) -> Iterator[MetricDifference]:
    for name in ("designs", "distinct_sequences"):
        yield _difference(
            ("selection", name),
            before["selection"][name],
            after["selection"][name],
            reason=reason,
        )
    for name in ("eligible_parts", "eligible_groups", "unused_parts", "unused_groups"):
        yield _difference(
            ("supply", name),
            before["supply"][name],
            after["supply"][name],
            reason=reason,
        )
    names = set(before["composition"]) | set(after["composition"])
    ordered = (
        *(name for name in COMPOSITION_METRICS if name in names),
        *sorted(names - set(COMPOSITION_METRICS)),
    )
    for name in ordered:
        a, b = before["composition"].get(name), after["composition"].get(name)
        missing = "metric_not_reported" if a is None or b is None else None
        for statistic in ("min", "max", "mean"):
            yield _difference(
                ("composition", name, statistic),
                None if a is None else a[statistic],
                None if b is None else b[statistic],
                None if a is None else a["count"],
                None if b is None else b["count"],
                reason=reason or missing,
            )
    a, b = before["concentration"], after["concentration"]
    yield _difference(
        ("concentration", "highest_part_occurrence_share"),
        a["highest_part_occurrence_share"],
        b["highest_part_occurrence_share"],
        a["occurrence_denominator"],
        b["occurrence_denominator"],
        reason=reason,
    )
    yield from _search_differences(before["search"], after["search"], reason)


def _search_differences(
    a: Mapping, b: Mapping, reason: str | None
) -> Iterator[MetricDifference]:
    availability = None
    if "not_included" in {a["availability"], b["availability"]}:
        availability = "search_history_not_included"
    elif "partial" in {a["availability"], b["availability"]}:
        availability = "partial_search_history"
    counts_a, counts_b = a["attempt_counts"], b["attempt_counts"]
    denominator_a = None if counts_a is None else counts_a["started"]
    denominator_b = None if counts_b is None else counts_b["started"]
    from dense_arrays.artifacts.records import OUTCOMES  # noqa: PLC0415

    for name in ("started", *OUTCOMES):
        yield _difference(
            ("search", "attempt_counts", name),
            None if counts_a is None else counts_a[name],
            None if counts_b is None else counts_b[name],
            denominator_a,
            denominator_b,
            reason=reason,
            unavailable=availability,
        )
    yield _difference(
        ("search", "active_seconds"),
        a["active_seconds"],
        b["active_seconds"],
        denominator_a,
        denominator_b,
        reason=reason,
        unavailable=availability,
    )


def _difference(  # noqa: PLR0913 - paired observations with independent denominators
    path: tuple[str, ...],
    before: float | None,
    after: float | None,
    before_denominator: int | None = None,
    after_denominator: int | None = None,
    *,
    reason: str | None = None,
    unavailable: str | None = None,
) -> MetricDifference:
    if reason:
        status = "incomparable"
    elif unavailable or before is None or after is None:
        status, reason = "unavailable", unavailable or "empty_population"
    else:
        status = "comparable"
    return MetricDifference(
        path, before, after, before_denominator, after_denominator, status, reason
    )
