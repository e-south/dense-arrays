"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/retention.py

Reconcile MMR admission with fixed or target-relative choice-pool limits.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, object_fields
from dense_arrays.parts.retention.pool import PoolSize

if TYPE_CHECKING:
    from collections.abc import Mapping

    from dense_arrays.parts.candidates import Candidate


def sizing_report(
    policy: PoolSize, *, target: int, available: int, admitted: int
) -> dict[str, object]:
    """Keep desired choice supply, hard effort cap and observed alternatives apart."""
    requested = policy.requested(target)
    limit = policy.resolve(target)
    return {
        **policy.to_dict(),
        "requested": requested,
        "limit": limit,
        "available": available,
        "capped": requested > limit,
        "shortfall": max(0, limit - admitted),
        "has_choice": target > 0 and admitted > target,
    }


def retention_summary(
    candidates: tuple[Candidate, ...], *, target: int, sizing: PoolSize | None
) -> dict[str, object]:
    """Count recorded pool admissions without rerunning selection."""
    counts = {
        name: sum(
            c.selection is not None and c.selection.pool_status == status
            for c in candidates
        )
        for name, status in (
            ("pool_size", "included"),
            ("below_score", "below_score"),
            ("beyond_limit", "beyond_limit"),
        )
    }
    return {
        "policy": "mmr",
        **counts,
        **(
            {
                "sizing": sizing_report(
                    sizing,
                    target=target,
                    available=counts["pool_size"] + counts["beyond_limit"],
                    admitted=counts["pool_size"],
                )
            }
            if sizing is not None
            else {}
        ),
    }


def validate_retention(
    value: object, *, counts: Mapping[str, int], target: int
) -> dict[str, object]:
    """Reconcile admission counts and any declared sizing evidence."""
    fields = {"policy", "pool_size", "below_score", "beyond_limit"}
    selection = object_fields(value, fields | {"sizing"}, "pool retention")
    if not fields <= selection.keys() or selection["policy"] != "mmr":
        msg = "unsupported or incomplete pool retention summary"
        raise ValueError(msg)
    for name in ("pool_size", "below_score", "beyond_limit"):
        integer(selection[name], field_name=f"retention.{name}", minimum=0)
    if (
        selection["pool_size"] + selection["below_score"] + selection["beyond_limit"]
        != counts["eligible_unique"]
        or counts["retained"] > selection["pool_size"]
    ):
        msg = "retention pool and eligible-unique counts disagree"
        raise ValueError(msg)
    if "sizing" in selection:
        _validate_sizing(selection, target)
    return selection


def _validate_sizing(selection: dict, target: int) -> None:
    declaration = {"policy", "per_retained", "maximum"}
    fields = declaration | {
        "requested",
        "limit",
        "available",
        "capped",
        "shortfall",
        "has_choice",
    }
    sizing = object_fields(selection["sizing"], fields, "retention sizing")
    if set(sizing) != fields:
        msg = "incomplete retention sizing"
        raise ValueError(msg)
    for name in ("requested", "limit", "available", "shortfall"):
        integer(sizing[name], field_name=f"retention sizing {name}", minimum=0)
    for name in ("capped", "has_choice"):
        if not isinstance(sizing[name], bool):
            msg = f"retention sizing {name} must be boolean"
            raise TypeError(msg)
    policy = PoolSize.from_dict({k: sizing[k] for k in declaration})
    expected = sizing_report(
        policy,
        target=target,
        available=selection["pool_size"] + selection["beyond_limit"],
        admitted=selection["pool_size"],
    )
    if sizing != expected or selection["pool_size"] != min(
        expected["available"], expected["limit"]
    ):
        msg = "retention sizing disagrees with declared target, cap or admission"
        raise ValueError(msg)
