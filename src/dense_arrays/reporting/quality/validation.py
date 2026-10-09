"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/quality/validation.py

Validate recorded quality aggregates without claiming source-artifact replay.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from collections.abc import Mapping

from dense_arrays._record_validation import (
    digest,
    integer,
    object_fields,
    required_text,
)
from dense_arrays.artifacts.provenance import Producer
from dense_arrays.artifacts.records import OUTCOMES
from dense_arrays.reporting.metrics import COMPOSITION_METRICS, validate_metric_value


def validate_quality(value: object, *, schema: str, policy: str) -> dict[str, object]:
    """Check the report envelope and all populations used by aggregate comparisons."""
    keys = {
        "schema",
        "policy",
        "run_id",
        "revision",
        "status",
        "population",
        "cell_ref",
        "attainment",
        "selection",
        "source_runs",
        "supply",
        "concentration",
        "composition",
        "padding_sides",
        "occupancy",
        "part_usage",
        "group_usage",
        "requirements",
        "cells",
        "search",
        "examined",
        "state_entries",
        "next_cursor",
        "usage_order",
        "usage_offset",
        "cost",
    }
    data = object_fields(value, keys, "quality report")
    if set(data) != keys or data["schema"] != schema or data["status"] != "exact":
        msg = "unsupported or incomplete quality report schema"
        raise ValueError(msg)
    required_text(data["policy"], field_name="quality.policy")
    required_text(data["population"], field_name="quality.population")
    selection = _required(
        data["selection"],
        {"sources", "query_id", "filter", "designs", "distinct_sequences"},
        "selection",
    )
    if not isinstance(selection, Mapping) or not isinstance(
        selection.get("sources"), list
    ):
        msg = "quality selection requires recorded sources"
        raise TypeError(msg)
    integer(selection["designs"], field_name="selection.designs", minimum=0)
    count = selection["designs"]
    integer(
        selection["distinct_sequences"],
        field_name="selection.distinct_sequences",
        minimum=0,
    )
    if selection["distinct_sequences"] > count or (
        count > 0 and selection["distinct_sequences"] == 0
    ):
        msg = "distinct sequences do not match the selected population"
        raise ValueError(msg)
    _origins(data["source_runs"], count)
    if data["policy"] == policy and set(data["composition"]) != set(
        COMPOSITION_METRICS
    ):
        msg = "composition fields do not match the declared metric policy"
        raise ValueError(msg)
    known_policy = data["policy"] == policy
    _distributions(data["composition"], count, known_policy=known_policy)
    _supply(data["supply"], data["concentration"])
    if (
        known_policy
        and sum(
            row["value"] * row["count"]
            for row in data["composition"]["placement_count"]["histogram"]
        )
        != data["concentration"]["occurrence_denominator"]
    ):
        msg = "placement_count does not match the occurrence denominator"
        raise ValueError(msg)
    _search(data["search"])
    _usage(data["part_usage"], count, data["concentration"]["occurrence_denominator"])
    _usage(
        data["group_usage"],
        count,
        data["concentration"]["occurrence_denominator"],
        groups=True,
    )
    return data


def _origins(value: object, selected: int) -> None:
    if not isinstance(value, list) or not value:
        msg = "quality requires source-run populations"
        raise TypeError(msg)
    seen, total = set(), 0
    for record in value:
        source = _required(
            record,
            {
                "run_id",
                "plan_id",
                "revision",
                "attainment",
                "selected_designs",
                "included_designs",
            },
            "source run",
        )
        required_text(source["run_id"], field_name="source.run_id")
        if "producer" in source:
            Producer.from_dict(source["producer"])
        digest(source["plan_id"], field_name="source.plan_id")
        integer(source["revision"], field_name="source.revision", minimum=0)
        if source["run_id"] in seen:
            msg = "quality report repeats a source run"
            raise ValueError(msg)
        seen.add(source["run_id"])
        counts = _required(
            source["attainment"], {"target", "accepted", "shortfall"}, "attainment"
        )
        for name in ("target", "accepted", "shortfall"):
            integer(counts[name], field_name=f"attainment.{name}", minimum=0)
        for name in ("included_designs", "selected_designs"):
            integer(source[name], field_name=name, minimum=0)
        if (
            counts["accepted"] + counts["shortfall"] != counts["target"]
            or not source["selected_designs"]
            <= source["included_designs"]
            <= counts["accepted"]
        ):
            msg = "quality source population does not reconcile"
            raise ValueError(msg)
        total += source["selected_designs"]
    if total != selected:
        msg = "quality selected population does not match source counts"
        raise ValueError(msg)


def _distributions(value: object, population: int, *, known_policy: bool) -> None:
    if not isinstance(value, Mapping) or not value:
        msg = "quality requires composition distributions"
        raise TypeError(msg)
    for name, record in value.items():
        required_text(name, field_name="composition metric")
        distribution = _required(
            record,
            {"count", "histogram", "min", "max", "mean", "denominator"},
            "distribution",
        )
        required_text(
            distribution["denominator"], field_name="distribution.denominator"
        )
        integer(distribution["count"], field_name="distribution.count", minimum=0)
        if distribution["count"] != population:
            msg = "quality metric count does not match selected population"
            raise ValueError(msg)
        bins = distribution["histogram"]
        values = _histogram_values(
            bins, name=name, population=population, known_policy=known_policy
        )
        if not population:
            valid = all(distribution[n] is None for n in ("min", "max", "mean"))
        else:
            for statistic in ("min", "max", "mean"):
                _number(distribution[statistic], statistic)
            mean = (
                math.fsum(item["value"] * item["count"] for item in bins) / population
            )
            valid = (
                distribution["min"] == values[0]
                and distribution["max"] == values[-1]
                and math.isclose(
                    distribution["mean"], mean, rel_tol=1e-12, abs_tol=1e-12
                )
            )
        if not valid:
            msg = "quality distribution does not match its recorded histogram"
            raise ValueError(msg)


def _histogram_values(
    bins: object, *, name: str, population: int, known_policy: bool
) -> list[int | float]:
    """Validate bin domains and population before checking derived statistics."""
    if not isinstance(bins, list):
        msg = "quality histogram must be an array"
        raise TypeError(msg)
    observed, values = 0, []
    for row in bins:
        item = _required(row, {"value", "count"}, "histogram bin")
        _number(item["value"], "histogram.value")
        if known_policy:
            validate_metric_value(name, item["value"])
        integer(item["count"], field_name="histogram.count", minimum=1)
        observed += item["count"]
        values.append(item["value"])
    if observed != population or values != sorted(set(values)):
        msg = "quality histogram does not reconcile with its population"
        raise ValueError(msg)
    return values


def _supply(supply: Mapping, concentration: Mapping) -> None:
    supply = _required(
        supply,
        {"eligible_parts", "eligible_groups", "unused_parts", "unused_groups"},
        "supply",
    )
    concentration = _required(
        concentration,
        {"occurrence_denominator", "highest_part_occurrence_share"},
        "concentration",
    )
    for name in ("eligible_parts", "eligible_groups", "unused_parts", "unused_groups"):
        integer(supply[name], field_name=name, minimum=0)
    if (
        supply["unused_parts"] > supply["eligible_parts"]
        or supply["unused_groups"] > supply["eligible_groups"]
    ):
        msg = "unused supply cannot exceed eligible supply"
        raise ValueError(msg)
    integer(
        concentration["occurrence_denominator"],
        field_name="occurrence_denominator",
        minimum=0,
    )
    share = concentration["highest_part_occurrence_share"]
    if share is not None:
        _number(share, "highest_part_occurrence_share")
    if (concentration["occurrence_denominator"] == 0) != (share is None) or (
        share is not None and not 0 <= share <= 1
    ):
        msg = "concentration does not match its occurrence denominator"
        raise ValueError(msg)


def _search(search: Mapping) -> None:
    search = _required(
        search, {"availability", "attempt_counts", "active_seconds"}, "search"
    )
    if search["availability"] not in {"complete", "partial", "not_included"}:
        msg = "unknown search availability"
        raise ValueError(msg)
    counts = search["attempt_counts"]
    if search["availability"] == "not_included":
        if counts is not None or search["active_seconds"] is not None:
            msg = "unavailable search history cannot supply attempt counts"
            raise ValueError(msg)
        return
    counts = object_fields(counts, {"started", *OUTCOMES}, "attempt_counts")
    if set(counts) != {"started", *OUTCOMES}:
        msg = "search history requires every attempt outcome"
        raise ValueError(msg)
    for name, count in counts.items():
        integer(count, field_name=name, minimum=0)
    if counts["started"] != sum(counts[name] for name in OUTCOMES):
        msg = "search counts do not reconcile"
        raise ValueError(msg)
    _number(search["active_seconds"], "active_seconds")


def _usage(
    value: object, designs: int, occurrences: int, *, groups: bool = False
) -> None:
    """Validate a displayed usage page without claiming it contains every row."""
    if not isinstance(value, list):
        msg = "quality usage must be an array"
        raise TypeError(msg)
    identity = "group" if groups else "part_id"
    keys = {
        identity,
        "occurrences",
        "designs",
        "occurrence_denominator",
        "design_denominator",
        "occurrence_share",
        "design_fraction",
    }
    for item in value:
        row = _required(item, keys, "usage row")
        required_text(row[identity], field_name=identity)
        for name in (
            "occurrences",
            "designs",
            "occurrence_denominator",
            "design_denominator",
        ):
            integer(row[name], field_name=f"usage.{name}", minimum=0)
        if (
            row["design_denominator"] != designs
            or row["occurrence_denominator"] != occurrences
            or row["designs"] > designs
            or row["occurrences"] > occurrences
        ):
            msg = "quality usage does not reconcile with its population"
            raise ValueError(msg)
        for numerator, denominator, fraction in (
            (row["occurrences"], occurrences, "occurrence_share"),
            (row["designs"], designs, "design_fraction"),
        ):
            expected = numerator / denominator if denominator else None
            if expected is not None:
                _number(row[fraction], fraction)
            if row[fraction] != expected:
                msg = f"quality {fraction} does not match its denominator"
                raise ValueError(msg)


def _number(value: object, name: str) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        msg = f"{name} must be a finite number"
        raise ValueError(msg)


def _required(value: object, keys: set[str], name: str) -> Mapping:
    """Require fields consumed by comparisons while preserving other recorded fields."""
    if not isinstance(value, Mapping):
        msg = f"{name} must be an object"
        raise TypeError(msg)
    if missing := keys - value.keys():
        msg = f"{name} is missing fields: {sorted(missing)}"
        raise ValueError(msg)
    return value


def report_entries(value: object) -> int:
    """Bound retained report structure without trusting a serialized work counter."""
    if isinstance(value, Mapping):
        return len(value) + sum(report_entries(v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return len(value) + sum(report_entries(v) for v in value)
    return 0
