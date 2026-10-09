"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/quality/library.py

Plot declared library-report values without computing acceptance or metrics.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from itertools import pairwise
from typing import TYPE_CHECKING

from dense_arrays.reporting.quality.models import QUALITY_POLICY, QUALITY_SCHEMA

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def quality_figure(report: Mapping[str, object]) -> Figure:
    """Show ranked part use, GC, density and mutually exclusive attempt outcomes."""
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    if report.get("schema") != QUALITY_SCHEMA or report.get("policy") != QUALITY_POLICY:
        msg = "unsupported quality report schema or metric policy"
        raise ValueError(msg)
    figure = Figure(figsize=(11, 8), layout="constrained", facecolor="white")
    usage, gc, density, outcomes = figure.subplots(2, 2).flat
    color = "#28556A"
    for axis, label in zip(
        (usage, gc, density, outcomes),
        ("part_usage", "gc_fraction", "density", "outcomes"),
        strict=True,
    ):
        axis.set_label(label)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=9)
    selection = report["selection"]
    figure.suptitle(
        f"Library quality · selected designs: {selection['designs']}\n"
        f"Distinct sequences: {selection['distinct_sequences']} · "
        f"Source runs: {len(report['source_runs'])}",
        fontsize=16,
        fontweight="normal",
    )
    rows = report["part_usage"][:12]
    usage.barh(
        range(len(rows)),
        [row["occurrences"] for row in rows],
        color=color,
    )
    usage.set_yticks(
        range(len(rows)),
        [
            row["part_id"]
            + (f" [{row['collection_id'][:8]}]" if "collection_id" in row else "")
            for row in rows
        ],
    )
    usage.invert_yaxis()
    usage.set_title(
        f"Part use · {len(rows)} shown · "
        f"eligible IDs: {report['supply']['eligible_parts']}",
        loc="left",
    )
    usage.set_xlabel(
        "Selected occurrences"
        + (" · [collection digest]" if rows and "collection_id" in rows[0] else "")
    )
    usage.xaxis.set_major_locator(MaxNLocator(integer=True))
    if not rows:
        usage.text(
            0.5, 0.5, "No eligible parts", transform=usage.transAxes, ha="center"
        )
    for axis, name, title, xlabel in (
        (gc, "gc_fraction", "GC composition", "GC fraction of final sequence"),
        (density, "density", "Packing density", "Covered final bases / final length"),
    ):
        histogram = report["composition"][name]["histogram"]
        values = [item["value"] for item in histogram]
        width = min(
            0.04,
            min((b - a for a, b in pairwise(values)), default=0.05) * 0.8,
        )
        _histogram(
            axis,
            values,
            [item["count"] for item in histogram],
            width=width,
            color=color,
        )
        axis.set(
            title=title, xlabel=xlabel, ylabel="Selected designs", xlim=(-0.05, 1.05)
        )
        axis.yaxis.set_major_locator(MaxNLocator(integer=True))
        if not histogram:
            axis.text(
                0.5, 0.5, "No selected designs", transform=axis.transAxes, ha="center"
            )
    outcomes.set_title("Source search outcomes")
    if report["search"]["availability"] == "not_included":
        outcomes.set_axis_off()
        outcomes.text(
            0.5,
            0.5,
            "Attempt records not included",
            transform=outcomes.transAxes,
            ha="center",
        )
    else:
        counts = {
            name: value
            for name, value in report["search"]["attempt_counts"].items()
            if name != "started" and value
        }
        outcomes.barh(list(counts), list(counts.values()), color=color)
        outcomes.invert_yaxis()
        outcomes.set_xlabel(
            "Available native attempts (mutually exclusive)"
            if report["search"]["availability"] != "complete"
            else "All source attempts (mutually exclusive)"
        )
        outcomes.xaxis.set_major_locator(MaxNLocator(integer=True))
    figure.supxlabel(
        f"{_source_caption(report)}\n"
        "Composition uses the selected designs. "
        "Density counts overlapping placements once.",
        fontsize=9,
    )
    return figure


def _histogram(
    axis: Axes,
    values: Sequence[float],
    counts: Sequence[int],
    *,
    width: float,
    color: str,
) -> None:
    """Draw exact bin rectangles with one artist, retaining every value and count."""
    from matplotlib.collections import PolyCollection

    vertices = [
        (
            (value - width / 2, 0),
            (value - width / 2, count),
            (value - width / 2 + width, count),
            (value - width / 2 + width, 0),
        )
        for value, count in zip(values, counts, strict=True)
    ]
    collection = PolyCollection(vertices, facecolors=color, edgecolors="none")
    collection.sticky_edges.y[:] = [0]
    axis.add_collection(collection)
    axis.autoscale_view()


def _source_caption(report: Mapping[str, object]) -> str:
    sources = report["source_runs"]
    if len(sources) == 1:
        source = sources[0]
        attainment = source["attainment"]
        return (
            f"Source: {attainment['accepted']} / {attainment['target']} accepted, "
            f"{attainment['shortfall']} shortfall ({source['state']}); "
            f"revision {source['revision']}"
        )
    caption_limit = 3
    shown = [
        f"{source['attainment']['accepted']} / {source['attainment']['target']} "
        f"({source['state']})"
        for source in sources[:caption_limit]
    ]
    remainder = (
        f"; {len(sources) - caption_limit} more in report"
        if len(sources) > caption_limit
        else ""
    )
    return "Source attainment in supplied order: " + "; ".join(shown) + remainder
