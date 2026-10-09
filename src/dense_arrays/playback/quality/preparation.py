"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/quality/preparation.py

Plot recipe-local preparation evidence without sampling or scoring.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from decimal import Decimal
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Mapping

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

_COLOR = "#28556A"
_MAX_RECIPES = 12
_MAX_BANDS = 24
_MAX_LABEL = 64


def preparation_figure(report: Mapping[str, object]) -> Figure:
    """Show exact saved counts and sequential MMR distances for each recipe."""
    from matplotlib.figure import Figure

    if report.get("schema") != "dense_arrays.pool_quality.v1":
        msg = "unsupported preparation quality report schema"
        raise ValueError(msg)
    recipes = report.get("recipes", [{"id": None, "accounting": report}])
    if len(recipes) > _MAX_RECIPES or any(
        len(recipe["accounting"].get("score_bands", {}).get("bands", ())) > _MAX_BANDS
        for recipe in recipes
    ):
        msg = (
            "preparation-quality PNG supports at most 12 recipes and 24 bands "
            "per recipe; export the full quality JSON for larger reports"
        )
        raise ValueError(msg)
    diversity = {entry["recipe_id"]: entry for entry in report.get("diversity", ())}
    figure = Figure(
        figsize=(13, 1.2 + 3.7 * len(recipes)), layout="constrained", facecolor="white"
    )
    rows = figure.subfigures(len(recipes), 1, squeeze=False)[:, 0]
    figure.suptitle(f"Prepared parts · {report['state']}", fontsize=16)
    for index, (recipe, row) in enumerate(zip(recipes, rows, strict=True)):
        accounting = recipe["accounting"]
        name = recipe["id"] or "Preparation"
        shown_name = name if len(name) <= _MAX_LABEL else name[: _MAX_LABEL - 1] + "…"
        row.suptitle(
            f"{shown_name} · {accounting['counts']['retained']} / "
            f"{accounting['requested_retention']} retained · "
            f"{accounting['stop_reason'].replace('_', ' ')}",
            fontsize=11,
        )
        axes = row.subplots(1, 3)
        for axis, label in zip(axes, ("yield", "diversity", "bands"), strict=True):
            axis.set_label(f"{label}_{index}")
            axis.spines[["top", "right"]].set_visible(False)
            axis.tick_params(labelsize=8)
        _yield(axes[0], accounting)
        _diversity(axes[1], diversity.get(recipe["id"]))
        _bands(axes[2], accounting.get("score_bands"))
    figure.supxlabel(
        "Each row uses its own recipe. MMR distances compare a choice with earlier "
        "selected cores; they are not all-pairs diversity.\n"
        "Score bands include boundary ties. Missing metrics are not zero.",
        fontsize=8,
    )
    return figure


def _yield(axis: Axes, report: Mapping[str, object]) -> None:
    from matplotlib.ticker import MaxNLocator

    counts = report["counts"]
    values = [
        counts[key] for key in ("processed", "eligible", "eligible_unique", "retained")
    ]
    axis.plot(values, range(4), "o", color=_COLOR)
    axis.set_yticks(range(4), ["Processed", "Eligible", "Representatives", "Retained"])
    axis.invert_yaxis()
    target = report["requested_retention"]
    axis.axvline(target, color="0.6", linestyle="--", linewidth=0.8)
    for index, value in enumerate(values):
        axis.annotate(
            str(value),
            (value, index),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            fontsize=8,
        )
    axis.set_xlim(-0.03 * max(1, *values, target), max(1, *values, target) * 1.25)
    axis.set_title("Candidate yield", loc="left", fontsize=10)
    target_evidence = report.get("mining_target")
    mining = (
        "\nMining target: "
        f"{target_evidence['eligible_unique']} representatives, "
        f"{target_evidence['minimum_candidates']} candidates "
        f"({'met' if target_evidence['met'] else 'unmet'})"
        if target_evidence is not None
        else ""
    )
    axis.set_xlabel(
        f"Candidates · dashed: retention target {target}\n"
        f"Effort cap: {report['candidate_budget']}{mining}",
        fontsize=8,
    )
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))


def _diversity(axis: Axes, report: Mapping[str, object] | None) -> None:
    from matplotlib.ticker import MaxNLocator

    axis.set_title("Selection distance", loc="left", fontsize=10)
    choices = [] if report is None else report["choices"][1:]
    if not choices:
        text = (
            "Distance evidence not recorded"
            if report is None
            else "No earlier-core comparisons\nFewer than two retained choices"
        )
        _unavailable(axis, text)
        return
    axis.plot(
        [c["rank"] for c in choices],
        [c["nearest_distance"] for c in choices],
        ".",
        color=_COLOR,
    )
    extent = max(1, *(c["nearest_distance"] for c in choices))
    axis.set_ylim(-0.02 * extent, 1.05 * extent)
    axis.set_xlabel("Selection rank · first choice excluded", fontsize=8)
    axis.set_ylabel("Weighted distance to nearest earlier core", fontsize=8)
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))


def _bands(axis: Axes, report: Mapping[str, object] | None) -> None:
    from matplotlib.ticker import MaxNLocator

    axis.set_title("Score bands", loc="left", fontsize=10)
    if report is None:
        _unavailable(axis, "Score bands not declared")
        return
    bands = report["bands"]
    locations = list(range(len(bands)))
    axis.barh(
        [n - 0.18 for n in locations],
        [b["count"] for b in bands],
        height=0.32,
        color="#C9D7DD",
        label="Eligible representatives",
    )
    axis.barh(
        [n + 0.18 for n in locations],
        [b["retained"] for b in bands],
        height=0.32,
        color=_COLOR,
        label="Retained subset",
    )
    lower = [0, *(b["upper_fraction"] for b in bands[:-1])]
    axis.set_yticks(
        locations,
        [
            f"{_percent(lo)} to {_percent(b['upper_fraction'])}"
            for lo, b in zip(lower, bands, strict=True)
        ],
    )
    axis.invert_yaxis()
    axis.set_xlabel(
        "Candidates · nominal upper-rank intervals\nActual sizes include score ties",
        fontsize=8,
    )
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))
    axis.legend(frameon=False, fontsize=7, loc="best")


def _unavailable(axis: Axes, text: str) -> None:
    axis.set_axis_off()
    axis.text(0, 0.5, text, transform=axis.transAxes, fontsize=9, color="0.35")


def _percent(value: float) -> str:
    """Keep declared tail boundaries distinct without binary rounding artifacts."""
    text = format(Decimal(str(value)) * 100, "f")
    return (text.rstrip("0").rstrip(".") if "." in text else text) + "%"
