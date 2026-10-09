"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/preparation/bands.py

Reconciled score-band summaries from persisted representative decisions.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from math import ceil
from statistics import median
from typing import TYPE_CHECKING

from dense_arrays._record_validation import digest, integer, object_fields
from dense_arrays.parts.retention.bands import BAND_POLICY, ScoreBands, score_cutoffs
from dense_arrays.parts.scoring.configuration import finite

if TYPE_CHECKING:
    from dense_arrays.parts.candidates import Candidate


def summarize_bands(
    candidates: tuple[Candidate, ...], policy: ScoreBands, *, scoring_id: str
) -> dict[str, object]:
    """Report actual band sizes and retained membership under the declared metric."""
    cutoffs = (*score_cutoffs(candidates, policy), None)
    groups: list[list[Candidate]] = [[] for _ in cutoffs]
    for candidate in candidates:
        if candidate.representative == candidate.index:
            if candidate.score_band is None or candidate.score_band > len(groups):
                msg = "score band is missing or outside its declared range"
                raise ValueError(msg)
            groups[candidate.score_band - 1].append(candidate)
    bands = []
    for index, (group, fraction, cutoff) in enumerate(
        zip(groups, (*policy.upper_fractions, 1.0), cutoffs, strict=True), 1
    ):
        scores = [c.score.raw for c in group]
        bands.append(
            {
                "band": index,
                "upper_fraction": fraction,
                "cutoff": cutoff,
                "count": len(group),
                "retained": sum(c.retained for c in group),
                "scores": {
                    "min": min(scores),
                    "median": median(scores),
                    "max": max(scores),
                }
                if scores
                else None,
            }
        )
    return {
        "schema": "dense_arrays.score_bands.v1",
        "policy": BAND_POLICY,
        "metric": "best_hit.raw",
        "units": "fimo_log2_odds",
        "population": "eligible_unique",
        "scoring_id": digest(scoring_id, field_name="score bands scoring_id"),
        "total": sum(len(g) for g in groups),
        "bands": bands,
    }


def band_report_size(value: object) -> int:
    """Bound nested recorded report entries before typed validation and copying."""
    if value is None:
        return 0
    if not isinstance(value, Mapping) or not isinstance(
        value.get("bands"), (list, tuple)
    ):
        msg = "score band report requires a mapping and ordered bands"
        raise TypeError(msg)
    size = len(value)
    for band in value["bands"]:
        if not isinstance(band, Mapping):
            msg = "score band entries must be objects"
            raise TypeError(msg)
        size += len(band)
        scores = band.get("scores")
        if scores is not None:
            if not isinstance(scores, Mapping):
                msg = "score band statistics must be objects"
                raise TypeError(msg)
            size += len(scores)
    return size


def accounting_bands_size(value: object) -> int:
    """Count score summaries in a single recipe or a complete preparation set."""
    if value is None:
        return 0
    if not isinstance(value, Mapping):
        msg = "preparation accounting must be a mapping"
        raise TypeError(msg)
    return band_report_size(value.get("score_bands")) + sum(
        band_report_size(recipe["accounting"].get("score_bands"))
        for recipe in value.get("recipes", ())
    )


def validate_band_report(value: object, *, total: int, retained: int) -> None:
    """Check the declared population, rank boundaries and complete disjoint counts."""
    fields = {
        "schema",
        "policy",
        "metric",
        "units",
        "population",
        "scoring_id",
        "total",
        "bands",
    }
    data = object_fields(value, fields, "score band report")
    expected = {
        "schema": "dense_arrays.score_bands.v1",
        "policy": BAND_POLICY,
        "metric": "best_hit.raw",
        "units": "fimo_log2_odds",
        "population": "eligible_unique",
    }
    if set(data) != fields or any(data[k] != v for k, v in expected.items()):
        msg = "unsupported or incomplete score band report"
        raise ValueError(msg)
    digest(data["scoring_id"], field_name="score bands scoring_id")
    integer(data["total"], field_name="score band total", minimum=0)
    bands = data["bands"]
    if not isinstance(bands, (list, tuple)) or len(bands) <= 1:
        msg = "score band report requires boundaries and a remainder"
        raise ValueError(msg)
    values = [
        _validate_band(b, i, total=total, last=i == len(bands))
        for i, b in enumerate(bands, 1)
    ]
    ScoreBands(tuple(b["upper_fraction"] for b in values[:-1]))
    if (
        values[-1]["upper_fraction"] != 1
        or data["total"] != total
        or sum(b["count"] for b in values) != total
        or sum(b["retained"] for b in values) != retained
    ):
        msg = "score band populations do not reconcile"
        raise ValueError(msg)
    previous = None
    for band in values:
        if band["scores"] is not None:
            if previous is not None and band["scores"]["max"] >= previous:
                msg = "score bands split ties or reverse score order"
                raise ValueError(msg)
            previous = band["scores"]["min"]
    _validate_boundaries(values, total)


def _validate_boundaries(bands: list[dict], total: int) -> None:
    cumulative = 0
    cutoff = None
    for band in bands[:-1]:
        previous_count = cumulative
        cumulative += band["count"]
        rank = ceil(Fraction(str(band["upper_fraction"])) * total)
        if rank > cumulative or (band["count"] > 0 and rank <= previous_count):
            msg = "score band size disagrees with its cumulative rank boundary"
            raise ValueError(msg)
        expected = band["scores"]["min"] if band["scores"] is not None else cutoff
        if band["cutoff"] != expected:
            msg = "score band cutoff disagrees with its observed boundary score"
            raise ValueError(msg)
        cutoff = expected


def _validate_band(
    value: object, index: int, *, total: int, last: bool
) -> dict[str, object]:
    fields = {"band", "upper_fraction", "cutoff", "count", "retained", "scores"}
    band = object_fields(value, fields, "score band")
    if set(band) != fields:
        msg = "incomplete score band"
        raise ValueError(msg)
    for name in ("band", "count", "retained"):
        integer(band[name], field_name=f"score band {name}", minimum=0)
    if band["band"] != index or band["retained"] > band["count"]:
        msg = "score band order or retained count is invalid"
        raise ValueError(msg)
    finite(band["upper_fraction"], "score band upper fraction")
    cutoff = band["cutoff"]
    if (cutoff is None) != (last or total == 0):
        msg = "score band cutoff must match its population and remainder"
        raise ValueError(msg)
    if cutoff is not None:
        finite(cutoff, "score band cutoff")
    scores = band["scores"]
    if (scores is None) != (band["count"] == 0):
        msg = "score band statistics must match its observed count"
        raise ValueError(msg)
    _validate_statistics(scores, cutoff)
    return band


def _validate_statistics(scores: object, cutoff: float | None) -> None:
    if scores is not None:
        scores = object_fields(
            scores, {"min", "median", "max"}, "score band statistics"
        )
        if set(scores) != {"min", "median", "max"}:
            msg = "incomplete score band statistics"
            raise ValueError(msg)
        for v in scores.values():
            finite(v, "score band statistic")
        if not scores["min"] <= scores["median"] <= scores["max"] or (
            cutoff is not None and scores["min"] < cutoff
        ):
            msg = "score band statistics disagree with score order or cutoff"
            raise ValueError(msg)
