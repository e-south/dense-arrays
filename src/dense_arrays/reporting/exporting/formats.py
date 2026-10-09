"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/exporting/formats.py

Supported projections and formats, checked before population reads.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header


def validate_format(view: str, format_name: str) -> None:
    """Reject unsupported combinations before source population reads."""
    if format_name not in {
        "json",
        "csv",
        "tsv",
        "fasta",
        "bundle",
        "selection",
        "jsonl",
    }:
        msg = (
            "supported export formats: json, jsonl, csv, tsv, fasta, bundle, selection"
        )
        raise ValueError(msg)
    if view not in {
        "designs",
        "arrays",
        "sequences",
        "placements",
        "attempts",
        "batches",
        "parts",
        "candidates",
        "summary",
        "request",
        "plan",
        "plans",
        "quality",
        "diagnostics",
        "comparison",
        "selection",
    }:
        msg = f"unsupported export view {view!r}"
        raise ValueError(msg)
    if format_name == "fasta" and view != "sequences":
        msg = "FASTA requires view='sequences'"
        raise ValueError(msg)
    if format_name == "bundle" and view not in {"designs", "arrays"}:
        msg = "bundle requires view='designs' or view='arrays'"
        raise ValueError(msg)
    if format_name == "selection" and view not in {"designs", "selection"}:
        msg = "selection format requires designs or selection view"
        raise ValueError(msg)
    if format_name in {"csv", "tsv"} and view not in {"sequences", "placements"}:
        msg = "CSV/TSV require the scalar sequences or placements view"
        raise ValueError(msg)

    if format_name == "jsonl" and view != "arrays":
        msg = "JSONL requires view='arrays'"
        raise ValueError(msg)
