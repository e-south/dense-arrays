"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/screening/sequence.py

Literal sequence checks over final DNA and declared placement intervals.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from collections.abc import Iterator

from dense_arrays.constraints import GC, Avoid
from dense_arrays.realized import RealizedArray
from dense_arrays.sequence import reverse_complement


def evaluate_screen(rule: Avoid | GC, realized: RealizedArray) -> dict[str, object]:
    """Report exact matches or GC counts, including zero-padding applicability."""
    if isinstance(rule, Avoid):
        violations = _violations(rule, realized.sequence, realized)
        return {"id": rule.id, "observed": violations, "passed": not violations}
    sequence = realized.sequence
    if rule.scope == "padding":
        assembly = realized.provenance["assembly"]
        start, length = assembly["packed_start"], assembly["packed_length"]
        sequence = sequence[:start] + sequence[start + length :]
    return _gc(rule, sequence)


def _gc(rule: GC, sequence: str) -> dict[str, object]:
    if not sequence:
        return {
            "id": rule.id,
            "observed": None,
            "passed": True,
            "status": "not_applicable",
        }
    count = sequence.count("G") + sequence.count("C")
    fraction = count / len(sequence)
    return {
        "id": rule.id,
        "observed": {"gc_bases": count, "bases": len(sequence), "fraction": fraction},
        "passed": rule.min <= fraction <= rule.max,
    }


def _violations(
    rule: Avoid, sequence: str, realized: RealizedArray | None = None
) -> list[dict[str, object]]:
    return sorted(
        _iter_violations(rule, sequence, realized),
        key=lambda m: (m["start"], m["end"], m["pattern"], m["strand"]),
    )


def _iter_violations(
    rule: Avoid, sequence: str, realized: RealizedArray | None = None
) -> Iterator[dict[str, object]]:
    """Yield overlapping matches under one strand and placement-exception policy."""
    exceptions = (
        []
        if realized is None
        else [p for p in realized.placements if p.feature_id in rule.except_placements]
    )
    for pattern in rule.patterns:
        orientations = [("forward", pattern)]
        reverse = reverse_complement(pattern)
        if rule.strands == "both" and reverse != pattern:
            orientations.append(("reverse", reverse))
        for strand, query in orientations:
            start = sequence.find(query)
            while start >= 0:
                end = start + len(query)
                if not any(p.start <= start and end <= p.end for p in exceptions):
                    yield {
                        "pattern": pattern,
                        "strand": strand,
                        "start": start,
                        "end": end,
                        "intersections": []
                        if realized is None
                        else _intersections(realized, start, end),
                    }
                start = sequence.find(query, start + 1)


def _intersections(
    realized: RealizedArray, start: int, end: int
) -> list[dict[str, object]]:
    """Retain full source intervals intersecting a zero-based half-open match."""
    intervals = [
        {
            "kind": "part",
            "part_id": p.feature_id,
            "placement_id": p.placement_id,
            "start": p.start,
            "end": p.end,
        }
        for p in realized.placements
    ]
    assembly = realized.provenance.get("assembly")
    if assembly is not None:
        packed_start = assembly["packed_start"]
        packed_end = packed_start + assembly["packed_length"]
        for side, a, b in (
            ("left", 0, packed_start),
            ("right", packed_end, len(realized.sequence)),
        ):
            if a < b:
                intervals.append(
                    {"kind": "padding", "side": side, "start": a, "end": b}
                )
    return [item for item in intervals if item["start"] < end and start < item["end"]]


def evaluate_sequence(rule: Avoid | GC, sequence: str) -> dict[str, object]:
    """Apply sequence-scoped screens without inventing placements or assembly."""
    _validate_sequence_scope(rule)
    if isinstance(rule, Avoid):
        violations = _violations(rule, sequence)
        return {"id": rule.id, "observed": violations, "passed": not violations}
    return _gc(rule, sequence)


def passes_sequence(rule: Avoid | GC, sequence: str) -> bool:
    """Check acceptance without collecting unused literal-match observations."""
    _validate_sequence_scope(rule)
    if isinstance(rule, Avoid):
        return next(_iter_violations(rule, sequence), None) is None
    return bool(_gc(rule, sequence)["passed"])


def _validate_sequence_scope(rule: Avoid | GC) -> None:
    """Keep detailed and boolean sequence checks on the same supported scope."""
    if isinstance(rule, Avoid):
        if rule.except_placements:
            msg = "sequence-only screening cannot resolve placement exceptions"
            raise ValueError(msg)
    elif rule.scope != "sequence":
        msg = "sequence-only GC screening requires sequence scope"
        raise ValueError(msg)
