"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/formats/jaspar.py

Read JASPAR frequency matrices with preserved counts and exact motif IDs.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
import re

from dense_arrays.parts.motifs.models import BASES, Motif


def _matrix(lines: list[str]) -> tuple[tuple[float, ...], ...]:
    if len(lines) != len(BASES):
        msg = "JASPAR requires exactly four base rows per motif"
        raise ValueError(msg)
    labeled = bool(re.match(r"[ACGT]\s*\[", lines[0]))
    counts = {}
    for index, line in enumerate(lines):
        if labeled:
            match = re.fullmatch(r"([ACGT])\s*\[\s*(.*?)\s*\]", line)
            if match is None or match[1] in counts:
                msg = "JASPAR requires each labeled ACGT row exactly once"
                raise ValueError(msg)
            base, text = match[1], match[2]
        else:
            base, text = "ACGT"[index], line
        try:
            row = tuple(float(x) for x in text.split())
        except ValueError as err:
            msg = "JASPAR rows must contain numbers with consistent labels"
            raise ValueError(msg) from err
        if not row or any(not math.isfinite(x) or x < 0 for x in row):
            msg = "JASPAR counts must be finite, nonnegative and nonempty"
            raise ValueError(msg)
        counts[base] = row
    if len({len(row) for row in counts.values()}) != 1:
        msg = "JASPAR base rows must have equal width"
        raise ValueError(msg)
    return tuple(zip(*(counts[b] for b in "ACGT"), strict=True))


def _motif(header: str, lines: list[str]) -> Motif:
    names = header[1:].strip().split(maxsplit=1)
    if not names:
        msg = "JASPAR motif requires an ID"
        raise ValueError(msg)
    counts = _matrix(lines)
    try:
        totals = tuple(math.fsum(row) for row in counts)
    except OverflowError as err:
        msg = "JASPAR count total exceeds supported numeric range"
        raise ValueError(msg) from err
    if any(x <= 0 for x in totals):
        msg = "JASPAR positions require a positive count total"
        raise ValueError(msg)
    probabilities = tuple(
        tuple(x / total for x in row) for row, total in zip(counts, totals, strict=True)
    )
    return Motif(
        names[0],
        probabilities,
        (0.25,) * 4,
        None,
        "dense_arrays.import",
        {
            "import_policy": "jaspar_counts.v1",
            "background_origin": "uniform_default",
            "source_counts": counts,
            "column_totals": totals,
            **({"name": names[1]} if len(names) > 1 else {}),
        },
    )


def read_jaspar(payload: bytes) -> tuple[Motif, ...]:
    """Read bracketed or ordered ACGT rows, without adding pseudocounts."""
    motifs, rows, seen = [], [], set()
    header = None
    for raw in [*payload.decode("utf-8-sig").splitlines(), ">"]:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if header is not None:
                motif = _motif(header, rows)
                if motif.motif_id in seen:
                    msg = f"duplicate JASPAR motif ID: {motif.motif_id}"
                    raise ValueError(msg)
                seen.add(motif.motif_id)
                motifs.append(motif)
            header, rows = line, []
        elif header is None:
            msg = "JASPAR input requires a motif header before count rows"
            raise ValueError(msg)
        else:
            rows.append(line)
    if not motifs:
        msg = "JASPAR input contains no motifs"
        raise ValueError(msg)
    return tuple(motifs)
