"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/scoring/fimo.py

Execute bounded FIMO scans and validate every reported hit against its source.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import csv
import io
import math
import tempfile
import time
from dataclasses import replace
from pathlib import Path

from dense_arrays.parts.motifs.models import BASES
from dense_arrays.problem import motif_library
from dense_arrays.sequence import reverse_complement

from .binding import FimoBinding
from .configuration import finite
from .process import ScoringError, invoke, remaining_seconds
from .records import FimoHit, FimoResult

HEADER = (
    "motif_id",
    "motif_alt_id",
    "sequence_name",
    "start",
    "stop",
    "strand",
    "score",
    "p-value",
    "q-value",
    "matched_sequence",
)


def _write_motif(binding: FimoBinding, path: Path) -> None:
    background = " ".join(
        f"{base} {value:.17g}"
        for base, value in zip(BASES, binding.effective_background, strict=True)
    )
    rows = [" ".join(f"{x:.17g}" for x in row) for row in binding.motif.probabilities]
    path.write_text(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
        f"Background letter frequencies:\n{background}\n\nMOTIF motif\n"
        f"letter-probability matrix: alength= 4 w= {binding.motif.width} nsites= 20\n"
        + "\n".join(rows)
        + "\n",
    )


def _command(
    binding: FimoBinding, directory: Path, fasta: Path, threshold: float
) -> list[str]:
    command = [
        str(binding.executable),
        "--text",
        "--verbosity",
        "1",
        "--thresh",
        str(threshold),
        "--bfile",
        "--motif--",
        "--motif-pseudo",
        str(binding.settings.pseudocount),
    ]
    if binding.settings.strands == "single":
        command.append("--norc")
    return [*command, str(directory / "motif.meme"), str(fasta)]


def _write_candidates(sequences: tuple[str, ...], path: Path) -> None:
    with path.open("w") as handle:
        for index, sequence in enumerate(sequences):
            handle.write(f">candidate_{index}\n{sequence}\n")


def _hit(
    row: dict[str, str], sequence: str, binding: FimoBinding, maximum: float | None
) -> FimoHit:
    if row["motif_id"] != "motif" or row["motif_alt_id"] or row["q-value"]:
        msg = "unexpected motif or q-value in text-mode evidence"
        raise ValueError(msg)
    start, end = int(row["start"]) - 1, int(row["stop"])
    if not 0 <= start < end <= len(sequence) or end - start != binding.motif.width:
        msg = "FIMO coordinates disagree with candidate or motif width"
        raise ValueError(msg)
    allowed = {"+"} if binding.settings.strands == "single" else {"+", "-"}
    if row["strand"] not in allowed:
        msg = "unexpected FIMO strand"
        raise ValueError(msg)
    core = sequence[start:end]
    if row["strand"] == "-":
        core = reverse_complement(core)
    if row["matched_sequence"] != core:
        msg = "FIMO matched sequence disagrees with oriented source interval"
        raise ValueError(msg)
    raw = finite(float(row["score"]), "FIMO score")
    return FimoHit(
        start,
        end,
        "forward" if row["strand"] == "+" else "reverse",
        core,
        raw,
        float(row["p-value"]),
        raw if maximum is None else maximum,
    )


def _parse(
    payload: bytes,
    sequences: tuple[str, ...],
    binding: FimoBinding,
    maximum: float | None,
) -> tuple[tuple[FimoHit | None, ...], int]:
    """Validate all rows while retaining only one deterministic hit per candidate."""
    try:
        return _parse_rows(payload, sequences, binding, maximum)
    except (KeyError, TypeError, ValueError, csv.Error) as err:
        reason = "malformed"
        raise ScoringError(reason, str(err)) from err


def _parse_rows(
    payload: bytes,
    sequences: tuple[str, ...],
    binding: FimoBinding,
    maximum: float | None,
) -> tuple[tuple[FimoHit | None, ...], int]:
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8")), delimiter="\t")
    if reader.fieldnames != list(HEADER):
        msg = "unsupported or missing FIMO TSV header"
        raise ValueError(msg)
    lookup = {f"candidate_{index}": index for index in range(len(sequences))}
    hits: list[FimoHit | None] = [None] * len(sequences)
    strands = 1 if binding.settings.strands == "single" else 2
    seen = [
        bytearray(((len(sequence) - binding.motif.width + 1) * strands + 7) // 8)
        for sequence in sequences
    ]
    total = 0
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            msg = "incomplete or extra FIMO TSV fields"
            raise ValueError(msg)
        index = lookup[row["sequence_name"]]
        hit = _hit(row, sequences[index], binding, maximum)
        position = hit.start * strands + (hit.strand == "reverse")
        offset, bit = divmod(position, 8)
        if seen[index][offset] & (1 << bit):
            msg = "duplicate FIMO hit for an oriented candidate window"
            raise ValueError(msg)
        seen[index][offset] |= 1 << bit
        previous = hits[index]
        if previous is None or _rank(hit) < _rank(previous):
            hits[index] = hit
        total += 1
    return tuple(hits), total


def _rank(hit: FimoHit) -> tuple[float, int, bool]:
    return -hit.raw, hit.start, hit.strand != "forward"


def _maximum_core(binding: FimoBinding) -> str:
    """Maximize each additive position before the backend's monotone quantization."""
    return "".join(
        BASES[
            max(
                range(len(BASES)),
                key=lambda index: (
                    math.log(row[index]) - math.log(binding.effective_background[index])
                    if row[index] > 0
                    else -math.inf
                ),
            )
        ]
        for row in binding.motif.probabilities
    )


def scan_fimo(
    binding: FimoBinding, sequences: tuple[str, ...] | list[str]
) -> FimoResult:
    """Score a bounded batch, calibrating its maximum with the same FIMO model.

    One maximizing core is scanned at threshold 1 to obtain the backend's own
    rounded maximum. Candidate thresholds remain inside FIMO, so rounded TSV
    p-values are never used to reclassify a borderline hit.
    """
    if not isinstance(binding, FimoBinding):
        msg = "scoring requires FimoBinding"
        raise TypeError(msg)
    if not isinstance(sequences, (tuple, list)) or not sequences:
        msg = "scoring requires a nonempty ordered candidate batch"
        raise ValueError(msg)
    sequences = tuple(sequences)
    motif_library(sequences)
    if any(len(sequence) < binding.motif.width for sequence in sequences):
        msg = "candidate length must be at least the motif width"
        raise ValueError(msg)
    strands = 1 if binding.settings.strands == "single" else 2
    windows = (
        sum(len(sequence) - binding.motif.width + 1 for sequence in sequences) * strands
    )
    if windows + strands > binding.settings.limits.windows:
        msg = "scoring window limit exceeded, including maximum calibration"
        raise ValueError(msg)
    binding.verify()
    deadline = time.monotonic() + binding.settings.limits.seconds
    with tempfile.TemporaryDirectory(prefix="dense-arrays-fimo-") as temporary:
        directory = Path(temporary)
        _write_motif(binding, directory / "motif.meme")
        reference = (_maximum_core(binding),)
        _write_candidates(reference, directory / "reference.fa")
        limits = binding.settings.limits
        calibration = invoke(
            _command(binding, directory, directory / "reference.fa", 1),
            cwd=directory,
            limits=replace(limits, seconds=remaining_seconds(deadline)),
        )
        reference_hits, _ = _parse(calibration.stdout, reference, binding, None)
        if reference_hits[0] is None:
            reason = "malformed"
            raise ScoringError(reason, "FIMO omitted the maximum calibration hit")
        maximum = reference_hits[0].raw
        remaining_seconds(deadline)
        available = limits.output_bytes - calibration.total_bytes
        if available <= 0:
            reason = "output_limit"
            raise ScoringError(reason, "no output budget remains after calibration")
        _write_candidates(sequences, directory / "candidates.fa")
        payload = invoke(
            _command(
                binding,
                directory,
                directory / "candidates.fa",
                binding.settings.hit_pvalue_max,
            ),
            cwd=directory,
            limits=replace(
                limits, seconds=remaining_seconds(deadline), output_bytes=available
            ),
        )
        hits, reported = _parse(payload.stdout, sequences, binding, maximum)
    binding.verify()
    remaining_seconds(deadline)
    return FimoResult(binding.binding_id, hits, windows, strands, reported, maximum)
