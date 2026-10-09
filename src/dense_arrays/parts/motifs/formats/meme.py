"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/formats/meme.py

Read minimal MEME DNA probabilities without rounding them to inferred counts.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
import re
from dataclasses import replace

from dense_arrays.parts.motifs.models import BASES, Motif, numeric_row


def _motif_start(line: str) -> bool:
    return re.match(r"MOTIF(?:\s|$)", line) is not None


def _statistics(text: str) -> dict[str, object]:
    fields = {}
    rest = text.removeprefix("letter-probability matrix:").strip()
    while rest:
        match = re.match(r"(alength|w|nsites|E|S)\s*=\s*(\S+)\s*", rest)
        if match is None or match[1] in fields:
            msg = "unsupported or repeated MEME matrix statistic"
            raise ValueError(msg)
        key, raw = match[1], match[2]
        value = float(raw)
        if not math.isfinite(value) or (key != "S" and value < 0):
            msg = "MEME matrix statistics require finite valid numbers"
            raise ValueError(msg)
        if key in {"w", "alength"}:
            if not re.fullmatch(r"[1-9]\d*", raw):
                msg = "MEME width and alphabet length must be positive integers"
                raise ValueError(msg)
            value = int(raw)
        if key == "nsites" and value <= 0:
            msg = "MEME nsites must be positive when supplied"
            raise ValueError(msg)
        fields[key] = value
        rest = rest[match.end() :]
    if fields.get("alength", len(BASES)) != len(BASES):
        msg = "MEME import supports the ACGT alphabet only"
        raise ValueError(msg)
    return fields


def _background(lines: list[str], index: int) -> tuple[tuple[float, ...], int]:
    values = {}
    while index < len(lines) and lines[index] and not _motif_start(lines[index]):
        tokens = lines[index].split()
        if len(tokens) % 2:
            msg = "MEME background requires ACGT frequency pairs"
            raise ValueError(msg)
        for base, value in zip(tokens[::2], tokens[1::2], strict=True):
            if base not in "ACGT" or len(base) != 1 or base in values:
                msg = "MEME background requires each ACGT base exactly once"
                raise ValueError(msg)
            values[base] = float(value)
        index += 1
        if len(values) == len(BASES):
            break
    if set(values) != set("ACGT"):
        msg = "incomplete MEME background frequencies"
        raise ValueError(msg)
    return numeric_row(
        tuple(values[b] for b in "ACGT"), probability=True, positive=True
    ), index


def _motif(
    lines: list[str], index: int, background: tuple, metadata: dict
) -> tuple[Motif, int]:
    names = lines[index].split()
    if len(names) not in {2, 3} or any("=" in name for name in names[1:]):
        msg = "MEME MOTIF requires an ID and optional alternate name"
        raise ValueError(msg)
    _, identifier, *alternate = names
    index += 1
    while index < len(lines) and not lines[index]:
        index += 1
    if index == len(lines) or not lines[index].startswith("letter-probability matrix:"):
        msg = "MEME motif requires a letter-probability matrix; use minimal MEME format"
        raise ValueError(msg)
    statistics = _statistics(lines[index])
    index += 1
    rows = []
    while (
        index < len(lines)
        and lines[index]
        and not (_motif_start(lines[index]) or lines[index].startswith("URL "))
    ):
        rows.append(tuple(float(x) for x in lines[index].split()))
        index += 1
        if len(rows) == statistics.get("w"):
            break
    if not rows or ("w" in statistics and len(rows) != statistics["w"]):
        msg = "MEME matrix rows do not match declared width"
        raise ValueError(msg)
    model = Motif(
        identifier,
        tuple(rows),
        background,
        None,
        "dense_arrays.import",
        {
            **metadata,
            **statistics,
            **({"alternate_name": alternate[0]} if alternate else {}),
        },
    )
    return model, index


def _header(lines: list[str]) -> tuple[int, tuple, dict]:
    index = next((i for i, line in enumerate(lines) if line), len(lines))
    if index == len(lines) or not re.fullmatch(r"MEME version \S+", lines[index]):
        msg = "MEME input requires a MEME version header"
        raise ValueError(msg)
    metadata = {
        "import_policy": "meme_probability.v1",
        "format_version": lines[index].split()[2],
        "background_origin": "uniform_default",
        "alphabet_origin": "default_acgt",
    }
    background, seen = (0.25,) * len(BASES), set()
    index += 1
    while index < len(lines) and not _motif_start(lines[index]):
        line = lines[index]
        if not line:
            index += 1
            continue
        if re.match(r"ALPHABET\s*=", line) and "alphabet" not in seen:
            if line.partition("=")[2].strip() != BASES:
                msg = "MEME import supports the ACGT alphabet only"
                raise ValueError(msg)
            seen.add("alphabet")
            metadata["alphabet_origin"] = "declared"
        elif line.startswith("strands:") and "strands" not in seen:
            value = " ".join(line.partition(":")[2].split())
            if value not in {"+", "+ -"}:
                msg = "MEME strands must be + or + -"
                raise ValueError(msg)
            metadata["source_strands"] = value
            seen.add("strands")
        elif (
            line.startswith("Background letter frequencies")
            and "background" not in seen
        ):
            background, index = _background(lines, index + 1)
            seen.add("background")
            metadata["background_origin"] = "declared"
            continue
        else:
            msg = f"unsupported or repeated MEME header at line {index + 1}"
            raise ValueError(msg)
        index += 1
    return index, background, metadata


def read_meme(payload: bytes) -> tuple[Motif, ...]:
    """Validate every record; retain source statistics separately from probabilities."""
    lines = [
        "" if line.lstrip().startswith("#") else line.strip()
        for line in payload.decode("utf-8-sig").splitlines()
    ]
    index, background, metadata = _header(lines)
    motifs, seen = [], set()
    while index < len(lines):
        line = lines[index]
        if not line:
            index += 1
            continue
        if _motif_start(line):
            motif, index = _motif(lines, index, background, metadata)
            if motif.motif_id in seen:
                msg = f"duplicate MEME motif ID: {motif.motif_id}"
                raise ValueError(msg)
            seen.add(motif.motif_id)
            motifs.append(motif)
            continue
        if line.startswith("URL ") and motifs and "url" not in motifs[-1].metadata:
            motifs[-1] = replace(
                motifs[-1], metadata={**motifs[-1].metadata, "url": line[4:].strip()}
            )
        else:
            msg = f"unsupported minimal MEME content at line {index + 1}: {line[:80]}"
            raise ValueError(msg)
        index += 1
    if not motifs:
        msg = "MEME input contains no motifs"
        raise ValueError(msg)
    return tuple(motifs)
