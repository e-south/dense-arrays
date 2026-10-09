"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/serialization.py

Stream portable supplied-array inputs with an explicit end-of-stream checksum.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
import json
from typing import TYPE_CHECKING, BinaryIO, TextIO

from dense_arrays._record_validation import canonical_json, mutable_json, object_fields
from dense_arrays.artifacts.reading import ReadBudget, ReadLimits
from dense_arrays.parts.serialization import part_from_dict, part_to_dict
from dense_arrays.playback.serialization import (
    realized_array_from_dict,
    realized_array_to_dict,
)

from .geometry import validate_array
from .models import ArrayCollection

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

    from dense_arrays.realized import RealizedArray

TRANSPORT_SCHEMA = "dense_arrays.array_collection_input.v1"
END_SCHEMA = "dense_arrays.array_collection_end.v1"
MAX_LINE_BYTES = 16 * 1024 * 1024


def _pairs(items: list[tuple[str, object]]) -> dict[str, object]:
    result = {}
    for key, value in items:
        if key in result:
            msg = f"duplicate JSON key {key!r}"
            raise ValueError(msg)
        result[key] = value
    return result


def _line(stream: BinaryIO) -> tuple[bytes, dict[str, object]]:
    line = stream.readline(MAX_LINE_BYTES + 1)
    if len(line) > MAX_LINE_BYTES:
        msg = "array input record exceeds the 16 MiB line limit"
        raise ValueError(msg)
    if not line:
        msg = "array input is truncated; missing completion record"
        raise ValueError(msg)
    value = json.loads(line, object_pairs_hook=_pairs)
    if not isinstance(value, dict) or canonical_json(value).encode() + b"\n" != line:
        msg = "array input requires canonical JSON objects, one per line"
        raise ValueError(msg)
    return line, value


def read_input(path: Path, limits: ReadLimits) -> ArrayCollection:
    """Bind a header; verify the stream while publication consumes it."""
    with path.open("rb") as stream:
        header, value = _line(stream)
    value = object_fields(value, {"schema", "parts", "provenance"}, "array input")
    if (
        set(value) != {"schema", "parts", "provenance"}
        or value["schema"] != TRANSPORT_SCHEMA
    ):
        msg = "unsupported array input schema"
        raise ValueError(msg)
    if not isinstance(value["parts"], list):
        msg = "array input parts must be a list"
        raise TypeError(msg)
    ReadBudget(limits).retain(len(value["parts"]))
    parts = tuple(part_from_dict(part) for part in value["parts"])
    return ArrayCollection(parts, _arrays(path, header), value["provenance"])


def _arrays(path: Path, header: bytes) -> Iterator[RealizedArray]:
    with path.open("rb") as stream:
        actual, _ = _line(stream)
        if actual != header:
            msg = "array input changed after binding"
            raise ValueError(msg)
        checksum = hashlib.sha256(header)
        count = 0
        while True:
            raw, value = _line(stream)
            if value.get("schema") == END_SCHEMA:
                if canonical_json(value) != canonical_json(
                    {
                        "schema": END_SCHEMA,
                        "arrays": count,
                        "sha256": checksum.hexdigest(),
                    }
                ) or stream.read(1):
                    msg = (
                        "array input completion count, checksum "
                        "or trailing content differs"
                    )
                    raise ValueError(msg)
                return
            array = realized_array_from_dict(value)
            if realized_array_to_dict(array) != value:
                msg = "array input contains noncanonical realized geometry"
                raise ValueError(msg)
            checksum.update(raw)
            count += 1
            yield array


def write_input(source: ArrayCollection, stream: TextIO, limits: ReadLimits) -> int:
    """Write a bounded stream suitable for the matching Python and CLI importer."""
    budget = ReadBudget(limits)
    budget.retain(len(source.parts))
    parts = {part.part_id: part for part in source.parts}
    header = {
        "schema": TRANSPORT_SCHEMA,
        "parts": [part_to_dict(p) for p in source.parts],
        "provenance": mutable_json(source.provenance),
    }
    checksum = hashlib.sha256()

    def write(value: dict[str, object]) -> None:
        line = canonical_json(value) + "\n"
        raw = line.encode()
        if len(raw) > MAX_LINE_BYTES:
            msg = "array input record exceeds the 16 MiB line limit"
            raise ValueError(msg)
        stream.write(line)
        checksum.update(raw)

    write(header)
    count = 0
    identities = set()
    for array in source.arrays:
        count += 1
        budget.examine()
        validate_array(array, parts)
        if array.source_id in identities:
            msg = f"duplicate array identity {array.source_id!r}"
            raise ValueError(msg)
        budget.retain()
        identities.add(array.source_id)
        write(realized_array_to_dict(array))
    stream.write(
        canonical_json(
            {"schema": END_SCHEMA, "arrays": count, "sha256": checksum.hexdigest()}
        )
        + "\n"
    )
    return count
