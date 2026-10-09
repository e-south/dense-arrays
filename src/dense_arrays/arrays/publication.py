"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/publication.py

Stream checked supplied arrays into an exclusively owned portable collection.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
import shutil
import sqlite3
from contextlib import closing
from importlib.metadata import version
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    canonical_json,
    mutable_json,
    semantic_digest,
)
from dense_arrays.artifacts.publication import write_new
from dense_arrays.artifacts.reading import ReadBudget, ReadLimits
from dense_arrays.artifacts.receipts import ExportReceipt
from dense_arrays.parts.serialization import part_to_dict
from dense_arrays.playback.serialization import realized_array_to_dict

from .geometry import validate_array
from .models import (
    BOUNDARY,
    DATABASE,
    MANIFEST,
    SCHEMA,
    ArrayCollection,
    sequence_identity,
)
from .storage import file_digest

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.parts import Part


def publish_collection(
    source: ArrayCollection, out: Path, limits: ReadLimits
) -> ExportReceipt:
    """Validate records; publish the manifest after the database commits."""
    budget = ReadBudget(limits)
    budget.retain(len(source.parts))
    parts = {part.part_id: part for part in source.parts}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.mkdir(exist_ok=False)
    owned = out.stat()
    marker = out / ".collection-pending"
    marker.touch(exist_ok=False)
    try:
        manifest = _store(source, out, parts, budget)
        manifest["collection_id"] = semantic_digest(manifest)
        write_new(out / MANIFEST, canonical_json(manifest) + "\n")
        marker.unlink()
        descriptor = os.open(out, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    except BaseException:
        if out.exists() and (out.stat().st_dev, out.stat().st_ino) == (
            owned.st_dev,
            owned.st_ino,
        ):
            shutil.rmtree(out)
        raise
    return ExportReceipt(
        str(out),
        "bundle",
        "arrays",
        manifest["arrays"],
        ({"collection_id": manifest["collection_id"]},),
        (),
        tuple(
            {
                "name": name,
                "bytes": (out / name).stat().st_size,
                "sha256": file_digest(out / name),
            }
            for name in (MANIFEST, DATABASE)
        ),
    )


def _store(
    source: ArrayCollection, out: Path, parts: dict[str, Part], budget: ReadBudget
) -> dict[str, object]:
    with closing(sqlite3.connect(out / DATABASE)) as connection:
        connection.execute("PRAGMA synchronous=FULL")
        connection.executescript("""
            CREATE TABLE parts(ordinal INTEGER PRIMARY KEY,
                part_id TEXT UNIQUE NOT NULL,
                payload TEXT NOT NULL, digest TEXT NOT NULL);
            CREATE TABLE arrays(ordinal INTEGER PRIMARY KEY,
                array_id TEXT UNIQUE NOT NULL,
                sequence_id TEXT NOT NULL, payload TEXT NOT NULL, digest TEXT NOT NULL);
        """)
        count, placements = 0, 0
        with connection:
            for number, part in enumerate(source.parts, 1):
                value = part_to_dict(part)
                payload = canonical_json(value)
                budget.examine(payload)
                connection.execute(
                    "INSERT INTO parts VALUES (?,?,?,?)",
                    (number, part.part_id, payload, semantic_digest(value)),
                )
            for count, array in enumerate(source.arrays, 1):
                budget.examine()
                validate_array(array, parts)
                value = realized_array_to_dict(array)
                try:
                    connection.execute(
                        "INSERT INTO arrays VALUES (?,?,?,?,?)",
                        (
                            count,
                            array.source_id,
                            sequence_identity(array.sequence),
                            canonical_json(value),
                            semantic_digest(value),
                        ),
                    )
                except sqlite3.IntegrityError as error:
                    msg = f"duplicate array identity {array.source_id!r}"
                    raise ValueError(msg) from error
                placements += len(array.placements)
        sequences = connection.execute(
            "SELECT COUNT(DISTINCT sequence_id) FROM arrays"
        ).fetchone()[0]
    return {
        "schema": SCHEMA,
        "evidence_boundary": BOUNDARY,
        "arrays": count,
        "placements": placements,
        "parts": len(parts),
        "sequences": sequences,
        "database": {
            "sha256": file_digest(out / DATABASE),
            "bytes": (out / DATABASE).stat().st_size,
        },
        "exporter": {"package": "dense-arrays", "version": version("dense-arrays")},
        "provenance": mutable_json(source.provenance),
    }
