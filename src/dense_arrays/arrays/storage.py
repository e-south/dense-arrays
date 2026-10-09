"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/storage.py

Read immutable collection manifests and checksum-bound SQLite records.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
import json
import sqlite3
from contextlib import closing, contextmanager
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    integer,
    object_fields,
    required_text,
    semantic_digest,
)
from dense_arrays.artifacts.errors import ArtifactIntegrityError, integrity_boundary
from dense_arrays.artifacts.store import checked_payload
from dense_arrays.parts.serialization import part_from_dict

from .models import (
    BOUNDARY,
    DATABASE,
    MANIFEST,
    SCHEMA,
    ArrayCollection,
    CollectionSummary,
)

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadBudget, ReadLimits
    from dense_arrays.parts import Part


def is_collection(value: object) -> bool:
    """Recognize typed sources or a committed/pending collection directory."""
    return isinstance(value, ArrayCollection) or (
        isinstance(value, (str, Path))
        and (
            (Path(value) / MANIFEST).is_file()
            or (Path(value) / ".collection-pending").exists()
        )
    )


def file_digest(path: Path) -> str:
    """Hash bytes incrementally, including collections larger than memory."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_summary(path: Path, limits: ReadLimits) -> CollectionSummary:
    """Reject incomplete, unknown or internally contradictory manifests."""
    with integrity_boundary(path):
        if (path / ".collection-pending").exists():
            msg = "array collection publication is incomplete"
            raise ValueError(msg)
        data = json.loads((path / MANIFEST).read_text())
        keys = {
            "schema",
            "collection_id",
            "arrays",
            "placements",
            "parts",
            "sequences",
            "database",
            "exporter",
            "provenance",
            "evidence_boundary",
        }
        data = object_fields(data, keys, "array collection")
        if (
            set(data) != keys
            or data["schema"] != SCHEMA
            or data["evidence_boundary"] != BOUNDARY
        ):
            msg = "unsupported or incomplete array collection schema"
            raise ValueError(msg)
        identity = data.pop("collection_id")
        if identity != semantic_digest(data):
            msg = "array collection manifest digest differs"
            raise ValueError(msg)
        for name in ("arrays", "placements", "parts", "sequences"):
            integer(data[name], field_name=name, minimum=0)
        if data["sequences"] > data["arrays"] or data["placements"] < data["arrays"]:
            msg = "array collection counts contradict their populations"
            raise ValueError(msg)
        exporter = object_fields(
            data["exporter"], {"package", "version"}, "collection exporter"
        )
        if set(exporter) != {"package", "version"}:
            msg = "collection exporter requires package and version"
            raise ValueError(msg)
        required_text(exporter["package"], field_name="exporter.package")
        required_text(exporter["version"], field_name="exporter.version")
        database = object_fields(
            data["database"], {"sha256", "bytes"}, "collection database"
        )
        digest(database.get("sha256"), field_name="database.sha256")
        integer(database.get("bytes"), field_name="database.bytes", minimum=1)
        if (path / DATABASE).is_symlink() or (
            path / DATABASE
        ).stat().st_size != database["bytes"]:
            msg = "collection database bytes differ"
            raise ValueError(msg)
        return CollectionSummary(
            {**data, "collection_id": identity}, read_limits=limits
        )


@contextmanager
def reader(path: Path) -> Iterator[sqlite3.Connection]:
    """Open a read-only connection and close it on exhaustion or early exit."""
    try:
        with closing(
            sqlite3.connect((path / DATABASE).resolve().as_uri() + "?mode=ro", uri=True)
        ) as connection:
            connection.execute("PRAGMA query_only=ON")
            yield connection
    except sqlite3.DatabaseError as error:
        raise ArtifactIntegrityError(str(error), artifact=path) from error


def read_parts(
    connection: sqlite3.Connection, summary: CollectionSummary, budget: ReadBudget
) -> dict[str, Part]:
    """Load the bounded part catalog once, preserving unused input identities."""
    budget.retain(summary.parts)
    parts = {}
    with integrity_boundary(None):
        for identity, payload, checksum in connection.execute(
            "SELECT part_id,payload,digest FROM parts ORDER BY ordinal"
        ):
            budget.examine(payload)
            part = part_from_dict(checked_payload((payload, checksum)))
            if part.part_id != identity or identity in parts:
                msg = "collection part identity differs"
                raise ValueError(msg)
            parts[identity] = part
        if len(parts) != summary.parts:
            msg = "collection part count differs"
            raise ValueError(msg)
    return parts
