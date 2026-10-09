"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/bound.py

Resolved part sources with immutable import evidence and explicit file bindings.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays._record_validation import (
    digest,
    object_fields,
    records,
    semantic_digest,
)

from .ingestion import validate_parts
from .provenance import ImportReport
from .serialization import part_from_dict, part_to_dict

if TYPE_CHECKING:
    from .models import Part

BOUND_SCHEMA = "dense_arrays.bound_parts.v1"


def encoded_bound_size(value: object) -> int:
    """Count embedded source state before constructing its typed records."""
    if not isinstance(value, dict) or not all(
        isinstance(value.get(key), list) for key in ("parts", "inputs")
    ):
        msg = "bound parts require part and input arrays"
        raise TypeError(msg)
    report = value.get("import_report")
    if not isinstance(report, dict) or not isinstance(
        report.get("transformations"), list
    ):
        msg = "bound parts require an import report with transformations"
        raise TypeError(msg)
    return len(value["parts"]) + len(value["inputs"]) + len(report["transformations"])


def collection_identity(
    parts: tuple[Part, ...], report: ImportReport, input_digests: tuple[str, ...]
) -> str:
    """Keep the part namespace independent of design rules and source locations."""
    return report.collection_id or semantic_digest(
        {
            "schema": "dense_arrays.part_collection.v1",
            "parts": [part_to_dict(p) for p in parts],
            "import_report": report.to_dict(),
            "input_digests": list(input_digests),
        }
    )


@dataclass(frozen=True, repr=False)
class BoundParts:
    """Resolved parts and their origins; locations require byte verification.

    Absent locations mean the parts are embedded evidence. Origin fingerprints
    are retained without claiming to read or verify the unavailable source files.
    Replace the whole source to deliberately start a new part collection.
    """

    parts: tuple[Part, ...]
    import_report: ImportReport
    input_digests: tuple[str, ...] = ()
    locations: tuple[Path, ...] | None = None
    snapshot_id: str = field(init=False)

    def __post_init__(self) -> None:
        """Freeze one complete source with matching provenance and ordered locators."""
        object.__setattr__(self, "parts", validate_parts(self.parts))
        if not isinstance(
            self.import_report, ImportReport
        ) or self.import_report.rows != len(self.parts):
            msg = "bound parts require an import report for the same records"
            raise ValueError(msg)
        fingerprints = records(self.input_digests, str, field_name="input_digests")
        for value in fingerprints:
            digest(value, field_name="input.sha256")
        object.__setattr__(self, "input_digests", fingerprints)
        if self.locations is not None:
            if not isinstance(self.locations, (tuple, list)) or len(
                self.locations
            ) != len(fingerprints):
                msg = "bound part locations must match every input fingerprint"
                raise ValueError(msg)
            object.__setattr__(
                self, "locations", tuple(Path(p).absolute() for p in self.locations)
            )
        object.__setattr__(self, "snapshot_id", semantic_digest(self._content()))

    @property
    def collection_id(self) -> str:
        """Identify the part namespace while preserving source-pool identity."""
        return collection_identity(self.parts, self.import_report, self.input_digests)

    @property
    def identities(self) -> int:
        """Count embedded records, source fingerprints and normalization evidence."""
        return (
            len(self.parts)
            + len(self.input_digests)
            + len(self.import_report.transformations)
        )

    def _content(self) -> dict[str, object]:
        return {
            "schema": BOUND_SCHEMA,
            "parts": [part_to_dict(p) for p in self.parts],
            "import_report": self.import_report.to_dict(),
            "input_digests": list(self.input_digests),
        }

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Encode locations relative to the request without changing part identity."""
        value = self._content()
        del value["input_digests"]
        return {
            **value,
            "snapshot_id": self.snapshot_id,
            "inputs": [
                {
                    "sha256": fingerprint,
                    "path": None
                    if self.locations is None
                    else (
                        str(self.locations[i])
                        if base is None
                        else os.path.relpath(self.locations[i], base)
                    ),
                }
                for i, fingerprint in enumerate(self.input_digests)
            ],
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> BoundParts:
        """Reject changed records and incomplete source binding declarations."""
        keys = {"schema", "parts", "import_report", "inputs", "snapshot_id"}
        data = object_fields(value, keys, "bound parts")
        if set(data) != keys or data["schema"] != BOUND_SCHEMA:
            msg = "unsupported or incomplete bound parts schema"
            raise ValueError(msg)
        if not isinstance(data["parts"], list) or not isinstance(data["inputs"], list):
            msg = "bound parts records and inputs must be arrays"
            raise TypeError(msg)
        inputs = [
            object_fields(i, {"path", "sha256"}, "bound input") for i in data["inputs"]
        ]
        if any(set(i) != {"path", "sha256"} for i in inputs):
            msg = "bound inputs require path and sha256"
            raise ValueError(msg)
        paths = [i["path"] for i in inputs]
        if any(p is None for p in paths) and any(p is not None for p in paths):
            msg = "bound parts require all source locations or embedded evidence"
            raise ValueError(msg)
        result = cls(
            tuple(part_from_dict(p) for p in data["parts"]),
            ImportReport.from_dict(data["import_report"]),
            tuple(i["sha256"] for i in inputs),
            None
            if not paths or paths[0] is None
            else tuple((base or Path.cwd()) / p for p in paths),
        )
        if result.snapshot_id != data["snapshot_id"]:
            msg = "bound parts identity does not match its records and provenance"
            raise ValueError(msg)
        return result

    def __repr__(self) -> str:
        """Keep notebook display bounded and disclose the source verification mode."""
        mode = "embedded" if self.locations is None else "bound"
        return f"BoundParts({self.snapshot_id[:12]}, {len(self.parts)} parts, {mode})"
