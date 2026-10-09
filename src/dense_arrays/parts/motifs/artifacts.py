"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/artifacts.py

Bind explicit motif input formats without importing preparation tools.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from dense_arrays._record_validation import digest, integer, required_text

from .models import BASES, Motif
from .windows import MotifWindow


@dataclass(frozen=True)
class PWMArtifact:
    """One motif file, explicit format and optional exact motif ID selection."""

    path: Path | str
    motif_ids: tuple[str, ...] = ()
    window: MotifWindow | None = None
    format: str = "artifact"

    def __post_init__(self) -> None:
        """Freeze the locator and reject ambiguous selection declarations."""
        object.__setattr__(self, "path", Path(self.path))
        if self.format not in {"artifact", "meme", "jaspar"}:
            msg = "supported motif input formats: artifact, meme, jaspar"
            raise ValueError(msg)
        if not isinstance(self.motif_ids, (tuple, list)):
            msg = "motif_ids must be an ordered array"
            raise TypeError(msg)
        for value in self.motif_ids:
            required_text(value, field_name="motif_ids")
        if len(set(self.motif_ids)) != len(self.motif_ids):
            msg = "motif_ids must not repeat identities"
            raise ValueError(msg)
        object.__setattr__(self, "motif_ids", tuple(self.motif_ids))
        if self.window is not None and not isinstance(self.window, MotifWindow):
            msg = "PWM window requires MotifWindow"
            raise TypeError(msg)


@dataclass(frozen=True)
class MotifInput:
    """Validated motif meaning plus the source bytes used to resolve it."""

    motif: Motif
    path: Path
    sha256: str
    format: str = "artifact"

    def __post_init__(self) -> None:
        """Require validated motif meaning and a complete source fingerprint."""
        if not isinstance(self.motif, Motif):
            msg = "motif input requires Motif"
            raise TypeError(msg)
        object.__setattr__(self, "path", Path(self.path))
        digest(self.sha256, field_name="motif input.sha256")
        if self.format not in {"artifact", "meme", "jaspar"}:
            msg = "unsupported motif input format"
            raise ValueError(msg)

    def verify(self) -> None:
        """Refuse changed input before preparation starts."""
        if hashlib.sha256(self.path.read_bytes()).hexdigest() != self.sha256:
            msg = "motif input changed since planning; create a new plan"
            raise ValueError(msg)


def _row(value: object) -> tuple:
    if not isinstance(value, dict) or set(value) != set(BASES):
        msg = "motif rows require exactly the ACGT keys"
        raise ValueError(msg)
    return tuple(value[base] for base in BASES)


def _rows(value: object) -> tuple:
    if not isinstance(value, list):
        msg = "motif matrices must be ordered arrays"
        raise TypeError(msg)
    return tuple(_row(row) for row in value)


def read_artifact(source: PWMArtifact) -> MotifInput:
    """Resolve one motif from JSON, minimal MEME or JASPAR and bind source bytes.

    Multi-record files require one exact motif ID. The complete file is
    validated before selection; labels and source statistics remain metadata.
    """
    if not isinstance(source, PWMArtifact):
        msg = "motif import requires PWMArtifact"
        raise TypeError(msg)
    payload = source.path.read_bytes()
    if source.format != "artifact":
        from .formats import read_jaspar, read_meme  # noqa: PLC0415

        motifs = (read_meme if source.format == "meme" else read_jaspar)(payload)
        selected = [
            m for m in motifs if not source.motif_ids or m.motif_id in source.motif_ids
        ]
        missing = set(source.motif_ids) - {m.motif_id for m in motifs}
        if missing:
            msg = f"missing motif IDs: {', '.join(sorted(missing))}"
            raise ValueError(msg)
        if len(selected) != 1:
            msg = "preparation requires one motif; select one exact motif ID"
            raise ValueError(msg)
        return MotifInput(
            selected[0],
            source.path.absolute(),
            hashlib.sha256(payload).hexdigest(),
            source.format,
        )
    value = json.loads(payload, object_pairs_hook=_unique_object)
    required = {
        "schema_version",
        "producer",
        "motif_id",
        "alphabet",
        "matrix_semantics",
        "background",
        "probabilities",
        "log_odds",
    }
    if not isinstance(value, dict) or not required <= value.keys():
        msg = "motif artifact is missing required fields"
        raise ValueError(msg)
    if (
        value["schema_version"] != "1.0"
        or value["alphabet"] != BASES
        or value["matrix_semantics"] != "probabilities"
    ):
        msg = "unsupported motif artifact schema, alphabet or matrix semantics"
        raise ValueError(msg)
    motif = Motif(
        value["motif_id"],
        _rows(value["probabilities"]),
        _row(value["background"]),
        _rows(value["log_odds"]),
        value["producer"],
        {k: v for k, v in value.items() if k not in required | {"length"}},
    )
    if "length" in value:
        integer(value["length"], field_name="motif.length", minimum=1)
        if value["length"] != motif.width:
            msg = "declared motif length does not match matrix rows"
            raise ValueError(msg)
    if source.motif_ids and source.motif_ids != (motif.motif_id,):
        msg = "requested motif identities do not match the single-motif artifact"
        raise ValueError(msg)
    return MotifInput(
        motif, source.path.absolute(), hashlib.sha256(payload).hexdigest()
    )


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate keys instead of accepting a silently overwritten matrix."""
    result = {}
    for key, value in pairs:
        if key in result:
            msg = f"duplicate motif artifact key: {key!r}"
            raise ValueError(msg)
        result[key] = value
    return result
