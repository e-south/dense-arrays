"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/scoring/binding.py

Bind motif, background and executable identities without scoring candidates.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import hashlib
import os
import re
import shutil
import tempfile
from dataclasses import asdict, dataclass, replace
from pathlib import Path

from dense_arrays._record_validation import digest, object_fields, semantic_digest
from dense_arrays.parts.motifs import Motif
from dense_arrays.parts.motifs.models import BASES, numeric_row

from .configuration import FimoScoring, ScoringLimits
from .process import ScoringError, invoke
from .records import FimoHit


def fingerprint(path: Path) -> str:
    """Hash exact source bytes without loading an entire executable into memory."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _background(path: Path) -> tuple[float, ...]:
    fields = []
    for line in path.read_text().splitlines():
        fields.extend(line.partition("#")[0].split())
    if len(fields) != 2 * len(BASES) or set(fields[::2]) != set(BASES):
        msg = "background requires exactly one probability for each ACGT base"
        raise ValueError(msg)
    try:
        values = dict(zip(fields[::2], map(float, fields[1::2]), strict=True))
        return numeric_row(
            tuple(values[x] for x in BASES), probability=True, positive=True
        )
    except ValueError as err:
        msg = f"invalid background: {err}"
        raise ValueError(msg) from err


@dataclass(frozen=True)
class FimoBinding:
    """Resolved local scorer, model and background with verified source bytes."""

    motif: Motif
    settings: FimoScoring
    executable: Path
    executable_sha256: str
    version: str
    background: tuple[float, ...]
    background_path: Path | None
    background_sha256: str | None

    def __post_init__(self) -> None:
        """Validate restored bindings without reading files or invoking a scorer."""
        if not isinstance(self.motif, Motif) or not isinstance(
            self.settings, FimoScoring
        ):
            msg = "FIMO binding requires Motif and FimoScoring"
            raise TypeError(msg)
        object.__setattr__(self, "executable", Path(self.executable))
        digest(self.executable_sha256, field_name="FIMO executable_sha256")
        if (
            not isinstance(self.version, str)
            or re.fullmatch(r"\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.]+)?", self.version)
            is None
        ):
            msg = "unsupported FIMO version"
            raise ValueError(msg)
        object.__setattr__(
            self,
            "background",
            numeric_row(self.background, probability=True, positive=True),
        )
        if self.background_path is None:
            if (
                self.background_sha256 is not None
                or self.settings.background is not None
                or self.background != self.motif.background
            ):
                msg = "FIMO background binding disagrees with motif background"
                raise ValueError(msg)
        else:
            object.__setattr__(self, "background_path", Path(self.background_path))
            digest(self.background_sha256, field_name="FIMO background_sha256")
            if self.settings.background != self.background_path:
                msg = "FIMO background locator disagrees with settings"
                raise ValueError(msg)

    @classmethod
    def from_dict(cls, value: object) -> FimoBinding:
        """Load normalized evidence without reapplying defaults or calling FIMO."""
        keys = {
            "schema",
            "motif",
            "executable",
            "executable_sha256",
            "version",
            "background",
            "effective_background",
            "background_path",
            "background_sha256",
            "background_policy",
            "strands",
            "hit_pvalue_max",
            "pseudocount",
            "site_count",
            "encoding",
            "tie_policy",
            "limits",
        }
        data = object_fields(value, keys, "FIMO binding")
        if set(data) != keys or data["schema"] != "dense_arrays.fimo_binding.v1":
            msg = "unsupported or incomplete FIMO binding"
            raise ValueError(msg)
        limits = object_fields(
            data["limits"], {"seconds", "windows", "output_bytes"}, "scoring limits"
        )
        if len(limits) != len({"seconds", "windows", "output_bytes"}):
            msg = "incomplete scoring limits"
            raise ValueError(msg)
        settings = FimoScoring(
            data["hit_pvalue_max"],
            data["background_path"],
            data["strands"],
            data["pseudocount"],
            data["executable"],
            ScoringLimits(**limits),
        )
        result = cls(
            Motif.from_dict(data["motif"]),
            settings,
            data["executable"],
            data["executable_sha256"],
            data["version"],
            data["background"],
            data["background_path"],
            data["background_sha256"],
        )
        if data != result.to_dict():
            msg = "FIMO binding fields disagree"
            raise ValueError(msg)
        return result

    @property
    def effective_background(self) -> tuple[float, ...]:
        """Use complement-averaged frequencies for explicitly double-strand scans."""
        if self.settings.strands == "single":
            return self.background
        return tuple(
            (a + b) / 2
            for a, b in zip(self.background, reversed(self.background), strict=True)
        )

    @property
    def binding_id(self) -> str:
        """Identify scoring semantics independently of local file locations."""
        data = self.to_dict()
        data.pop("executable")
        data.pop("background_path")
        return semantic_digest(data)

    def to_dict(self) -> dict[str, object]:
        """Expose the scoring interpretation and its source fingerprints."""
        return {
            "schema": "dense_arrays.fimo_binding.v1",
            "motif": self.motif.to_dict(),
            "executable": str(self.executable),
            "executable_sha256": self.executable_sha256,
            "version": self.version,
            "background": list(self.background),
            "effective_background": list(self.effective_background),
            "background_path": str(self.background_path)
            if self.background_path
            else None,
            "background_sha256": self.background_sha256,
            "background_policy": "complement_average.v1"
            if self.settings.strands == "double"
            else "as_supplied.v1",
            "strands": self.settings.strands,
            "hit_pvalue_max": self.settings.hit_pvalue_max,
            "pseudocount": self.settings.pseudocount,
            "site_count": 20,
            "encoding": "meme_17_digits.v1",
            "tie_policy": "score_start_forward.v1",
            "limits": asdict(self.settings.limits),
        }

    def verify_hit(self, hit: FimoHit | None) -> None:
        """Check saved hit scope without rescoring or reclassifying rounded p-values."""
        if hit is None:
            return
        if not isinstance(hit, FimoHit):
            msg = "bound FIMO observations require FimoHit or None"
            raise TypeError(msg)
        if len(hit.core) != self.motif.width:
            msg = "FIMO hit width disagrees with its bound motif"
            raise ValueError(msg)
        if self.settings.strands == "single" and hit.strand != "forward":
            msg = "reverse FIMO hit is outside the single-strand scoring policy"
            raise ValueError(msg)

    def verify(self) -> None:
        """Refuse changed executable or background bytes before scoring."""
        if fingerprint(self.executable) != self.executable_sha256:
            msg = "FIMO executable changed since planning"
            raise ValueError(msg)
        if (
            self.background_path is not None
            and fingerprint(self.background_path) != self.background_sha256
        ):
            msg = "FIMO background changed since planning"
            raise ValueError(msg)


def bind_fimo(motif: Motif, settings: FimoScoring) -> FimoBinding:
    """Resolve sources and query the version; never score or sample a candidate."""
    if not isinstance(motif, Motif) or not isinstance(settings, FimoScoring):
        msg = "FIMO binding requires Motif and FimoScoring"
        raise TypeError(msg)
    background_path = (
        settings.background.absolute() if settings.background is not None else None
    )
    background_digest = fingerprint(background_path) if background_path else None
    background = _background(background_path) if background_path else motif.background
    found = (
        shutil.which("fimo")
        if settings.executable is None
        else str(settings.executable.absolute())
    )
    if found is None or not Path(found).is_file() or not os.access(found, os.X_OK):
        reason = "unavailable"
        raise ScoringError(reason, "install FIMO or provide its executable path")
    executable = Path(found).resolve()
    executable_digest = fingerprint(executable)
    with tempfile.TemporaryDirectory(prefix="dense-arrays-fimo-") as directory:
        output = invoke(
            [str(executable), "--version"], cwd=Path(directory), limits=settings.limits
        )
    version = output.stdout.decode("utf-8", errors="replace").strip()
    if re.fullmatch(r"\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.]+)?", version) is None:
        reason = "malformed"
        raise ScoringError(reason, "unrecognized FIMO version response")
    settings = replace(settings, executable=executable, background=background_path)
    result = FimoBinding(
        motif,
        settings,
        executable,
        executable_digest,
        version,
        background,
        background_path,
        background_digest,
    )
    result.verify()
    return result
