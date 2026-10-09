"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/preparation/motifs.py

Frozen preparation motif/scorer bindings shared by sources and screens.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import cached_property
from typing import TYPE_CHECKING

from dense_arrays._record_validation import integer, object_fields
from dense_arrays.parts.motifs import Motif, MotifInput, read_artifact
from dense_arrays.parts.motifs.windows import (
    MotifWindow,
    WindowSelection,
    select_window,
)
from dense_arrays.parts.scoring import FimoBinding, bind_fimo

from .requests import _absolute, _locator

if TYPE_CHECKING:
    from pathlib import Path

    from dense_arrays.parts import FimoScoring, PWMArtifact


@dataclass(frozen=True)
class MotifSource:
    """One input model and its frozen scorer and background identities."""

    input: MotifInput
    scoring: FimoBinding
    window: MotifWindow | None = None

    def __post_init__(self) -> None:
        """Require the source and scorer to describe the same model."""
        if not isinstance(self.input, MotifInput) or not isinstance(
            self.scoring, FimoBinding
        ):
            msg = "resolved PWM source and scoring motif must agree"
            raise TypeError(msg)
        if self.window is not None and not isinstance(self.window, MotifWindow):
            msg = "resolved PWM window requires MotifWindow"
            raise TypeError(msg)
        expected = self.selection.motif if self.selection else self.input.motif
        if expected != self.scoring.motif:
            msg = "resolved PWM source and scoring motif must agree"
            raise ValueError(msg)

    @cached_property
    def selection(self) -> WindowSelection | None:
        """Resolve a declared window from the preserved original model."""
        return select_window(self.input.motif, self.window) if self.window else None

    @property
    def motif(self) -> Motif:
        """Use one effective model for proposals, scoring and core diversity."""
        return self.scoring.motif

    def verify(self) -> None:
        """Recheck the source and scorer bytes before execution."""
        self.input.verify()
        self.scoring.verify()

    def window_bound(self, length: int, count: int) -> int:
        """Bound oriented candidate windows and one maximum-score calibration."""
        strands = 1 if self.scoring.settings.strands == "single" else 2
        return (count * max(0, length - self.motif.width + 1) + 1) * strands

    def to_dict(self, *, base: Path | None = None) -> dict[str, object]:
        """Serialize complete evidence using destination-relative locators."""
        scoring = self.scoring.to_dict()
        for key in ("executable", "background_path"):
            scoring[key] = _locator(scoring[key], base)
        return {
            "kind": "pwm_artifact",
            "path": _locator(self.input.path, base),
            "sha256": self.input.sha256,
            "motif": self.input.motif.to_dict(),
            "scoring": scoring,
            **(
                {"format": self.input.format} if self.input.format != "artifact" else {}
            ),
            **({"window": self.selection.to_dict()} if self.selection else {}),
        }

    @classmethod
    def from_dict(cls, value: object, *, base: Path | None = None) -> MotifSource:
        """Restore bindings without reading inputs or executing external tools."""
        keys = {"kind", "path", "sha256", "motif", "scoring"}
        data = object_fields(
            value, keys | {"window", "format"}, "resolved motif source"
        )
        window = data.pop("window", None)
        input_format = data.pop("format", "artifact")
        if set(data) != keys or data["kind"] != "pwm_artifact":
            msg = "incomplete or unsupported resolved motif source"
            raise ValueError(msg)
        scoring = dict(data["scoring"])
        for key in ("executable", "background_path"):
            scoring[key] = (
                str(_absolute(scoring[key], base)) if scoring[key] is not None else None
            )
        declared = None
        if window is not None:
            if not isinstance(window, dict):
                msg = "resolved motif window must be an object"
                raise ValueError(msg)
            integer(window.get("start"), field_name="motif window.start", minimum=0)
            integer(window.get("end"), field_name="motif window.end", minimum=1)
            background = window.get("background_source")
            declared = MotifWindow(
                window["end"] - window["start"],
                background=window.get("background")
                if background == "explicit"
                else background,
            )
        result = cls(
            MotifInput(
                Motif.from_dict(data["motif"]),
                _absolute(data["path"], base),
                data["sha256"],
                input_format,
            ),
            FimoBinding.from_dict(scoring),
            declared,
        )
        if window is not None and result.selection.to_dict() != window:
            msg = "resolved motif window evidence disagrees with source model"
            raise ValueError(msg)
        return result


def resolve_motif(
    source: PWMArtifact, scoring: FimoScoring
) -> tuple[PWMArtifact, MotifSource]:
    """Read a motif and resolve its scorer without sampling or scoring."""
    loaded = read_artifact(source)
    selected = select_window(loaded.motif, source.window) if source.window else None
    motif = selected.motif if selected else loaded.motif
    bound = MotifSource(loaded, bind_fimo(motif, scoring), source.window)
    return replace(source, path=loaded.path), bound
