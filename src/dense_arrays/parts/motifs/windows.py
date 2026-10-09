"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/motifs/windows.py

Contiguous motif windows selected by relative entropy against a declared null.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from dense_arrays._record_validation import integer, object_fields

from .models import Motif, numeric_row

WINDOW_POLICY = "max_relative_entropy.v1"


@dataclass(frozen=True)
class MotifWindow:
    """Select a fixed-width window by relative entropy; ties choose the earliest.

    Background is 'motif' (the default), 'uniform', or four positive ACGT
    probabilities. It defines selection only, independently of sampling/scoring.
    """

    length: int
    background: str | tuple[float, ...] = "motif"

    def __post_init__(self) -> None:
        """Require an explicit positive width independent of sequence length."""
        integer(self.length, field_name="motif window.length", minimum=1)
        if isinstance(self.background, str):
            if self.background not in {"uniform", "motif"}:
                msg = (
                    "motif window background must be motif, uniform "
                    "or ACGT probabilities"
                )
                raise ValueError(msg)
        else:
            object.__setattr__(
                self,
                "background",
                numeric_row(self.background, probability=True, positive=True),
            )

    def to_dict(self) -> dict[str, object]:
        """Encode the declared width; the resolved plan binds the algorithm."""
        return {
            "length": self.length,
            "background": self.background
            if isinstance(self.background, str)
            else list(self.background),
        }

    @classmethod
    def from_dict(cls, value: object) -> MotifWindow:
        """Reject undeclared strategy or coordinate overrides."""
        return cls(**object_fields(value, {"length", "background"}, "motif window"))


@dataclass(frozen=True)
class WindowSelection:
    """Selected model and zero-based, half-open coordinates in its source model."""

    source_model_id: str
    motif: Motif
    start: int
    end: int
    information_bits: float
    source_information_bits: float
    background: tuple[float, ...]
    background_source: str

    def to_dict(self) -> dict[str, object]:
        """Expose coordinates, information retained and both computational IDs."""
        return {
            "schema": "dense_arrays.motif_window.v1",
            "policy": WINDOW_POLICY,
            "source_model_id": self.source_model_id,
            "model_id": self.motif.model_id,
            "start": self.start,
            "end": self.end,
            "information_bits": self.information_bits,
            "source_information_bits": self.source_information_bits,
            "discarded_information_bits": self.source_information_bits
            - self.information_bits,
            "background": list(self.background),
            "background_source": self.background_source,
            "retained_information_fraction": (
                self.information_bits / self.source_information_bits
                if self.source_information_bits > 0
                else None
            ),
        }


def select_window(motif: Motif, window: MotifWindow) -> WindowSelection:
    """Maximize summed D(P_i || B), preserving aligned supplied score rows.

    Zero-probability contributions are zero. No inferred counts or pseudocounts
    modify the model. Accurate sums give equal windows equal scores. The uniform
    background recovers 2-H; this quantity does not establish retained affinity.
    """
    if not isinstance(motif, Motif) or not isinstance(window, MotifWindow):
        msg = "window selection requires Motif and MotifWindow"
        raise TypeError(msg)
    if window.length > motif.width:
        msg = "motif window length exceeds source motif width"
        raise ValueError(msg)
    background_source = (
        window.background if isinstance(window.background, str) else "explicit"
    )
    background = (
        motif.background
        if background_source == "motif"
        else (0.25,) * 4
        if background_source == "uniform"
        else window.background
    )
    columns = tuple(
        max(
            0.0,
            math.fsum(
                p * (math.log2(p) - math.log2(b))
                for p, b in zip(row, background, strict=True)
                if p > 0
            ),
        )
        for row in motif.probabilities
    )
    scores = (
        math.fsum(columns[start : start + window.length])
        for start in range(motif.width - window.length + 1)
    )
    start, information = max(enumerate(scores), key=lambda item: item[1])
    end = start + window.length
    selected = replace(
        motif,
        probabilities=motif.probabilities[start:end],
        log_odds=motif.log_odds[start:end] if motif.log_odds is not None else None,
    )
    return WindowSelection(
        motif.model_id,
        selected,
        start,
        end,
        information,
        math.fsum(columns),
        background,
        background_source,
    )
