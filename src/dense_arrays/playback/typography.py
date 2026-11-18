"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/typography.py

Publication typography shared by dense-array playback projections.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, fields


@dataclass(frozen=True, slots=True)
class NucleotideTypography:
    """Backend-calibrated typography for one rendered nucleotide scale."""

    family: str
    weight: str
    graph_font_size_pt: float
    duplex_font_size_pt: float
    target_cap_height_px: float

    def __post_init__(self) -> None:
        """Validate nonempty font names and positive sizing values."""
        if not self.family.strip():
            msg = "nucleotide font family must not be empty"
            raise ValueError(msg)
        if not self.weight.strip():
            msg = "nucleotide font weight must not be empty"
            raise ValueError(msg)
        for field_name in (
            "graph_font_size_pt",
            "duplex_font_size_pt",
            "target_cap_height_px",
        ):
            if float(getattr(self, field_name)) <= 0:
                msg = f"{field_name} must be positive"
                raise ValueError(msg)


PUBLICATION_NUCLEOTIDE_TYPOGRAPHY = NucleotideTypography(
    family="Arial",
    weight="normal",
    graph_font_size_pt=13.2,
    duplex_font_size_pt=30.0,
    target_cap_height_px=18.0,
)

PUBLICATION_LABEL_FONT_SIZE_PT = 11.5


@dataclass(frozen=True, slots=True)
class PlaybackTypography:
    """Linked native type sizes; explicit nucleotide sizes are never shrunk."""

    nucleotide_font_size_pt: float | None = None
    annotation_scale: float = (
        PUBLICATION_LABEL_FONT_SIZE_PT
        / PUBLICATION_NUCLEOTIDE_TYPOGRAPHY.graph_font_size_pt
    )

    def __post_init__(self) -> None:
        """Reject invalid dimensions before measuring or allocating media."""
        for name in ("nucleotide_font_size_pt", "annotation_scale"):
            value = getattr(self, name)
            if name == "nucleotide_font_size_pt" and value is None:
                continue
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                msg = f"{name} must be numeric"
                raise TypeError(msg)
            if not math.isfinite(value) or value <= 0:
                msg = f"{name} must be finite and positive"
                raise ValueError(msg)

    @classmethod
    def from_mapping(cls, raw: Mapping[str, object]) -> PlaybackTypography:
        """Read the same strict contract from API mappings and CLI JSON."""
        if not isinstance(raw, Mapping):
            msg = "typography must be a mapping"
            raise TypeError(msg)
        unknown = set(raw) - {item.name for item in fields(cls)}
        if unknown:
            msg = f"unknown typography fields: {sorted(unknown)}"
            raise ValueError(msg)
        return cls(**raw)
