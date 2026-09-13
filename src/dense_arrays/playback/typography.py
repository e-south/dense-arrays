"""Publication typography shared by dense-array playback projections.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import dataclass


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
