"""Fit a compact nucleotide grid to a native duplex viewport.

Module Author(s): Eric J. South
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from .timeline import complement_sequence
from .typography import (
    PUBLICATION_LABEL_FONT_SIZE_PT,
    PUBLICATION_NUCLEOTIDE_TYPOGRAPHY,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.font_manager import FontProperties
    from matplotlib.text import Text

FEATURE_HEIGHT = 0.56


@dataclass(frozen=True)
class DuplexGeometry:
    """One physical nucleotide scale shared by every strand and motif row."""

    font: FontProperties
    cap_height_px: float
    label_font_size_pt: float


def fit_duplex_grid(
    axis: Axes,
    sequence: str,
    vertical_bounds: tuple[float, float],
    *,
    bottom_padding_pt: float = 0,
) -> DuplexGeometry:
    """Center measured cells; shrink uniformly only when the viewport requires it."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    typography = PUBLICATION_NUCLEOTIDE_TYPOGRAPHY
    font = FontProperties(
        family=typography.family,
        weight=typography.weight,
        size=typography.graph_font_size_pt,
    )
    outlines = [TextPath((0, 0), base, prop=font).get_extents() for base in "ACGT"]
    cap_height = max(bounds.height for bounds in outlines)
    symbols = set(sequence + complement_sequence(sequence))
    widest_glyph = max(
        TextPath((0, 0), base, prop=font).get_extents().width for base in symbols
    )
    cell_width = widest_glyph + 0.2 * cap_height
    length = len(sequence)
    # A motif bar contains the cap height and 15% padding per side.
    vertical_unit = cap_height * 1.3 / FEATURE_HEIGHT
    viewport = axis.get_window_extent(axis.figure.canvas.get_renderer())
    pixels_per_point = axis.figure.dpi / 72
    width_pt, height_pt = (
        viewport.width / pixels_per_point,
        viewport.height / pixels_per_point,
    )
    bottom, top = vertical_bounds
    scale = min(
        1.0,
        width_pt / ((length + 7) * cell_width),
        (height_pt - bottom_padding_pt) / ((top - bottom) * vertical_unit),
    )
    font.set_size(typography.graph_font_size_pt * scale)
    horizontal_span = width_pt / (cell_width * scale)
    vertical_span = height_pt / (vertical_unit * scale)
    axis.set_xlim(length / 2 - horizontal_span / 2, length / 2 + horizontal_span / 2)
    middle = (bottom + top) / 2 - bottom_padding_pt / (2 * vertical_unit * scale)
    axis.set_ylim(middle - vertical_span / 2, middle + vertical_span / 2)
    return DuplexGeometry(
        font,
        cap_height * scale * pixels_per_point,
        PUBLICATION_LABEL_FONT_SIZE_PT * scale,
    )


def draw_nucleotide(
    axis: Axes,
    position: tuple[float, float],
    base: str,
    color: str,
    geometry: DuplexGeometry,
) -> Text:
    """Center the actual glyph outline in its coordinate cell and motif box."""
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import ScaledTranslation

    bounds = TextPath((0, 0), base, prop=geometry.font).get_extents()
    offset = ScaledTranslation(
        -(bounds.x0 + bounds.width / 2) / 72,
        -(bounds.y0 + bounds.height / 2) / 72,
        axis.figure.dpi_scale_trans,
    )
    return axis.text(
        *position,
        base,
        ha="left",
        va="baseline",
        color=color,
        fontproperties=geometry.font,
        transform=axis.transData + offset,
    )
