"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/graph/arrowheads.py

Measure tapered heads and draw shafts separately from unstroked fills.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.path import Path

    from ..graph_style import GraphArrowDimensions


def arrow_head(
    path: Path, *, length: float, width: float, shape: str = "triangle"
) -> Path | None:
    """Return a filled-head path with its tip at the existing route endpoint."""
    from matplotlib.path import Path

    tip = path.vertices[-1]
    tangent = tip - path.vertices[-2]
    norm = float(np.hypot(*tangent))
    if norm == 0:
        return None
    along = tangent / norm
    across = np.array([-along[1], along[0]])
    if shape == "triangle":
        points = np.array([(-1, 0.5), (0, 0), (-1, -0.5), (-1, 0.5)]) * (length, width)
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    elif shape == "notched":
        points = np.array(
            [(-11, 4), (0, 0), (-11, -4), (-9, -1.67), (-9, 1.66), (-11, 4), (-11, 4)]
        ) * (length / 11, width / 8)
        codes = [
            Path.MOVETO,
            Path.LINETO,
            Path.LINETO,
            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,
            Path.CLOSEPOLY,
        ]
    else:
        msg = "head shape must be triangle or notched"
        raise ValueError(msg)
    vertices = tip + points[:, :1] * along + points[:, 1:] * across
    return Path(vertices, codes)


def draw_arrow(
    axis: Axes,
    path: Path,
    dimensions: GraphArrowDimensions,
    *,
    color: str,
    zorder: float,
) -> None:
    """Draw a thin shaft and independently filled head; never stroke the head."""
    from matplotlib.patches import FancyArrowPatch, PathPatch

    shaft, length, width = (
        dimensions.shaft_width_pt,
        dimensions.head_length_pt,
        dimensions.head_width_pt,
    )
    opacity = dimensions.opacity
    axis.add_patch(
        FancyArrowPatch(
            path=path,
            arrowstyle="-",
            linewidth=shaft,
            color=color,
            alpha=opacity,
            capstyle="butt",
            zorder=zorder,
            gid="graph-arrow-shaft",
        )
    )
    head = arrow_head(path, length=length, width=width, shape=dimensions.head_shape)
    if head is not None:
        axis.add_patch(
            PathPatch(
                head,
                facecolor=color,
                alpha=opacity,
                edgecolor="none",
                linewidth=0,
                zorder=zorder,
                gid="graph-arrow-head",
            )
        )
