"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/graph/arrowheads.py

Draw filled, tapered arrowheads along an existing route's terminal tangent.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.path import Path


def tapered_arrow(
    path: Path, mutation_size: float, linewidth: float, aspect_ratio: float = 1.0
) -> tuple[list[Path], list[bool]]:
    """Use Matplotlib's public arrow callable contract; preserve the curved shaft.

    The head is eleven units long and eight units wide with a shallow curved
    notch. Its length follows ``mutation_size``, independent of emphasis width.
    The routed endpoint remains the tip. A zero terminal tangent has no head.
    """
    from matplotlib.path import Path

    del linewidth, aspect_ratio  # Presentation does not distort head geometry.
    tip = path.vertices[-1]
    tangent = tip - path.vertices[-2]
    norm = float(np.hypot(*tangent))
    if norm == 0:
        return [path], [False]
    along = tangent / norm
    across = np.array([-along[1], along[0]])
    shape = np.array(
        [(-11, 4), (0, 0), (-11, -4), (-9, -1.67), (-9, 1.66), (-11, 4), (-11, 4)]
    ) * (mutation_size / 11)
    vertices = tip + shape[:, :1] * along + shape[:, 1:] * across
    head = Path(
        vertices,
        [
            Path.MOVETO,
            Path.LINETO,
            Path.LINETO,
            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,
            Path.CLOSEPOLY,
        ],
    )
    return [path, head], [False, True]
