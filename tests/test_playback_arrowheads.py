"""Keep the tapered head directional without changing the routed shaft.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import numpy as np
import pytest
from matplotlib.path import Path

from dense_arrays.playback.graph.arrowheads import arrow_head


@pytest.mark.parametrize("end", [(20, 0), (0, 20), (-20, 0), (0, -20)])
def test_tapered_head_follows_the_terminal_tangent_without_route_deformation(
    end: tuple[int, int],
) -> None:
    route = Path(
        [(0, 0), (end[0] / 2, end[1] / 2), end], [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    )
    before = route.vertices.copy()
    head = arrow_head(route, length=7.8, width=5.6)
    assert np.array_equal(route.vertices, before)
    assert np.array_equal(head.vertices[1], route.vertices[-1])
    assert list(head.codes).count(Path.CURVE4) == 3
    # A filled, notched head has no vertices beyond the intended endpoint.
    tangent = np.asarray(end) / 20
    assert ((head.vertices[:-1] - end) @ tangent <= 0).all()


def test_zero_length_terminal_tangent_does_not_invent_a_direction():
    path = Path([(0, 0), (20, 0), (20, 0)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    assert arrow_head(path, length=7.8, width=5.6) is None
