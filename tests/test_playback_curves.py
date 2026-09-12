"""Independent numerical checks for playback curve geometry.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import pytest

from dense_arrays.playback.graph.model import QuadraticCurve
from dense_arrays.playback.graph.routing import (
    quadratic_arc_length,
    quadratic_arc_t,
    quadratic_point,
    quadratic_segment,
)


def test_quadratic_midpoint_uses_control_point() -> None:
    assert quadratic_point((0, 0), (2, 4), (4, 0), 0.5) == (2, 2)


def test_straight_curve_distance_and_fraction_are_exact() -> None:
    curve = QuadraticCurve((0, 0), (2, 0), (4, 0), (0, 0), (4, 0))
    assert quadratic_arc_length(curve) == pytest.approx(4)
    assert quadratic_arc_t(curve, 0.25) == pytest.approx(0.25)
    assert quadratic_segment(curve, 0.25, 0.75) == ((1, 0), (2, 0), (3, 0))
