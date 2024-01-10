"""Tests of the optimize.py module."""

import dense_arrays as da
import pytest


@pytest.mark.parametrize(
    "motifa, motifb, shift",
    [
        ("ATGCATTA", "CATTATG", 3),
        ("ATGCATTA", "TATGA", 6),
        ("ATGTTAACT", "TTAA", 8),
        ("AGC", "CAG", 2),
        ("CAG", "AGC", 1),
        ("CGT", "CAG", 3),
        ("", "CGT", 0),
    ],
)
def test_shift_metric(motifa: str, motifb: str, shift: int):
    assert da.optimize.shift_metric(motifa, motifb) == shift
