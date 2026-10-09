"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/planning/matrices/expansion.py

Resolve bounded named combinations without opening sources or constructing models.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

import itertools
import math
from collections.abc import Mapping

from .requests import MatrixSpec


def cell_id(choices: Mapping[str, str]) -> str:
    """Stable readable identity independent of axis enumeration order."""
    return ",".join(f"{name}={choices[name]}" for name in sorted(choices))


def expansion_size(request: MatrixSpec) -> int:
    """Validate cardinality without allocating combinations or reading inputs."""
    axes = request.axes
    if request.pairing == "cross_product":
        count = math.prod(len(options) for options in axes.values())
    elif request.pairing == "zip":
        names = tuple(next(iter(axes.values())))
        if any(set(options) != set(names) for options in axes.values()):
            msg = "zip pairing requires identical named choices in every axis"
            raise ValueError(msg)
        count = len(names)
    else:
        count = len(request.pairs)
    if count > request.max_cells:
        msg = (
            f"matrix expansion has {count} cells, "
            f"exceeding max_cells={request.max_cells}"
        )
        raise ValueError(msg)
    return count


def expand(request: MatrixSpec) -> tuple[dict[str, str], ...]:
    """Bound cardinality before allocating the canonical combinations."""
    expansion_size(request)
    axes = request.axes
    if request.pairing == "cross_product":
        combinations = tuple(
            dict(zip(axes, values, strict=True))
            for values in itertools.product(*axes.values())
        )
    elif request.pairing == "zip":
        names = tuple(next(iter(axes.values())))
        combinations = tuple(dict.fromkeys(axes, name) for name in names)
    else:
        combinations = tuple(dict(pair) for pair in request.pairs)
    for choices in combinations:
        if set(choices) != set(axes) or any(
            choice not in axes[name] for name, choice in choices.items()
        ):
            msg = "explicit pair must select one known choice for every axis"
            raise ValueError(msg)
    if len({cell_id(choices) for choices in combinations}) != len(combinations):
        msg = "matrix pairs repeat a cell"
        raise ValueError(msg)
    return combinations
