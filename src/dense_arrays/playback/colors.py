"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/colors.py

Validate and blend portable RGB colors without renderer dependencies.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import re


def validate_color(color: str) -> None:
    """Require a portable opaque RGB color for raster output."""
    if not isinstance(color, str) or re.fullmatch(r"#[0-9a-fA-F]{6}", color) is None:
        msg = "colors must use #RRGGBB notation"
        raise ValueError(msg)


def blend_color(start: str, end: str, progress: float) -> str:
    """Interpolate opaque colors without changing artist geometry."""
    channels = (
        round(
            int(start[index : index + 2], 16) * (1 - progress)
            + int(end[index : index + 2], 16) * progress
        )
        for index in (1, 3, 5)
    )
    return "#" + "".join(f"{channel:02x}" for channel in channels)
