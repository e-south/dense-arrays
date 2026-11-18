"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/playback/json_io.py

Reject duplicate keys consistently at playback JSON input boundaries.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import json


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result = {}
    for key, value in pairs:
        if key in result:
            msg = f"duplicate JSON object key: {key!r}"
            raise ValueError(msg)
        result[key] = value
    return result


def loads_json(payload: str) -> object:
    """Decode JSON without silently replacing duplicate object fields."""
    return json.loads(payload, object_pairs_hook=_unique_object)
