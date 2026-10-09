"""Selection requests specify allocation and randomness without hidden defaults.

Author: Eric J. South.
"""

import pytest

from dense_arrays import reporting


def test_selection_request_round_trip_and_strict_allocations():
    allocations = {"run/a": 6, "run/b": 0}
    request = reporting.LibrarySelection(
        filter=reporting.DesignFilter(groups=("A",)),
        take=reporting.Take(per_cell=allocations, policy="random", seed=23),
    )
    allocations["run/a"] = 999
    value = request.to_dict()
    assert value == {
        "schema": "dense_arrays.library-selection.v1",
        "filter": {
            "design_ids": [],
            "cells": [],
            "part_ids": [],
            "groups": ["A"],
            "metrics": {},
        },
        "take": {
            "per_cell": {"run/a": 6, "run/b": 0},
            "policy": "random",
            "seed": 23,
            "shortfall": "error",
        },
    }
    assert reporting.LibrarySelection.from_dict(value) == request
    assert reporting.Take(count=0).count == 0
    assert reporting.LibrarySelection().take is None
    for invalid in (
        {},
        {"count": 1, "per_cell": {}},
        {"count": True},
        {"count": -1},
        {"per_cell": {"a": 1.5}},
        {"per_cell": {"": 1}},
        {"count": 2, "policy": "random"},
        {"count": 2, "seed": 23},
        {"count": 2, "policy": "weighted"},
        {"count": 2, "shortfall": "redistribute"},
    ):
        with pytest.raises((ValueError, TypeError)):
            reporting.Take(**invalid)
    with pytest.raises(ValueError, match="unknown"):
        reporting.LibrarySelection.from_dict({**value, "limit": 6})
    with pytest.raises(ValueError, match="schema"):
        reporting.LibrarySelection.from_dict({**value, "schema": "unknown"})
