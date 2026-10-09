"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/reporting/collections/filtering.py

Resolve unqualified selectors across a bounded set of source namespaces.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.errors import InvalidQueryError
from dense_arrays.reporting.design_filters import DesignFilter

from .sources import aliases, annotations, cells

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.artifacts.reading import ReadBudget

    from .sources import Annotations, SourceView


def resolve_filter(
    sources: tuple[SourceView, ...], selected: DesignFilter | None, budget: ReadBudget
) -> DesignFilter | None:
    """Reject missing and ambiguous labels before emitting the first record."""
    if selected is None:
        return None
    names = ("design_ids", "cells", "part_ids", "groups")
    offered = {name: {v: set() for v in getattr(selected, name)} for name in names}
    budget.retain(selected.identities)

    def offer(name: str, local: str, full: str) -> None:
        for label in {local, full} & offered[name].keys():
            matches = offered[name][label]
            if full not in matches:
                budget.retain()
                matches.add(full)

    for source in sources:
        for cell in cells(source):
            offer("cells", cell.split("/")[1], cell)
        if selected.design_ids:
            for local, full in aliases(source, selected.design_ids):
                offer("design_ids", local, full)
        with annotations(
            source, budget, needed=bool(selected.part_ids or selected.groups)
        ) as lookup:
            for name, local, full in _part_labels(lookup):
                offer(name, local, full)
    resolved = _unambiguous(offered)
    # The temporary alias index is replaced by the resolved predicate.
    budget.identities -= selected.identities + sum(
        len(matches) for labels in offered.values() for matches in labels.values()
    )
    result = DesignFilter(**resolved, metrics=selected.metrics)
    budget.retain(result.identities)
    return result


def _part_labels(lookup: Annotations) -> Iterator[tuple[str, str, str]]:
    """Keep plan-scoped part identities distinct from shared group labels."""
    for parts, collection in lookup.values():
        for part in parts.values():
            yield "part_ids", part.part_id, f"{collection}/{part.part_id}"
            if part.group is not None:
                yield "groups", part.group, part.group


def _unambiguous(offered: dict[str, dict[str, set[str]]]) -> dict[str, tuple[str, ...]]:
    """Resolve each requested label to exactly one fully qualified identity."""
    resolved = {}
    for name, labels in offered.items():
        values = set()
        for label, matches in labels.items():
            if not matches:
                msg = f"unknown {name}: {label!r}"
                raise InvalidQueryError(msg)
            if len(matches) > 1:
                msg = f"ambiguous {name}: {label!r}; use a full reference"
                raise InvalidQueryError(msg)
            values.update(matches)
        resolved[name] = tuple(sorted(values))
    return resolved
