"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/inspection.py

Resolve collection queries and independently verify their supplied geometry.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import TYPE_CHECKING

from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.parts import PartFilter

from .models import DATABASE, ArrayFilter, CollectionSummary
from .reading import CollectionView
from .storage import file_digest, read_summary, reader

if TYPE_CHECKING:
    from dense_arrays.artifacts.cursors import Cursor
    from dense_arrays.artifacts.reading import ReadLimits


def inspect_collection(  # noqa: PLR0913 - resolved public query options
    path: Path | str,
    *,
    view: str,
    verify: bool,
    limit: int | None,
    selected: object,
    all_rows: bool,
    read_limits: ReadLimits,
    after: Cursor | None,
    compare: object,
) -> CollectionSummary | CollectionView:
    """Read supplied geometry and metadata through bounded collection views."""
    if compare is not None:
        msg = "supplied-array collections do not contain comparable generation plans"
        raise ValueError(msg)
    if view not in {"summary", "arrays", "parts", "sequences", "placements"}:
        msg = (
            "array collection views: summary, arrays, parts, sequences, "
            "placements; execution evidence is unavailable"
        )
        raise ValueError(msg)
    if selected is not None and view == "summary":
        msg = "filters apply to array records, not collection summary"
        raise ValueError(msg)
    if isinstance(selected, PartFilter) and view == "parts":
        if selected.metrics:
            msg = "collection parts support identity/group filters, not metric filters"
            raise ValueError(msg)
        selected = ArrayFilter(part_ids=selected.part_ids, groups=selected.groups)
    if selected is not None and not isinstance(selected, ArrayFilter):
        msg = "array collection record queries require ArrayFilter"
        raise TypeError(msg)
    if view == "parts" and selected is not None and selected.array_ids:
        msg = (
            "parts view covers the catalog; array filters require an array record view"
        )
        raise ValueError(msg)
    path = Path(path).absolute()
    summary = read_summary(path, read_limits)
    if verify:
        summary = _verify(path, summary, read_limits)
    if view == "summary":
        return summary
    return CollectionView(
        path,
        summary,
        view,
        None if all_rows else (100 if limit is None else limit),
        selected or ArrayFilter(),
        read_limits,
        after,
    )


def _verify(
    path: Path, summary: CollectionSummary, limits: ReadLimits
) -> CollectionSummary:
    checksum = summary.manifest["database"]["sha256"]
    if file_digest(path / DATABASE) != checksum:
        msg = "collection database checksum differs"
        raise ArtifactIntegrityError(msg, artifact=path)
    count, placements = 0, 0
    with CollectionView(
        path, summary, "arrays", limit=None, read_limits=limits
    ).records() as rows:
        for record in rows:
            count += 1
            placements += len(record.realized.placements)
    with reader(path) as connection:
        sequences = connection.execute(
            "SELECT COUNT(DISTINCT sequence_id) FROM arrays"
        ).fetchone()[0]
    if (count, placements, sequences) != (
        summary.arrays,
        summary.placements,
        summary.manifest["sequences"],
    ):
        msg = "array collection counts differ from stored records"
        raise ArtifactIntegrityError(msg, artifact=path)
    if file_digest(path / DATABASE) != checksum:
        msg = "collection database changed during verification"
        raise ArtifactIntegrityError(msg, artifact=path)
    return replace(summary, verified=True)
