"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/arrays/rendering.py

Render supplied placement geometry through the shared playback presentation.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

from dense_arrays.artifacts.receipts import ExportReceipt
from dense_arrays.playback.presentation import PlaybackDocument
from dense_arrays.playback.reconstruction import reconstruct_playback
from dense_arrays.reporting.rendering import publish_document

if TYPE_CHECKING:
    from pathlib import Path

    from .reading import CollectionView


def render_array(query: CollectionView, out: Path) -> ExportReceipt:
    """Require exactly one array and show only its supplied placement evidence."""
    if out.exists() or out.is_symlink():
        msg = f"output destination already exists: {out}"
        raise FileExistsError(msg)
    if out.suffix.lower() != ".png":
        msg = "array rendering requires a .png output"
        raise ValueError(msg)
    with query.records() as rows:
        record = next(rows, None)
        second = next(rows, None)
    if record is None or second is not None:
        msg = "array rendering requires exactly one selected array; use an array ID"
        raise ValueError(msg)
    document = PlaybackDocument(
        reconstruct_playback(record.realized),
        title=f"Dense array {record.array_id}",
        label_overrides={
            p.placement_id: p.feature_id for p in record.realized.placements
        },
    )
    rendered = publish_document(document, out)
    return ExportReceipt(str(out), "png", "array", 1, query.sources, (), rendered)
