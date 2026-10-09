"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/tables/parquet.py

Optional Parquet projection with bounded decoded batches.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import io
from contextlib import contextmanager
from typing import TYPE_CHECKING

from . import TableRows

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.parts.models import PartTable


@contextmanager
def open_parquet(source: PartTable, data: bytes) -> Iterator[TableRows]:
    """Read only mapped top-level columns, preserving row-group order and nulls."""
    try:
        import pyarrow.parquet as pq  # noqa: PLC0415 - optional input dependency
    except ImportError as err:
        msg = (
            f"{source.table}: Parquet input requires pyarrow; "
            "install 'dense-arrays[tables]'"
        )
        raise ValueError(msg) from err
    with io.BytesIO(data) as stream, pq.ParquetFile(stream) as reader:

        def rows(columns: tuple[str, ...]) -> Iterator[dict[str, object]]:
            for batch in reader.iter_batches(
                batch_size=1024, columns=list(columns), use_threads=False
            ):
                yield from batch.to_pylist()

        yield TableRows(tuple(reader.schema_arrow.names), rows)
