"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/tables/delimited.py

Strict CSV/TSV rows from captured bytes.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import csv
import io
from contextlib import contextmanager
from typing import TYPE_CHECKING

from . import TableRows

if TYPE_CHECKING:
    from collections.abc import Iterator

    from dense_arrays.parts.models import PartTable


@contextmanager
def open_delimited(source: PartTable, data: bytes) -> Iterator[TableRows]:
    """Retain text values and reject ragged rows, including ignored columns."""
    with io.StringIO(data.decode("utf-8-sig")) as stream:
        reader = csv.DictReader(
            stream, delimiter="," if source.format == "csv" else "\t", strict=True
        )

        def rows(columns: tuple[str, ...]) -> Iterator[dict[str, object]]:
            for number, row in enumerate(reader, 1):
                if None in row or any(value is None for value in row.values()):
                    msg = (
                        f"{source.table}: row {number}: "
                        "expected one value per header column"
                    )
                    raise ValueError(msg)
                yield {name: row[name] for name in columns}

        try:
            yield TableRows(tuple(reader.fieldnames or ()), rows)
        except csv.Error as err:
            msg = f"{source.table}: invalid {source.format} table: {err}"
            raise ValueError(msg) from err
