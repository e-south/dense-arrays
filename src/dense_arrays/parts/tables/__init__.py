"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/tables/__init__.py

Format readers expose headers and typed rows to the shared part validator.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

    from dense_arrays.parts.models import PartTable


@dataclass(frozen=True)
class TableRows:
    """One open snapshot, with projection after shared column resolution."""

    header: tuple[object, ...]
    rows: Callable[[tuple[str, ...]], Iterator[dict[str, object]]]


@contextmanager
def open_table(source: PartTable, data: bytes) -> Iterator[TableRows]:
    """Own parser lifetimes; optional packages are imported only on their route."""
    if source.format in {"csv", "tsv"}:
        from .delimited import open_delimited  # noqa: PLC0415

        reader = open_delimited
    elif source.format == "parquet":
        from .parquet import open_parquet  # noqa: PLC0415

        reader = open_parquet
    else:
        from .excel import open_excel  # noqa: PLC0415

        reader = open_excel
    with reader(source, data) as table:
        yield table
