"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/recovery.py

Exclusive ownership of an existing local run; no lock replacement or repair.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import sqlite3
from contextlib import closing, contextmanager
from typing import TYPE_CHECKING

from dense_arrays.artifacts.store import DATABASE

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path


class RecoveryError(ValueError):
    """A saved run cannot continue under its unchanged execution contract."""

    def __init__(self, code: str, message: str, *, artifact: Path) -> None:
        """Keep machine-readable recovery causes distinct from display text."""
        super().__init__(message)
        self.code = code
        self.artifact = artifact


@contextmanager
def own_run(path: Path) -> Iterator[sqlite3.Connection]:
    """Hold the original lock inode and open only an existing run database."""
    import fcntl  # noqa: PLC0415 - native local Unix lock, acquired before writes

    with (path / ".writer.lock").open("r+b") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as err:
            code = "writer_busy"
            message = "run already has an active writer; retry after it exits"
            raise RecoveryError(code, message, artifact=path) from err
        with closing(
            sqlite3.connect((path / DATABASE).as_uri() + "?mode=rw", uri=True)
        ) as connection:
            connection.execute("PRAGMA synchronous=FULL")
            yield connection
