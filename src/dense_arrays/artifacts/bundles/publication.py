"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/bundles/publication.py

Exclusive bundle publication with a manifest as the final commit marker.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import os
import shutil
import sqlite3
from contextlib import closing, contextmanager
from typing import TYPE_CHECKING

from dense_arrays.artifacts.bundles.models import BUNDLE_DATABASE

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path


@contextmanager
def destination(path: Path) -> Iterator[sqlite3.Connection]:
    """Reserve a new directory; an interrupted writer has no valid bundle commit."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.mkdir(exist_ok=False)
    owned = path.stat()
    marker = path / ".bundle-pending"
    marker.touch(exist_ok=False)
    try:
        with closing(sqlite3.connect(path / BUNDLE_DATABASE)) as connection:
            connection.execute("PRAGMA synchronous=FULL")
            connection.executescript("""
                CREATE TABLE plans(plan_id TEXT PRIMARY KEY, payload TEXT NOT NULL,
                    digest TEXT NOT NULL);
                CREATE TABLE batches(run_id TEXT NOT NULL, cell_id TEXT NOT NULL,
                    batch_index INTEGER NOT NULL, batch_id TEXT NOT NULL,
                    payload TEXT NOT NULL, digest TEXT NOT NULL,
                    PRIMARY KEY(run_id,cell_id,batch_index));
                CREATE TABLE designs(ordinal INTEGER PRIMARY KEY,
                    design_ref TEXT UNIQUE NOT NULL, run_id TEXT NOT NULL,
                    cell_id TEXT NOT NULL, local_id TEXT NOT NULL,
                    plan_id TEXT NOT NULL, payload TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE INDEX design_alias ON designs(local_id);
            """)
            with connection:
                yield connection
        marker.unlink()
        with _directory_descriptor(path) as descriptor:
            os.fsync(descriptor)
    except BaseException:
        if path.exists() and (path.stat().st_dev, path.stat().st_ino) == (
            owned.st_dev,
            owned.st_ino,
        ):
            shutil.rmtree(path)
        raise


@contextmanager
def _directory_descriptor(path: Path) -> Iterator[int]:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        yield descriptor
    finally:
        os.close(descriptor)
