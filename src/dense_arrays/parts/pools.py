"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/parts/pools.py

Explicit reusable pool sources, separate from mutable input tables.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from dataclasses import dataclass
from pathlib import Path

from dense_arrays._record_validation import digest
from dense_arrays.parts.filters import PartFilter


@dataclass(frozen=True)
class PoolHandle:
    """A location and immutable collection identity, with no open reader."""

    path: Path
    pool_id: str

    def __post_init__(self) -> None:
        """Normalize the location without opening or scanning the pool."""
        object.__setattr__(self, "path", Path(self.path).absolute())
        digest(self.pool_id, field_name="pool_id")


@dataclass(frozen=True)
class PoolSource:
    """Reference an immutable prepared collection and optional part predicate."""

    pool: str | Path | PoolHandle
    select: PartFilter | None = None

    def __post_init__(self) -> None:
        """Freeze explicit paths without opening the collection."""
        if not isinstance(self.pool, PoolHandle):
            object.__setattr__(self, "pool", Path(self.pool))
        if self.select is not None and not isinstance(self.select, PartFilter):
            msg = "pool selection must be PartFilter"
            raise TypeError(msg)

    @property
    def path(self) -> Path:
        """Resolve the caller-supplied location without reading it."""
        return self.pool.path if isinstance(self.pool, PoolHandle) else self.pool
