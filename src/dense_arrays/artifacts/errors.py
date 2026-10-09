"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/errors.py

Typed failures at persisted-evidence boundaries.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from contextlib import contextmanager
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path


class ArtifactIntegrityError(ValueError):
    """Stored evidence is invalid; retrying generation cannot repair it."""

    def __init__(self, message: str, *, artifact: Path | None = None) -> None:
        """Retain the affected artifact independently of human diagnostic text."""
        super().__init__(message)
        self.artifact = artifact


class InvalidQueryError(ValueError):
    """The caller's predicate cannot be resolved against the chosen snapshot."""


@contextmanager
def integrity_boundary(path: Path | None) -> Iterator[None]:
    """Classify invalid persisted data while preserving the original cause."""
    try:
        yield
    except (ArtifactIntegrityError, InvalidQueryError):
        raise
    except (ValueError, TypeError, KeyError) as err:
        raise ArtifactIntegrityError(str(err), artifact=path) from err
