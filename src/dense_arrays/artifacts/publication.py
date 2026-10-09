"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/artifacts/publication.py

Atomic create-only publication for explicit standalone artifact files.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

import os
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO


def write_new(path: Path, text: str) -> None:
    """Publish a complete file without replacing an existing destination."""
    with new_text_file(path) as stream:
        stream.write(text)


@contextmanager
def new_text_file(path: Path) -> Iterator[TextIO]:
    """Stream to owned staging and publish only after the entire write succeeds."""
    if path.exists() or path.is_symlink():
        msg = f"output destination already exists: {path}"
        raise FileExistsError(msg)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
        ) as stream:
            temporary = Path(stream.name)
            yield stream
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink()
