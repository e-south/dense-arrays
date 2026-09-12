"""Stage requested exports before publishing files to distinct output paths.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import os
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping


def _same_file(left: Path, right: Path) -> bool:
    return left == right or (left.exists() and right.exists() and left.samefile(right))


def _validate_outputs(
    source: Path, outputs: Mapping[str, Path], *, replace: bool
) -> dict[str, Path]:
    targets: dict[str, Path] = {}
    for name, path in outputs.items():
        if path.is_symlink():
            msg = f"output path must not be a symlink: {path}"
            raise ValueError(msg)
        target = path.resolve()
        if _same_file(source, target):
            msg = f"output path aliases the input: {path}"
            raise ValueError(msg)
        if any(_same_file(target, prior) for prior in targets.values()):
            msg = f"output paths must be distinct: {path}"
            raise ValueError(msg)
        if target.exists():
            if not target.is_file():
                msg = f"output path is not a regular file: {path}"
                raise ValueError(msg)
            if not replace:
                msg = f"output already exists: {path}; pass --replace to overwrite"
                raise FileExistsError(msg)
        for parent in target.parents:
            if parent.exists() and not parent.is_dir():
                msg = f"output parent is not a directory: {parent}"
                raise ValueError(msg)
        targets[name] = target
    return targets


def publish_exports(
    source: Path,
    outputs: Mapping[str, Path],
    render: Callable[[Mapping[str, Path]], None],
    *,
    replace: bool,
) -> tuple[Path, ...]:
    """Render the full requested set before publishing individual files.

    Rendering failures leave existing files untouched. Publication is atomic
    per file; a filesystem failure reports exactly which files were published.
    """
    targets = _validate_outputs(source.resolve(), outputs, replace=replace)
    with ExitStack() as stack:
        staged: dict[str, Path] = {}
        for name, target in targets.items():
            target.parent.mkdir(parents=True, exist_ok=True)
            folder = stack.enter_context(
                TemporaryDirectory(prefix=".dense-arrays-", dir=target.parent)
            )
            staged[name] = Path(folder) / name
        render(staged)
        if any(not path.is_file() for path in staged.values()):
            msg = "renderer did not produce every requested export"
            raise RuntimeError(msg)
        published: list[Path] = []
        try:
            for name, target in targets.items():
                if replace:
                    staged[name].replace(target)
                else:
                    os.link(staged[name], target)
                published.append(target)
        except OSError as exc:
            completed = ", ".join(str(path) for path in published) or "none"
            msg = f"export publication failed: {exc}; published files: {completed}"
            raise OSError(msg) from exc
    return tuple(published)
