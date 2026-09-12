"""Render saved placements with explicit input errors and export destinations.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import json
import shutil
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Annotated

import typer

from ..realized import REALIZED_ARRAY_SCHEMA_VERSION
from .models import PLAYBACK_PLAN_SCHEMA_VERSION, PlaybackPlan
from .output import publish_exports
from .presentation import PlaybackDocument
from .reconstruction import reconstruct_playback
from .serialization import loads_playback_plan, loads_realized_array

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_plan(input_path: Path) -> PlaybackPlan:
    payload = input_path.read_text(encoding="utf-8")
    root = json.loads(payload)
    if not isinstance(root, dict):
        msg = "playback input must be a JSON object"
        raise TypeError(msg)
    schema = root.get("schema_version")
    if schema == REALIZED_ARRAY_SCHEMA_VERSION:
        return reconstruct_playback(loads_realized_array(payload))
    if schema == PLAYBACK_PLAN_SCHEMA_VERSION:
        return loads_playback_plan(payload)
    msg = f"unsupported playback input schema: {schema!r}"
    raise ValueError(msg)


def _render_exports(document: PlaybackDocument, paths: Mapping[str, Path]) -> None:
    if "playback.mp4" in paths and shutil.which("ffmpeg") is None:
        msg = "MP4 export requires an FFmpeg executable on PATH"
        raise ValueError(msg)
    try:
        from .matplotlib_renderer import (
            render_collection_gif,
            render_collection_mp4,
            render_collection_poster_png,
        )

        if "poster.png" in paths:
            render_collection_poster_png((document,), paths["poster.png"])
        if "playback.mp4" in paths:
            render_collection_mp4((document,), paths["playback.mp4"])
        if "playback.gif" in paths:
            render_collection_gif((document,), paths["playback.gif"])
    except ImportError as exc:
        msg = (
            "Media export requires playback dependencies; "
            "run `uv sync --frozen --extra playback` in the checkout"
        )
        raise ValueError(msg) from exc


@app.command()
def render(
    input_path: Annotated[
        Path, typer.Argument(help="RealizedArray or PlaybackPlan JSON")
    ],
    title: Annotated[
        str, typer.Option(help="Publication title")
    ] = "Dense-array solution playback",
    subtitle: Annotated[str, typer.Option(help="Publication subtitle")] = "",
    poster_out: Annotated[
        Path | None, typer.Option("--poster", help="Optional poster PNG")
    ] = None,
    mp4_out: Annotated[
        Path | None, typer.Option("--mp4", help="Optional MP4 output")
    ] = None,
    gif_out: Annotated[
        Path | None, typer.Option("--gif", help="Optional animated GIF output")
    ] = None,
    *,
    replace: Annotated[
        bool, typer.Option("--replace", help="Replace existing output files")
    ] = False,
) -> None:
    """Validate saved placements, render all requested formats, and publish them."""
    outputs = {
        name: path
        for name, path in (
            ("poster.png", poster_out),
            ("playback.mp4", mp4_out),
            ("playback.gif", gif_out),
        )
        if path is not None
    }
    if not outputs:
        typer.echo(
            "Error: Choose at least one export: --poster, --mp4, or --gif", err=True
        )
        raise typer.Exit(code=1)
    try:
        plan = _load_plan(input_path)
        document = PlaybackDocument(plan=plan, title=title, subtitle=subtitle)
        published = publish_exports(
            input_path,
            outputs,
            lambda paths: _render_exports(document, paths),
            replace=replace,
        )
    except (
        OSError,
        ValueError,
        TypeError,
        RuntimeError,
        subprocess.SubprocessError,
    ) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    for path in published:
        typer.echo(f"Wrote {path}")


if __name__ == "__main__":
    app()
