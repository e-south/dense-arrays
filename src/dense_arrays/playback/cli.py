"""Standalone CLI for rendering serialized dense-array playback artifacts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer

from .html import PlaybackDocument, render_playback_html
from .reconstruction import reconstruct_playback
from .serialization import loads_playback_plan, loads_realized_array

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.command()
def render(
    input_path: Annotated[
        Path, typer.Argument(help="RealizedArray or PlaybackPlan JSON")
    ],
    html_out: Annotated[
        Path, typer.Option("--html", help="Self-contained HTML output")
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
) -> None:
    """Render a strict serialized contract without importing optimizer internals."""
    payload = input_path.read_text(encoding="utf-8")
    schema = json.loads(payload).get("schema_version")
    if schema == "dense_arrays.realized_array.v1":
        plan = reconstruct_playback(loads_realized_array(payload))
    elif schema == "dense_arrays.playback_plan.v1":
        plan = loads_playback_plan(payload)
    else:
        msg = f"unsupported playback input schema: {schema!r}"
        raise typer.BadParameter(msg, param_hint="input_path")
    document = PlaybackDocument(plan=plan, title=title, subtitle=subtitle)
    html_out.parent.mkdir(parents=True, exist_ok=True)
    html_out.write_text(
        render_playback_html(plan, title=title, subtitle=subtitle),
        encoding="utf-8",
    )
    if poster_out is not None or mp4_out is not None or gif_out is not None:
        try:
            from .matplotlib_renderer import (
                render_collection_gif,
                render_collection_mp4,
                render_collection_poster_png,
            )
        except ImportError as exc:
            msg = (
                "Raster playback export requires the optional playback dependencies. "
                "Install them with `pip install 'dense-arrays[playback]'`."
            )
            raise RuntimeError(msg) from exc
    if poster_out is not None:
        render_collection_poster_png((document,), poster_out)
    if mp4_out is not None:
        render_collection_mp4((document,), mp4_out)
    if gif_out is not None:
        render_collection_gif((document,), gif_out)


if __name__ == "__main__":
    app()
