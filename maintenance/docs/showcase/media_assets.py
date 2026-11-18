"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/maintenance/docs/showcase/media_assets.py

Encode README delivery media from the exact rendered packing movie.

Module Author(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET

from svg_assets import (
    CARD_PADDING,
    CARD_WIDTH,
    PACKING_HEIGHT,
    PACKING_WIDTH,
    SVG,
    card,
)


def _png(executable: str, source: Path, target: Path, width: int) -> None:
    subprocess.run(  # noqa: S603 - resolved executable and generated artwork
        [
            executable,
            str(source),
            "--export-width=" + str(width),
            "--export-filename=" + str(target),
        ],
        check=True,
        capture_output=True,
        timeout=60,
    )


def readme_media(output: Path, *, inkscape: str, ffmpeg: str, fps: int) -> None:
    """Keep rounded artwork, loop timing and the source movie's complete content."""
    _png(inkscape, output / "packing.svg", output / "packing-card.png", 2880)
    with TemporaryDirectory(prefix="dense-arrays-readme-") as temporary:
        stage = Path(temporary)
        background = ET.Element(
            f"{{{SVG}}}svg", {"viewBox": f"0 0 {PACKING_WIDTH} {PACKING_HEIGHT}"}
        )
        (stage / "card.svg").write_bytes(
            card(background, PACKING_WIDTH, PACKING_HEIGHT)
        )
        _png(inkscape, stage / "card.svg", stage / "card.png", CARD_WIDTH)
        inner = CARD_WIDTH - 2 * CARD_PADDING
        filters = (
            f"[1:v]scale={inner}:-1:flags=lanczos[clip];"
            f"[0:v][clip]overlay={CARD_PADDING}:{CARD_PADDING}:shortest=1:format=auto[loop]"
        )
        subprocess.run(  # noqa: S603 - resolved executable and generated inputs
            [
                ffmpeg,
                "-v",
                "error",
                "-loop",
                "1",
                "-framerate",
                str(fps),
                "-i",
                str(stage / "card.png"),
                "-i",
                str(output / "packing.mp4"),
                "-filter_complex",
                filters,
                "-map",
                "[loop]",
                "-c:v",
                "libwebp_anim",
                "-quality",
                "90",
                "-compression_level",
                "4",
                "-loop",
                "0",
                str(output / "packing-loop.webp"),
            ],
            check=True,
            capture_output=True,
            timeout=120,
        )
        if (output / "packing-loop.webp").stat().st_size > 1024 * 1024:
            msg = "README animation exceeds the 1 MiB delivery budget"
            raise ValueError(msg)
