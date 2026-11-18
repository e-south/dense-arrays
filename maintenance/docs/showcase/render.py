"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/maintenance/docs/showcase/render.py

Rebuild the published sequence demonstrations through public producer interfaces.

Module Author(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import shutil
from dataclasses import dataclass, field
from importlib.metadata import distribution
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
from dnadesign import baserender
from PIL import Image, ImageChops
from svg_assets import player_css, poster_card, workflow_card

from dense_arrays.playback import (
    PlaybackDocument,
    PlaybackPresentation,
    loads_realized_array,
    reconstruct_playback,
)
from dense_arrays.playback.matplotlib_renderer import (
    render_collection_mp4,
    render_collection_poster_png,
)

ROOT = Path(__file__).resolve().parent
MAX_NATIVE_FRAMES = 32


def digest(path: Path) -> str:
    """Identify selected inputs and delivered bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def producer_revision(name: str) -> str:
    """Record the immutable Git installation used by the locked environment."""
    value = json.loads(distribution(name).read_text("direct_url.json") or "{}")
    revision = value.get("vcs_info", {}).get("commit_id")
    if not revision:
        msg = "Run the recipe with its locked Git-backed environment"
        raise ValueError(msg)
    return revision


def load_inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Admit the source identities before preparing any derived frame."""
    config = json.loads((ROOT / "rendering.json").read_bytes())
    records = json.loads((ROOT / "inputs/records.json").read_bytes())
    arrays = {
        r["source_id"]: loads_realized_array(json.dumps(r)) for r in records["arrays"]
    }
    drawings = {r["id"]: r for r in records["drawing_records"]}
    selected = [key for group in config["groups"].values() for key in group]
    if set(selected) != arrays.keys() or arrays.keys() != drawings.keys():
        msg = "Realized and drawing records must have the same explicit identities"
        raise ValueError(msg)
    for key, array in arrays.items():
        drawing = drawings[key]
        expected = {p.placement_id: p for p in array.placements}
        if (
            drawing["sequence"] != array.sequence
            or {f["id"] for f in drawing["features"]} != expected.keys()
        ):
            msg = "Drawing sequence and features differ from their source record"
            raise ValueError(msg)
        for feature in drawing["features"]:
            source = expected[feature["id"]]
            if (
                feature["span"]["start"],
                feature["span"]["end"],
                feature["span"]["strand"],
            ) != (source.start, source.end, source.orientation.value):
                msg = "Drawing coordinates differ from their source placement"
                raise ValueError(msg)
    return config, arrays, drawings


def native_job(
    config: dict[str, Any], kind: str, *, svg: bool = False
) -> dict[str, Any]:
    """Create ordinary BaseRender jobs without interpreting scientific data."""
    bundle = kind + ("-vectors" if svg else "")
    job = {
        "version": 4,
        "contract": {"kind": "sequence_rows_render_v3"},
        "bundle": {"path": bundle},
        "input": {
            "kind": "json",
            "path": kind + "-records.json",
            "alphabet": "DNA",
            "adapter": {
                "kind": "generic_features",
                "columns": {
                    "id": "id",
                    "sequence": "sequence",
                    "features": "features",
                    "effects": "effects",
                },
                "policies": {},
            },
        },
        "render": {
            "renderer": "sequence_rows",
            "style": {"preset": None, "overrides": config["style"]},
            "options": {
                "layout": config["layouts"]["100" if kind == "packing" else "60"]
            },
            "native_svg": config["native_svg"],
        },
        "outputs": [{"kind": "images", "fmt": "svg" if svg else "png"}],
        "run": {"strict": True, "fail_on_skips": True},
    }
    if kind == "library" and not svg:
        job["outputs"].append(
            {
                "kind": "video",
                "path": "library.mp4",
                "fps": config["timing"]["fps"],
                "frames_per_record": round(
                    config["timing"]["fps"] * config["timing"]["library_hold_seconds"]
                ),
                "background": "#FFFFFF",
                "encode_timeout_s": 120,
            }
        )
    return job


def prepare_native(
    root: Path, config: dict[str, Any], arrays: dict[str, Any], drawings: dict[str, Any]
) -> dict[str, list[str]]:
    """Derive color-only states, retaining every coordinate and source identity."""
    packing = []
    orders = {}
    for key in config["groups"]["packing"]:
        order = [step.placement_id for step in reconstruct_playback(arrays[key]).steps]
        orders[key] = []
        for count in range(len(order) + 1):
            record = copy.deepcopy(drawings[key])
            record["id"] = key + "--" + str(count)
            orders[key].append(record["id"])
            for feature in record["features"]:
                if feature["id"] not in order[:count]:
                    feature["tags"] = ["pending"]
            packing.append(record)
    if sum(map(len, orders.values())) > MAX_NATIVE_FRAMES:
        msg = "This documentation recipe admits at most 32 native frames"
        raise ValueError(msg)
    (root / "packing-records.json").write_text(json.dumps(packing))
    (root / "library-records.json").write_text(
        json.dumps([drawings[k] for k in config["groups"]["library"]])
    )
    for kind, svg in (("packing", False), ("library", False), ("library", True)):
        baserender.run_job(native_job(config, kind, svg=svg), caller_root=root)
    return orders


@dataclass
class Frames:
    """Provide one bounded batch of producer-owned sequence images."""

    root: Path
    orders: dict[str, list[str]]
    native_nucleotide_cap_height_px: float
    preferred_figure_height_inches: float = 2.4
    renders_distance_brackets: bool = True
    images: dict[str, np.ndarray] = field(default_factory=dict)

    def frame(self, document: PlaybackDocument, step: int | None) -> np.ndarray:
        """Read each selected image once while graph playback reuses it."""
        name = self.orders[document.plan.source_id][0 if step is None else step + 1]
        if name not in self.images:
            with Image.open(self.root / (name + ".png")) as image:
                self.images[name] = np.asarray(image.convert("RGBA")).copy()
                self.images[name].setflags(write=False)
        return self.images[name]


def margin_clip(root: Path) -> str:
    """Retain all source ink with extra encoder clearance, clipping only margins."""
    boxes = []
    for path in root.glob("*.png"):
        raw = Image.open(path).convert("RGBA")
        paper = Image.new("RGBA", raw.size, "white")
        paper.alpha_composite(raw)
        boxes.append(
            ImageChops.difference(
                paper.convert("RGB"), Image.new("RGB", raw.size, "white")
            ).getbbox()
        )
    left, top = (
        max(0, min(b[0] for b in boxes) - 8),
        max(0, min(b[1] for b in boxes) - 8),
    )
    right, bottom = (
        min(raw.width, max(b[2] for b in boxes) + 8),
        min(raw.height, max(b[3] for b in boxes) + 8),
    )
    values = (
        top / raw.height * 100,
        (raw.width - right) / raw.width * 100,
        (raw.height - bottom) / raw.height * 100,
        left / raw.width * 100,
    )
    return "inset(" + " ".join(f"{value:.10f}%" for value in values) + ")"


def render(output: Path) -> None:
    """Render a fresh candidate and retain compact input/output provenance."""
    if output.exists():
        msg = "Use a fresh output directory"
        raise ValueError(msg)
    config, arrays, drawings = load_inputs()
    inputs = {
        str(p.relative_to(ROOT)): digest(p)
        for p in (ROOT / "inputs").iterdir()
        if p.is_file()
    }
    inputs["rendering.json"] = digest(ROOT / "rendering.json")
    inputs.update(
        {
            name: digest(ROOT / name)
            for name in ("render.py", "svg_assets.py", "pyproject.toml", "uv.lock")
        }
    )
    runtime = {
        name: producer_revision(name) for name in ("dense-arrays", "dnadesign-tools")
    }
    with TemporaryDirectory(prefix="dense-arrays-docs-") as temporary:
        stage = Path(temporary)
        orders = prepare_native(stage, config, arrays, drawings)
        dpi = config["style"]["dpi"]
        cap = (
            baserender.measure_text(
                "ACGT",
                family=config["style"]["font_mono"],
                size_pt=config["style"]["font_size_seq"],
            ).ink_height_pt
            * dpi
            / 72
        )
        provider = Frames(stage / "packing/images", orders, cap)
        settings = PlaybackPresentation.from_mapping(config["presentation"])
        documents = tuple(
            PlaybackDocument(
                reconstruct_playback(arrays[key]),
                "Recorded promoter realization",
                label_overrides={
                    p.placement_id: p.label for p in arrays[key].placements
                },
                color_overrides={
                    p.placement_id: config["style"]["palette"][p.label]
                    for p in arrays[key].placements
                },
                presentation=settings,
            )
            for key in config["groups"]["packing"]
        )
        output.mkdir(parents=True)
        render_collection_poster_png(
            documents,
            output / "packing.png",
            dpi=dpi,
            duplex_frame_renderer=provider.frame,
        )
        render_collection_mp4(
            documents,
            output / "packing.mp4",
            duplex_frame_renderer=provider.frame,
            **{
                key: config["timing"][key]
                for key in (
                    "fps",
                    "seconds_per_step",
                    "lead_seconds",
                    "hold_seconds",
                    "scene_transition_seconds",
                )
            },
        )
        (output / "packing.svg").write_bytes(
            poster_card((output / "packing.png").read_bytes())
        )
        context = (ROOT / "inputs/workflow-context.svg").read_bytes()
        native = (
            stage / "library-vectors/images" / (config["groups"]["library"][0] + ".svg")
        ).read_bytes()
        (output / "library.svg").write_bytes(
            workflow_card(context, native, config["workflow"]["array_bounds_css_px"])
        )
        shutil.copyfile(stage / "library/library.mp4", output / "library.mp4")
        (output / "player.css").write_text(
            player_css(
                context, config["workflow"], margin_clip(stage / "library/images")
            )
        )
        if any(digest(ROOT / name) != value for name, value in inputs.items()):
            msg = "Documentation inputs changed while rendering"
            raise ValueError(msg)
        receipt = {
            "schema": "dense-arrays.documentation-showcase/v1",
            "inputs": inputs,
            "producers": runtime,
            "record_ids": config["groups"],
            "authority": "placement_reconstructed",
            "timing": config["timing"],
            "outputs": {p.name: digest(p) for p in output.iterdir() if p.is_file()},
        }
        (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


def main() -> None:
    """Expose one reproduction command with an explicit fresh destination."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.output)
    print(args.output / "receipt.json")  # noqa: T201 - CLI result


if __name__ == "__main__":
    main()
