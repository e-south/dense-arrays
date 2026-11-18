"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/maintenance/docs/showcase/svg_assets.py

Compose closed vector artwork and measured media windows for documentation.

Module Author(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import base64
import copy
import itertools
import re
from typing import TYPE_CHECKING, Any
from xml.etree import ElementTree as ET

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"
ET.register_namespace("", SVG)


def _editor_metadata(root: ET.Element) -> None:
    for node in root.iter():
        if node.tag == f"{{{SVG}}}script":
            msg = "Showcase SVG must contain no scripts"
            raise ValueError(msg)
        for key in list(node.attrib):
            if key.startswith(("{http://www.inkscape.org/", "{http://sodipodi.")):
                del node.attrib[key]


def _inline_images(root: ET.Element, visit: Callable[[bytes], ET.Element]) -> None:
    for parent in list(root.iter()):
        for child in list(parent):
            ref = child.get("href", child.get(f"{{{XLINK}}}href", ""))
            if ref.startswith("data:image/svg+xml;base64,"):
                inner = visit(base64.b64decode(ref.split(",", 1)[1]))
                for key in ("x", "y", "width", "height", "transform"):
                    if key in child.attrib:
                        inner.set(key, child.attrib[key])
                index = list(parent).index(child)
                parent.remove(child)
                parent.insert(index, inner)
            elif ref and not ref.startswith(("#", "data:image/png;base64,")):
                msg = "Showcase artwork must be closed, without external resources"
                raise ValueError(msg)


def _rename_ids(root: ET.Element, prefix: str, counter: Iterator[int]) -> None:
    mapping = {}
    for node in root.iter():
        old = node.get("id")
        if old:
            names = {
                "aligned-heading": "array-heading",
                "array-heading": "array-heading",
                "aligned-takeaway": "array-takeaway",
                "array-takeaway": "array-takeaway",
            }
            mapping[old] = names.get(old, f"{prefix}{next(counter)}")

    def references(value: str) -> str:
        return re.sub(
            r"url\(#([\w.:-]+)\)",
            lambda m: "url(#" + mapping.get(m[1], m[1]) + ")",
            value,
        )

    for node in root.iter():
        old = node.get("id")
        if old:
            node.set("id", mapping[old])
        for key, value in list(node.attrib.items()):
            if key in ("href", f"{{{XLINK}}}href") and value.startswith("#"):
                node.set(key, "#" + mapping.get(value[1:], value[1:]))
            else:
                node.set(key, references(value))
        if node.tag == f"{{{SVG}}}style" and node.text:
            node.text = references(node.text)


def closed_svg(data: bytes, prefix: str) -> ET.Element:
    """Inline SVG images, namespace IDs and discard editor-only path metadata."""
    counter = itertools.count()

    def visit(payload: bytes) -> ET.Element:
        root = ET.fromstring(payload)  # noqa: S314 - controlled checked-in SVG
        _editor_metadata(root)
        _inline_images(root, visit)
        _rename_ids(root, prefix, counter)
        return root

    return visit(data)


def insert_svg(parent: ET.Element, graphic: ET.Element, bounds: list[float]) -> None:
    """Place a graphic without changing its coordinate system or aspect ratio."""
    item = copy.deepcopy(graphic)
    for name, value in zip(("x", "y", "width", "height"), bounds, strict=True):
        item.set(name, str(value))
    parent.append(item)


def card(graphic: ET.Element, width: float, height: float) -> bytes:
    """Add the reviewed rounded border and uniform padding."""
    output_width, padding = 1440, 32
    inner = output_width - 2 * padding
    output_height = height / width * inner + 2 * padding
    root = ET.Element(
        f"{{{SVG}}}svg",
        {
            "width": str(output_width),
            "height": str(output_height),
            "viewBox": f"0 0 {output_width} {output_height}",
        },
    )
    ET.SubElement(
        root,
        f"{{{SVG}}}rect",
        {
            "x": "1",
            "y": "1",
            "width": str(output_width - 2),
            "height": str(output_height - 2),
            "rx": "24",
            "fill": "white",
            "stroke": "#d8e0e7",
            "stroke-width": "2",
        },
    )
    insert_svg(root, graphic, [padding, padding, inner, output_height - 2 * padding])
    return ET.tostring(root, encoding="utf-8") + b"\n"


def poster_card(png: bytes) -> bytes:
    """Retain the complete graph/sequence raster at its measured 16 by 2.4 ratio."""
    root = ET.Element(f"{{{SVG}}}svg", {"viewBox": "0 0 2880 432"})
    ET.SubElement(
        root,
        f"{{{SVG}}}image",
        {
            "width": "2880",
            "height": "432",
            "href": "data:image/png;base64," + base64.b64encode(png).decode(),
        },
    )
    result = ET.fromstring(card(root, 2880, 432))  # noqa: S314 - generated SVG
    ET.SubElement(
        result,
        f"{{{SVG}}}rect",
        {
            "x": "1372",
            "y": "7",
            "width": "56",
            "height": "25",
            "rx": "6",
            "fill": "white",
            "stroke": "#b6c5d0",
        },
    )
    text = ET.SubElement(
        result,
        f"{{{SVG}}}text",
        {
            "x": "1400",
            "y": "24",
            "font-family": "Arial",
            "font-size": "13",
            "text-anchor": "middle",
            "fill": "#33475b",
        },
    )
    text.text = "Play"
    return ET.tostring(result, encoding="utf-8") + b"\n"


def workflow_card(context: bytes, native: bytes, bounds: list[float]) -> bytes:
    """Insert a newly rendered library sequence into the retained workflow."""
    root = closed_svg(context, "c")
    insert_svg(root, closed_svg(native, "s"), bounds)
    _, _, width, height = map(float, root.get("viewBox", "").split())
    return card(root, width, height)


def media_geometry(context: bytes, bounds: list[float]) -> dict[str, float]:
    """Use the same coordinate frame for the video and static vector fallback."""
    root = ET.fromstring(context)  # noqa: S314 - checked-in closed SVG
    x, y, width, height = map(float, root.get("viewBox", "").split())
    inner, padding = 1376, 32
    total_height = height / width * inner + 2 * padding
    left, top, span, rise = bounds
    return {
        "left": (padding + (left - x) / width * inner) / 1440 * 100,
        "top": (padding + (top - y) / width * inner) / total_height * 100,
        "width": span / width * inner / 1440 * 100,
        "height": rise / width * inner / total_height * 100,
    }


def player_css(context: bytes, config: dict[str, Any], clip: str) -> str:
    """Emit presentation geometry from one layout record rather than hand-fit CSS."""
    native = media_geometry(context, config["array_bounds_css_px"])
    packing = {
        "left": 32 / 1440 * 100,
        "top": 32 / (64 + 1376 * 0.15) * 100,
        "width": 1376 / 1440 * 100,
        "height": 1376 * 0.15 / (64 + 1376 * 0.15) * 100,
    }
    lines = [
        ".sequence-playback{position:relative;line-height:0;margin:1rem 0}",
        ".sequence-playback img{display:block;width:100%;height:auto}",
        ".sequence-playback video{position:absolute;object-fit:fill}",
        ".sequence-playback video[hidden]{display:none}",
        (
            ".sequence-playback button{position:absolute;right:6px;top:4px;z-index:2;"
            "font:11px Arial;padding:5px 10px;background:white;color:#33475b;"
            "border:1px solid #b6c5d0;border-radius:6px;cursor:pointer}"
        ),
    ]
    for kind, geometry in (("packing", packing), ("library", native)):
        body = "".join(f"{name}:{value:.10f}%;" for name, value in geometry.items())
        if kind == "library":
            body += f"clip-path:{clip};"
        lines.append(f".sequence-playback--{kind} video{{{body}}}")
    return "\n".join(lines) + "\n"
