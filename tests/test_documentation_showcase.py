"""Verify the published media and their independent reproduction inputs.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "docs/assets/showcase"
RECIPE = ROOT / "maintenance/docs/showcase"
SVG = "http://www.w3.org/2000/svg"


def test_readme_uses_a_pinned_png_for_branch_and_package_rendering():
    readme = (ROOT / "README.md").read_text()
    assert re.search(
        r"https://raw\.githubusercontent\.com/e-south/dense-arrays/"
        r"[a-f0-9]{40}/docs/assets/showcase/packing-card\.png",
        readme,
    )
    with Image.open(ASSETS / "packing-card.png") as image:
        assert image.format == "PNG"
        assert image.width == 2880
        image.verify()


def test_published_media_match_reproduction_receipt():
    receipt = json.loads((ASSETS / "receipt.json").read_text())
    for base, files in ((RECIPE, receipt["inputs"]), (ASSETS, receipt["outputs"])):
        for name, expected in files.items():
            assert hashlib.sha256((base / name).read_bytes()).hexdigest() == expected

    records = json.loads((RECIPE / "inputs/records.json").read_text())
    arrays = {r["source_id"]: r for r in records["arrays"]}
    drawings = {r["id"]: r for r in records["drawing_records"]}
    for group, length in (("packing", 100), ("library", 60)):
        for key in receipt["record_ids"][group]:
            assert len(arrays[key]["sequence"]) == length
            assert arrays[key]["sequence"] == drawings[key]["sequence"]
    assert receipt["authority"] == "placement_reconstructed"


@pytest.mark.parametrize("name", ["packing.svg", "library.svg"])
def test_published_svg_is_closed_and_has_resolvable_unique_ids(name: str):
    data = (ASSETS / name).read_bytes()
    root = ET.fromstring(data)  # noqa: S314 - checked-in publication artifact
    ids = [n.attrib["id"] for n in root.iter() if "id" in n.attrib]
    assert len(ids) == len(set(ids))
    references = set()
    for node in root.iter():
        assert node.tag != f"{{{SVG}}}script"
        for attribute, value in node.attrib.items():
            assert "/Users/" not in value
            if attribute.endswith("href"):
                assert value.startswith(("#", "data:image/png;base64,"))
                if value.startswith("#"):
                    references.add(value[1:])
            references.update(re.findall(r"url\(#([\w.:-]+)\)", value))
    assert references <= set(ids)


@pytest.mark.parametrize(
    "content",
    [
        '<image href="file:///private/artwork.svg"/>',
        '<image href="https://example.org/artwork.png"/>',
        "<script>alert(1)</script>",
    ],
)
def test_showcase_composer_rejects_external_or_executable_artwork(content: str):
    spec = importlib.util.spec_from_file_location(
        "showcase_svg", RECIPE / "svg_assets.py"
    )
    assert spec
    assert spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with pytest.raises(ValueError, match=r"scripts|external resources"):
        module.closed_svg(f'<svg xmlns="{SVG}">{content}</svg>'.encode(), "test")
