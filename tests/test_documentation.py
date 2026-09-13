"""Run reader examples and check links in the generated documentation.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest
from mkdocs.config import load_config

ROOT = Path(__file__).resolve().parents[1]


def _environment(tmp_path: Path) -> dict[str, str]:
    environment = dict(os.environ, TMPDIR=str(tmp_path))
    environment.pop("__PYVENV_LAUNCHER__", None)
    return environment


@pytest.mark.parametrize("page", ["quickstart.md", "constraints.md", "playback.md"])
def test_guide_python_examples(page: str, tmp_path: Path):
    source = (ROOT / "docs" / page).read_text(encoding="utf-8")
    examples = re.findall(r"^```python\n(.*?)^```", source, re.MULTILINE | re.DOTALL)
    assert examples, f"No maintained Python examples found in {page}"
    result = subprocess.run(  # noqa: S603
        [sys.executable, "-c", "\n\n".join(examples)],
        cwd=tmp_path,
        env=_environment(tmp_path),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, f"{page}:\n{result.stdout}\n{result.stderr}"


class _PageLinks(HTMLParser):
    """Collect browser navigation, media resources, and addressable anchors."""

    def __init__(self) -> None:
        super().__init__()
        self.anchors: set[str] = set()
        self.targets: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        """Read links, assets, and both modern and legacy named anchors."""
        attributes = dict(attrs)
        identifier = attributes.get("id")
        if identifier:
            self.anchors.add(identifier)
        if tag == "a" and attributes.get("name"):
            self.anchors.add(str(attributes["name"]))
        for name in ("href", "src", "poster"):
            target = attributes.get(name)
            if target:
                self.targets.append(target)


def _local_target(path: str, source: Path, site: Path, prefix: str) -> Path:
    if path.startswith("/"):
        target = (site / unquote(path.removeprefix(prefix)).lstrip("/")).resolve()
    elif path:
        target = (source.parent / unquote(path)).resolve()
    else:
        target = source
    return target / "index.html" if target.is_dir() else target


def test_built_documentation_links(tmp_path: Path):
    site = tmp_path / "site"
    prefix = urlsplit(load_config(str(ROOT / "mkdocs.yml")).site_url).path.rstrip("/")
    result = subprocess.run(  # noqa: S603
        [sys.executable, "-m", "mkdocs", "build", "--strict", "--site-dir", str(site)],
        cwd=ROOT,
        env=_environment(tmp_path),
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    pages: dict[Path, _PageLinks] = {}
    for path in site.rglob("*.html"):
        parser = _PageLinks()
        parser.feed(path.read_text(encoding="utf-8"))
        pages[path] = parser
    assert pages
    failures: list[str] = []
    for source, page in pages.items():
        for href in page.targets:
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            target = _local_target(url.path, source, site, prefix)
            description = f"{source.relative_to(site)}: {href}"
            if (
                not target.is_relative_to(site)
                or not target.is_file()
                or (
                    url.fragment
                    and target in pages
                    and unquote(url.fragment) not in pages[target].anchors
                )
            ):
                failures.append(description)
    assert not failures, "Broken built documentation links:\n" + "\n".join(failures)
