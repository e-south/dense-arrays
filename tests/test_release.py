"""A release cannot borrow successful checks from another commit or workflow."""

from __future__ import annotations

import re
import runpy
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _run(**overrides: object) -> dict:
    return {
        "id": 1,
        "head_sha": "abc123",
        "event": "push",
        "head_branch": "main",
        "path": ".github/workflows/ci.yml",
        "head_repository": {"full_name": "e-south/dense-arrays"},
        "status": "completed",
        "conclusion": "success",
        **overrides,
    }


def _check(runs: list[dict]) -> dict:
    module = runpy.run_path(str(ROOT / "maintenance/release.py"))
    return module["require_successful_ci"]({"workflow_runs": runs}, "abc123")


def test_exact_main_success_qualifies() -> None:
    assert _check([_run()])["id"] == 1


@pytest.mark.parametrize(
    "changes",
    [
        {"head_sha": "different"},
        {"event": "pull_request"},
        {"head_branch": "feature"},
        {"path": ".github/workflows/unrelated.yml"},
        {"head_repository": {"full_name": "another/dense-arrays"}},
    ],
)
def test_other_run_cannot_qualify(changes: dict) -> None:
    with pytest.raises(ValueError, match="No canonical"):
        _check([_run(**changes)])


@pytest.mark.parametrize(
    "changes", [{"conclusion": "failure"}, {"status": "in_progress"}]
)
def test_earlier_success_cannot_mask_latest_failure(changes: dict) -> None:
    with pytest.raises(ValueError, match="latest"):
        _check([_run(), _run(id=2, **changes)])


def test_package_readme_uses_durable_images_and_absolute_links() -> None:
    root = ROOT
    project = tomllib.loads((root / "pyproject.toml").read_text())["project"]
    assert project["readme"] == "README.md"
    text = (root / "README.md").read_text()
    targets = re.findall(r"\]\(([^)]+)\)", text)
    assert all(target.startswith(("https://", "#")) for target in targets)
    images = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
    banner = images[0]
    prefix = "https://raw.githubusercontent.com/e-south/dense-arrays/"
    assert banner.startswith(prefix)
    revision, relative = banner.removeprefix(prefix).split("/", 1)
    assert revision == "v" + project["version"] or re.fullmatch(
        r"[0-9a-f]{40}", revision
    )
    assert relative.endswith(".png")
    data = (root / relative).read_bytes()
    assert data.startswith(b"\x89PNG\r\n\x1a\n")
    assert len(data) < 10_000_000
