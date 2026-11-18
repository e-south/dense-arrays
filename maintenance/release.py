"""Require successful CI for the exact main commit before publishing a release.

Module Author(s): Eric J. South
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def package_banner_url(readme: str, version: str) -> str:
    """Admit the current release tag or an immutable source image.

    Returns
    -------
    str
        The absolute PNG URL selected by the package description.

    Raises
    ------
    ValueError
        If the banner uses a mutable reference or a different release version.
    """
    images = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", readme)
    match = re.fullmatch(
        r"https://raw\.githubusercontent\.com/e-south/dense-arrays/([^/]+)/"
        r"docs/assets/dense-arrays-banner\.png",
        images[0] if images else "",
    )
    if not match or not (
        match[1] == "v" + version or re.fullmatch(r"[a-f0-9]{40}", match[1])
    ):
        message = "Banner must select the release version or an immutable source commit"
        raise ValueError(message)
    return match[0]


def require_successful_ci(payload: dict, revision: str) -> dict:
    """Return the latest main-push CI run, rejecting missing or failed evidence.

    Raises
    ------
    ValueError
        If the exact revision lacks a successful canonical CI run.

    Returns
    -------
    dict
        The matching GitHub Actions run retained in release provenance.
    """
    runs = [
        run
        for run in payload["workflow_runs"]
        if run.get("head_sha") == revision
        and run.get("event") == "push"
        and run.get("head_branch") == "main"
        and run.get("path") == ".github/workflows/ci.yml"
        and (run.get("head_repository") or {}).get("full_name")
        == "e-south/dense-arrays"
    ]
    if not runs:
        message = "No canonical main-push CI run for the release commit"
        raise ValueError(message)
    latest = max(runs, key=lambda run: int(run["id"]))
    if latest.get("status") != "completed" or latest.get("conclusion") != "success":
        message = "The latest main-push CI run has not completed successfully"
        raise ValueError(message)
    return latest


def main() -> None:
    """Validate downloaded CI evidence and print its qualifying run."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checks", type=Path, required=True)
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            require_successful_ci(json.loads(args.checks.read_text()), args.revision),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
