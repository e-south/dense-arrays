"""Public help teaches each operation without probing tools or running work.

Author: Eric J. South.
"""

import re
import shutil
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from dense_arrays.cli import app
from dense_arrays.workflow import operations


@pytest.mark.parametrize(
    "command,input_word,output_word,action",
    [
        ("plan", "request", "plan", "fix"),
        ("prepare", "request", "pool", "inspect"),
        ("run", "request", "run", "inspect"),
        ("inspect", "saved", "stdout", "limit"),
        ("export", "saved", "file", "--all"),
        ("render", "saved", "PNG", "playback"),
    ],
)
def test_operation_help_explains_first_use(
    command: str, input_word: str, output_word: str, action: str
):
    result = CliRunner().invoke(
        app, [command, "--help"], terminal_width=80, env={"NO_COLOR": "1"}
    )
    assert result.exit_code == 0, result.output
    text = " ".join(re.sub(r"[│╭╮╰╯─]", " ", result.stdout).split())
    assert f"Example: dense-arrays {command} " in text
    assert "Inputs:" in text
    assert input_word in text
    assert "Outputs:" in text
    assert output_word in text
    assert "On failure:" in text
    assert action in text
    assert text.index("Example:") < text.index("Inputs:") < text.index("Outputs:")
    assert "library-workflow/curated-example/" in text
    assert "Versioned example files" in text
    assert max(map(len, result.stdout.splitlines())) <= 80


def test_root_help_distinguishes_requirements_without_probing_tools(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    def forbidden(*_args: object, **_kwargs: object) -> None:
        pytest.fail("help must not execute operations or discover external tools")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(shutil, "which", forbidden)
    monkeypatch.setattr(subprocess, "run", forbidden)
    monkeypatch.setattr(subprocess, "Popen", forbidden)
    for operation in ("plan", "prepare", "run", "inspect", "export", "render"):
        monkeypatch.setattr(operations, operation, forbidden)
    for command in (
        [],
        ["plan"],
        ["prepare"],
        ["run"],
        ["inspect"],
        ["export"],
        ["render"],
    ):
        result = CliRunner().invoke(
            app, [*command, "--help"], terminal_width=80, env={"NO_COLOR": "1"}
        )
        assert result.exit_code == 0, result.output
        text = " ".join(re.sub(r"[│╭╮╰╯─]", " ", result.stdout).split())
        for phrase in (
            "Example:",
            "Inputs:",
            "Outputs:",
            "On failure:",
            "Supported in base",
            "Optional dependency",
            "Not implemented",
            "CSV/TSV",
            "tables",
            "playback",
            "FIMO",
            "solver trace",
            "https://dunloplab.gitlab.io/dense-arrays/",
        ):
            assert phrase in text
        assert max(map(len, result.stdout.splitlines())) <= 80
    assert not list(tmp_path.iterdir())
