"""CLI smoke tests."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from click import unstyle
from typer.testing import CliRunner

from dense_arrays import DenseArray, Optimizer, SolverBackendError
from dense_arrays.cli import app

runner = CliRunner()

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path


@pytest.mark.parametrize("command", ["optimize", "solutions"])
def test_backend_failure_is_a_concise_cli_error(command: str) -> None:
    result = runner.invoke(
        app,
        [command, "--motif", "ACGT", "--length", "8", "--solver", "INVALID"],
    )
    assert result.exit_code == 1
    assert "Error:" in result.stderr
    assert "INVALID" in result.stderr
    assert "Traceback" not in result.output


@pytest.mark.parametrize("command", ["optimize", "solutions"])
def test_invalid_motif_file_encoding_is_reported(command: str, tmp_path: Path) -> None:
    path = tmp_path / "motifs.txt"
    path.write_bytes(b"\xff\xfe")
    result = runner.invoke(app, [command, "--motifs-file", str(path), "--length", "8"])
    assert result.exit_code == 1
    assert "Error:" in result.stderr
    assert "Traceback" not in result.output


def test_backend_failure_after_a_result_is_not_success(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def interrupted_solutions(self: Optimizer, **_: object) -> Iterator[DenseArray]:
        yield DenseArray(self.library, self.sequence_length, [0], [None])
        message = "Backend stopped after the first result"
        raise SolverBackendError(message)

    monkeypatch.setattr(Optimizer, "solutions", interrupted_solutions)
    result = runner.invoke(app, ["solutions", "--motif", "ACGT", "--length", "8"])
    assert "Solution 1" in result.stdout
    assert result.exit_code == 1
    assert "Backend stopped" in result.stderr


def _write_motifs(tmp_path: Path) -> Path:
    motifs = ["ATGC", "CGT", "ATTA", "TTATTA"]
    path = tmp_path / "motifs.txt"
    path.write_text("\n".join(motifs) + "\n", encoding="utf-8")
    return path


def test_cli_optimize(tmp_path: Path) -> None:
    motifs_path = _write_motifs(tmp_path)
    result = runner.invoke(
        app,
        [
            "optimize",
            "--motifs-file",
            str(motifs_path),
            "--length",
            "8",
            "--strands",
            "single",
        ],
        env={"RICH_DISABLE": "1"},
    )
    assert result.exit_code == 0
    assert "Optimal solution" in result.stdout


def test_cli_rejects_both_inputs(tmp_path: Path) -> None:
    motifs_path = _write_motifs(tmp_path)
    result = runner.invoke(
        app,
        [
            "optimize",
            "--motifs-file",
            str(motifs_path),
            "--motif",
            "ATGC",
            "--length",
            "8",
        ],
        env={"RICH_DISABLE": "1"},
    )
    assert result.exit_code != 0
    combined = result.stdout + (result.stderr or "")
    plain_output = " ".join(unstyle(combined).split())
    assert "either --motif or --motifs-file" in plain_output
