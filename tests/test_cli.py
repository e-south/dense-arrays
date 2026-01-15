"""CLI smoke tests."""

from __future__ import annotations

from typing import TYPE_CHECKING

from typer.testing import CliRunner

from dense_arrays.cli import app

runner = CliRunner()

if TYPE_CHECKING:
    from pathlib import Path


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
    assert "either --motif or --motifs-file" in combined
