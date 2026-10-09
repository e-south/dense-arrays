"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/cli_errors.py

Consistent CLI diagnostics and exit codes around shared operations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from typing import TYPE_CHECKING

import typer

from dense_arrays._record_validation import canonical_json, mutable_json
from dense_arrays.artifacts.errors import ArtifactIntegrityError
from dense_arrays.artifacts.recovery import RecoveryError
from dense_arrays.errors import OptimizationError
from dense_arrays.parts.scoring import ScoringError
from dense_arrays.reporting import ReadLimitError
from dense_arrays.reporting.selections import SelectionShortfall
from dense_arrays.workflow.execution import RunExecutionError

if TYPE_CHECKING:
    from collections.abc import Iterator


@contextmanager
def diagnostics(*, json_output: bool = False) -> Iterator[None]:
    """Translate failures without contaminating an already streaming data output."""
    try:
        yield
    except KeyboardInterrupt as err:
        typer.echo("Interrupted.", err=True)
        for note in getattr(err, "__notes__", ()):
            typer.echo(note, err=True)
        raise typer.Exit(130) from err
    except (
        ValueError,
        TypeError,
        OSError,
        sqlite3.Error,
        OptimizationError,
        ReadLimitError,
        ScoringError,
        RunExecutionError,
    ) as err:
        execution_error = isinstance(
            err,
            (
                sqlite3.Error,
                OptimizationError,
                ArtifactIntegrityError,
                ReadLimitError,
                RecoveryError,
                ScoringError,
                RunExecutionError,
            ),
        )
        code = (
            err.code
            if isinstance(err, RecoveryError)
            else "selection_shortfall"
            if isinstance(err, SelectionShortfall)
            else "read_limit"
            if isinstance(err, ReadLimitError)
            else "artifact_integrity"
            if isinstance(err, ArtifactIntegrityError)
            else "scoring_error"
            if isinstance(err, ScoringError)
            else "execution_error"
            if execution_error
            else "invalid_input"
        )
        exit_code = 4 if execution_error else 2
        if json_output:
            typer.echo(
                canonical_json(
                    {
                        "schema": "dense_arrays.error.v1",
                        "code": code,
                        "message": str(err),
                        "exit_code": exit_code,
                        **(
                            {"reason": err.reason}
                            if isinstance(err, ScoringError)
                            else {}
                        ),
                        **(
                            {"counts": mutable_json(err.counts)}
                            if isinstance(err, SelectionShortfall)
                            else {}
                        ),
                        "artifact": (
                            str(err.artifact)
                            if isinstance(
                                err,
                                (
                                    ArtifactIntegrityError,
                                    RecoveryError,
                                    RunExecutionError,
                                ),
                            )
                            and err.artifact is not None
                            else None
                        ),
                    }
                )
            )
        typer.echo(f"{code}: {err}", err=True)
        if isinstance(err, RunExecutionError):
            typer.echo(
                f"Committed run: {err.artifact}. Use inspect to review its evidence.",
                err=True,
            )
        raise typer.Exit(exit_code) from err
