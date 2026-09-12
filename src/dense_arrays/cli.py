"""Design motif arrays and report solver outcomes from the command line.

Module Author(s): Virgile Andreani, Eric J. South
Dunlop Lab
"""

from __future__ import annotations

import itertools as it
from enum import StrEnum
from pathlib import Path  # noqa: TC003
from typing import TYPE_CHECKING, Annotated

import typer
from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.text import Text

from .errors import OptimizationError
from .optimizer import Optimizer

if TYPE_CHECKING:
    from .solution import DenseArray

app = typer.Typer(
    add_completion=False,
    help="Design densely packed DNA arrays from motif libraries.",
)

console = Console()


class Strands(StrEnum):
    """Strand options for optimization."""

    single = "single"
    double = "double"


def _load_motifs(motif: list[str] | None, motifs_file: Path | None) -> list[str]:
    if motif and motifs_file is not None:
        msg = "Provide either --motif or --motifs-file, not both."
        raise typer.BadParameter(msg)
    motif = motif or []
    if motifs_file is None and not motif:
        msg = "Provide at least one motif via --motif or --motifs-file."
        raise typer.BadParameter(msg)

    if motifs_file is not None:
        lines = motifs_file.read_text(encoding="utf-8").splitlines()
        motifs = [
            line.strip()
            for line in lines
            if line.strip() and not line.lstrip().startswith("#")
        ]
    else:
        motifs = [m.strip() for m in motif if m.strip()]

    if not motifs:
        msg = "Motif list is empty after parsing."
        raise typer.BadParameter(msg)
    return motifs


def _print_solution(title: str, solution: DenseArray) -> None:
    header = Text(
        f"{title} | score {solution.nb_motifs} | length {len(solution.sequence)} / "
        f"{solution.sequence_length} limit | "
        f"compression {solution.compression_ratio:.3f}",
        style="bold green",
    )
    console.print(header)
    console.print(Panel.fit(str(solution), border_style="cyan"))


@app.command()
def optimize(
    length: Annotated[
        int,
        typer.Option("--length", min=1, help="Target sequence length."),
    ],
    strands: Annotated[
        Strands,
        typer.Option(
            "--strands",
            case_sensitive=False,
            help="single or double stranded optimization.",
        ),
    ] = Strands.double,
    solver: Annotated[
        str,
        typer.Option("--solver", help="OR-Tools solver backend."),
    ] = "CBC",
    motif: Annotated[
        list[str] | None,
        typer.Option("--motif", help="Motif sequence (repeatable)."),
    ] = None,
    motifs_file: Annotated[
        Path | None,
        typer.Option(
            "--motifs-file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            help="File with one motif per line (blank lines and # comments ignored).",
        ),
    ] = None,
) -> None:
    """Compute a single optimal dense array.

    Raises
    ------
    typer.Exit
        If input cannot be read or validated, or optimization fails.
    """
    try:
        motifs = _load_motifs(motif, motifs_file)
        opt = Optimizer(motifs, sequence_length=length, strands=strands.value)
        best = opt.optimal(solver=solver)
    except (OSError, ValueError, OptimizationError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    _print_solution("Optimal solution", best)


@app.command()
def solutions(  # noqa: PLR0913, PLR0917
    length: Annotated[
        int,
        typer.Option("--length", min=1, help="Target sequence length."),
    ],
    strands: Annotated[
        Strands,
        typer.Option(
            "--strands",
            case_sensitive=False,
            help="single or double stranded optimization.",
        ),
    ] = Strands.double,
    solver: Annotated[
        str,
        typer.Option("--solver", help="OR-Tools solver backend."),
    ] = "CBC",
    max_solutions: Annotated[
        int,
        typer.Option("--max-solutions", min=1, help="Number of solutions to display."),
    ] = 10,
    diverse: Annotated[
        bool,
        typer.Option("--diverse", help="Return solutions in diversity-biased order."),
    ] = False,
    motif: Annotated[
        list[str] | None,
        typer.Option("--motif", help="Motif sequence (repeatable)."),
    ] = None,
    motifs_file: Annotated[
        Path | None,
        typer.Option(
            "--motifs-file",
            exists=True,
            file_okay=True,
            dir_okay=False,
            readable=True,
            help="File with one motif per line (blank lines and # comments ignored).",
        ),
    ] = None,
) -> None:
    """List multiple solutions in decreasing score order.

    Raises
    ------
    typer.Exit
        If input cannot be read or validated, or optimization fails.
    """
    try:
        motifs = _load_motifs(motif, motifs_file)
        opt = Optimizer(motifs, sequence_length=length, strands=strands.value)
        iterator = (
            opt.solutions_diverse(solver=solver)
            if diverse
            else opt.solutions(solver=solver)
        )
        solutions_iter = it.islice(iterator, max_solutions)
        any_solution = False
        for idx, sol in enumerate(solutions_iter, start=1):
            if idx > 1:
                console.print(Rule(style="dim"))
            _print_solution(f"Solution {idx}", sol)
            any_solution = True
        if not any_solution:
            typer.echo("No feasible solution was found.", err=True)
            raise typer.Exit(code=1)
    except (OSError, ValueError, OptimizationError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc


if __name__ == "__main__":
    app()
