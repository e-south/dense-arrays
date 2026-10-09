"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/export_cli.py

Portable file and stdout handoffs for the shared export operation.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

import sys
from pathlib import Path
from typing import Annotated

import typer

from dense_arrays._record_validation import canonical_json
from dense_arrays.reporting import ReadLimits
from dense_arrays.reporting.exporting import validate_format
from dense_arrays.reporting.selections import LibrarySelection
from dense_arrays.workflow.cli_errors import diagnostics
from dense_arrays.workflow.exporting import (
    export_cost,
    publish_export,
    resolve_export,
)
from dense_arrays.workflow.queries import query_filter


def export_command(  # noqa: PLR0913 - explicit paired query/output options
    artifact: Annotated[
        list[Path],
        typer.Argument(help="Run/bundle paths, one pool, or a request/plan file."),
    ],
    *,
    out: Annotated[
        str, typer.Option(help="New file or bundle directory; '-' streams file data.")
    ],
    view: Annotated[
        str,
        typer.Option(
            help=(
                "designs, sequences, placements, attempts, batches, parts, "
                "candidates; summary, request, plan, quality, diagnostics; "
                "bundled plans."
            )
        ),
    ] = "designs",
    format: Annotated[  # noqa: A002 - paired public output option
        str,
        typer.Option(help="json; scalar csv/tsv; sequences fasta; bundle; selection."),
    ] = "json",
    all_rows: Annotated[
        bool,
        typer.Option(
            "--all",
            help="Complete filtered records; omit for bounded selections or documents.",
        ),
    ] = False,
    limit: Annotated[
        int | None,
        typer.Option(help="Usage/diagnostic display rows in a report export."),
    ] = None,
    compare: Annotated[
        Path | None,
        typer.Option(help="Compare --view plan rules or --view quality metrics."),
    ] = None,
    selection: Annotated[
        Path | None,
        typer.Option(help="Declared filter, allocation policy or saved snapshot."),
    ] = None,
    design_id: Annotated[list[str] | None, typer.Option("--design-id")] = None,
    plan_id: Annotated[list[str] | None, typer.Option("--plan-id")] = None,
    cell: Annotated[list[str] | None, typer.Option("--cell")] = None,
    part_id: Annotated[list[str] | None, typer.Option("--part-id")] = None,
    group: Annotated[list[str] | None, typer.Option("--group")] = None,
    recipe_id: Annotated[
        list[str] | None,
        typer.Option(
            "--recipe-id",
            help="Preparation recipe ID; repeatable, candidates view only.",
        ),
    ] = None,
    candidate_index: Annotated[
        list[int] | None,
        typer.Option(
            "--candidate-index", help="Preparation candidate index; repeatable."
        ),
    ] = None,
    reason: Annotated[
        list[str] | None,
        typer.Option("--reason", help="Preparation rejection reason; repeatable."),
    ] = None,
    attempt_id: Annotated[list[int] | None, typer.Option("--attempt-id")] = None,
    outcome: Annotated[list[str] | None, typer.Option("--outcome")] = None,
    max_read_records: Annotated[int | None, typer.Option()] = None,
    max_pairs: Annotated[int | None, typer.Option()] = None,
    max_identity_entries: Annotated[int | None, typer.Option()] = None,
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Versioned receipt; stderr when data uses stdout."),
    ] = False,
) -> None:
    """Export records, editable requests, resolved plans or scoped JSON reports.

    Example: dense-arrays export run --all --out designs.json

    Inputs: saved runs/bundles, one pool, or a request/plan file.

    Outputs: a new file or bundle; --out - streams supported text data.

    On failure: use --all for complete record exports, or supply a selection.
    Match --view to --format and choose a new destination if it exists.
    """
    with diagnostics(json_output=json_output and out != "-"):
        validate_format(view, format)
        selected = query_filter(
            view,
            selection=selection,
            design_id=design_id,
            cell=cell,
            part_id=part_id,
            group=group,
            attempt_id=attempt_id,
            candidate_index=candidate_index,
            recipe_id=recipe_id,
            reason=reason,
            outcome=outcome,
            plan_id=plan_id,
        )
        caps = {
            key: value
            for key, value in (
                ("records", max_read_records),
                ("pairs", max_pairs),
                ("identities", max_identity_entries),
            )
            if value is not None
        }
        if (
            isinstance(selected, LibrarySelection)
            or format == "selection"
            or view == "selection"
        ):
            from dense_arrays.workflow.selections import selection_cost  # noqa: PLC0415

            cost = selection_cost(
                artifact[0] if len(artifact) == 1 else artifact,
                selected,
                ReadLimits(**caps),
            )
            typer.echo("Read cost: " + canonical_json(cost.to_dict()), err=True)
        query = resolve_export(
            artifact[0] if len(artifact) == 1 else artifact,
            view=view,
            selected=selected,
            all_rows=all_rows,
            format_name=format,
            limit=limit,
            read_limits=ReadLimits(**caps),
            compare=compare,
        )
        typer.echo(
            "Read cost: " + canonical_json(export_cost(query, format).to_dict()),
            err=True,
        )
        receipt = publish_export(
            query,
            format_name=format,
            out=sys.stdout if out == "-" else out,
        )
        if json_output:
            typer.echo(canonical_json(receipt.to_dict()), err=out == "-")
        else:
            typer.echo(
                f"Exported {receipt.records} {receipt.view}: {receipt.destination}",
                err=True,
            )
        if receipt.selection is not None and receipt.selection["status"] == "partial":
            typer.echo(
                "Selection shortfall: "
                f"{receipt.selection['shortfall']} unfilled designs.",
                err=True,
            )
            raise typer.Exit(3)
