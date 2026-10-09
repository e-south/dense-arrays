"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/cli.py

Task-oriented CLI translation; domain behavior lives in shared operations.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from pathlib import Path  # noqa: TC003 - Typer evaluates command annotations
from typing import Annotated

import typer

from dense_arrays._record_validation import canonical_json
from dense_arrays.arrays import CollectionView
from dense_arrays.artifacts import RunHandle
from dense_arrays.parts import Part, PreparationSet, PreparationSpec
from dense_arrays.planning import (
    BatchSampling,
    DesignSpec,
    GenerationPlan,
    Length,
    MatrixPlan,
    PlanEvidence,
    PreparationPlan,
    Target,
)
from dense_arrays.reporting import (
    BundleView,
    DiagnosticReport,
    LibrarySelection,
    LibraryView,
    PlanComparison,
    QualityComparison,
    QualityReport,
    QualitySnapshot,
    ReadLimits,
    RecordView,
    RequestReport,
    RunSummary,
)
from dense_arrays.reporting.pools import PoolQualityReport, PoolQualitySnapshot
from dense_arrays.reporting.selections.views import SelectionView
from dense_arrays.workflow import operations
from dense_arrays.workflow.cli_errors import diagnostics
from dense_arrays.workflow.inputs import read_source
from dense_arrays.workflow.presentation import (
    display_diagnostics,
    display_generation_plan,
    display_matrix_plan,
    display_plan,
    display_preparation_plan,
    display_quality,
    display_quality_comparison,
    display_records,
    display_summary,
)
from dense_arrays.workflow.queries import query_filter
from dense_arrays.workflow.selections import selection_cost

HELP_EPILOG = """Capabilities | Requirements

Supported in base | CSV/TSV parts, packing, inspection and text exports.

Optional dependency | Excel/Parquet: tables extra; PNG: playback extra;
FIMO scoring: MEME Suite executable.

Not implemented | Playback of an exact solver trace.

Documentation: https://dunloplab.gitlab.io/dense-arrays/library-workflow/

Versioned example files (v1 requests; save the three inputs together):
https://dunloplab.gitlab.io/dense-arrays/library-workflow/curated-example/
"""


def plan_command(
    source: Annotated[
        Path,
        typer.Argument(
            help=(
                "Preparation, set, design, matrix or extension request; "
                "or saved plan (YAML/JSON)."
            )
        ),
    ],
    *,
    out: Annotated[
        Path | None, typer.Option(help="Create a resolved plan file.")
    ] = None,
    json_output: Annotated[
        bool, typer.Option("--json", help="Versioned JSON on stdout.")
    ] = False,
    max_read_records: Annotated[
        int | None, typer.Option(help="Maximum parent evidence records examined.")
    ] = None,
    max_identity_entries: Annotated[
        int | None, typer.Option(help="Maximum parent identity entries retained.")
    ] = None,
) -> None:
    """Preview inputs, requirements and effort without running a solver.

    Example: dense-arrays plan design.yaml --out plan.json

    Inputs: a YAML/JSON request or saved plan and its declared source files.

    Outputs: a preview on stdout; --out creates a resolved plan file.

    On failure: fix the named input or requirement, then plan again.
    If the output exists, choose a new path.
    """
    with diagnostics(json_output=json_output):
        if out is not None and out.exists():
            msg = f"output destination already exists: {out}"
            raise FileExistsError(msg)
        request = read_source(source)
        caps = {
            key: value
            for key, value in (
                ("records", max_read_records),
                ("identities", max_identity_entries),
            )
            if value is not None
        }
        if caps and isinstance(request, (GenerationPlan, PreparationPlan, MatrixPlan)):
            msg = "read caps apply to resolving a parent, not an already resolved plan"
            raise ValueError(msg)
        result = (
            request
            if isinstance(request, (GenerationPlan, PreparationPlan, MatrixPlan))
            else operations.plan(
                request, read_limits=ReadLimits(**caps) if caps else None
            )
        )
        if out is not None:
            result.write(out)
        if json_output:
            typer.echo(canonical_json(result.to_dict()))
        elif isinstance(result, MatrixPlan):
            display_matrix_plan(result)
        elif isinstance(result, PreparationPlan):
            display_preparation_plan(result)
        else:
            display_generation_plan(result)
        if out is not None and not json_output:
            typer.echo(f"Plan saved: {out}", err=True)


def prepare_command(  # noqa: PLR0913 - shared preparation and batch options
    source: Annotated[
        Path,
        typer.Argument(
            help="Preparation request/plan or generation plan for batch sampling."
        ),
    ],
    *,
    out: Annotated[
        Path, typer.Option(help="New pool directory or prepared batch plan file.")
    ],
    batch_size: Annotated[
        int | None, typer.Option(help="Number of offered parts per active cell.")
    ] = None,
    batch_strategy: Annotated[
        str | None, typer.Option(help="uniform or group_balanced.")
    ] = None,
    batch_seed: Annotated[
        int | None, typer.Option(help="Explicit batch sampling seed (default 0).")
    ] = None,
    unique_sequences: Annotated[
        bool,
        typer.Option("--unique-sequences", help="Offer distinct supplied sequences."),
    ] = False,
    unique_cores: Annotated[
        bool,
        typer.Option(
            "--unique-cores", help="Offer distinct oriented cores within each group."
        ),
    ] = False,
    max_per_group: Annotated[
        int | None, typer.Option(help="Hard cap on offered parts per group.")
    ] = None,
    batch_count: Annotated[
        int | None, typer.Option(help="Prepared batches per active cell (default 1).")
    ] = None,
    attempts_per_batch: Annotated[
        int | None, typer.Option(help="Attempt cap for each scheduled batch.")
    ] = None,
    accepted_per_batch: Annotated[
        int | None, typer.Option(help="Accepted-design cap for each scheduled batch.")
    ] = None,
    json_output: Annotated[
        bool, typer.Option("--json", help="Versioned JSON receipt on stdout.")
    ] = False,
) -> None:
    """Prepare curated, PWM or background parts, or a candidate-batch plan.

    Example: dense-arrays prepare prepare.yaml --out pool

    Inputs: a preparation request/plan; batch sampling takes a generation plan.

    Outputs: a reusable pool directory, or a batch plan with --batch-size.

    On failure: inspect an incomplete pool with --view quality before retrying.
    Choose a new destination; install FIMO only for recipes that require it.
    """
    with diagnostics(json_output=json_output):
        if out.exists() or out.is_symlink():
            msg = f"output destination already exists: {out}"
            raise FileExistsError(msg)
        if batch_size is None and (
            batch_strategy is not None
            or batch_seed is not None
            or batch_count is not None
            or attempts_per_batch is not None
            or accepted_per_batch is not None
            or unique_sequences
            or unique_cores
            or max_per_group is not None
        ):
            msg = "batch options require --batch-size"
            raise ValueError(msg)
        policy = (
            None
            if batch_size is None
            else BatchSampling(
                batch_size,
                batch_strategy or "uniform",
                0 if batch_seed is None else batch_seed,
                unique_sequences=unique_sequences,
                unique_cores=unique_cores,
                max_per_group=max_per_group,
            )
        )
        result = operations.prepare(
            read_source(source),
            out=out,
            sampling=policy,
            batch_count=1 if batch_count is None else batch_count,
            attempts_per_batch=attempts_per_batch,
            accepted_per_batch=accepted_per_batch,
        )
        if isinstance(result, (GenerationPlan, MatrixPlan)):
            if json_output:
                typer.echo(canonical_json({**result.to_dict(), "path": str(out)}))
            else:
                typer.echo(
                    f"Prepared batch plan: {out}\n"
                    f"Run: dense-arrays run {out} --out <run-directory>",
                    err=True,
                )
            return
        summary = operations.inspect(result)
        if json_output:
            typer.echo(canonical_json({**summary.to_dict(), "path": str(result.path)}))
        elif summary.preparation is not None:
            typer.echo(
                f"{summary.state.capitalize()}: {summary.retained_parts} / "
                f"{summary.preparation.requested_retention} requested parts retained; "
                f"{summary.source_parts} candidates processed.\nPool: {result.path}\n"
                f"Details: dense-arrays inspect {result.path} --view quality",
                err=True,
            )
        else:
            typer.echo(
                f"Prepared: {summary.retained_parts} / {summary.source_parts} "
                f"supplied parts retained.\nPool: {result.path}\n"
                f"Details: dense-arrays inspect {result.path} --view parts",
                err=True,
            )
        if summary.state == "incomplete":
            raise typer.Exit(4 if summary.preparation.counts["execution_error"] else 3)


def run_command(  # noqa: PLR0913 - explicit CLI options share one request
    source: Annotated[
        Path | None,
        typer.Argument(
            help="Design, matrix or extension request; or saved generation plan."
        ),
    ] = None,
    *,
    out: Annotated[
        Path | None, typer.Option(help="New run directory; required unless resuming.")
    ] = None,
    resume: Annotated[
        Path | None,
        typer.Option(
            help="Resume an unchanged interrupted run, or verify a completed run."
        ),
    ] = None,
    motif: Annotated[
        list[str] | None,
        typer.Option(
            "--motif", help="Supplied motif; repeat for distinct occurrences."
        ),
    ] = None,
    length: Annotated[
        int | None, typer.Option(help="Maximum packed sequence length.")
    ] = None,
    count: Annotated[
        int | None, typer.Option(help="Accepted unique designs requested.")
    ] = None,
    seed: Annotated[int | None, typer.Option(help="Recorded root seed.")] = None,
    strands: Annotated[
        str | None, typer.Option(help="single or double; default double.")
    ] = None,
    json_output: Annotated[
        bool, typer.Option("--json", help="Versioned JSON receipt on stdout.")
    ] = False,
) -> None:
    """Generate a bounded collection and preserve results for inspection.

    Example: dense-arrays run design.yaml --out run

    Inputs: a design/matrix/extension request, saved plan, or inline motifs.

    Outputs: a native run directory with accepted designs and attempt evidence.

    On failure: inspect run --view diagnostics to explain a shortfall.
    Use run --resume run only for unchanged interrupted work; new requests
    need a new output directory.
    """
    with diagnostics(json_output=json_output):
        if resume is not None:
            if any(
                value is not None
                for value in (source, out, motif, length, count, seed, strands)
            ):
                msg = "resume is exclusive with source, output and design options"
                raise ValueError(msg)
            result = operations.run(resume=resume)
        elif source is not None:
            if any(
                value is not None for value in (motif, length, count, seed, strands)
            ):
                msg = "use a source request or inline options, not both"
                raise ValueError(msg)
            request = read_source(source)
        else:
            if not motif or length is None:
                msg = "provide SOURCE or --motif with --length"
                raise ValueError(msg)
            request = DesignSpec(
                parts=tuple(
                    Part(f"row:{i}", sequence) for i, sequence in enumerate(motif, 1)
                ),
                length=Length(maximum=length),
                target=Target(count=1 if count is None else count),
                seed=0 if seed is None else seed,
                strands="double" if strands is None else strands,
            )
        if resume is None:
            result = operations.run(request, out=out)
        summary = operations.inspect(result)
        if json_output:
            typer.echo(canonical_json({**summary.to_dict(), "path": str(result.path)}))
        else:
            typer.echo(
                f"{summary.state.capitalize()}: {summary.accepted} / {summary.target} "
                f"designs accepted; {summary.termination_reason}.\n"
                f"Results: {result.path}\n"
                f"Details: dense-arrays inspect {result.path} --verify",
                err=True,
            )
        if summary.state != "completed":
            raise typer.Exit(4 if summary.state == "failed" else 3)


def inspect_command(  # noqa: PLR0913 - explicit CLI options share one query
    artifact: Annotated[
        list[Path],
        typer.Argument(
            help="Run, bundle, array collection, pool, or saved plan paths."
        ),
    ],
    *,
    array_id: Annotated[
        list[str] | None,
        typer.Option("--array-id", help="Supplied array ID; repeatable."),
    ] = None,
    view: Annotated[
        str,
        typer.Option(
            help=(
                "summary; runs/bundles: designs, sequences, placements, batches, "
                "quality, selection; bundles: plans, plan; "
                "runs: attempts, diagnostics, plan, request; "
                "pools: parts, candidates, quality; "
                "array collections: arrays, parts, sequences, placements."
            )
        ),
    ] = "summary",
    verify: Annotated[
        bool, typer.Option(help="Independently verify committed evidence.")
    ] = False,
    compare: Annotated[
        Path | None,
        typer.Option(
            help="Compare --view plan rules or --view quality library metrics."
        ),
    ] = None,
    limit: Annotated[
        int | None, typer.Option(help="Display limit: 100 records or 20 diagnostics.")
    ] = None,
    all_rows: Annotated[
        bool, typer.Option("--all", help="Stream the complete declared record scope.")
    ] = False,
    group: Annotated[
        list[str] | None,
        typer.Option(
            "--group", help="Pool part group or placed group in designs; repeatable."
        ),
    ] = None,
    part_id: Annotated[
        list[str] | None,
        typer.Option(
            "--part-id", help="Pool part ID or placed part in designs; repeatable."
        ),
    ] = None,
    selection: Annotated[
        Path | None,
        typer.Option(help="Declared filter file; exclusive with filter flags."),
    ] = None,
    outcome: Annotated[
        list[str] | None,
        typer.Option(
            "--outcome", help="Attempt or preparation candidate outcome; repeatable."
        ),
    ] = None,
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
    attempt_id: Annotated[
        list[int] | None,
        typer.Option("--attempt-id", help="Attempt ordinal; repeatable."),
    ] = None,
    cell: Annotated[
        list[str] | None,
        typer.Option("--cell", help="Cell ID or run/cell reference; repeatable."),
    ] = None,
    design_id: Annotated[
        list[str] | None,
        typer.Option("--design-id", help="Design ID or full reference; repeatable."),
    ] = None,
    plan_id: Annotated[
        list[str] | None,
        typer.Option(
            "--plan-id", help="Full included plan ID; repeatable for --view plans."
        ),
    ] = None,
    max_read_records: Annotated[
        int | None, typer.Option(help="Maximum data records examined, not page size.")
    ] = None,
    max_pairs: Annotated[
        int | None,
        typer.Option(help="Maximum pair evaluations (record views use none)."),
    ] = None,
    max_identity_entries: Annotated[
        int | None,
        typer.Option(help="Maximum identity entries retained during reading."),
    ] = None,
    after: Annotated[
        str | None,
        typer.Option(help="Continue a record query at its original snapshot."),
    ] = None,
    json_output: Annotated[
        bool, typer.Option("--json", help="Versioned JSON report on stdout.")
    ] = False,
) -> None:
    """Read saved evidence without generating, exporting or repairing.

    Example: dense-arrays inspect run --verify

    Inputs: saved runs, bundles, pools, requests or plans; views vary by input.

    Outputs: a bounded report on stdout; --json selects versioned JSON.

    On failure: choose a matching view/filter; narrow the selection or raise
    an explicit read limit only after reviewing the reported cost.
    Use export to save data.
    """
    with diagnostics(json_output=json_output):
        source = artifact[0] if len(artifact) == 1 else artifact
        selected = query_filter(
            view,
            source=artifact[0] if len(artifact) == 1 else artifact,
            array_id=array_id,
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
        if verify:
            preview = operations.inspect(source, read_limits=ReadLimits(**caps))
            typer.echo(
                "Read cost: " + canonical_json(preview.verification_cost.to_dict()),
                err=True,
            )
            if isinstance(preview, RunSummary):
                source = RunHandle(source, preview.run_id, revision=preview.revision)
        if view == "selection" or isinstance(selected, LibrarySelection):
            cost = selection_cost(source, selected, ReadLimits(**caps))
            if cost is not None:
                typer.echo("Read cost: " + canonical_json(cost.to_dict()), err=True)
        result = operations.inspect(
            source,
            view=view,
            verify=verify,
            limit=limit,
            all=all_rows,
            select=selected,
            read_limits=ReadLimits(**caps) if caps else None,
            after=after,
            compare=compare,
        )
        if isinstance(
            result,
            (
                RequestReport,
                DesignSpec,
                PreparationSpec,
                PreparationSet,
                GenerationPlan,
                MatrixPlan,
                PlanEvidence,
                PreparationPlan,
                PlanComparison,
            ),
        ):
            display_plan(result, json_output=json_output)
        elif isinstance(
            result, (RecordView, LibraryView, BundleView, SelectionView, CollectionView)
        ):
            # A streamed page may already contain bytes; do not append an error object.
            with diagnostics():
                display_records(result, json_output=json_output)
        elif isinstance(result, DiagnosticReport):
            display_diagnostics(result, json_output=json_output)
        elif isinstance(result, QualityComparison):
            display_quality_comparison(result, json_output=json_output)
        elif isinstance(
            result,
            (QualityReport, QualitySnapshot, PoolQualityReport, PoolQualitySnapshot),
        ):
            display_quality(result, json_output=json_output)
        else:
            display_summary(result, json_output=json_output)


def register(app: typer.Typer) -> None:
    """Register the shared workflow operations on the command-line application."""
    from dense_arrays.workflow.export_cli import export_command  # noqa: PLC0415

    app.command("plan", epilog=HELP_EPILOG)(plan_command)
    app.command("prepare", epilog=HELP_EPILOG)(prepare_command)
    app.command("run", epilog=HELP_EPILOG)(run_command)
    app.command("inspect", epilog=HELP_EPILOG)(inspect_command)
    app.command("render", epilog=HELP_EPILOG)(render_command)
    app.command("export", epilog=HELP_EPILOG)(export_command)


def render_command(  # noqa: PLR0913 - shared rendering options
    artifact: Annotated[
        list[Path],
        typer.Argument(
            help="Run, bundle, array collection, pool or quality report paths."
        ),
    ],
    *,
    out: Annotated[Path, typer.Option(help="Create a PNG file.")],
    view: Annotated[
        str, typer.Option(help="design, array, library-quality or preparation-quality.")
    ] = "design",
    selection: Annotated[
        Path | None,
        typer.Option(
            help="Filter or saved selection; design requires exactly one match."
        ),
    ] = None,
    array_id: Annotated[
        list[str] | None,
        typer.Option("--array-id", help="Supplied array ID; repeatable."),
    ] = None,
    design_id: Annotated[list[str] | None, typer.Option("--design-id")] = None,
    cell: Annotated[list[str] | None, typer.Option("--cell")] = None,
    part_id: Annotated[list[str] | None, typer.Option("--part-id")] = None,
    group: Annotated[list[str] | None, typer.Option("--group")] = None,
    max_read_records: Annotated[
        int | None, typer.Option(help="Maximum report records examined.")
    ] = None,
    max_pairs: Annotated[
        int | None,
        typer.Option(help="Maximum pair evaluations during report verification."),
    ] = None,
    max_identity_entries: Annotated[
        int | None, typer.Option(help="Maximum retained report lookup entries.")
    ] = None,
    json_output: Annotated[
        bool, typer.Option("--json", help="Versioned receipt on stdout.")
    ] = False,
) -> None:
    """Render stored placements and requirement evidence without generating again.

    Example: dense-arrays render run --view library-quality --out quality.png

    Inputs: saved run, bundle, sampled pool or quality report evidence.

    Outputs: one PNG at a new path; no designs are generated.

    On failure: install the playback extra for missing rendering dependencies.
    For a design view, use inspect --view designs to select exactly one ID.
    """
    with diagnostics(json_output=json_output):
        caps = {
            key: value
            for key, value in (
                ("records", max_read_records),
                ("pairs", max_pairs),
                ("identities", max_identity_entries),
            )
            if value is not None
        }
        limits = ReadLimits(**caps) if caps else None
        source = artifact[0] if len(artifact) == 1 else artifact
        selected = query_filter(
            "quality" if view == "library-quality" else view,
            source=artifact[0] if len(artifact) == 1 else artifact,
            array_id=array_id,
            selection=selection,
            design_id=design_id,
            cell=cell,
            part_id=part_id,
            group=group,
        )
        if view in {"library-quality", "preparation-quality"}:
            if isinstance(selected, LibrarySelection):
                cost = selection_cost(source, selected, limits or ReadLimits())
                typer.echo("Read cost: " + canonical_json(cost.to_dict()), err=True)
            source = operations.inspect(
                source, view="quality", read_limits=limits, select=selected
            )
            typer.echo("Read cost: " + canonical_json(source.cost.to_dict()), err=True)
        receipt = operations.render(
            source,
            out=out,
            view=view,
            read_limits=None
            if isinstance(
                source,
                (
                    QualityReport,
                    QualitySnapshot,
                    PoolQualityReport,
                    PoolQualitySnapshot,
                ),
            )
            else limits,
            select=selected if view in {"design", "array"} else None,
        )
        if json_output:
            typer.echo(canonical_json(receipt.to_dict()))
        else:
            population = (
                "retained parts"
                if view == "preparation-quality"
                else "arrays"
                if view == "array"
                else "designs"
            )
            typer.echo(
                f"Rendered {receipt.view} for {receipt.records} {population}: "
                f"{receipt.destination}",
                err=True,
            )
