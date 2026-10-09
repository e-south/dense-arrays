"""
--------------------------------------------------------------------------------
Dense Arrays
dense-arrays/src/dense_arrays/workflow/presentation.py

Human and JSON projections of shared workflow reports and record pages.

Module Author(s): Eric J. South
Maintainer(s): Eric J. South
--------------------------------------------------------------------------------
"""  # noqa: D205, D400 - structured module header

from __future__ import annotations

from typing import TYPE_CHECKING

import typer

from dense_arrays._record_validation import canonical_json, mutable_json
from dense_arrays.artifacts.preparation.candidates import PoolCandidate
from dense_arrays.parts import BoundParts, PreparationSet, PreparationSpec
from dense_arrays.planning import (
    DesignSpec,
    GenerationPlan,
    MatrixPlan,
    MatrixSpec,
    PlanEvidence,
    PreparationPlan,
)
from dense_arrays.planning.preparation import preparation_to_dict
from dense_arrays.planning.serialization import request_to_dict
from dense_arrays.reporting import (
    BundleView,
    DiagnosticReport,
    LibraryView,
    PlanComparison,
    QualityComparison,
    QualityReport,
    RecordView,
    RequestReport,
    RunSummary,
)
from dense_arrays.reporting.pools import PoolQualityReport, PoolQualitySnapshot

if TYPE_CHECKING:
    from dense_arrays.artifacts.provenance import Producer


def batch_summary(request: DesignSpec) -> str | None:
    """Describe the offered search and its local allowance in plain output."""
    policy = request.schedule or request.resampling
    continuation = (
        " Unproven search: next batch."
        if policy is not None and policy.on_unproven == "next_batch"
        else ""
    )
    if request.schedule is not None:
        schedule = request.schedule
        maximum = max(len(batch.part_ids) for batch in schedule.batches)
        accepted = (
            f"; accepted per batch: {schedule.accepted_per_batch}"
            if schedule.accepted_per_batch is not None
            else ""
        )
        return (
            f"Batches: {len(schedule.batches)}; maximum offered parts: {maximum}; "
            f"attempts per batch: {schedule.attempts_per_batch}{accepted}."
            f"{continuation}"
        )
    if request.resampling is not None:
        policy = request.resampling
        return (
            f"Runtime batches: at most {policy.max_batches}; "
            f"offered parts: {policy.sampling.size}; "
            f"attempts per batch: {policy.attempts_per_batch}; "
            f"feedback: {'enabled' if policy.feedback else 'disabled'}. {continuation}"
        )
    if request.batch is not None:
        return (
            f"Batch: {len(request.batch.part_ids)} offered / "
            f"{len(request.parts)} eligible parts."
        )
    return None


def display_generation_plan(result: GenerationPlan) -> None:
    """Show resolved search scope and bounded effort without running a solver."""
    preview = result.preview
    length = result.request.length
    length_text = (
        f"maximum {length.maximum}"
        if length.maximum is not None
        else f"exact {length.exact}"
    )
    variable_bound = "at most " if result.request.schedule is not None else ""
    search = (
        "Search: greedy; optimality unproven; one proposal per offered batch."
        if result.request.search == "greedy"
        else (
            "Search: exact CBC; feasibility unknown; "
            f"{variable_bound}{preview['path_variables']} path variables."
        )
    )
    typer.echo(
        f"Plan: {preview['parts']} parts; one design combination; "
        f"target {preview['target']}; "
        f"length {length_text}.\n"
        f"{search}\n"
        f"Limits: {result.request.limits.attempts} attempts, "
        f"{result.request.limits.active_seconds:g} active seconds."
    )
    if batch := batch_summary(result.request):
        typer.echo(batch)
    _display_search_policy(result.request, include_method=False)
    if result.parent is not None:
        typer.echo(
            f"Additional target: {result.request.target.count}; "
            f"{len(result.parent.exclusions)} parent/ancestor "
            "sequences excluded.\n"
            f"Parent stays {result.parent.state}: "
            f"{result.parent.accepted} / {result.parent.target}."
        )
    if result.request.exclude is not None:
        scope = result.request.exclude
        mapping = ", ".join(f"{a} -> {b}" for a, b in scope.cell_mapping.items())
        typer.echo(
            f"{len(result.exclusions)} excluded sequences; "
            f"cells {mapping}; uniqueness {scope.uniqueness}."
        )


def _display_search_policy(
    request: DesignSpec, *, indent: str = "", include_method: bool = True
) -> None:
    """Name the search method and optional preference beside the recipe."""
    if include_method and request.search == "greedy":
        typer.echo(
            f"{indent}Search: greedy; optimality unproven; "
            "one proposal per offered batch."
        )
    if request.packing_preference:
        typer.echo(
            f"{indent}Packing preference: underused parts; "
            "proposed usage within each offered batch."
        )


def display_matrix_plan(result: MatrixPlan, *, limit: int = 10) -> None:
    """Show allocation with a bounded cell preview and execution availability."""
    typer.echo(
        f"Matrix plan {result.plan_id}: {len(result.cells)} cells "
        "(design combinations); "
        f"{result.preview['active_cells']} active; target {result.total}."
    )
    for cell in result.cells[:limit]:
        status = "active" if cell.active else "inactive"
        excluded = (
            f"; {len(cell.plan.exclusions)} excluded" if cell.plan.exclusions else ""
        )
        typer.echo(
            f"{cell.cell_id}: target {cell.target} ({status}); "
            f"{len(cell.plan.request.parts)} eligible parts{excluded}"
        )
        if result.request.sources:
            typer.echo(f"  Part collection {cell.plan.collection_id[:12]}")
        if batch := batch_summary(cell.plan.request):
            typer.echo(f"  {batch}")
        _display_search_policy(cell.plan.request, indent="  ")
    if len(result.cells) > limit:
        typer.echo(f"{len(result.cells) - limit} more cells; use --json for all cells.")
    typer.echo("Effort is shared; exact sequence uniqueness applies within each cell.")


def display_plan(
    result: RequestReport
    | DesignSpec
    | MatrixSpec
    | PreparationSpec
    | PreparationSet
    | GenerationPlan
    | MatrixPlan
    | PlanEvidence
    | PreparationPlan
    | PlanComparison,
    *,
    json_output: bool,
) -> None:
    """Display resolved choices or their semantic differences without executing."""
    if isinstance(result, RequestReport) and not json_output:
        _display_request(result)
        return
    if json_output:
        typer.echo(
            canonical_json(
                request_to_dict(result)
                if isinstance(result, DesignSpec)
                else preparation_to_dict(result)
                if isinstance(result, (PreparationSpec, PreparationSet))
                else result.to_dict()
            )
        )
    elif isinstance(result, MatrixPlan):
        display_matrix_plan(result)
    elif isinstance(result, MatrixSpec):
        typer.echo(
            f"Matrix request: {len(result.axes)} axes; up to {result.max_cells} cells."
        )
    elif isinstance(result, DesignSpec):
        _display_design_request(result)
    elif isinstance(result, (PreparationSpec, PreparationSet)):
        _display_preparation_request(result)
    elif isinstance(result, PreparationPlan):
        display_preparation_plan(result)
    elif isinstance(result, PlanComparison):
        _display_comparison(result)
    else:
        label = "Plan evidence" if isinstance(result, PlanEvidence) else "Plan"
        typer.echo(
            f"{label} {result.plan_id}; target {result.request.target.count}; "
            f"{len(result.request.parts)} parts."
        )
        if batch := batch_summary(result.request):
            typer.echo(batch)
        _display_search_policy(result.request)


def _display_preparation_request(result: PreparationSpec | PreparationSet) -> None:
    """Summarize preparation input families without opening source files."""
    if isinstance(result, PreparationSet):
        typer.echo(f"Preparation set: {len(result.recipes)} independent recipes.")
    else:
        typer.echo(f"Preparation request: {type(result.source).__name__}.")


def _display_comparison(result: PlanComparison) -> None:
    """Bound human comparison output without omitting its continuation route."""
    typer.echo("Changed: " + ", ".join(result.changed_fields))
    typer.echo("Unchanged: " + ", ".join(result.unchanged_fields))
    shown = result.changes[:20]
    for change in shown:
        before = _change_value(change.before)
        after = _change_value(change.after)
        typer.echo(f"{change.kind} {change.pointer}: {before} -> {after}")
    if len(shown) < len(result.changes):
        typer.echo(
            f"{len(result.changes) - len(shown)} more changes; "
            "use --json for all values."
        )


def _display_design_request(request: DesignSpec) -> None:
    """Describe the declared source without reading it or claiming execution."""
    source = request.parts
    if isinstance(source, BoundParts):
        mode = (
            "embedded parts; original fingerprints retained"
            if source.locations is None
            else "bound parts; source files checked during planning and execution"
        )
        description = f"{len(source.parts)} {mode}"
    elif isinstance(source, tuple):
        description = f"{len(source)} resolved parts"
    else:
        description = f"unresolved {type(source).__name__} source"
    typer.echo(f"Design request: {description}; target {request.target.count}.")
    _display_search_policy(request)


def _display_request(result: RequestReport) -> None:
    """Show editable rules with origin and exclusion scope beside the request."""
    display_plan(result.request, json_output=False)
    if isinstance(result.request, DesignSpec) and result.request.lineage is not None:
        parent = result.request.lineage.parent
        typer.echo(f"Parent: {parent.run_id} at revision {parent.revision}.")
    if isinstance(result.request, DesignSpec) and result.request.exclude is not None:
        scope = result.request.exclude
        mapping = ", ".join(f"{a} -> {b}" for a, b in scope.cell_mapping.items())
        typer.echo(
            f"Library exclusion: cells {mapping}; uniqueness {scope.uniqueness}."
        )


def _change_value(value: object) -> str:
    """Bound human examples while marking omitted content explicitly."""
    text = canonical_json(mutable_json(value))
    width = 160
    return text if len(text) <= width else text[:width] + "... [full value: --json]"


def display_quality(result: QualityReport, *, json_output: bool) -> None:
    """Present exact scope and denominators before the paginated usage tables."""
    if isinstance(result, (PoolQualityReport, PoolQualitySnapshot)):
        display_pool_quality(result, json_output=json_output)
        return
    typer.echo("Read cost: " + canonical_json(result.cost.to_dict()), err=True)
    value = result.to_dict()
    if json_output:
        typer.echo(canonical_json(value))
        return
    typer.echo(
        f"{value['selection']['designs']} selected designs; "
        f"{value['selection']['distinct_sequences']} distinct sequences.\n"
        f"Eligible supply: {canonical_json(value['supply'])}"
    )
    search = value["search"]
    if search["availability"] == "not_included":
        typer.echo("Search history: not included in the supplied artifacts.")
    else:
        scope = (
            "all sources"
            if search["availability"] == "complete"
            else "available native sources; partial history"
        )
        typer.echo(
            f"Source search outcomes ({scope}): "
            f"{canonical_json(search['attempt_counts'])}"
        )
    for source in value["source_runs"]:
        attainment = source["attainment"]
        typer.echo(
            f"Source {source['run_id']} ({source['state']}): "
            f"{attainment['accepted']} / {attainment['target']} accepted; "
            f"{attainment['shortfall']} shortfall; "
            f"{source['included_designs']} included; "
            f"{source['selected_designs']} selected."
        )
    for name, metric in value["composition"].items():
        typer.echo(
            f"{name}: mean={metric['mean']}; n={metric['count']} selected designs"
        )
    for part in value["part_usage"]:
        typer.echo(
            f"{part.get('part_ref', part['part_id'])}: "
            f"{part['occurrences']} occurrences; "
            f"{part['designs']} / {part['design_denominator']} designs"
        )
    if value["next_cursor"]:
        typer.echo(
            f"Continue usage tables with --after {value['next_cursor']}", err=True
        )


def display_quality_comparison(result: QualityComparison, *, json_output: bool) -> None:
    """Show population sizes and only changed or unavailable aggregate metrics."""
    typer.echo("Read cost: " + canonical_json(result.cost.to_dict()), err=True)
    value = result.to_dict()
    if json_output:
        typer.echo(canonical_json(value))
        return
    typer.echo(
        f"Selected designs: {value['before']['selected_designs']} -> "
        f"{value['after']['selected_designs']}. Differences are after minus before."
    )
    for metric in result.metrics:
        if metric.status != "comparable" or metric.delta:
            typer.echo(
                f"{'.'.join(metric.path)}: {metric.before} -> {metric.after}; "
                f"delta={metric.delta}; n={metric.before_denominator} -> "
                f"{metric.after_denominator}; {metric.reason or metric.status}"
            )


def display_diagnostics(result: DiagnosticReport, *, json_output: bool) -> None:
    """Show a shared evidence report after disclosing the scan cost."""
    typer.echo("Read cost: " + canonical_json(result.cost.to_dict()), err=True)
    value = result.to_dict()
    if json_output:
        typer.echo(canonical_json(value))
        return
    typer.echo(
        f"{value['accepted']} / {value['target']} accepted; "
        f"shortfall: {value['shortfall']}.\n"
        f"Attempt outcomes: {canonical_json(value['attempt_counts'])}\n"
        f"Reason totals (may overlap): {canonical_json(value['reason_counts'])}"
    )
    for diagnostic in result.diagnostics:
        rule = (
            ""
            if diagnostic.requirement_id is None
            else f" [{diagnostic.requirement_id}]"
        )
        typer.echo(
            f"{diagnostic.severity}: {diagnostic.code}{rule}\n"
            f"  Observed: {canonical_json(mutable_json(diagnostic.observed))}\n"
            f"  Expected: {canonical_json(mutable_json(diagnostic.expected))}\n"
            f"  Evidence: {', '.join(diagnostic.evidence_refs)}\n"
            f"  Proof scope: {diagnostic.proof_scope or 'unknown'}\n"
            f"  Next: {diagnostic.next_action}"
        )
    if result.omitted:
        typer.echo(
            f"{result.omitted} additional diagnostics omitted; increase --limit."
        )


def display_records(result: RecordView | LibraryView, *, json_output: bool) -> None:
    """Stream rows without materializing an unbounded CLI result list."""
    typer.echo("Read cost: " + canonical_json(result.cost.to_dict()), err=True)
    if json_output:
        header = canonical_json(
            {
                "schema": "dense_arrays.record_page.v1",
                "revision": result.revision,
                "view": result.view,
                "cost": result.cost.to_dict(),
                **(
                    {"sources": list(result.sources)}
                    if isinstance(result, (LibraryView, BundleView))
                    else {}
                ),
            }
        )
        typer.echo(header[:-1] + ',"records":[', nl=False)
    with result.records() as records:
        for index, record in enumerate(records):
            if isinstance(record, PoolCandidate) and not json_output:
                display_candidate(record)
                continue
            if isinstance(record, PlanEvidence) and not json_output:
                display_plan(record, json_output=False)
                continue
            if json_output and index:
                typer.echo(",", nl=False)
            typer.echo(canonical_json(record.to_dict()), nl=not json_output)
    if json_output:
        tail = canonical_json(
            {
                "examined": records.examined,
                "returned": records.returned,
                "next_cursor": records.next_cursor,
            }
        )
        typer.echo("]," + tail[1:])
    elif records.next_cursor is not None:
        typer.echo(f"Continue with --after {records.next_cursor}", err=True)


def display_summary(result: object, *, json_output: bool) -> None:
    """Present allocation and source status without exposing record membership."""
    from dense_arrays.artifacts.bundles import BundleSummary  # noqa: PLC0415
    from dense_arrays.artifacts.pool_records import PoolSummary  # noqa: PLC0415
    from dense_arrays.reporting.selections import SelectionSnapshot  # noqa: PLC0415

    if json_output:
        typer.echo(canonical_json(result.to_dict()))
    elif isinstance(result, SelectionSnapshot):
        typer.echo(
            f"Selection: {result.selected} / {result.requested} requested; "
            f"{result.available} available; {result.shortfall} shortfall "
            f"({result.status})."
        )
    elif isinstance(result, BundleSummary):
        typer.echo(
            f"Selected collection: {result.designs} designs; "
            f"verified: {result.verified}."
        )
    elif isinstance(result, PoolSummary) and result.preparation is not None:
        typer.echo(
            f"{result.state.capitalize()}: {result.retained_parts} / "
            f"{result.preparation.requested_retention} requested parts retained; "
            f"{result.source_parts} candidates processed; verified: {result.verified}."
        )
    elif isinstance(result, PoolSummary):
        typer.echo(
            f"{result.state.capitalize()}: {result.retained_parts} / "
            f"{result.source_parts} supplied parts retained; "
            f"verified: {result.verified}."
        )
    else:
        typer.echo(
            f"{result.state.capitalize()}: {result.accepted} / {result.target} "
            f"designs; {result.termination_reason}.\n"
            f"Attempts: {result.counts['started']}; "
            f"verified: {result.verified}; resumable: {result.resumable}."
        )
    if not json_output and isinstance(result, (PoolSummary, RunSummary)):
        _display_producer(result.producer)
    if not json_output and isinstance(result, RunSummary) and result.cells:
        _display_cells(result)


def _display_cells(summary: RunSummary, *, limit: int = 10) -> None:
    """Expose attainment and inactive cells without printing an unbounded matrix."""
    for cell in tuple(summary.cells.values())[:limit]:
        typer.echo(
            f"{cell.cell_id}: {cell.accepted} / {cell.target}; "
            f"{cell.state}; {cell.termination_reason or 'in progress'}"
        )
    if len(summary.cells) > limit:
        typer.echo(
            f"{len(summary.cells) - limit} more cells; use --json for all cells."
        )


def _display_producer(producer: Producer | None) -> None:
    """Keep runtime context concise; complete recorded versions remain in JSON."""
    if producer is None:
        typer.echo("Producer: not recorded.")
        return
    solver = (
        "solver not recorded"
        if producer.solver is None
        else f"{producer.solver.name}: {producer.solver.version}"
    )
    typer.echo(
        f"Producer: Dense Arrays {producer.package_version or 'version unavailable'}; "
        f"{producer.python_implementation} {producer.python_version}; {solver}."
    )


def display_preparation_plan(result: PreparationPlan) -> None:
    """Show requested retention separately from unknown sampled yield."""
    preview = result.preview
    if result.sampled:
        typer.echo(
            f"Preparation: at most {preview['candidate_budget']} candidates; "
            f"requested retention {preview['requested_retention']}; retained "
            f"count unknown.\n"
            f"Required external tools: "
            f"{', '.join(preview['required_tools']) or 'none'}."
        )
        if "sampled_length" in preview:
            length = preview["sampled_length"]
            distribution = (
                "uniform length prior conditioned on constraints"
                if length["distribution"] == "uniform_prior_conditioned_on_constraints"
                else "uniform integer draws"
            )
            typer.echo(
                f"Candidate lengths: {length['minimum']}..{length['maximum']} "
                f"bases, inclusive; {distribution}."
            )
        if "mining_target" in preview:
            _display_mining_target(preview["mining_target"])
        _display_pool_limit(preview)
        if "motif_import" in preview:
            imported = preview["motif_import"]
            typer.echo(
                f"Motif input: {imported['format']}; supplied score matrix: "
                f"{imported['score_matrix']}."
            )
        if "recipes" in preview:
            _display_preparation_recipes(preview["recipes"])
            typer.echo(
                "Across recipes: sequence collisions "
                f"{preview['sequence_collisions']}; observed core collisions "
                f"{preview.get('core_collisions', 'preserve')} "
                "(parts without a core are excluded)."
            )
        if "motif_window" in preview:
            _display_motif_window(preview["motif_window"])
        if "screening_windows" in preview:
            typer.echo(
                f"Windowed exclusion motifs: {len(preview['screening_windows'])}; "
                "use --json for coordinates and information."
            )
        if "proposal" in preview:
            _display_proposal(preview["proposal"])
    else:
        typer.echo(
            f"Preparation: {preview['retained_parts']} / "
            f"{preview['source_parts']} supplied parts retained.\n"
            "Retention count: exact; required external tools: none."
        )


def _display_proposal(proposal: dict[str, object]) -> None:
    """Explain declared proposal support and bounded construction before execution."""
    typer.echo(
        f"Proposal: {proposal['strategy']}; motif placement: "
        f"{proposal['motif_placement']}.\n"
        f"Sampling ACGT probabilities: {proposal['base_probabilities']} "
        f"from {proposal['background_source']}."
    )
    if "construction_limits" in proposal:
        limits = proposal["construction_limits"]
        typer.echo(
            "Background probabilities conditioned on constraints; "
            f"{len(proposal['constraint_ids'])} sequence rules.\n"
            f"Construction limits: {limits['states']} states; "
            f"{limits['automaton_states']} automaton states; "
            f"{limits['mass_bits']} stored mass bits; "
            f"{limits['seconds']:g}s cooperative allowance."
        )


def _display_motif_window(window: dict[str, object], *, indent: str = "") -> None:
    """Show the selected source interval and its declared information reference."""
    typer.echo(
        f"{indent}Motif window: [{window['start']}, {window['end']}) in source; "
        f"{window['information_bits']:.6g} / "
        f"{window['source_information_bits']:.6g} information bits "
        f"relative to {window['background_source']} background."
    )


def display_pool_quality(result: PoolQualityReport, *, json_output: bool) -> None:
    """Show attainment, every preparation stage and the exact stopping reason."""
    typer.echo("Read cost: " + canonical_json(result.cost.to_dict()), err=True)
    value = result.to_dict()
    if json_output:
        typer.echo(canonical_json(value))
        return
    if isinstance(result, PoolQualitySnapshot):
        typer.echo("Recorded pool report; candidate evidence is not included.")
    counts = value["counts"]
    typer.echo(
        f"{value['state'].capitalize()}: {counts['retained']} / "
        f"{value['requested_retention']} requested parts retained.\n"
        f"Processed: {counts['processed']}; rejected: "
        f"{counts['eligibility_rejected']}; execution errors: "
        f"{counts['execution_error']}.\n"
        f"Eligible: {counts['eligible']}; duplicate discarded: "
        f"{counts['duplicate_discarded']}; eligible unique: "
        f"{counts['eligible_unique']}; not selected: {counts['not_selected']}.\n"
        f"Stopped: {value['stop_reason']}."
    )
    if "mining_target" in value:
        _display_mining_outcome(value)
    _display_score_bands(value)
    _display_construction(value)
    if value["rejections"]:
        reasons = sorted(value["rejections"].items())
        displayed = reasons[:20]
        typer.echo(
            "Rejection reasons: "
            + ", ".join(f"{name}={count}" for name, count in displayed)
        )
        if len(reasons) > len(displayed):
            typer.echo(
                f"{len(reasons) - len(displayed)} more reasons; "
                "use --json for all counts."
            )
    recipe_limit = 10
    for recipe in value.get("recipes", [])[:recipe_limit]:
        account = recipe["accounting"]
        typer.echo(
            f"{recipe['id']}: {account['counts']['retained']} / "
            f"{account['requested_retention']} retained; "
            f"stopped: {account['stop_reason']}."
        )
        if "mining_target" in account:
            _display_mining_outcome(account)
        _display_score_bands(account)
        _display_construction(account)
        _display_retention(account)
    if len(value.get("recipes", [])) > recipe_limit:
        typer.echo("Use --json for all per-recipe accounting.")
    _display_retention(value)


def _display_construction(account: dict[str, object]) -> None:
    """Distinguish a recorded zero-support proof from unfinished counting."""
    report = account.get("construction")
    if report is None:
        return
    outcomes = {
        "feasible": "valid sampling support counted",
        "infeasible": (
            "no sequence in the declared sampling support meets the compiled rules"
        ),
        "limited": f"feasibility unknown; reached {report['reason']} limit",
    }
    typer.echo(
        f"Conditional construction: {outcomes[report['status']]}. "
        f"Work: {report['states']} states, {report['mass_bits']} stored mass bits."
    )


def _display_retention(account: dict[str, object]) -> None:
    """Show pool admission and target-relative choice supply separately."""
    if account.get("retention") is not None:
        selection = account["retention"]
        typer.echo(
            f"MMR pool: {selection['pool_size']}; below retention score cutoff: "
            f"{selection['below_score']}; beyond pool limit: "
            f"{selection['beyond_limit']}."
        )
        if sizing := selection.get("sizing"):
            typer.echo(
                f"MMR choice pool: {selection['pool_size']} / {sizing['limit']}; "
                f"requested before cap: {sizing['requested']}; "
                f"available above cutoff: {sizing['available']}; "
                f"has alternatives: {'yes' if sizing['has_choice'] else 'no'}."
            )


def _display_pool_limit(preview: dict[str, object]) -> None:
    """State the resolved selection bound before candidate supply is known."""
    retention = preview.get("retention", {})
    if "pool_limit" in retention:
        typer.echo(
            f"MMR choice-pool limit: {retention['pool_limit']}; "
            f"distance work bound: {retention['distance_work_bound']} core bases."
        )


def display_candidate(record: PoolCandidate) -> None:
    """Summarize one decision without dumping sequences or scoring observations."""
    value = record.candidate
    details = (
        (
            f"recipe={value.recipe_id}; local index={value.recipe_index}; "
            if value.recipe_id is not None
            else ""
        )
        + f"representative={value.representative}; rank={value.rank}; "
        + (f"score band={value.score_band}; " if value.score_band is not None else "")
        + f"reasons={','.join(value.reasons[:3]) or '-'}; error={value.error or '-'}"
    )
    width = 180
    if len(details) > width:
        details = details[:width] + "..."
    typer.echo(
        f"Candidate {value.index}: {record.outcome}; "
        f"{len(value.part.sequence)} bases.\n"
        f"  {details}\n  Full sequence and decision evidence: --json"
    )


def _display_score_bands(account: dict[str, object]) -> None:
    """Show recipe-local empirical bands with bounded output and explicit units."""
    report = account.get("score_bands")
    if report is None:
        return
    typer.echo(
        f"Score bands: {report['total']} eligible unique; FIMO log2-odds; "
        f"scoring {report['scoring_id'][:12]}. Boundary ties stay together."
    )
    limit = 10
    for band in report["bands"][:limit]:
        scores = band["scores"]
        span = "empty" if scores is None else f"{scores['min']:g} to {scores['max']:g}"
        typer.echo(
            f"  Band {band['band']}: {band['count']} eligible; "
            f"{band['retained']} retained; score {span}."
        )
    if len(report["bands"]) > limit:
        typer.echo("Use --json for all score bands.")


def _display_mining_outcome(account: dict[str, object]) -> None:
    """Distinguish retained attainment from the declared eligible-supply target."""
    target = account["mining_target"]
    typer.echo(
        f"Mining target: {account['counts']['eligible_unique']} / "
        f"{target['eligible_unique']} eligible unique; "
        f"{'met' if target['met'] else 'unmet'}; "
        f"minimum candidates {target['minimum_candidates']}."
    )


def _display_mining_target(target: dict[str, int]) -> None:
    """State the resolved supply goal and batch-level stopping rule."""
    typer.echo(
        f"Mining target: {target['eligible_unique']} eligible unique; "
        f"minimum {target['minimum_candidates']} processed candidates; "
        "checked after each scored batch."
    )


def _display_preparation_recipes(recipes: tuple[dict[str, object], ...]) -> None:
    """Keep multi-recipe previews bounded while showing each supply target."""
    recipe_limit = 10
    for recipe in recipes[:recipe_limit]:
        typer.echo(
            f"{recipe['id']}: at most {recipe['candidate_budget']} candidates; "
            f"requested retention {recipe['requested_retention']}."
        )
        if "mining_target" in recipe:
            _display_mining_target(recipe["mining_target"])
        _display_pool_limit(recipe)
        if "motif_window" in recipe:
            _display_motif_window(recipe["motif_window"], indent="  ")
    if len(recipes) > recipe_limit:
        typer.echo("Use --json for the complete recipe preview.")
