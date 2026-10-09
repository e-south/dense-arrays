---
title: Code map
description: Find the module, contract, and tests that own a change.
---

# Code map

Dense Arrays separates optimization, persisted placement contracts, and
presentation. Find the owner of the behavior before editing:

| Surface | Owns |
| --- | --- |
| `problem`, `constraints` | Immutable packing inputs and requirement validation |
| `model`, `optimizer`, `errors` | Model construction, solving, enumeration, and distinct solver outcomes |
| `greedy` | Greedy selection with explicit entry occurrences and orientations |
| `solution` | Immutable `DenseArray` results and terminal representation |
| `parts/`, `planning/` | Part identities, input parsing, preparation policies, typed requirements and immutable previews |
| `generation/` | Packing translation, bounded assembly, versioned random streams and independent final checks |
| `artifacts/` | Native records, exclusive output ownership, transactional revisions and publication |
| `reporting/` | Snapshot queries, diagnostics, quality metrics, integrity verification and native-to-playback projection |
| `workflow/` | Shared operations, execution sequencing and thin CLI translation |
| `diagnostics` | Immutable diagnostic records shared by validation and inspection |
| `realized`, `_record_validation` | Producer-neutral placements, JSON snapshots, and sequence alignment |
| `playback.models`, `playback.validation` | Supported plan fields and cross-record geometry/evidence checks |
| `playback.reconstruction`, `playback.serialization` | Deterministic compilation and strict JSON boundaries |
| `playback.presentation`, renderers | Document choices, visible evidence, NetworkX graph layout, stills, and video |

Optimization returns a `DenseArray`; it does not automatically serialize a
`RealizedArray`. Producer adapters own that translation. Study identities,
biological labels, selected records, and publication interpretation remain
caller-owned.

## Find the files for a change

Paths below are relative to `src/dense_arrays/` and `tests/`, respectively.
Public entrypoints stay small; helper modules own the listed decisions.

| Task | Implementation under `src/dense_arrays/` | Tests under `tests/` |
| --- | --- | --- |
| Packing inputs, counts, ranges, or mutation checks | `problem.py`, `constraints.py`, `optimizer.py` | `test_core_contracts.py`, `test_optimize.py` |
| Exact model construction | `model.py` | `test_optimize.py`, `test_core_contracts.py` |
| Solver status or enumeration failure | `optimizer.py`, `errors.py` | `test_solver_outcomes.py`, `test_cli.py` |
| Greedy selected entries and occurrences | `greedy.py` | `test_greedy.py` |
| Workflow search methods, heuristic outcomes and replay | `generation/heuristic.py`, `generation/packing.py`; `artifacts/search.py`; `reporting/search.py`; `workflow/execution.py` | `workflow/test_greedy_search.py` |
| Empirical score bands, reconciled summaries and recipe-local queries | `parts/retention/bands.py`; `artifacts/preparation/bands.py`; `reporting/pools/filters.py` | `workflow/preparation/test_score_bands.py` |
| Target-relative MMR choice pools and admission accounting | `parts/retention/pool.py`, `parts/retention/mmr.py`; `artifacts/preparation/retention.py` | `workflow/preparation/test_pool_sizing.py` |
| Bounded named-window expansion into ordinary preparation sets | `parts/preparation/sets.py`; `planning/preparation/requests.py`; `workflow/inputs.py` | `workflow/preparation/test_window_sets.py` |
| Cross-recipe retained sequence and observed-core collisions | `parts/retention/collisions.py`; `parts/preparation/sets.py`; `planning/preparation/sets.py`; shared preparation execution and verification | `workflow/preparation/test_core_collisions.py`, `test_sets.py` |
| Sequence overlap calculations | `sequence.py` | `test_optimize.py` |
| Result construction or terminal layout | `solution.py` | `test_core_contracts.py`, `test_optimize.py`, `test_cli.py` |
| Optimizer command inputs or errors | `cli.py` | `test_cli.py` |
| Curated input, workflow planning and persisted generation | `parts/`, `planning/`, `generation/`, `workflow/` | `workflow/` |
| Table row errors and impossible occurrence counts | `parts/tables/diagnostics.py`, `parts/ingestion.py`; `planning/diagnostics.py`, `validation.py`; `diagnostics.py` | `workflow/inputs/test_import_diagnostics.py`, `workflow/test_planning_diagnostics.py` |
| Command help, examples and feature requirements | `cli.py`; `workflow/cli.py`, `export_cli.py` | `workflow/test_help.py` |
| Bounded matrix expansion, source selections, substitutions, requirement additions and allocation | `planning/matrices/`, `resolution.py`, `serialization.py`; `parts/bound.py`; `workflow/inputs.py`, `presentation.py` | `workflow/test_matrix_planning.py`, `workflow/test_matrix_sources.py`, `workflow/test_legacy_fixture.py`, `workflow/test_curated_example.py` |
| Matrix execution, cell state and shared effort | `workflow/matrices.py`, `execution.py`; `artifacts/run_state.py`, `run_plans.py`, `store.py`; `reporting/verification.py` | `workflow/test_matrix_execution.py` |
| Candidate-batch policies, compact models and executable replay | `planning/batches/`; `generation/batches/`; `workflow/batches.py`; `artifacts/store.py` | `workflow/test_batches.py` |
| Joint batch uniqueness, core orientation and group caps | `planning/batches/validation.py`; `generation/batches/constrained.py`; `parts/models.py` | `workflow/test_batch_eligibility.py` |
| Runtime resampling, feedback, saved decisions and recovery | `planning/batches/resampling.py`; `planning/batches/feedback.py`; `workflow/batch_cursor.py`; `artifacts/batches/`; `reporting/runtime_batches.py` | `workflow/resampling/` |
| Ordered batch schedules, local limits and replay frontiers | `planning/batches/schedules.py`; `workflow/schedules.py`; `reporting/batch_accounting.py` | `workflow/test_batch_schedules.py`, `workflow/test_batch_limits.py` |
| Curated preparation, pool reuse and import evidence | `parts/preparation/`, `provenance.py`, `pools.py`; `planning/preparation/`; `artifacts/pools.py` | `workflow/test_preparation.py`, `workflow/test_provenance.py` |
| CSV/TSV, Parquet and Excel part inputs | `parts/ingestion.py`, `tables/`, `models.py`, `serialization.py` | `workflow/inputs/test_table_formats.py` |
| Independent recipe sets and combined pool accounting | `parts/preparation/sets.py`; `planning/preparation/sets.py`; `artifacts/preparation/sets.py`; `workflow/preparation.py` | `workflow/preparation/test_sets.py` |
| Sampled recipes, mining targets, pool decisions and stage accounting | `parts/sampling.py`, `mining.py`, `eligibility.py`, `screening/`, `retention/`; `planning/preparation/`; `workflow/preparation.py`; `artifacts/preparation/`; `reporting/pools/` | `workflow/preparation/test_sampled.py`, `test_mining_targets.py`, `test_mmr.py`, `test_pwm_exclusion.py`, `test_reading.py`, `test_proposals.py`, `test_lengths.py` |
| Conditional background support, exact sampling and resource outcomes | `parts/background/`; `parts/sampling.py`; `planning/preparation/sampled.py`; `workflow/preparation.py`; `artifacts/preparation/` | `workflow/preparation/test_conditional_counts.py`, `test_conditional_workflow.py` |
| Motif inputs, format parsing, score units and optional FIMO execution | `parts/motifs/`; `parts/motifs/formats/`; `parts/scoring/` | `workflow/preparation/test_motifs.py`, `workflow/preparation/test_motif_formats.py`, `workflow/preparation/test_fimo.py`, `workflow/preparation/test_windows.py` |
| Bounded record reads, snapshots and cursors | `artifacts/reading.py`, `cursors.py`; `reporting/readers.py`, `filters.py` | `workflow/test_read_costs.py`, `workflow/test_pagination.py`, `workflow/test_attempt_filters.py` |
| Recorded producer and solver versions | `artifacts/provenance.py`, `store.py`, `pool_records.py`; `solver.py`, `optimizer.py`; `reporting/summary.py` | `workflow/test_execution_provenance.py` |
| Saved candidates, rejection recount and packing restoration | `artifacts/candidates.py`, `records.py`; `generation/packing.py`, `assembly.py`; `reporting/candidates.py`, `verification.py` | `workflow/test_candidate_evidence.py` |
| Part-usage packing preference and objective replay | `model.py:part_usage_weights`; `generation/objectives.py`; `artifacts/objectives.py`; `reporting/objectives.py` | `workflow/test_packing_preference.py` |
| Design predicates, sequence/placement joins and portable text exports | `reporting/design_filters.py`, `design_queries.py`, `projections.py`, `exporting/records.py`, `exporting/source_bindings.py`; `artifacts/publication.py` | `workflow/test_design_queries.py`, `workflow/test_exports.py`, `workflow/test_export_bindings.py` |
| Editable request, resolved plan and scoped report exports | `reporting/exporting/documents.py`, `formats.py`; `workflow/exporting.py`; `planning/serialization.py`, `preparation/` | `workflow/test_document_exports.py` |
| Immutable request edits, bound parts, lineage and library exclusions | `parts/bound.py`; `planning/models.py`, `resolution.py`, `lineage.py`, `libraries.py`; `reporting/plans/requests.py`, `reading.py`; `workflow/extensions.py` | `workflow/test_request_revision.py` |
| Portable selected collections and contained verification | `planning/evidence.py`; `artifacts/bundles/`; `reporting/bundles/`, `exporting/bundles.py`; `workflow/bundle_inspection.py` | `workflow/test_bundles.py` |
| Total/cell quotas, saved membership and selected projections | `reporting/selections/`; `workflow/selections.py` | `workflow/test_selection_contracts.py`, `test_selections.py` |
| Combined native/bundled libraries, namespace resolution and bounded identity unions | `reporting/collections/`; shared `readers.py` and `projections.py` | `workflow/test_collections.py`, `workflow/test_portable_collections.py`, `workflow/test_extension.py` |
| Query flag translation and pure data stdout | `workflow/queries.py`, `export_cli.py`, `cli_errors.py` | `workflow/test_exports.py`, `workflow/test_cli.py` |
| Single-cell/matrix additional targets, parent exclusions and plan comparison | `planning/extension.py`; `workflow/extensions.py`, `plans.py`; `artifacts/store.py`; `reporting/plans/`, `verification.py` | `workflow/test_extension.py`, `workflow/test_matrix_extension.py`, `workflow/test_plan_comparison.py`, `workflow/test_matrix_comparison.py` |
| Included plan inventories and full-ID selection | `reporting/plans/bundles.py`, `filters.py`; `reporting/bundles/`; `workflow/bundle_inspection.py` | `workflow/test_bundle_plans.py` |
| Explain search outcomes, available histories and selected or combined quality | `reporting/diagnostics.py`, `accounting.py`, `quality/`, `metrics.py` | `workflow/test_diagnostics.py`, `workflow/test_quality.py`, `workflow/test_selected_quality.py`, `workflow/test_portable_collections.py` |
| Compare quality or reuse saved metrics | `reporting/quality/comparison.py`, `differences.py`, `snapshots.py`, `validation.py`; `workflow/quality.py` | `workflow/test_quality_comparison.py` |
| Display reports without changing policy | `workflow/presentation.py`; `playback/quality/` | `workflow/test_quality.py`, `workflow/test_diagnostics.py` |
| Plot sampled-pool yield, saved selection distances and score bands | `reporting/pools/quality.py`, `diversity.py`, `snapshots.py`; `playback/quality/preparation.py`; `reporting/rendering.py` | `workflow/preparation/test_rendering.py`, `test_diversity_reports.py` |
| Identity-safe fixed geometry and exact packing length | `constraints.py`, `problem.py`, `model.py`, `optimizer.py` | `packing/test_fixed_occurrences.py`, `packing/test_length_coordinates.py` |
| Final assembly, padding and sequence screens | `generation/assembly.py`, `geometry.py`, `randomness.py`, `screening.py` | `workflow/test_assembly.py`, `workflow/test_screening.py`, `workflow/test_geometry.py` |
| Shared literal matching, detailed observations and bounded acceptance predicates | `parts/screening/sequence.py`; preparation eligibility and verification | `workflow/preparation/test_sequence_predicate.py`, `workflow/test_screening.py` |
| Render one selected design from a run, bundle or saved panel | `reporting/rendering.py`; shared record, collection and selection readers; existing playback owners | `workflow/test_render.py`, `workflow/test_render_selection.py` |
| Admit candidate bases and offered packing-model dimensions before allocation | `parts/sampling.py`; `planning/models.py`, `batches/bindings.py`; preparation and generation execution boundaries | `workflow/preparation/test_resource_admission.py`, `workflow/test_model_admission.py` |
| Transaction failures, corruption or competing writers | `artifacts/`, `reporting/` | `workflow/test_integrity.py`, `workflow/test_process_ownership.py` |
| Resume admission, cell replay, stable locking and consumed budgets | `artifacts/recovery.py`, `store.py`; `workflow/recovery.py`, `execution.py`, `matrices.py` | `workflow/test_recovery.py`, `workflow/test_matrix_recovery.py`, `workflow/test_recovery_processes.py` |
| Playback CLI preflight and publication | `playback/cli.py`, `playback/output.py` | `test_playback_cli.py` |
| Realized fields, nested provenance, or placement alignment | `realized.py`, `_record_validation.py` | `test_playback_contracts.py` |
| Plan JSON and evidence/geometry validation | `playback/models.py`, `playback/validation.py`, `playback/serialization.py` | `test_playback_contracts.py` |
| Coordinate reconstruction and caller notices | `playback/reconstruction.py` | `test_playback_contracts.py`, `test_playback.py` |
| Document labels, colors, and visible evidence | `playback/presentation.py`, `playback/theme.py` | `test_playback_presentation.py` |
| Graph projection, selected relations, layout, routing | `playback/graph/`, `playback/graph_drawing.py` | `test_playback_graph.py`, plus rendered stills |
| Raster scene and sequence frames | `playback/scene_drawing.py`, `playback/duplex_drawing.py`, `playback/duplex_frames.py` | `test_playback.py`, `test_playback_presentation.py` |
| Native nucleotide cells, glyph centering, and shared typography | `playback/duplex_geometry.py`, `playback/typography.py` | `test_playback_typography.py`, `test_playback_resting.py` |
| Frame timing, writers, and figure cleanup | `playback/frame_schedule.py`, `playback/export.py`, `playback/matplotlib_renderer.py` | `test_playback_exports.py` |
| Dependency import boundaries | `__init__.py`, `playback/__init__.py`, `playback/graph/__init__.py` | `test_optional_playback_imports.py` |
| Runnable documentation and built links | `docs/`, `mkdocs.yml` | `test_documentation.py` |

Use `rg -n 'def |class ' <file>` to locate symbols. Run focused checks from the
repository root with
`uv run --no-sync pytest -q tests/<file>`, then finish with the
[full development gate](../development.md#local-verification).

## Change boundaries

Input validation is independent of solver execution. The optimizer holds an
immutable problem and delegates model construction and greedy realization.
A `DenseArray` can describe a valid contiguous layout more broadly than an
exact path; `forbid()` checks the narrower exact-path requirement separately.

Plan validation and evidence projection are independent of drawing. Graph
projection is independent of raster geometry. Default and injected graph
layout engines receive the same selected topology. Raster scene drawing,
producer-frame adaptation, timing, and writer lifecycle have separate owners;
public export functions compose them.

Keep these boundaries when adding behavior. A useful module owns a decision
that can change without editing unrelated solver or renderer code.

## Continue by task

- [Saved-library guide](../library-workflow.md): prepare parts, generate libraries and inspect saved evidence.
- [Library workflow architecture](library-workflow/index.md): domain, operation and native evidence contracts.
- [Playback contract](solution-playback.md): coordinates, authority, and producer handoffs.
- [Product brief](animation-product-spec.md): visual and publication goals that need output review.
- [Playback guide](../playback.md): runnable entrypoints.
- [API reference](../api.md): inputs, results, and failures.
- [Caller migrations](../migration.md): deliberate input and integration changes.
- [Development](../development.md): repository verification and release checks.
