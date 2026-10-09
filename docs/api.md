---
title: Python API
description: Choose operations and typed records for preparing parts, generating libraries, inspecting evidence and rendering placements.
---

# Python API

Use the six operations below for persisted libraries. Their Python and CLI
interfaces share requests, validation and Dense Arrays records. Start with the
[saved-library guide](library-workflow.md) for a complete example, or the
[first-array tutorial](quickstart.md) for direct optimization in memory.

## Library operations

```python
import dense_arrays as da
from dense_arrays import parts, planning

# Two synthetic 16-base sites share eight bases and fit within 24 bases.
request = planning.DesignSpec(
    parts=(
        parts.Part("upstream", "ACGTTGCAAGTCCTGA"),
        parts.Part("downstream", "AGTCCTGATCGTACCG"),
    ),
    length=planning.Length(maximum=24),
)
# Planning validates the typed request without generating sequences.
resolved = da.plan(request)
# A new destination retains accepted designs and the attempt history.
run = da.run(resolved, out="runs/library")
# Inspection reads that saved evidence; it does not repeat generation.
quality = da.inspect(run, view="quality")
print(quality.to_dict()["attainment"])  # Read accepted and requested counts.
```

Run this example in a new working directory with the
[library workflow installed](installation.md#use-the-library-workflow).
It creates `runs/library`. For file inputs and reusable pools, follow
[the curated-parts guide](library-workflow/preparation/curated.md).

| Operation | Main typed inputs | Result and next task |
| --- | --- | --- |
| `prepare` | `parts.PreparationSpec`, `parts.PreparationSet`, or a resolved preparation plan | A reusable pool; [prepare parts](library-workflow/preparation.md). A generation or matrix plan with `planning.BatchSampling` instead freezes candidate batches. |
| `plan` | Preparation, design, extension or matrix request | `PreparationPlan`, `GenerationPlan` or `MatrixPlan`; [preview requests](library-workflow/handoffs.md). |
| `run` | `planning.DesignSpec`, `planning.ExtensionSpec`, `planning.MatrixSpec`, or a resolved generation plan | `artifacts.RunHandle`; [generation and assembly](library-workflow/generation/assembly.md). `run(resume=...)` continues an unchanged eligible run. |
| `inspect` | Saved artifact, typed report, or ordered source collection | A summary, record view, comparison or selection snapshot; [inspect results](library-workflow/results/quality.md). |
| `export` | Saved evidence and explicit output scope | `artifacts.ExportReceipt`; [export records](library-workflow/results/export.md) or [editable documents](library-workflow/handoffs.md). |
| `render` | Saved design evidence or quality report | `artifacts.ExportReceipt`; [choose a figure](reference/outputs.md). |

The [operation contract](architecture/library-workflow/operations.md) lists
signatures, effects, stream ownership and failure behavior. Specialized types
live in `parts`, `planning`, `artifacts` and `reporting`.

## Requests and generation

| Task | Types and route |
| --- | --- |
| Set length, requirements, padding and effort | `planning.DesignSpec`, `Length`, `Limits` and typed requirements; [assembly](library-workflow/generation/assembly.md), [constraints](constraints.md). |
| Choose exact or greedy search and part-use preferences | `DesignSpec.search`, `packing_preference`; [search methods](library-workflow/search.md). |
| Expand combinations and allocate targets | `MatrixSpec`, `Variant`, `Allocation`, `MatrixPlan`; [matrices](library-workflow/matrices.md). |
| Freeze offered parts or resample during generation | `BatchSampling`, `CandidateBatch`, `BatchSchedule`, `Resampling`, `FeedbackPolicy`; [candidate batches](library-workflow/batches.md). |
| Change a request or compare resolved effects | `reporting.RequestReport`, `planning.Lineage`, `parts.BoundParts`, `reporting.PlanComparison`; [document handoffs](library-workflow/handoffs.md). |
| Add designs with unchanged requirements | `ExtensionSpec`, `ParentRun`; [extensions](library-workflow/extension.md). |
| Exclude a declared accepted library | `planning.LibraryExclusion`; [exclusions and lineage](library-workflow/handoffs.md). |
| Continue an interrupted run | `run(resume=...)`, `artifacts.recovery.RecoveryError`; [recovery](library-workflow/recovery.md). |

`dense_arrays.RunExecutionError` carries an inspectable `.run` and its `.artifact`
when execution fails after saving the run state. Its cause preserves the
original exception. Input errors before run creation retain their original types;
keyboard interruption remains `KeyboardInterrupt`.

When an occurrence minimum exceeds the available parts, `planning.PlanningError`
carries a `.diagnostic` with the requirement ID, requested minimum, available
count and input references. Direct table inputs include one-based data rows.
The CLI returns the same record under `diagnostic` with `--json` and exits with
status 2. Add matching parts or revise the minimum before generating.

## Inspect, select and share

Record views (`RecordView`, `LibraryView`, `BundleView` and `SelectionView`)
expose `.records()`, `.cost` and `.sources`. Iterate inside a context manager
so early termination releases resources. Their filters and cursors remain tied
to the selected source revisions. `reporting.ReadLimits` bounds records, pair
work and identity entries; a displayed page size does not bound report scope.

| Task | Types and route |
| --- | --- |
| Explain attainment and search loss | `DiagnosticReport`, `QualityReport`, `AttemptFilter`; [shortfalls and quality](library-workflow/results/quality.md). |
| Read accepted sequences and placements | `DesignFilter`, `SequenceRecord`, `PlacementRecord`; [record selection and exports](library-workflow/results/export.md). |
| Select a seeded panel or per-cell quota | `LibrarySelection`, `Take`, `SelectionSnapshot`; [saved panels](library-workflow/selection.md). |
| Compare quality with explicit populations | `QualityComparison`, `MetricDifference`, `QualitySnapshot`; [quality comparisons](library-workflow/handoffs.md#compare-library-quality). |
| Move a library with its required evidence | `BundleSummary`, `BundleView`, `planning.PlanEvidence`, `PlanFilter`; [portable libraries](library-workflow/bundles.md). |
| Check schemas, producers and committed evidence | `RunSummary`, `PoolSummary`, `artifacts.provenance.Producer`; [native evidence](architecture/library-workflow/artifacts.md). |

A saved `SelectionSnapshot` can be reused by inspection, export and rendering
without drawing again. Library-quality reports distinguish the selected design
population from each source run's original attainment and effort. Per-design
rendering requires exactly one selected design, including when its source is a
bundle. [Outputs](reference/outputs.md) lists figures and useful summaries.

Attempt JSON can include `artifacts.Attempt.candidate`: a `CandidateEvidence`
record with packed and final realized arrays. A missing final array means
assembly stopped before final evaluation; absent candidate evidence remains
unknown. Accepted-design exports retain their own population.

## Part preparation

| Task | Types and route |
| --- | --- |
| Import curated parts and annotations | `parts.Part`, `PartTable`, `PartFilter`; [table inputs](library-workflow/tables.md). CSV/TSV use the base install; Parquet/XLSX require the `tables` extra. |
| Diagnose invalid input rows | `parts.TableImportError`, `RowDiagnostic`; [table diagnostics](library-workflow/tables.md#correct-invalid-rows). |
| Sample motif candidates and retain eligible parts | `PWMArtifact`, `PreparationSpec`, `Sampling`, `CandidateBudget`, `FimoScoring`, `Eligibility`, `Uniqueness`, `Retention`; [sampled pools](library-workflow/preparation.md). |
| Declare proposal lengths and strategies | `Sampling`, `LengthRange`, `planning.Length`; [proposal strategies](library-workflow/preparation/motifs.md#choose-a-pwm-proposal-strategy) and [candidate lengths](library-workflow/preparation/windows.md#vary-candidate-length). |
| Generate constrained backgrounds | `Background`, `ConditionalLimits`; [background sampling](library-workflow/background.md). |
| Stop after enough eligible candidates | `MiningTarget`; [mining targets](library-workflow/preparation/effort.md#stop-after-enough-eligible-candidates). Candidate effort remains independently capped. |
| Balance score and core distance | `MMR`, `PoolSize`; [retention](library-workflow/preparation/retention.md#retain-score-and-core-diversity). |
| Describe score distributions | `ScoreBands`; [score-band reports](library-workflow/preparation/retention.md#describe-the-eligible-score-distribution). |
| Combine independently budgeted recipes | `PreparationSet`; [multiple recipes](library-workflow/preparation/sets.md#prepare-several-recipes-together). |
| Expand named motif windows | `PreparationSet.from_windows`, `MotifWindow`; [named windows](library-workflow/preparation/windows.md#prepare-named-motif-windows). |
| Inspect recorded candidate decisions | `reporting.CandidateFilter`, `PoolQualityReport`, `PoolQualitySnapshot`; [candidate evidence](library-workflow/preparation/inspection.md#inspect-candidate-decisions) and [preparation figures](library-workflow/preparation-quality.md). |

`plan` previews a preparation request without sampling or scoring. Requested
retention and unknown yield remain separate. `prepare` publishes the pool;
`inspect(pool, view="quality")` reports reconciled stage counts. Saved preparation
quality JSON can be read and rendered without its original pool.

## Optimization

[Optimizer](reference/optimizer.md) supports direct in-memory packing: construct
a problem, add requirements, solve with CBC, enumerate arrangements, and handle
distinct failure outcomes. Use the workflow operations above when the task needs
persisted run accounting, recovery or a reusable library.

## Dense-array results

[DenseArray](reference/results.md) exposes the sequence, offsets, motif count and
compression ratio. `Optimizer` and `DenseArray` are exported by `dense_arrays`.

## Motif inputs and scoring

[Motif inputs and scoring](reference/motif-scoring.md) covers native motif JSON,
minimal MEME probabilities, JASPAR counts, score matrices and optional FIMO
scoring, with explicit geometry, units, background and resource limits.

## Sequence utilities

[Sequence utilities](reference/sequence.md) covers complements, pairwise overlaps
and sequence-display helpers.

## Realized arrays and playback

- [Realized arrays](reference/realized.md): describe a sequence and its placements.
- [Playback](reference/playback.md): reconstruct, serialize and render those
  placements through `dense_arrays.playback` and its rendering module.

For command-line use, see [CLI options and failures](reference/cli.md).
