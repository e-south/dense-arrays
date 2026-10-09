---
title: Choose a result view
description: Match a design question to a native summary, figure or portable export.
---

# Choose a result view

Use `inspect` to read evidence, `render` to make a figure and `export` to save
records for another tool. These operations reuse saved results. They do not
generate a replacement library or invoke a motif scorer.

For direct packing, [optimize motifs or enumerate arrangements](cli.md#optimize-motifs)
in the terminal, or use the [Python optimizer](optimizer.md).
Its `DenseArray` results provide sequence, selected-entry count, compression
ratio and strand-specific offsets without a persisted library workspace.

## Figures

| Question | View | What it shows |
| --- | --- | --- |
| How is one sequence assembled? | `render --view design` | Selected parts, strand, final coordinates, padding and supported requirement evidence. Select exactly one design. |
| What did the selected library contain? | `render --view library-quality` | Ranked part use, GC composition, packing density and available search outcomes. Design filters and saved selections define the population. |
| Where were preparation candidates lost? | `render --view preparation-quality` | Per-recipe yield, saved MMR distances to earlier selected cores and declared score bands. Missing metrics remain unavailable. |

All three views write PNG files with the optional `playback` dependencies.
[Design figures](../library-workflow/generation/assembly.md#assemble-and-render-an-exact-length-design)
preserve placement identities. [Library reports](../library-workflow/results/quality.md#explain-shortfalls-and-assess-a-library)
keep selected composition separate from source-run attainment.
[Preparation figures](../library-workflow/preparation-quality.md) keep scoring
models and recipe populations separate. Quality figures embed their report and
digest; exported quality JSON provides the same metrics for reuse.

For sequence-level playback stills and animations from realized placements, use
the [playback guide](../playback.md). Those views explain recorded geometry and
an order reconstructed from placement coordinates. They do not add missing
construction evidence.

## Summaries and records

| Question | Inspection view | Scope and next step |
| --- | --- | --- |
| What inputs, requirements and limits were used? | `request`, `plan`, `plans` | Read effective requests and bound plans; compare plans before revising or extending a library. |
| Did the work reach its target? | `summary` | Read attainment, stopping state and available recovery information. A matrix retains per-cell targets. |
| Why did search stop or reject candidates? | `diagnostics`, `attempts` | Inspect effort, outcome categories and saved rejection evidence. Exhausted effort does not establish infeasibility. |
| Which parts were offered to each search batch? | `batches` | Read recorded batch membership and sampling decisions from a run or bundle. |
| Which prepared candidates were retained? | `candidates`, `parts` | Follow eligibility, representatives, scores, core coordinates and recorded selection decisions. Candidate filters select recipes, outcomes or declared score bands. |
| What does one result contain? | `quality` | Summarize preparation yield and retention, or library composition and available search outcomes. |
| How do libraries differ? | `quality --compare` | Compare supported library metrics with their populations and denominators. This is separate from plan comparison. |
| Which sequences and placements should another tool receive? | `designs`, `sequences`, `placements` | Export native identities and coordinate joins; choose a movable bundle when the recipient needs complete verification evidence. |
| Which panel was selected? | Saved selection | Reuse its pinned source revisions and ordered design references for inspection, export and rendering. |

See [saved reports and comparisons](../library-workflow/handoffs.md),
[portable exports](../library-workflow/bundles.md) and
[reproducible selection](../library-workflow/selection.md) for paired Python and
CLI examples. `dense-arrays inspect --help`, `export --help` and `render --help`
list supported options and read limits. Live preparation-quality verification
may repeat MMR comparisons within its pair-work allowance; detached quality
reports read recorded metrics only.
