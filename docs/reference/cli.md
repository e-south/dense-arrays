---
title: Command-line reference
description: Choose saved-library, terminal optimization or playback commands and interpret their outputs.
---

# Command-line reference

Use the [saved-library workflow](../library-workflow.md) for persisted designs,
reusable requests and bounded generation. Use `optimize` or `solutions` for a
terminal result, and `dense-arrays-playback` to render a saved placement file.
Follow [installation](../installation.md) for the commands needed by your task.

## Saved-library workflow

The six operations prepare parts, resolve requests, generate designs, inspect
evidence, export results and render selected designs or library reports:

```bash
dense-arrays prepare --help
dense-arrays plan --help
dense-arrays run --help
dense-arrays inspect --help
dense-arrays export --help
dense-arrays render --help
```

Sampled preparation requests can declare `mining_target` with `eligible_unique`
or `max_retained_fraction`, plus optional `minimum_candidates`. `plan` shows the
resolved supply goal; `prepare` stops at the first attained batch or an effort
limit. `inspect --view quality` reports attainment. An unmet target exits **3**
even when the retained count is complete. See [mining targets](../library-workflow/preparation/effort.md#stop-after-enough-eligible-candidates).

`plan MATRIX --out PLAN --json` saves a [design matrix](../library-workflow/matrices.md)
with explicit cell targets. `run PLAN --out RUN --json` generates one native run
with per-cell results and a shared effort budget. `inspect RUN --verify --json`
checks the saved evidence and reports each cell's attainment. A matrix cell is
one named combination of axis choices; its target counts accepted designs.
The matrix `sources` mapping can select different tables or pool subsets for
each cell. See [matrix sources](../library-workflow/matrices.md#select-parts-for-each-combination).

`inspect BEFORE --view plan --compare AFTER` compares two saved matrix plans
or native runs by named combination. Add `--json` for all changes, or use
`export` with the same arguments and `--out changes.json` to save the report.
See [matrix comparison](../library-workflow/matrices.md#compare-planned-combinations).

`prepare PLAN --batch-size 16 --batch-strategy group_balanced --batch-seed 7
--out BATCH_PLAN` saves an executable [candidate batch](../library-workflow/batches.md).
Pass that plan to `run` to reuse its exact offered parts.
Combine `--unique-sequences`, `--unique-cores` and `--max-per-group N` to restrict
each offered batch. Core uniqueness requires complete core annotations and group
labels; group caps require group labels. Incompatible policies fail before output.
Add `--batch-count 3 --attempts-per-batch 10` to prepare an ordered search
schedule. Every batch shares the cell target and the run's global limits.
Use `--accepted-per-batch 2` with the attempt cap to advance after two accepted
designs; duplicates and rejections consume only the attempt allowance.

For editable inputs and JSON reports, see
[document exports](../library-workflow/handoffs.md). `export --view request|plan|summary|quality|diagnostics`
omits `--all`; record exports without a bounded selection require it.
For seeded or per-cell panels, use [saved selections](../library-workflow/selection.md).
For measured interruptions, use [run recovery](../library-workflow/recovery.md):
`run --resume RUN` preserves the original request and accepts no design overrides.
`export --view plan --compare OTHER` writes a generation-plan comparison;
`export --view quality --compare OTHER` writes aggregate metric differences with
both populations and denominators. `--out -` keeps document data on stdout and
diagnostics/receipts on stderr.

Runtime sampling is declared in a generation request's `resampling` field.
`plan` resolves the same policy used by Python; `run` saves each selected batch.
Use `inspect RUN --view batches --limit 10 --json` or
`export RUN --view batches --all --out batches.json` to read recorded decisions.
The [batch guide](../library-workflow/batches.md#resample-during-generation)
provides paired examples and feedback semantics.

For preparation evidence, `inspect POOL --view candidates` shows saved decisions.
Use repeatable `--candidate-index`, `--outcome`, `--reason` or set `--recipe-id`
to narrow the page;
`--json` includes full sequences, scores and representative links.
`export POOL --view candidates --all --out candidates.json` saves the selected
records. `inspect REPORT.json --view quality` reopens an exported pool report.
See [candidate inspection](../library-workflow/preparation/inspection.md#inspect-candidate-decisions)
for supported outcomes, limits and verification scope.

`render RUN --design-id REF --out design.png` renders one accepted design using
its persisted placements and cell plan. Use a full reference from
`inspect RUN --view designs --json`. Normal design filters and `--selection`
files are supported; the selection must contain exactly one design.
`render RUN --view library-quality --out quality.png` renders an aggregate
report. See [design rendering](../library-workflow/generation/assembly.md#assemble-and-render-an-exact-length-design)
for source types, dependencies and selection errors.
`render POOL_OR_REPORT --view preparation-quality --out preparation.png`
plots sampled-pool yield, recorded MMR distances and score bands. Use
`--max-pairs` to bound MMR verification work on a live pool. See
[preparation figures](../library-workflow/preparation-quality.md) for portable
reports and metric scope.

## Exit codes and machine output

Use `--json` for machine-readable plans, reports and errors. Plan and inspection
results go to stdout; human-readable receipts and diagnostics go to stderr.
`export --out -` writes the exported data to stdout.

| Exit code | Meaning |
| --- | --- |
| `0` | The operation succeeded. Inspecting a valid stopped run also succeeds; its report retains the stopped state. |
| `2` | Invalid input or a selection shortfall rejected by the requested policy. |
| `3` | Generation or preparation ended below its target, or a partial selection export was explicitly allowed. |
| `4` | Execution, integrity or read-limit failure. |
| `130` | Interrupted execution. |

Before streaming begins, JSON domain errors use `dense_arrays.error.v1`, with
`code`, `message`, `exit_code` and the affected artifact when known. An execution
failure after a run was saved includes its path. Python exposes the same saved
run through `RunExecutionError.run`, with the original exception as its cause.

After streaming begins, errors go to stderr and leave the data prefix on stdout.
Check the exit code before treating streamed output as complete.

## Optimize motifs

`optimize` and `solutions` return results directly to the terminal:

```bash
dense-arrays optimize --help
dense-arrays solutions --help
```

| Option | Meaning |
| --- | --- |
| `--motif` | One motif; repeat for each library entry |
| `--motifs-file` | One motif per line; replaces `--motif` |
| `--length` | Required positive integer sequence-length limit |
| `--strands` | `single` or `double`; defaults to `double` |
| `--solver` | Backend name passed to OR-Tools; defaults to `CBC` |
| `--solver-seconds` | Cooperative time limit for each solve; unset by default |
| `--solver-threads` | Explicit backend thread count; supported for SCIP only |
| `--max-solutions` | Maximum displayed results for `solutions`; defaults to 10 |
| `--diverse` | Bias `solutions` toward less represented motif entries |

Positional and regulator constraints use the [Python API](../constraints.md).
`--solver-seconds` applies to each solve; it is cooperative and can be exceeded
by the backend. It does not bound the total duration of enumeration.
`--max-solutions` limits result count only. Saved-library `run` requests also
support an accumulated active-time limit.
Output is a terminal display. It is not the persisted placement JSON expected
by playback.

Bad options, unreadable motif files, malformed inputs, no feasible first
result, and solver failures produce errors on stderr and a nonzero exit.
If a failure follows an already printed result, the command still exits
nonzero; preceding output does not imply that enumeration completed.
See [solver outcomes](optimizer.md#solver-outcomes) for the Python exception types.

## Render saved placements

`dense-arrays-playback` renders a supplied placement or playback file:

```bash
dense-arrays-playback --help
```

Supply a realized-array or playback-plan JSON file and at least one output:
`--poster` for PNG, `--mp4`, or `--gif`. Multiple formats can be requested in
one command. `--title` and `--subtitle` supply artifact metadata; the
[evidence reference](playback-presentation.md#read-the-evidence) explains where
each format stores it. All formats require the playback extra, and MP4 also
requires FFmpeg. See the
[export guide](../playback.md#export-a-still-or-video) for working-directory
instructions and the [media presentation reference](playback-presentation.md)
for Python settings and producer frame callbacks.

Inputs must pass schema and semantic validation before export. The command
rejects input/output aliases, colliding destinations, symlink output paths,
and existing files unless `--replace` is given. Invalid input, missing media
dependencies, and export failures produce concise errors on stderr and exit
nonzero.

Every requested format is rendered to temporary files before any destination
is published. A rendering failure leaves existing destination files untouched.
Publication then occurs atomically **per file**, not as one filesystem
transaction across all formats. If publication fails partway through, the
error lists the files already published. Successful commands print each
written path.
