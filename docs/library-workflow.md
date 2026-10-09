---
title: Generate a saved library
description: Prepare pools, generate bounded designs, inspect outcomes, and export sequences with joinable placements.
author: Eric J. South
---

# Generate a saved library

Generate designs into an explicit directory, then inspect their sequences,
placements and search outcomes without rerunning the solver. See
[resource limits](library-workflow/resources.md) to bound candidate and model sizes.
Start with explicit
parts or [prepare a sampled pool](library-workflow/preparation.md) from a motif
or a declared background distribution.

Follow [installation](installation.md#use-the-library-workflow) for the library
workflow. Run each example in a new working directory; output destinations must not already exist. The base installation
includes CBC and YAML/JSON request parsing. Rendering uses the optional
`playback` extra.

For a complete table-based example with two combinations, bounded resampling
and portable results, use [Build a curated library](library-workflow/curated-example.md).

## Generate and inspect one design

```bash
# Generate into a new directory with explicit bounds.
dense-arrays run --motif ACGTTGCAAGTCCTGA --motif AGTCCTGATCGTACCG \
  --motif TCGTACCGATGCTTAG --motif ATGCTTAGGACGTTCA \
  --length 40 --count 1 --seed 7 --out runs/first
# Read or verify saved evidence without generating again.
dense-arrays inspect runs/first --verify
# Read or verify saved evidence without generating again.
dense-arrays inspect runs/first --view designs --limit 1 --json
```

`--length` is a maximum. Use an explicit assembly request below for exact final
length. A successful receipt reports accepted/target designs and
the destination. An incomplete search preserves its accepted prefix and reports
why it stopped. Repeating the command does not overwrite an existing run.
Use an explicit [packing preference](library-workflow/search.md) to favor underused
parts while preserving occurrence count as the primary objective.

The equivalent Python operations return typed values:

```python
from pathlib import Path

import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

# Declare inputs and bounds before running the solver.
request = planning.DesignSpec(
    parts=(
        parts.Part("row:1", "ACGTTGCAAGTCCTGA"),
        parts.Part("row:2", "AGTCCTGATCGTACCG"),
        parts.Part("row:3", "TCGTACCGATGCTTAG"),
        parts.Part("row:4", "ATGCTTAGGACGTTCA"),
    ),
    length=planning.Length(maximum=40),  # Permit shorter accepted sequences.
    target=planning.Target(count=1),  # Request one unique final sequence.
    seed=7,  # Fix seeded sampling and padding streams.
)
preview = da.plan(request)  # Validate and bind inputs without solving.
print(dict(preview.preview))  # Read resolved counts, constraints and effort limits.
result = da.run(preview, out=Path("runs/python-first"))  # Create a new run directory.
summary = da.inspect(result, verify=True)  # Recount and verify saved records.
assert summary.accepted == summary.target == 1  # Confirm the target was reached.
assert summary.producer.solver.name == "CBC"  # Read the solver used during execution.
print(summary.producer.to_dict())  # Show recorded software and platform versions.
with da.inspect(result, view="designs", limit=1).records() as records:
    design = next(
        records
    )  # Take the one requested design; context exit closes the reader.
    print(design.reference, design.realized.sequence)  # Keep identity beside DNA.
    for placement in design.realized.placements:  # Check every recorded interval.
        assert (
            design.realized.sequence[placement.start : placement.end]
            == placement.sequence
        )
```

A record view holds no open file. Each `records()` call creates an independent
iterator at the same committed revision. Use its context manager when stopping
early; exhaustion also closes its reader.

Run and pool summaries record the Dense Arrays, Python and OR-Tools versions,
plus operating-system family and machine architecture. A run records its solver
name and reported version after building the model. Preparation and failures
before model creation have no solver identity. Inspection shows these recorded
values without querying the current solver; `--json` includes the complete
producer record. Versions help diagnose differences between executions; they
do not establish deterministic replay or identify unpublished source edits.

## Choose the next task

| Task | Guide |
| --- | --- |
| Import sites and prepare a reusable pool | [Curated parts](library-workflow/preparation/curated.md) |
| Add fixed positions, final length and padding | [Assembly](library-workflow/generation/assembly.md) |
| Select sequences and export annotated placements | [Record exports](library-workflow/results/export.md) |
| Page through records with explicit work limits | [Inspection](library-workflow/results/inspection.md) |
| Explain shortfalls and plot composition | [Diagnostics and quality](library-workflow/results/quality.md) |
| Combine design requirements or sample part batches | [Matrices](library-workflow/matrices.md) and [batches](library-workflow/batches.md) |
| Add designs or continue interrupted execution | [Extension](library-workflow/extension.md) and [recovery](library-workflow/recovery.md) |

## Interpret completion and failures

The default search enumerates exact packing paths with proven optimality for
each offered model. [Greedy search](library-workflow/search.md#generate-a-greedy-proposal)
provides one unproven packing per offered batch. Neither method enumerates
arbitrary gaps or every DNA sequence.
Different paths can produce the same final DNA; within a cell, the run accepts
only the first occurrence of an exact final sequence.

Defaults are one requested design, double-strand eligibility, seed zero,
1,000 attempts, 300 accumulated active seconds and 30 seconds per solve. These limits bound effort; they do not guarantee completion. Time limits are cooperative; model construction and
cleanup are not hard wall-clock deadlines. Exact packing enumeration does not
draw randomly; recording a seed does not promise solver tie order. Padding uses
a versioned SHAKE-256 stream bound to seed, cell, batch, attempt and trial. The
assembly record retains the coordinate transform, stream identity and trial.
Changing a policy deliberately changes its version; historical runs are not
reinterpreted as the new policy.

A solver attempt can make several padding proposals. Their count is separate
from solver attempts. A failed final screen records `screening_rejection`; a
bounded padding search records `padding_trials_exhausted`. Neither establishes
that all possible assembled sequences are infeasible. Accepted sequences are
unique after assembly; the policy excludes a packing path after its
accepted, duplicate or rejected candidate and does not enumerate every padding.

| Observation | Meaning |
| --- | --- |
| `completed` / `target_attained` | The original accepted-design target was reached. |
| `stopped` / `attempt_limit` or `active_time_limit` | Effort ended with a visible shortfall. |
| `batch_infeasible` | CBC proved no path for the offered initial model. |
| `batch_exhausted` | No further path remains after this batch's exclusions. |
| `batch_schedule_exhausted` | The declared schedule ended below target; this does not prove global infeasibility. |
| `solver_unproven` or `solver_unknown` | No qualifying optimal result; an unknown cause stays unknown. |
| `failed` | Execution or backend failure; inspect any committed prefix. |

`inspect --verify` checks record checksums, coordinates, part identities,
requirements, design/attempt joins and count reconciliation. Plain summary
inspection reads bounded metadata instead of scanning the library. Record
views default to 100 rows. A clean interruption can
[resume its original request](library-workflow/recovery.md) when measured time,
remaining budgets and saved packing evidence permit continuation. Abrupt exits
with unknown active time require a new linked request.

Verification accepts the same read caps and exposes `verification_cost` before
execution. Its returned `verification` names the checked plan/attempt/design or
preparation/part boundary, record count and checked UTF-8 JSON bytes. This does
not count physical SQLite/index bytes or reproduce generation. In Python,
`artifacts.RunHandle(path, run_id, revision=N)` pins reads, exports and rendering
to a committed prefix while a writer continues.

Machine output uses `--json`; receipts and diagnostics otherwise use stderr,
while plan/inspection reports use stdout. Workflow exit codes are `0` for a
successful operation, `2` for invalid input or a rejected selection shortfall,
`3` for a generation shortfall or an explicitly allowed partial selection export,
`4` for execution failure and `130` for interruption. Inspecting a valid stopped
run succeeds; its status still says stopped. Integrity and read-limit failures
also exit `4`. Before streaming begins, `--json` domain failures emit a
`dense_arrays.error.v1` envelope with code, message, exit code and affected
artifact when known. If generation raises after saving run state, the CLI reports
`execution_error` with that run's path; Python raises `RunExecutionError` with
an inspectable `run` handle and the original exception as its cause.
After streaming begins, diagnostics remain on stderr and
the incomplete data prefix is not followed by an unrelated error object.

For low-level transient solving, use [Optimizer](reference/optimizer.md).
For parts from motif models, use [sampled preparation](library-workflow/preparation.md).
For figures and reusable reports, choose a [result view](reference/outputs.md).
