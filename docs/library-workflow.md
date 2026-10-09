---
title: Generate a saved library
description: Prepare pools, generate bounded designs, inspect outcomes, and export sequences with joinable placements.
author: Eric J. South
---

# Generate a saved library

Generate a library, inspect its sequences and binding-site placements, and save
results for later analysis. The same requests work in Python and the CLI.

Follow [installation](installation.md#use-the-library-workflow), then run these
examples from a new working directory. Each output path is created by the
command that uses it. The base installation includes CBC and request-file
support; figures use the optional `playback` extra.

## Generate and inspect one design

```bash
# Pack four 16-base motifs into a saved design of at most 40 bases.
dense-arrays run --motif ACGTTGCAAGTCCTGA --motif AGTCCTGATCGTACCG \
  --motif TCGTACCGATGCTTAG --motif ATGCTTAGGACGTTCA \
  --length 40 --count 1 --seed 7 --out runs/first
# Verify the saved sequence, placements and attempt counts.
dense-arrays inspect runs/first --verify
# Show the accepted design and its recorded placements.
dense-arrays inspect runs/first --view designs --limit 1 --json
```

`--length` sets a maximum. For a fixed final length and optional padding, use an
[assembly request](library-workflow/generation/assembly.md). The receipt shows
accepted and requested counts, the run directory, and why generation stopped.
Use a new destination for another run.

The equivalent Python operations return typed values:

```python
from pathlib import Path

import dense_arrays as da

# Describe the same part collection and length limit in Python.
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
with da.inspect(result, view="designs", limit=1).records() as records:
    design = next(records)  # Read one design; the context manager closes the reader.
    print(design.reference, design.realized.sequence)  # Keep identity beside DNA.
    for placement in design.realized.placements:  # Check every recorded interval.
        assert (
            design.realized.sequence[placement.start : placement.end]
            == placement.sequence
        )
```

The `with` block closes the record reader when finished, including when you stop
early. See [record inspection](library-workflow/results/inspection.md) for
pagination, work limits and recorded software versions.

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

Check the saved summary before choosing the next step:

| Result | Next step |
| --- | --- |
| `completed` / `target_attained` | Inspect, select or export the accepted designs. |
| `stopped` | Read the stopping reason and accepted count, then [inspect the diagnostics](library-workflow/results/quality.md). |
| `batch_infeasible` or `batch_exhausted` | Review the offered parts and constraints; [search outcomes](reference/optimizer.md#solver-outcomes) distinguish infeasibility from exhausted alternatives. |
| `failed` | Inspect the saved run for the execution error and any committed designs. |

The default request asks for one design, considers both strands, and permits
1,000 attempts, 300 active seconds and 30 seconds per solve. See
[resource limits](library-workflow/resources.md) to size larger requests. Limits
can stop a run below its target; time allowances are cooperative.

Accepted sequences are unique within each design combination. Exact search
maximizes placed occurrences within the offered packing model. Use the
[search guide](library-workflow/search.md) for greedy proposals, part-use
preferences and the interpretation of solver outcomes.

`inspect --verify` checks saved identities, coordinates, requirements and
accounting. Summary inspection reads the recorded totals. For an interrupted
run, the [recovery guide](library-workflow/recovery.md) explains when to resume
and when to start a linked request. The [CLI reference](reference/cli.md#exit-codes-and-machine-output)
explains exit codes and JSON errors.
