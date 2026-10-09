---
title: Build a curated library
description: Generate two design combinations from a reusable part pool, inspect their evidence and share the results.
author: Eric J. South
---

# Build a curated library

Prepare 750 demonstration sequences in three named groups, then request 50
100-base designs for each of two combinations. The baseline uses those sequences
alone. The fixed-pair combination also requires two named anchor parts with a
16–18-base spacer. Group labels are input annotations; the example makes no
claim about measured binding or expression.

Follow [library-workflow installation](../installation.md#use-the-library-workflow).
Save these three files together in a new working directory:

- [parts.csv](../examples/curated-library/parts.csv): 750 sequence records and two anchors, each with an explicit ID.
- [prepare.yaml](../examples/curated-library/prepare.yaml): curate the table into a reusable pool.
- [design.yaml](../examples/curated-library/design.yaml): declare combinations, targets, sampling, requirements and effort limits.

Run the following commands from that directory. Destinations must be new.

## Prepare, preview and generate

```bash
# Prepare the declared pool or offered batch.
dense-arrays prepare prepare.yaml --out pool --json
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan design.yaml --out design.plan.json --json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run design.plan.json --out run --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect run --verify --json
```

The preview lists `architecture=baseline` and `architecture=fixed_pair`, each
with a target of 50 accepted designs. The baseline selects the three sequence
groups from the pool; the fixed-pair combination includes the anchors as well.
Its `add_requirements` declaration adds the two fixed occurrences and their
spacing rule to the shared padding check.

Each offered batch contains ten parts, including any required anchors. Sampling
balances the declared groups and requires distinct full sequences. These inputs
have no annotated motif cores, so core uniqueness is disabled. Group balancing
applies to offered parts; it does not require every group in each accepted design.

The run shares a maximum of 2,000 attempts and 600 active seconds. Each solver
call has a cooperative ten-second allowance. Each combination can sample up to
200 batches, with at most 30 attempts or ten accepted designs per batch. A time
limit can stop the run before its targets are met. `on_unproven: next_batch`
records unresolved exact-search outcomes and moves to another batch; it does
not accept unproven designs. Backend errors and invalid results still stop.
CLI exit status `3` means incomplete output; inspect its recorded reason before
deciding on [recovery](recovery.md) or a new request.

Every accepted sequence is exactly 100 bases. Added bases go on the left and
must have a GC fraction from 0.4 to 0.6, inclusive. This is a **padding** check,
not a GC bound on the complete sequence. Empty padding is not applicable. A
one-base pad cannot meet this range; its candidate is rejected after the declared
trial limit. The range is never relaxed. Final coordinates include the padding
shift, so the spacer remains the downstream start minus the upstream end.

## Inspect and share

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect run --view quality --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect run --view diagnostics --json
# Write the declared selection or document to a new destination.
dense-arrays export run --view sequences --all --format fasta --out library.fasta
# Write the declared selection or document to a new destination.
dense-arrays export run --all --format bundle --out library
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect library --verify --json
```

Quality reports retain separate targets and attainment for both combinations.
The portable bundle retains sequences, placements, cell identities and the
evidence needed to check their requirements after the original run and pool
are unavailable. FASTA carries sequence records; use the bundle for complete
design evidence. Equal sequences from different cells retain distinct design
references.

With the `playback` extra installed, render the same quality report:

```bash
# Render a figure from the selected saved records.
dense-arrays render run --view library-quality --out quality.png
```

## Use the Python interface

Start in another fresh directory containing the same three inputs. The public
request reader validates settings without opening their part sources. Planning
then binds the input records; generation uses that saved plan.

```python
import dense_arrays as da

preparation = da.inspect("prepare.yaml", view="request").request
# Prepare the declared parts or batch and save its identities for reuse.
pool = da.prepare(preparation, out="pool")
request = da.inspect("design.yaml", view="request").request
# Resolve the request and bind its input records before execution.
plan = da.plan(request)
# Save the resolved plan with its input bindings; keep the destination new.
plan.write("design.plan.json")
# Generate under the declared bounds into a new output directory.
run = da.run(plan, out="run")
# Read saved run state and attainment. Recount stored evidence before returning.
summary = da.inspect(run, verify=True)
for cell in summary.cells.values():
    print(cell.cell_id, cell.accepted, cell.target, cell.termination_reason)
# Publish the selected final sequences to a new destination.
da.export(run, view="sequences", all=True, format="fasta", out="library.fasta")
# Publish the declared records to a new destination.
da.export(run, all=True, format="bundle", out="library")
```

To change the requested number, edit `allocation.per_cell` in YAML or use
`request.with_changes(allocation=planning.Allocation(per_cell=10))` after
`from dense_arrays import planning`. This creates a new request; it does not
change an existing run. Use the same [six operations](../library-workflow.md)
when changing sources, constraints or targets.
