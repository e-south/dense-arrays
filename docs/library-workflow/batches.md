---
title: Choose candidate batches
description: Prepare fixed selections or resample bounded batches using committed feedback.
---

# Choose candidate batches

A candidate batch limits which eligible parts one packing search can use. Prepare
it once, inspect its identities, and run the saved plan as often as needed in new
output directories. Each plan retains the full eligible collection and the
ordered offered subset. Its packing model contains only the offered parts.

## Select and save parts for search

Start with a resolved generation plan. Choose a batch size and an explicit seed:

```python
from pathlib import Path
import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

# Resolve the request and bind its input records before execution.
source = da.plan(
    planning.DesignSpec(
        parts=[
            parts.Part("a", "ACGTTGCAAGTCCTGA", group="A"),
            parts.Part("b", "GATCAGTACCTAGGTC", group="A"),
            parts.Part("c", "TTGACCGATAGCTACG", group="B"),
            parts.Part("d", "CAGTTCGATGACCTAG", group="B"),
        ],
        length=planning.Length(maximum=16),
        strands="single",
        target=planning.Target(3),
    )
)
# Save the resolved plan with its input bindings; keep the destination new.
source.write("source.plan.json")
# Prepare the declared parts or batch and save its identities for reuse.
prepared = da.prepare(
    source,
    sampling=planning.BatchSampling(size=2, strategy="group_balanced", seed=7),
    out="batch.plan.json",
)
assert prepared.preview["parts"] == 4
assert prepared.preview["offered_parts"] == 2
assert len(prepared.request.batch.part_ids) == 2
# Generate under the declared bounds into a new output directory.
result = da.run(prepared, out="batch-run")
# Read saved run state and attainment. Recount stored evidence before returning.
summary = da.inspect(result, verify=True)
assert summary.accepted == 2
assert summary.termination_reason == "batch_exhausted"
```

The CLI accepts the same saved source plan and sampling settings. Run these
commands in a directory containing `source.plan.json`; each output must be new:

```bash
# Prepare the declared pool or offered batch.
dense-arrays prepare source.plan.json --batch-size 2 \
  --batch-strategy group_balanced --batch-seed 7 --out batch.cli.plan.json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect batch.cli.plan.json --view plan --json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run batch.cli.plan.json --out batch-cli-run
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect batch-cli-run --view attempts --all --json
```

Preparation verifies source files before sampling. The result embeds the parts
and original input fingerprints, so moving the plan or removing the source table
does not prevent replay. Running, inspecting, exporting or resuming this plan
never samples again. Existing output files are never replaced.

## Choose a policy

| Setting | Meaning |
| --- | --- |
| `size` | Exact number of distinct supplied part identities offered to search. |
| `strategy="uniform"` | Select by seeded part priorities after including fixed occurrences. |
| `strategy="group_balanced"` | Repeatedly select from the least represented available group, counting fixed occurrences already selected. Every part must have a group. |
| `seed=0` | Default seed for batch selection, independent of the generation seed used for assembly. |

Default eligibility uses `part_priority_sha256.v1`. Priorities depend on the seed,
logical cell and part identity. Reordering input rows does not change the selected
identities. Matrix preparation creates one independent batch per active cell;
zero-target cells do not sample. Changing a cell's parts requires a batch bound
to that cell's updated collection.

Fixed occurrences always occupy slots in a sampled batch. A size below their
count, above the eligible count, or below one fails before writing output.
Sampling selects identities; equal DNA sequences remain distinct parts unless
uniqueness is requested below. Packing requirements, accepted-design diversity
and final-sequence checks retain their own contracts.

## Bound repetition within each batch

Combine hard membership policies when each batch should contain distinct
sequences or limit the contribution from each group:

```python
# Prepare the declared parts or batch and save its identities for reuse.
unique = da.prepare(
    source,
    sampling=planning.BatchSampling(
        size=2,
        seed=7,
        unique_sequences=True,
        max_per_group=1,
    ),
    out="unique.plan.json",
)
offered = [
    p for p in source.request.parts if p.part_id in unique.request.batch.part_ids
]
assert len({p.sequence for p in offered}) == 2
assert len({p.group for p in offered}) == 2
```

```bash
# Prepare the declared pool or offered batch.
dense-arrays prepare source.plan.json --batch-size 2 --batch-seed 7 \
  --unique-sequences --max-per-group 1 --out unique.cli.plan.json
```

| Setting | Membership rule |
| --- | --- |
| `unique_sequences=True` / `--unique-sequences` | At most one part with each exact supplied sequence across the batch. Reverse complements remain distinct. |
| `unique_cores=True` / `--unique-cores` | At most one part with each oriented core within each group. Every eligible part must have a group and complete core annotations. |
| `max_per_group=N` / `--max-per-group N` | At most N offered parts from any one group. Every eligible part must have a group. |

Core comparison uses the annotated half-open interval and `core_orientation`.
For example, an annotated `TT` interval in reverse orientation has core sequence
`AA`. The same core may occur in different groups. Missing annotations remain
unknown and fail a core-uniqueness request before publication.

These settings use `part_flow_priority_sha256.v1`: select a feasible full-size
batch, then prefer seeded part ranks. With `group_balanced`, first minimize the
sum of squared group counts among feasible selections; seeded ranks break
objective ties. Fixed parts contribute to all counts and remain mandatory.
Input-row order does not choose representatives. Tied choices may depend on the
installed backend version; saved batch membership is the replay authority.

All selected policies apply together. A fixed-part conflict or an impossible
requested size fails before writing output. Group caps never relax automatically.
Revise the explicit size or policy to change that search. These rules restrict
the offered batch; they do not prove that its packing constraints are feasible.

## Supply a known batch

An explicit ordered list can reproduce a selection obtained elsewhere:

```python
batch = planning.CandidateBatch(
    part_ids=("b", "c"),
    collection_id=source.collection_id,
)
# Resolve the request and bind its input records before execution.
explicit = da.plan(source.request.with_changes(batch=batch))
# Save the resolved plan with its input bindings; keep the destination new.
explicit.write("explicit.plan.json")
assert explicit.request.batch.sampling is None
```

The collection binding covers part identities, sequences and import provenance.
Unknown IDs, repeated IDs, mismatched collections and unsupported policy versions
fail before execution. A batch digest binds its ordered membership and selection
metadata. For a matrix, declare explicit batches through
`MatrixSpec.batches`, keyed by the complete cell identity; omitted cells offer
their full eligible collection. A matrix base cannot silently broadcast a batch.

## Interpret a shortfall

The example requests three designs from a two-part batch with room for one part
per design. It accepts two and stops with `batch_exhausted`. Another batch may
produce additional sequences. Likewise, `batch_infeasible` applies to the offered
subset and constraints, not to every possible selection from the eligible pool.

Every attempt records `evidence.batch_id`, including interrupted attempts. The
saved plan contains the membership; accepted placements retain their original
part and placement identities. [Resume](recovery.md) restores that same batch and
its committed packing exclusions.

To prepare a different selection, revise the source plan explicitly; `prepare`
rejects plans that already contain a batch or schedule. Use
[plan comparison](handoffs.md) to inspect changes.

## Search an ordered schedule

Prepare several selections and give each a finite attempt allowance:

```python
# Prepare the declared parts or batch and save its identities for reuse.
scheduled = da.prepare(
    source,
    sampling=planning.BatchSampling(size=2, strategy="group_balanced", seed=7),
    batch_count=3,
    attempts_per_batch=2,
    accepted_per_batch=1,
    out="schedule.plan.json",
)
assert scheduled.preview["batches"] == 3
assert scheduled.preview["max_offered_parts"] == 2
# Generate under the declared bounds into a new output directory.
scheduled_run = da.run(scheduled, out="scheduled-run")
# Read saved run state and attainment. Recount stored evidence before returning.
scheduled_summary = da.inspect(scheduled_run, verify=True)
assert scheduled_summary.counts["started"] <= 6
assert scheduled_summary.accepted <= 3
```

```bash
# Prepare the declared pool or offered batch.
dense-arrays prepare source.plan.json --batch-size 2 \
  --batch-strategy group_balanced --batch-seed 7 --batch-count 3 \
  --attempts-per-batch 2 --accepted-per-batch 1 --out schedule.cli.plan.json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run schedule.cli.plan.json --out scheduled-cli-run
```

Each cell searches batches in the saved order. It advances when the offered
model is infeasible, its enumeration is exhausted, or a local cap is reached.
`accepted_per_batch` optionally limits each batch's contribution to the cell.
Only accepted, distinct final designs consume that cap. Rejections, duplicates
and interrupted reservations still consume attempts. An accepted cap requires
`attempts_per_batch`, so a batch with no accepted designs remains bounded.
Omit `accepted_per_batch` to search until exhaustion or the attempt limit.
Backend errors and invalid results stop the affected cell. Unresolved solver
outcomes also stop by default. Set `on_unproven="next_batch"` on a
`BatchSchedule` or `Resampling` policy to advance after an `unknown` or
`unproven` exact-search outcome. The failed search consumes an attempt and keeps
its original evidence; no unproven design is accepted. Unknown causes remain
unknown, even with a configured solver time limit. A measured interruption
after that recorded boundary resumes at the next batch under the same limits.

Targets and sequence uniqueness apply across the entire cell. Global time and
attempt limits also remain in force across all batches and matrix cells. A
cell below target after its final batch stops with `batch_schedule_exhausted`.
This records the end of the declared search allowance, not global infeasibility.
The CLI returns exit status 3 for such an incomplete run.

`planning.BatchSchedule(batches=(...), attempts_per_batch=...)` accepts explicit
`CandidateBatch` selections. Assign it to `DesignSpec.schedule`, or to a cell's
entry in `MatrixSpec.batches`. A single-cell request chooses one of `batch`, `schedule` or `resampling`.

Attempt records include one-based `batch_index` and `batch_attempt` coordinates
alongside `batch_id`. Scheduled designs retain `batch_id`, including in portable
bundles, so their placements can be checked against the offered membership.
Resume recounts the committed batch positions and restores packing exclusions
only for batches that still have work. It does not restart consumed batches or
replenish any allowance. The accepted count resumes within its original cap.

## Resample during generation

Use `Resampling` when later batches should respond to accepted designs or failed
searches. The plan declares a finite number of batches and attempts per batch.
Selections are created during execution and saved before search begins.

This example continues from the four-part `source` plan above:

```python
# Resolve the request and bind its input records before execution.
adaptive = da.plan(
    source.request.with_changes(
        resampling=planning.Resampling(
            sampling=planning.BatchSampling(size=1, seed=7),
            max_batches=20,
            attempts_per_batch=1,
            accepted_per_batch=1,
            feedback=planning.FeedbackPolicy(coverage_alpha=1e6, coverage_power=20),
        ),
    )
)
# Save the resolved plan with its input bindings; keep the destination new.
adaptive.write("adaptive.plan.json")
# Generate under the declared bounds into a new output directory.
result = da.run(adaptive, out="adaptive-run")
assert da.inspect(result, verify=True).accepted == 3
with da.inspect(result, view="batches", all=True).records() as records:
    decisions = list(records)
assert len(decisions) == 3
assert sum(decisions[-1].batch.feedback.used.values()) == 2
```

Use the same saved policy through the CLI:

```bash
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run adaptive.plan.json --out adaptive-cli-run
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect adaptive-cli-run --view batches --limit 10 --json
# Write the declared selection or document to a new destination.
dense-arrays export adaptive-cli-run --view batches --all --out batches.json
```

Each decision records the cell, one-based batch index, preceding global attempt
count, exact membership, sampling policy and feedback counts. `inspect` and
`export` support paginated batch records in runs and portable bundles. Reader
limits include membership and feedback entries. Verification of a native run
recounts feedback from committed attempts; bundles verify saved membership and
policy but do not include the attempt history needed to recount feedback.

A matrix can share `base.resampling` or override it per cell through
`MatrixSpec.batches`. Each cell keeps its own observations and named random
stream. Zero-target cells make no selections. A cell that reaches `max_batches`
below target stops with `batch_limit`; global attempts and active-time limits
still apply. This result does not establish global infeasibility.

### Interpret feedback

`FeedbackPolicy` defaults to `coverage_alpha=1`, `coverage_power=1`,
`failure_alpha=0`, and `failure_power=1`. Omit `feedback` for sampling independent
of outcomes. Observations aggregate by group and exact supplied sequence, so
aliases share a weight:

```text
weight = (1 + coverage_alpha / (1 + used)^coverage_power)
         / (1 + failure_alpha * failed)^failure_power
```

Accepted placements increment `used`. A final-screening rejection or an initially
infeasible offered model increments `failed` for every offered part. This is a
selection heuristic; it does not identify which part caused rejection.
Duplicates, exhausted enumeration, interruptions, and unresolved or failed
backend calls do not increment either count. Backend failures and unknown solver
outcomes stop the affected cell instead of trying another batch.

Weights use versioned seeded priorities. Equal weights preserve unweighted
priority order. Group balancing remains the primary policy when requested;
weights order alternatives within that policy. Hard sequence/core uniqueness
and group caps still apply. Invalid or numerically unusable weights fail
explicitly.

A saved batch and its first attempt reservation commit together. Resume reuses
saved selections, recounts observations, and rebuilds exclusions for the current
batch. It draws a new batch only after the current batch ends. Limits and accepted
counts retain their original meaning after interruption.
