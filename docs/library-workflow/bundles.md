---
title: Share a portable library
description: Package selected designs and their verification evidence in a movable directory.
---

# Share a portable library

Export a bundle to carry selected sequences, placements, design rules and source
context in one directory. Move or copy the entire directory, then inspect its
records and verify their geometry and requirements using the contained evidence.

From the [saved-library example](../library-workflow.md), run:

```bash
# Write the declared selection or document to a new destination.
dense-arrays export runs/first --view designs --all \
  --format bundle --out handoff/library
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect handoff/library --verify
# Write the declared selection or document to a new destination.
dense-arrays export handoff/library --view sequences --all \
  --format fasta --out handoff/library.fasta
```

The bundle destination must be a new directory. `--all` declares the complete
selection; add `--design-id`, `--cell`, `--part-id`, `--group` or a declared filter
file to choose its members. A displayed inspection page cannot become a bundle
implicitly. Native runs and existing bundles can contribute to the same bundle.

## Export a selected collection

Run this Python example in a new directory after [installing the library workflow](../installation.md#use-the-library-workflow):

```python
from pathlib import Path

import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning, reporting

# Generate under the declared bounds into a new output directory.
run = da.run(
    planning.DesignSpec(
        parts=[
            parts.Part("a", "ACGTTGCAAGTCCTGA", group="A"),
            parts.Part("b", "GATCAGTACCTAGGTC", group="B"),
        ],
        length=planning.Length(maximum=32),
        strands="single",
        target=planning.Target(count=2),
    ),
    out="runs/source",
)
with da.inspect(run, view="designs", limit=1).records() as records:
    chosen = next(records).reference
# Choose designs by recorded identities and final-sequence properties.
selected = reporting.DesignFilter(design_ids=(chosen,))
# Publish the declared records to a new destination.
receipt = da.export(
    run,
    select=selected,
    all=True,
    format="bundle",
    out="handoff/library",
)
assert receipt.records == 1

Path("handoff/library").rename("handoff/moved-library")
# Read saved run state and attainment. Recount stored evidence before returning.
summary = da.inspect("handoff/moved-library", verify=True)
assert summary.designs == 1
assert summary.verified
assert summary.to_dict()["scope"] == "selected_collection"
assert summary.to_dict()["source_runs"][0]["accepted"] == 2
```

The selected collection contains one design. Its original source summary still
reports two accepted designs. Bundle counts describe included records; source
summaries preserve the original targets, effort and completion state.

Full design references survive the handoff. Repeated references with identical
records count once when combining artifacts. Equal DNA with different design
references retains each construction history. Conflicting records for the same
reference fail export.

## Read, filter and export bundled records

`inspect` supports `summary`, `designs`, `sequences`, `placements`, `selection`,
`quality`, `plans` and `plan` for a bundle.
Record views share the normal filters, work limits, closeable iterators and
continuation cursors. Continue a placement page inside a design without dropping
or repeating rows.

```python
bundle = "handoff/moved-library"
with da.inspect(bundle, view="placements", limit=1).records() as rows:
    first = list(rows)
    cursor = rows.next_cursor
with da.inspect(bundle, view="placements", all=True, after=cursor).records() as rows:
    remaining = list(rows)
assert len(first + remaining) == 2

# Publish the selected final sequences to a new destination.
da.export(
    bundle, view="sequences", all=True, format="fasta", out="handoff/sequences.fasta"
)
# Publish the selected coordinate annotations to a new destination.
da.export(
    bundle, view="placements", all=True, format="tsv", out="handoff/placements.tsv"
)
# Publish the declared records to a new destination.
da.export(bundle, view="summary", out="handoff/summary.json")
```

A bundle can also be filtered into a new bundle with `export(..., format="bundle")`.
The new collection retains its source context and resolved plan evidence. An
empty selection is valid and reports zero included designs.

## Combine artifacts and assess included designs

Pass native runs and bundles in the desired source order. Full design identities
are deduplicated before filters and projections; a repeated source never adds
another copy. Filtering uses the original run/cell and part collection identities.
Bare labels that match multiple identities fail with a request for a full reference.

```python
# Read saved accepted designs with placements.
combined = da.inspect([bundle, run], view="designs", all=True)
with combined.records() as records:
    assert len(list(records)) == 2

quality = da.inspect(bundle, view="quality").to_dict()
assert quality["selection"]["designs"] == 1
assert quality["source_runs"][0]["included_designs"] == 1
assert quality["source_runs"][0]["attainment"]["accepted"] == 2
assert quality["search"]["availability"] == "not_included"
assert quality["search"]["attempt_counts"] is None

combined_quality = da.inspect([bundle, run], view="quality").to_dict()
assert combined_quality["selection"]["designs"] == 2
assert combined_quality["search"]["availability"] == "complete"
```

Quality reports measure the included design population, with optional filters or
a [saved panel](selection.md). Original target and acceptance counts remain in
`source_runs[].attainment`; `included_designs` and `selected_designs` describe
the union before and after selection. Supplying a native run provides its checked
attempt history. A bundle alone supplies no attempt records: absent histories
are `null`, not zero. `search.availability` is `complete`, `partial` or
`not_included`, and `unavailable_source_refs` identifies missing histories.
Quality requires one consistent original revision per run across all inputs.

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect handoff/moved-library --view quality --json
# Render a figure from the selected saved records.
dense-arrays render handoff/moved-library --view library-quality \
  --out handoff/quality.png
```

Rendering uses the optional playback extra and preserves these evidence limits
in the figure and embedded report. Complete record scans check contained counts;
use `inspect --verify` for the full geometry, requirement and checksum checks.

## Inspect the included design rules

Use `plans` to list the included plan evidence. Plans shared by several source
runs appear once, by semantic ID. The view uses the same closeable readers,
pagination and work limits as design records. The human CLI lists short plan
summaries; `--json` includes the full evidence records.

```python
with da.inspect(bundle, view="plans", all=True).records() as records:
    included_plans = list(records)
assert len(included_plans) == 1
plan_filter = reporting.PlanFilter(plan_ids=(included_plans[0].plan_id,))
# Read saved resolved rules and input bindings.
evidence = da.inspect(bundle, view="plan", select=plan_filter)
assert evidence == included_plans[0]
# Publish resolved plan evidence to a new destination.
da.export(
    bundle,
    view="plan",
    select=plan_filter,
    out="handoff/design-rules.json",
)
assert da.inspect("handoff/design-rules.json", view="plan") == evidence
```

`plan` returns one `PlanEvidence`. If a bundle contains several plans, provide a
full plan ID; ambiguous requests fail. An unknown ID also fails, even when a
different requested ID is valid. Empty filters include the whole inventory.

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect handoff/moved-library --view plans --all --json
# For this single-plan bundle, the ID is unambiguous.
dense-arrays export handoff/moved-library --view plan \
  --out handoff/cli-design-rules.json
```

For multiple plans, add `--plan-id PLAN_ID` to `inspect` or `export --view plan`.
`--plan-id` is repeatable with `--view plans`; `PlanFilter` provides the matching
Python predicate and selection-file schema. `export --view plans --all` writes a
record collection; `export --view plan` writes one native evidence document.

Evidence records contain normalized parts, rules, policy versions, fingerprints
and any frozen parent exclusions. They support
[semantic comparison](extension.md#compare-design-rules) and verification.
Execution requires a design request or an executable generation plan with its
input bindings. The evidence document is not accepted by `run`.

## What travels with the directory

`bundle.json` identifies the collection, selection, origin snapshots, plan IDs,
metric policy and database checksum. `bundle.sqlite3` holds the selected design
records and their resolved plan evidence. Placements and requirement results are
part of each design record; part annotations and design rules come from its bound
plan. Keep both files together.

Resolved plan evidence contains normalized parts, rules, policy versions, import
provenance and input fingerprints. It preserves the generation plan's semantic
identity while omitting executable source locations. Prepared-pool identity and
annotations remain bound to the resolved parts. Parent exclusions travel as
sequence identities and design references.

Origin run manifests also retain their recorded producer and solver versions.
Summary and quality inspection expose those values after the original runs are
removed. Missing provenance in older source artifacts stays unavailable.

## Interpret verification

`inspect --verify` checks the database checksum, plan identities, record joins,
sequence uniqueness, placement geometry, parent exclusions and independently
recounted final requirements. Its receipt names the included evidence boundary
and reports checked record bytes separately from physical database bytes.

Original input tables, complete attempt histories and ancestor design records
are outside that boundary. A bundle supports verification of its included
designs. It does not reproduce the original search or verify source-wide
accounting. Use a [resolved executable plan](handoffs.md) and its bound input
files when running generation again. Attempt inspection uses a run’s recorded history. Individual-design rendering
also accepts bundles; select exactly one design.

Export checks selected design evidence before committing the bundle. Concurrent
writers cannot replace a destination owned by another export. An interrupted
write retains an incomplete marker and fails inspection explicitly; choose a
new output directory for a retry. Work limits apply to source metadata, plan
evidence and record scans, including the identity state retained for receipts.
These are explicit counters, not hard memory or time limits.
