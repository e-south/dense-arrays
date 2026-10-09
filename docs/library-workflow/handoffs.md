---
title: Save requests, plans and reports
description: Export editable inputs and scoped JSON reports through the same Python and CLI operations.
---

# Save requests, plans and reports

Use `export` to save an editable request, preserve a resolved plan, or share a
report with its original scope. JSON documents retain their native schema.
Every destination must be new.

Start in an installed Dense Arrays checkout with the
[saved-library workflow](../library-workflow.md). These commands use its
`runs/first` example:

```bash
# Write the declared selection or document to a new destination.
dense-arrays export runs/first --view request \
  --out handoff/request.json
# Write the declared selection or document to a new destination.
dense-arrays export runs/first --view plan \
  --out handoff/plan.json
# Write the declared selection or document to a new destination.
dense-arrays export runs/first --view quality \
  --limit 20 --out handoff/quality.json
# Write the declared selection or document to a new destination.
dense-arrays export runs/first --view diagnostics \
  --out handoff/diagnostics.json
```

Omit `--all` for documents. It is required only for
[record exports](results/export.md#select-and-export-sequences) without a
bounded [selection](selection.md). A quality or
diagnostic report keeps its declared table limit; exporting it does not silently
expand a page into every row. Quality aggregates still cover the entire selected
population.

## Choose the document

| View | Input | Result |
| --- | --- | --- |
| `request` | Request file or typed request | Editable inputs; no input table read or solver call |
| `request` | Native run or generation plan | Editable rules, bound parts and any explicit exclusions; native runs add parent-revision lineage |
| `request` | Portable bundle | One included plan's editable rules; choose `--plan-id` when more than one plan is included |
| `request` | Pool or preparation plan | Preparation rules with their table locator |
| `plan` | Run, pool or resolved plan | Frozen rules, input fingerprints and preparation/generation evidence |
| `summary` | One run or pool | Its original manifest scope |
| `quality` | One or more runs | Selected composition, original source attainment and search effort |
| `diagnostics` | One run | Bounded explanations from its attempt history |
| `plan` with `--compare` | Two generation plans or runs | Semantic changes and unchanged categories |
| `quality` with `--compare` | Runs, bundles, bound reports or saved quality JSON | Aggregate metric differences with each population and denominator |

Reports use JSON. Sequence and placement records additionally support CSV/TSV,
and sequences support FASTA. Unsupported combinations fail before source scans.

## Save and revise design inputs

The following Python example runs independently in a new working directory:

```python
from pathlib import Path

import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

table = Path("parts.csv")
# Write the example input or request so it can also be used from the CLI.
table.write_text("part_id,sequence,group\na,ACGTTGCAAGTCCTGA,A\nb,GATCAGTACCTAGGTC,B\n")
# Declare the part collection, sequence bounds and generation policy.
request = planning.DesignSpec(
    parts=parts.PartTable(table, "csv"),
    length=planning.Length(maximum=32),
    strands="single",
)
# Publish the declared records to a new destination.
da.export(request, out="handoff/source-request.json")
# Generate under the declared bounds into a new output directory.
run = da.run(request, out="runs/original")

effective = da.inspect(run, view="request").request
assert isinstance(effective.parts, parts.BoundParts)
revised = effective.with_changes(target=planning.Target(count=2), seed=19)
# Publish the declared records to a new destination.
da.export(revised, out="handoff/revised-request.json")
# Read saved resolved rules and input bindings.
original_plan = da.inspect(run, view="plan")
assert da.plan(revised).collection_id == original_plan.collection_id
# Save the resolved plan with its input bindings; keep the destination new.
original_plan.write("handoff/original-plan.json")
# Read saved resolved rules and input bindings.
comparison = da.inspect(run, view="plan", compare=da.plan(revised))
# Publish the declared records to a new destination.
da.export(comparison, out="handoff/changes.json")
```

Run the saved request with
`dense-arrays run handoff/revised-request.json --out runs/revised`.
Use `inspect runs/original --view request --json` to display the same effective
request. A request recovered from a run retains its resolved parts, import
transformations and input fingerprints in `parts.BoundParts`. Its collection
identity stays the same when you change design rules, target or seed.
`RequestReport.request` is immutable, and
`with_changes(...)` validates field types without reading files or running a
solver. Planning validates the combined requirements against the supplied parts.
`lineage.parent` records the source run, plan and committed revision; it does not
exclude that run's sequences.

For a native table or pool input, recovered requests also retain the source
locations. Both `plan` and `run` check those files against the saved fingerprints;
missing or changed files fail before a run destination is created. Inspection
and request export do not reopen them. Exported locations are relative to the
new request file, so moving a folder containing both request and sources works.

A request recovered from a portable bundle retains the same part identities,
fingerprints and import history, with no source-file locations. Planning uses
its embedded parts; it does not claim to verify unavailable original files.
The plan preview declares `input_verification: embedded_parts`. Frozen extensions
likewise retain their verified part evidence without requiring the original
input files. To replace the parts deliberately, use
`effective.with_changes(parts=(parts.Part("new", "GCTACGTTCAGATCGA"),))`; the replacement
starts its own collection. Editing embedded records in exported JSON without
updating their evidence fails an identity check.
To change parts in JSON, replace the `parts` object with an array of new part
records instead of editing the saved snapshot.

For direct comparison export, both interfaces accept the same two inputs:
`da.export(run, view="plan", compare=da.plan(revised), out="handoff/changes-2.json")`
and `dense-arrays export runs/original --view plan --compare PLAN.json --out changes.json`.
`PLAN.json` must be a resolved generation plan, not an unresolved request.

An extension request can be exported before its parent is available. Exporting
a resolved extension as an editable request retains its parent and ancestor
exclusions as an explicit `exclude` policy. Changing a rule preserves that policy;
removing it requires an explicit edit such as `with_changes(exclude=None)`.
Saved plans and bundle evidence preserve their declared lineage; only native
run inspection adds a new parent reference. A shared plan in a bundle does not
imply one particular origin run.

To exclude an accepted library when revising an ordinary request, continue the
example with an explicit source, cell mapping and uniqueness policy:

```python
without_reuse = revised.with_changes(
    exclude=planning.LibraryExclusion(
        source=planning.ParentRun(run.path),
        cell_mapping={"default": "default"},
        uniqueness="exact_sequence_per_cell.v1",
    )
)
# Resolve the request and bind its input records before execution.
revised_plan = da.plan(without_reuse)
assert revised_plan.preview["excluded_sequences"] == 1
# Save the resolved plan with its input bindings; keep the destination new.
revised_plan.write("handoff/revised.plan.json")
# Publish editable settings to a new destination.
da.export(revised_plan, view="request", out="handoff/revised-with-exclusions.json")
```

Planning verifies a terminal source with no active writer, then embeds its
accepted-sequence exclusions and inherited exclusions. The resolved request can
be exported, edited and planned again after the source run is unavailable.
Only the explicit `default` → `default` cell mapping and exact final-sequence
identity are currently supported. A changed matrix mapping or uniqueness policy
fails before execution. `plan(..., read_limits=...)` bounds source-evidence reads.
The CLI plan preview displays the excluded count, mapping and uniqueness policy.

A `PoolSource` containing a path can be saved as an editable request. If it
contains a `PoolHandle`, export the resolved plan to retain the checked pool
identity instead of reducing it to an unchecked path.

## Keep report scope visible

Continue the Python example:

```python
# Use the typed requests and operations needed by this example.
from dense_arrays import reporting

# Read saved composition and search metrics.
report = da.inspect(
    run,
    view="quality",
    select=reporting.DesignFilter(groups=("A",)),
    limit=2,
)
# Publish the declared records to a new destination.
receipt = da.export(report, out="handoff/selected-quality.json")
assert receipt.view == "quality"
assert receipt.records == 1
# Publish the declared records to a new destination.
da.export(run, view="diagnostics", out="handoff/diagnostics.json")
```

The matching CLI query is
`dense-arrays export runs/original --view quality --group A --limit 2 --out selected-quality.json`.
To retain an inspection cursor or custom diagnostic limit in Python, export the
typed report returned by `inspect`. A bound report cannot be given a different
filter or read limit at export time. Source attainment and search effort retain
their original population even when composition uses a design filter.

`--out -` emits only document JSON on stdout. Cost diagnostics and receipts go
to stderr, including with `--json`. For a file destination, `--json` emits the
receipt on stdout. Python returns an `ExportReceipt` separately; writable text
streams remain open. Receipts include the document schema, checksum and source
identities. They are not inserted into editable input JSON.

## Compare library quality

Use `view="quality"` and `compare=` to measure differences between two libraries.
Continue the example with a second run:

```python
# Generate under the declared bounds into a new output directory.
second = da.run(revised, out="runs/revised")
# Read saved composition and search metrics.
difference = da.inspect(run, view="quality", compare=second)
assert isinstance(difference, reporting.QualityComparison)
counts = next(
    item for item in difference.metrics if item.path == ("selection", "designs")
)
assert (counts.before, counts.after, counts.delta) == (1, 2, 1)
# Publish the declared records to a new destination.
da.export(difference, out="handoff/quality-comparison.json")
```

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/original --view quality --compare runs/revised --json
# Write the declared selection or document to a new destination.
dense-arrays export runs/original --view quality --compare runs/revised \
  --out handoff/quality-comparison-cli.json
```

Each `MetricDifference` reports before/after values, their denominators, and an
after-minus-before `delta`. Comparisons cover selected counts, eligible/unused
supply, composition minima/maxima/means, concentration, and search counts/time.
Both sides retain their source attainment and selection scope. They can be
single runs, combined sources, portable bundles, or separately filtered reports.
For different filters, pass two `QualityReport` values from `inspect`.

Composition uses selected designs. Search metrics describe available source-wide
attempt history; filtering designs does not filter the search denominator.
Missing or partial search history produces an explicit reason and no numeric
delta. Empty populations retain null metrics. Different or unsupported metric
policy versions produce `incomparable` records with no delta.

`scope_changed` describes a change in the declared selection.
`population_changed` is null when differently scoped, equally sized populations
cannot be distinguished from aggregate evidence. A comparison does not infer
causality or statistical significance.

Report usage-table limits do not limit aggregate comparisons. Do not pass `limit`,
`after` or `all` to a comparison. Its shared `ReadLimits` cover both scans and
retained state while preserving stricter limits attached to either report. To
compare sampled panels, first create saved selections or bind their quality
reports; comparison does not draw a new panel.

## Reuse saved metrics

Quality JSON retains the reported aggregates and usage page. Read it after
relocating or removing the source artifacts:

```python
saved = reporting.QualitySnapshot.from_report(da.inspect(second, view="quality"))
# Publish the declared records to a new destination.
da.export(saved, out="handoff/revised-quality.json")
# Read saved composition and search metrics.
reopened = da.inspect("handoff/revised-quality.json", view="quality")
assert reopened.to_dict() == saved.to_dict()
# Read saved composition and search metrics.
replayed_difference = da.inspect(saved, view="quality", compare=reopened)
assert all(
    item.delta == 0
    for item in replayed_difference.metrics
    if item.status == "comparable"
)
# Draw the selected saved evidence without generating new sequences.
da.render(reopened, view="library-quality", out="handoff/revised-quality.png")
```

`QualitySnapshot` checks histogram populations, source attainment and available
attempt totals. Completed sources must have reached their targets; aggregate
search counts and active time must agree with the included source histories.
It records metrics without rereading the original design or attempt records.
Comparisons label each side as `saved_report` or
`artifact_records`. Read limits cap the retained document structure separately
from historical work counters stored in the report. Rendering uses the saved
metrics and requires a supported metric policy and the playback extra.

## Reuse preparation settings

```python
# Prepare the declared parts or batch and save its identities for reuse.
pool = da.prepare(
    parts.PreparationSpec(parts.PartTable(table, "csv")),
    out="pools/curated",
)
# Read saved resolved rules and input bindings.
saved_preparation = da.inspect(pool, view="plan")
assert isinstance(saved_preparation, planning.PreparationPlan)
# Publish editable settings to a new destination.
da.export(pool, view="request", out="handoff/preparation-request.json")
# Publish the declared records to a new destination.
da.export(saved_preparation, out="handoff/preparation-plan.json")
preparation_request = da.inspect(pool, view="request").request
restricted = preparation_request.with_changes(
    retain=parts.Retention(select=parts.PartFilter(groups=("A",)))
)
# Save the resolved plan with its input bindings; keep the destination new.
da.plan(restricted).write("handoff/restricted-preparation.plan.json")
```

Use `dense-arrays prepare handoff/preparation-plan.json --out pools/replayed`
to check the bound input bytes and prepare a new pool. Inspecting or exporting
the stored preparation plan does not require the original table to be present.
Executing its preparation does.

## Understand file locations

Exported input locators are relative to the destination document. Stream exports
use the current working directory as their base. Relative paths preserve plan
identity; changing bound source bytes does not. Unknown metadata strings are
left unchanged.

Both resolved plan types support `write(path)`. It uses the same create-only
publication and relative locators as `export` and the CLI's `plan --out`.

A JSON plan is not a self-contained library bundle: it preserves evidence and
locators but does not copy external input files. Ordinary generation plans and
preparation plans recheck those inputs before execution. Frozen extension plans
use their embedded resolved evidence. For selected records with contained
verification evidence, export a [portable library bundle](bundles.md).

Exports honor `ReadLimits` before publishing. Report scans must complete before
the output file is created, and plan exports check their embedded identity state.
These are work counters, not hard byte or wall-clock limits. See
[reader limits](results/inspection.md#bound-record-inspection) for the
common inspection contract.
