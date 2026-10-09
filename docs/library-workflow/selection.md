---
title: Select and save a library panel
description: Allocate designs by total count or cell quota, then reuse the same identities for reports and exports.
---

# Select and save a library panel

Use a `LibrarySelection` to choose a fixed number of eligible designs. Save its
`SelectionSnapshot` to reuse the exact ordered identities, source revisions and
record content in inspection, export and quality figures.

`first` keeps source order. `random` requires an explicit seed and samples
without replacement. A total count and per-cell quotas are separate choices;
neither policy redistributes an unfilled quota.

## Choose designs in Python

Run this example in a new working directory with Dense Arrays installed. It
creates two small native libraries and selects one design from each:

```python
import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning, reporting

# Declare the part collection, sequence bounds and generation policy.
request = planning.DesignSpec(
    parts=[
        parts.Part("a", "ACGTTGCAAGTCCTGA", group="A"),
        parts.Part("b", "GATCAGTACCTAGGTC", group="B"),
    ],
    length=planning.Length(maximum=32),
    strands="single",
    target=planning.Target(count=2),
)
# Generate under the declared bounds into a new output directory.
left = da.run(request, out="runs/left")
# Generate under the declared bounds into a new output directory.
right = da.run(request, out="runs/right")
sources = [left, right]
# Declare the panel size and reproducible membership policy.
selection = reporting.LibrarySelection(
    filter=reporting.DesignFilter(groups=("A",)),
    take=reporting.Take(
        per_cell={f"{left.run_id}/default": 1, f"{right.run_id}/default": 1},
        policy="random",
        seed=23,
    ),
)
# Read saved run state and attainment.
panel = da.inspect(sources, view="selection", select=selection)
assert panel.selected == panel.requested == 2
assert panel.available == 4
assert panel.status == "complete"
assert len(list(panel.references())) == 2

# Publish the declared records to a new destination.
da.export(panel, format="selection", out="handoff/panel.selection.json")
# Publish the selected final sequences to a new destination.
receipt = da.export(
    sources,
    view="sequences",
    select=panel,
    format="fasta",
    out="handoff/panel.fasta",
)
assert receipt.design_refs == tuple(panel.references())
```

The snapshot's representation shows totals, status and source count. Its
`references()` iterator reads already saved membership and performs no source
scan. `to_dict()` explicitly serializes the full membership document.

For a total of two designs, replace the allocation with `Take(count=2)`.
This chooses the first two eligible designs across the ordered source union.
Use `Take(count=2, policy="random", seed=23)` for a total random sample.
`count` and `per_cell` are mutually exclusive. Zero counts are valid, including
an explicit empty per-cell mapping; cells omitted from that mapping are outside
its allocation scope.

## Use the same policy from the CLI

Save the policy as a declared YAML or JSON document. For a total sample:

```yaml
schema: dense_arrays.library-selection.v1  # Request type and wire-format version.
filter:
  groups: [A]
take:
  count: 2
  policy: random
  seed: 23
  shortfall: error
```

With this document saved as `panel.yaml`, use the two libraries above:

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/left runs/right --view selection \
  --selection panel.yaml --json
# Write the declared selection or document to a new destination.
dense-arrays export runs/left runs/right --selection panel.yaml \
  --format selection --out handoff/cli.selection.json
# Write the declared selection or document to a new destination.
dense-arrays export runs/left runs/right --selection handoff/cli.selection.json \
  --view sequences --format fasta --out handoff/cli.fasta
```

For quotas, replace `count` with a `per_cell` mapping whose keys are full
`run_id/cell_id` references. The current single-cell generation workflow uses
`default` as its cell ID. A bare cell ID is accepted only when it resolves to
one cell across all sources. Unknown and ambiguous IDs fail even for zero quotas.

A selection file and convenience filter flags are mutually exclusive. Omit
`--all` when using a bounded `take` or a saved snapshot: combining them is an
error. A filter without a `take` still requires `--all` for a complete record
export. `--limit` bounds a displayed page; it never changes the sampling policy.

## Reuse the saved panel

Continue the Python example:

```python
# Read saved composition and search metrics.
quality = da.inspect(sources, view="quality", select=panel)
assert quality.to_dict()["selection"]["designs"] == 2
assert quality.to_dict()["source_runs"][0]["attainment"]["accepted"] == 2
# Publish the declared records to a new destination.
da.export(quality, out="handoff/panel-quality.json")
# Publish the declared records to a new destination.
da.export(sources, select=panel, format="bundle", out="handoff/panel-library")
```

With the optional playback extra, render the same scope:

```bash
# Render a figure from the selected saved records.
dense-arrays render runs/left runs/right \
  --selection handoff/panel.selection.json --view library-quality \
  --out handoff/panel-quality.png
```

Quality reports retain each source's original attainment. Search denominators
cover the available native histories and explicitly identify missing histories.
Composition describes the selected designs. The figure and export receipt
identify the saved selection.

Later commits do not advance a snapshot. Missing revisions, changed committed
manifests, missing members or changed selected records fail. Explicit source
paths may locate moved artifacts, but must preserve the saved source order and
identities. Snapshot files save relative source locators; Python can reload one
with `SelectionSnapshot.from_dict(value, base=path.parent)`.

A snapshot stores references, not the design evidence itself. Share a
[portable bundle](bundles.md) when the original sources may be unavailable.
Selection supports a native run, a bundle, or an ordered list containing either.
The same saved membership applies to quality reports and rendering from these
sources. Bundle reports describe included designs and label unavailable search
history explicitly.

## Handle shortfalls and work limits

By default, insufficient eligible designs raise `SelectionShortfall`; a file
export publishes nothing. Its `counts` identify requested, available, selected
and missing designs per allocation. Excess candidates in another cell never
fill a deficit.

Set `shortfall="allow_partial"` only when a smaller panel is acceptable. The
snapshot and export receipt retain `status="partial"` and the exact shortfall.
CLI inspection succeeds with exit `0`; an export of the qualified partial
selection publishes its output and exits `3`. The source run's status does not
change. Default shortfall errors exit `2` with code `selection_shortfall` and
per-allocation counts in the JSON diagnostic.

Materialization scans the entire eligible source scope to count availability,
even for `first` and zero allocations. It creates no spill files. `ReadLimits`
and the matching CLI read/state caps bound the scan and retained identity
entries; exceeding a cap fails, including with `allow_partial`. Limits count
records and lookup entries, not physical memory or elapsed time. CLI cost
descriptors precede population scans.

Random selection uses `sha256_priority.v1`: priorities hash the declared seed
and full design reference under `dense_arrays.selection_priority.v1`. The
lowest priorities win; source ordinal breaks ties. Chosen members are returned
in source order. Equal DNA with different construction identities remains
eligible as distinct designs. Repeated copies of one full design identity count
once; conflicting copies fail. The policy uses no global random state.
