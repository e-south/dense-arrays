---
title: Library selection contracts
description: Select accepted designs with typed filters, explicit quotas and reusable source-bound snapshots.
author: Eric J. South
---

# Library selection contracts

Selection chooses accepted records from declared source revisions. It preserves
run attainment and effort as source-scoped evidence. The [saved-panel guide](../../library-workflow/selection.md)
provides a complete workflow; [exports](exports.md) defines publication.
The filtering example uses `runs/curated` from the
[generation recipe](generation.md#constrained-design). The quota example assumes
a saved matrix run with cell IDs `a`, `b`, `c` and `d`, each with six eligible designs.

## Result selection

Filters are typed for their records; there is no universal `Selection` whose
fields are mostly inapplicable.

| Type | Fields and meaning |
| --- | --- |
| `PartFilter` | Scoped part IDs, groups and available named part metrics; applies to pool reuse, curated retention and part inspection. |
| `CandidateFilter` | Saved preparation candidate indices, outcomes, reasons, recipe IDs and recipe-local score bands. |
| `DesignFilter` | Full/bare design IDs, cell references, selected part/group identities and available design metrics. |
| `AttemptFilter` | Attempt IDs, cell references and outcome codes; applies to attempts/diagnostics. |
| `PlanFilter` | Complete semantic plan IDs in an explicit plan inventory. |
| `LibrarySelection` | A `DesignFilter` plus optional bounded `take` policy over accepted designs. |
| `SelectionSnapshot` | Materialized ordered design references, source revisions, policy/version, requested/available/selected counts, and shortfall status. |

Multiple values within a filter field mean OR; different fields mean AND.
Part/group matching uses placements, not incidental substrings. Unknown IDs,
metrics and mismatched view/filter types are errors. Bare references must resolve
unambiguously; full design references are `(run_id, cell_id, design_id)`, full
cell references are `(run_id, cell_id)`, and parts include collection identity.
Selectors do not infer equivalence across cells of different runs.

CLI convenience flags `--design-id`, `--cell`, `--part-id`, `--group`,
and `--outcome` lower into the matching filter type. They are exclusive with
`--selection FILE`. The file has a declared part-filter, design-filter,
attempt-filter, library-selection or selection-snapshot schema; reject a type
unsupported by the view. There is no executable expression language or hidden
merge precedence. Range endpoints are inclusive.

Views include summary, request, plan, quality, parts, candidates, designs,
sequences, placements, attempts, diagnostics and selection. Pools support
parts/quality and saved preparation candidate evidence; design queries
apply to accepted designs and their projections; attempt queries apply to attempts
and their diagnostics. Request/plan views reject record filters. Help lists the
supported combinations. Comparisons consume like report types and label
incomparable metrics; they never imply causality.

Multiple explicit run/bundle sources form a union deduplicated by full design
reference, never by DNA string. A repeated reference with conflicting content
is an integrity failure. Order is source argument order, cell order and accepted
ordinal, keeping the first identical reference. Repeated sources do not multiply
records. Sources with missing/incompatible evidence fail before export.

Record pages use `--limit`/`--after`; `--all` / Python `all=True` streams
the declared scope. Quality aggregates cover that scope, not the page size.
A limit is presentation pagination; it is never a library selection policy.

Design metric filters support length, GC fraction, placement count and
packing density; part filters support length and available named scores.
Unavailable metrics fail without rescoring. Null values do not match a range;
reports count exclusions due to null observations.

```yaml
schema: dense_arrays.library-selection.v1
filter:
  groups: [A]
  metrics:
    gc_fraction: {min: 0.2, max: 0.8}
```

```python
import dense_arrays as da
from dense_arrays import reporting

selected = reporting.LibrarySelection(
    filter=reporting.DesignFilter(
        groups=("A",),
        metrics={"gc_fraction": reporting.Range(min=0.2, max=0.8)},
    ),
)
rows = da.inspect("runs/curated", view="sequences", select=selected, all=True)
```

```text
dense-arrays inspect runs/curated --view sequences --selection selection.yaml --all
```

### Bounded selection

Use optional `take` to request exactly one total `count`
or explicit `per_cell` allocations, with policy `first` or `random`.
Counts are nonnegative integers; zero allocations are explicit. A total never
implies per-cell balancing. Unlisted cells in a per-cell request are outside
the selection. There is no automatic redistribution.

`first` preserves the documented candidate order. `random` samples without
replacement and requires a seed. The snapshot records the algorithm version,
seed, exact input revisions/order, candidate scope and selected full references.
Stable tie/order rules and the specified seed derivation make the same request
over the same ordered snapshots repeatable. Do not reuse preparation MMR or
generation weighting as an unstated result-selection policy.

```yaml
schema: dense_arrays.library-selection.v1
filter: {cells: [a, b, c, d]}
take:
  per_cell: {a: 6, b: 6, c: 6, d: 6}
  policy: random
  seed: 23
  shortfall: error
```

```python
import dense_arrays as da
from dense_arrays import reporting

panel = reporting.LibrarySelection(
    filter=reporting.DesignFilter(cells=("a", "b", "c", "d")),
    take=reporting.Take(
        per_cell={"a": 6, "b": 6, "c": 6, "d": 6},
        policy="random",
        seed=23,
        shortfall="error",
    ),
)
selection = da.inspect("runs/matrix", view="selection", select=panel)
# Reuse the pinned panel; this export does not draw again.
da.export("runs/matrix", select=selection, format="bundle", out="selected-library")
```

```text
dense-arrays export runs/matrix --selection panel.yaml --format selection --out panel.selection.json
dense-arrays export runs/matrix --selection panel.selection.json --format bundle --out selected-library
dense-arrays render runs/matrix --selection panel.selection.json --view library-quality --out panel.png
```

In this single-run example cell IDs resolve unambiguously. Multi-run requests
use full cell references when necessary. Scientific quotas, grouping and criteria
remain caller-supplied; the tool does not decide which designs deserve experiments.

Default `shortfall: error` returns requested/available counts and publishes
nothing when a requested partition is short. Explicit `allow_partial` retains
the shortfall in the snapshot and export receipt; it never borrows from another
cell or labels the original run complete. An inspection that successfully
reports a partial selection exits 0; exporting an explicitly allowed partial
selection publishes the qualified result and exits 3. Unknown IDs remain errors,
not zero-availability shortfalls.

A snapshot can be reused by inspect/export/render without drawing again.
`format=selection` writes a native selection snapshot, not copies of source
evidence; moving it alone does not make it a portable library.
`format=bundle` supplies that evidence. Reading a snapshot resolves its pinned
sources or reports missing evidence, never substitutes a newer revision.
A bounded take or snapshot is a complete selection declaration; otherwise
record exports require `all=True` / `--all` so a default page cannot become a library.
Contradictory flags, such as `--all` with bounded take, fail.
Summary/quality reports and editable request exports already declare their scope
and do not require record-pagination flags. A filtered empty set is valid;
export preserves zero records and headers/manifest scope. A positive take over
that set follows the explicit shortfall rule.
