---
title: Generate a design matrix
description: Generate named part and requirement combinations with explicit targets and per-cell results.
author: Eric J. South
---

# Generate a design matrix

Use a matrix to generate named alternatives in one saved run. `plan` expands the
combinations, validates each concrete request and saves its resolved parts and
policies. `run` shares the effort budget across cells and records each cell's
attainment, shortfall and termination reason.

A **matrix cell** is one design combination, identified by its named choices.
For example, `parts=pool_a,spacer=short` selects one choice from each of
two axes. Each cell has a resolved recipe and a target number of accepted designs.

| Term | Meaning |
| --- | --- |
| Axis | A named dimension of alternatives, such as `parts` or `spacer`. |
| Choice | One named alternative on an axis, with explicit part or rule substitutions. |
| Matrix cell | One combination of choices, with its own source, rules and target. |
| Candidate batch | Parts offered to one packing search within a cell. |
| Attempt | One recorded search outcome within a cell and candidate batch. |
| Design | One accepted sequence with placements and construction evidence. |

Use a single `DesignSpec` when there is only one recipe. Pools supply reusable
parts; a matrix organizes generation requests. The same pool can supply several
cells, and several designs can belong to one cell.

A target of 20 requests 20 accepted designs. It does not declare biological or
technical replicates. Axis labels describe caller-selected alternatives;
experimental conditions, cell types and replicate relationships require an
explicit experimental design. The tool does not infer them from names.

## Declare variants and targets

This example combines two upstream and two downstream alternatives. Each
variant replaces a complete part with the same identity. An empty `Variant()`
retains the base. Requirements can also be replaced by identity using
`Variant(requirements=(...))`.

Use `Variant(add_requirements=(...))` to introduce constraints only in selected
choices. Added IDs must be absent from the base and every other selected variant;
collisions fail during planning. Replacements still require an existing ID.
Additions follow the base requirements in declared axis order, and the combined
rules are validated together. The [curated library example](curated-example.md)
uses this to compare a baseline with a fixed pair and spacer.

Run in a new directory after [installing the library workflow](../installation.md#use-the-library-workflow):

```python
import json
from pathlib import Path

import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

# Declare the part collection, sequence bounds and generation policy.
base = planning.DesignSpec(
    parts=[
        parts.Part("up", "ACGTTGCAAGTCCTGA"),
        parts.Part("down", "GATCAGTACCTAGGTC"),
    ],
    length=planning.Length(maximum=32),
    strands="single",
    requirements=[
        planning.Fixed("up-required", "up", "forward"),
        planning.Fixed("down-required", "down", "forward"),
    ],
)
# Combine named choices and assign an explicit target to each combination.
request = planning.MatrixSpec(
    base=base,
    axes={
        "up": {
            "original": planning.Variant(),
            "alternate": planning.Variant(parts=[parts.Part("up", "TTGACCGATAGCTACG")]),
        },
        "down": {
            "original": planning.Variant(),
            "alternate": planning.Variant(
                parts=[parts.Part("down", "CAGTTCGATGACCTAG")]
            ),
        },
    },
    allocation=planning.Allocation(total=6, policy="balanced"),
    max_cells=4,
)
# Resolve the request and bind its input records before execution.
resolved = da.plan(request)
assert [cell.target for cell in resolved.cells] == [2, 2, 1, 1]
for cell in resolved.cells:
    print(cell.cell_id, cell.target, cell.plan.plan_id)
# Save the resolved plan with its input bindings; keep the destination new.
resolved.write("matrix-plan.json")
# Write the example input or request so it can also be used from the CLI.
Path("matrix.json").write_text(json.dumps(request.to_dict(), indent=2))
```

The same request can be previewed and saved through the CLI:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan matrix.json
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan matrix.json --out matrix-plan-cli.json --json
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan matrix-plan.json --json
```

Human output shows up to ten cells, their eligible-part counts and inactive status. JSON includes
every cell, resolved request, input fingerprint and allocation. Save destinations
must be new. Reading a saved plan uses its frozen records;
`resolved.verify_inputs()` explicitly checks whether its source files still match.
Relative source locations resolve against the saved document's directory.

## Select parts for each combination

Cells use `base.parts` by default. Set `sources` to replace the entire eligible
collection for named cells. Each source accepts the same `PartTable`,
`PoolSource`, `BoundParts` or inline `Part` records as a single design request.
This example selects two groups from one prepared pool:

```python
# Write the example input or request so it can also be used from the CLI.
Path("binding-sites.csv").write_text(
    "part_id,sequence,group\nsite_a,ACGTTGCAAGTCCTGA,A\nsite_b,GATCAGTACCTAGGTC,B\n"
)
# Prepare the declared parts or batch and save its identities for reuse.
pool = da.prepare(
    parts.PreparationSpec(parts.PartTable("binding-sites.csv", "csv")),
    out="binding-sites",
)
# Combine named choices and assign an explicit target to each combination.
source_request = planning.MatrixSpec(
    base=planning.DesignSpec(
        parts.PoolSource(pool), planning.Length(maximum=16), strands="single"
    ),
    axes={"parts": {"A": planning.Variant(), "B": planning.Variant()}},
    sources={
        "parts=A": parts.PoolSource(pool, parts.PartFilter(groups=("A",))),
        "parts=B": parts.PoolSource(pool, parts.PartFilter(groups=("B",))),
    },
    allocation=planning.Allocation(per_cell=1),
    max_cells=2,
)
# Resolve the request and bind its input records before execution.
source_plan = da.plan(source_request)
assert [len(cell.plan.request.parts) for cell in source_plan.cells] == [1, 1]
# Save the resolved plan with its input bindings; keep the destination new.
source_plan.write("source-plan.json")
# Generate under the declared bounds into a new output directory.
source_run = da.run(source_plan, out="source-run")
assert da.inspect(source_run, verify=True).accepted == 2
```

Use the saved plan through the CLI:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan source-plan.json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run source-plan.json --out source-run-cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect source-run-cli --verify --json
```

Source keys use complete cell IDs, including every axis. Unknown cells and empty
selections fail before generation. Choice labels do not select sources implicitly.
Each source replacement applies before the axis variants; part substitutions
must reference IDs in that selected source. The combined parts and rules are
then validated together. Sources are never implicitly joined.

Saved plans preserve each selection's pool identity, filter, part records and
input fingerprints. Normal execution verifies the source files before creating
output. [Prepared batches](batches.md) and [extensions](extension.md) embed the
verified parts and their origin fingerprints so they can run after the sources
have moved or been removed.

## Generate and inspect the library

Continue from the plan above:

```python
# Generate under the declared bounds into a new output directory.
run = da.run(resolved, out="matrix-run")
# Read saved run state and attainment. Recount stored evidence before returning.
summary = da.inspect(run, verify=True)
assert summary.accepted == summary.target == 6
for cell in summary.cells.values():
    print(cell.cell_id, cell.accepted, cell.target, cell.termination_reason)
# Publish the declared records to a new destination.
da.export(run, all=True, format="bundle", out="matrix-library")
```

The CLI uses the same saved plan and records:

```bash
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run matrix-plan.json --out matrix-run-cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect matrix-run-cli --verify --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect matrix-run-cli --view quality --json
# Write the declared selection or document to a new destination.
dense-arrays export matrix-run-cli --all --format bundle --out matrix-library-cli
```

Each active cell receives one solver attempt per round, in expansion order.
The base request's attempt and active-time limits cover the whole run; solver
time applies to each attempt. Model construction consumes active time too.
Completed or exhausted cells leave the rotation. Their unused targets are never
transferred to another cell. Zero-target cells remain visible as `inactive` and
never build a model.

Exact-sequence uniqueness applies within each cell. Equal DNA from different
cells retains distinct design references and construction evidence. Attempt
records carry both a global ordinal and `evidence.cell_attempt`, the per-cell
ordinal used for padding randomness. Filters, saved selections, quality reports,
exports and rendering retain the design's cell and resolved part annotations.
See [saved panels](selection.md) for explicit per-cell selection quotas.

`verify=True` independently checks placements, acceptance, global and per-cell
accounting, and every cell's plan binding. A portable bundle retains cell targets
and selected designs; its quality report distinguishes original attainment from
the included population and labels absent search history.

After a measured interruption, `run --resume matrix-run` continues eligible cells
under the original shared budget. It preserves the committed prefix and resumes
the rotation after the last recorded attempt. A completed run is verified and
returned unchanged. See [recovery](recovery.md) for eligibility and crash limits.

To request more designs with the same rules, [extend the matrix](extension.md#extend-a-design-matrix)
with explicit additional targets for every cell. The new run excludes earlier
accepted sequences within each corresponding cell. Use [prepared batches](batches.md)
to save offered parts per cell before generation.

## Compare planned combinations

Compare resolved plans before choosing which one to run. The report matches cells
by their full IDs and parts and requirements by their IDs within each cell.
This example reallocates the same total to different combinations:

```python
revised_request = request.with_changes(
    allocation=planning.Allocation(
        counts={
            cell.cell_id: 3 if cell.cell_id == "down=original,up=original" else 1
            for cell in resolved.cells
        }
    )
)
# Resolve the request and bind its input records before execution.
revised_plan = da.plan(revised_request)
# Save the resolved plan with its input bindings; keep the destination new.
revised_plan.write("revised-matrix-plan.json")
# Read saved resolved rules and input bindings.
comparison = da.inspect(resolved, view="plan", compare=revised_plan)
changed = {change.path: change for change in comparison.changes}
assert changed[("cells", "down=original,up=original", "target", "count")].after == 3
assert "base" in comparison.unchanged_fields
# Publish the declared records to a new destination.
da.export(comparison, out="matrix-changes.json")
```

The CLI uses the same report:

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect matrix-plan.json --view plan --compare revised-matrix-plan.json
# Write the declared selection or document to a new destination.
dense-arrays export matrix-plan.json --view plan --compare revised-matrix-plan.json --out matrix-changes-cli.json
```

`cells` reports the resolved effects; `matrix` reports changes to axes, source
selections, pairing, allocation and sampling declarations. `base` identifies
shared recipe edits. Added and removed cells carry their complete records.
`cell_order`, `axis_order` and `choice_order` expose reordering separately from
content changes. With balanced allocation, a reorder can also change cell targets.

Source locations appear separately under `locations`, with base and per-cell
bindings. Moving unchanged inputs produces no semantic change. Comparison reads
saved evidence without reopening those inputs or running a solver. Run paths
are also accepted, using their persisted matrix plans. Compare two matrix plans;
to compare individual recipes in Python, select their `cell.plan` values.

Human output shows up to twenty changes; `--json` or export retains all values.
`--max-read-records` and `--max-identity-entries` bound the two retained plans.
The Python equivalent is `read_limits=reporting.ReadLimits(...)`.

## Choose combinations and allocation

`pairing="cross_product"` expands all choices in declared axis and choice order.
`pairing="zip"` matches identical choice names across axes, in the first axis's
order. `pairing="explicit"` uses `pairs=({"up": "original", "down": "alternate"}, ...)`
in the supplied order. Repeated pairs, missing choices and unknown identities
fail. Two axes cannot replace the same part or requirement in a cell.
`max_cells` is required and checked before reading source parts.

Choose exactly one allocation:

| Declaration | Result |
| --- | --- |
| `Allocation(per_cell=2)` | Every expanded cell receives two designs. |
| `Allocation(counts={cell_id: count, ...})` | Every cell has an explicit count, including zeros. Missing and unknown cells fail. |
| `Allocation(total=6, policy="balanced")` | Divide the total among active cells; distribute the remainder to the first cells in expansion order. |

For a total smaller than the expanded cell count, name inactive cells with
`zero_cells=(cell_id, ...)`. Every remaining cell must receive at least one
design. Zero-target cells stay in the plan. Targets belong to the allocation;
omit `target` from the base recipe. Per-cell exclusions belong in the matrix's
`exclude` mapping; the base recipe cannot broadcast an exclusion across cells.

Cell IDs sort axis names, for example `down=alternate,up=original`. A cell ID is
a readable address within its matrix, not a global content fingerprint. Labels use
letters or digits initially, then letters, digits, underscores, dots or hyphens.
Reordering axes or choices changes expansion order but preserves the identity
and seed of each unchanged cell. Changing targets, parts or rules changes the
corresponding cell plan fingerprint, even when its cell ID stays the same.
Moving unchanged source files preserves those fingerprints. Native JSON stores axes and choices as ordered arrays
so object-key sorting cannot change allocations; compact input files also accept
named mappings.

Allocation uses `ordered_balanced.v1` or `explicit_counts.v1`; cell seeds use
`matrix_cell_sha256.v1`. Feasibility and search yield remain unknown until
execution. A seed does not promise a particular CBC tie order. The preview names
matrix-wide effort limits and exact-sequence uniqueness within each cell;
it does not multiply the base limits into a total estimate.
