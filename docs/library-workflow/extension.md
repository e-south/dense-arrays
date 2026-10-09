---
title: Add designs to a saved library
description: Generate an additional collection with inherited rules and frozen parent exclusions, without changing earlier runs.
author: Eric J. South
---

# Add designs to a saved library

Create a child run when you need more designs under the same requirements.
The child has its own target, seed and effort budget. Earlier runs keep their
original status and results, and their accepted sequences are excluded from the
new collection.

## Extend an existing run from the CLI

Start with a terminal run from the [saved-library guide](../library-workflow.md).
Save the following as `extend.yaml`; its parent path is relative to that file:

```yaml
schema: dense_arrays.extension.v1  # Request type and wire-format version.
parent: {run: runs/first}
additional: 1
limits: {attempts: 1000, active_seconds: 300, solver_seconds: 30}  # Global and per-solver work allowances.
seed: 19  # Seed for versioned candidate streams.
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan extend.yaml --out extension.plan.json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/first --view plan --compare extension.plan.json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run extension.plan.json --out runs/additional
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/additional --verify
# Write the declared selection or document to a new destination.
dense-arrays export runs/additional --view sequences --all \
  --format fasta --out additional.fasta
```

The preview reports the additional target and number of excluded sequences.
For single-cell plans, comparison separates changed target, effort, seed, exclusions and lineage from
unchanged parts, requirements, length, assembly and policies. Input file locations
are reported separately from semantic differences. Use `--json` to save the same
comparison fields for another program.

## Extend a stopped library in Python

This self-contained example starts from eight accepted designs toward a target
of twelve. Run it in a new working directory after [installing the library workflow](../installation.md#use-the-library-workflow).

```python
import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

# Synthetic 16-base sites make input identities distinct from their DNA.
sites = {
    "A": "ACGTTGCAAGTCCTGA",
    "C": "GATCAGTACCTAGGTC",
    "G": "TTGACCGATAGCTACG",
    "T": "CAGTTCGATGACCTAG",
}
# Generate under the declared bounds into a new output directory.
parent = da.run(
    planning.DesignSpec(
        parts=[parts.Part(name, seq, group=name) for name, seq in sites.items()],
        length=planning.Length(maximum=64),
        strands="single",
        target=planning.Target(count=12),
        limits=planning.Limits(attempts=8),
    ),
    out="runs/parent",
)
# Read saved run state and attainment.
parent_summary = da.inspect(parent)
assert (parent_summary.accepted, parent_summary.target) == (8, 12)
assert parent_summary.state == "stopped"

# Request additional designs while excluding the saved parent library.
extension = planning.ExtensionSpec(
    parent=planning.ParentRun(parent.path),
    additional=4,
    limits=planning.Limits(attempts=100, active_seconds=300, solver_seconds=30),
    seed=19,
)
# Resolve the request and bind its input records before execution.
extension_plan = da.plan(extension)
# Read saved resolved rules and input bindings.
comparison = da.inspect(parent, view="plan", compare=extension_plan)
assert "requirements" in comparison.unchanged_fields
assert "exclusions" in comparison.changed_fields

# Generate under the declared bounds into a new output directory.
child = da.run(extension_plan, out="runs/child")
# Read saved run state and attainment. Recount stored evidence before returning.
child_summary = da.inspect(child, verify=True)
assert (child_summary.accepted, child_summary.target) == (4, 4)
assert child_summary.state == "completed"
assert da.inspect(parent).state == "stopped"
```

## Compare design rules

`inspect(..., view="plan", compare=...)` accepts saved generation plans, native
runs, and individual [plan evidence records](bundles.md#inspect-the-included-design-rules).
For matrices, [compare named combinations](matrices.md#compare-planned-combinations)
through the same operation. The comparison reads persisted rules and never solves. `changed_fields` and
`unchanged_fields` summarize categories; `changes` contains typed `PlanChange`
values with `path`, `kind`, `before` and `after`.

```python
target_change = next(
    change for change in comparison.changes if change.path == ("target", "count")
)
assert (target_change.before, target_change.after) == (12, 4)
excluded = [change for change in comparison.changes if change.path[0] == "exclusions"]
assert len(excluded) == parent_summary.accepted
assert all(change.kind == "added" for change in excluded)
```

Part and requirement changes use their declared IDs, so changing a rule's minimum
count identifies that exact field. Added or removed records carry the complete
record on the applicable side. `kind` distinguishes an absent value from an
explicit null. Ordered part, requirement and input lists have separate order
fields; a reorder is not reported as changed record content. Exclusions use
sequence identities and retain their design references and record digests.

`dense_arrays.plan_comparison.v2` keeps input locations in `locations`, separately
from semantic changes. A location-free evidence record reports `null` there; an
executable plan with no file inputs reports an empty list. Comparisons preserve
import provenance and policy versions, and reject unsupported plan encodings.
They do not claim that a changed field caused a difference in generated designs.

The CLI displays bounded change examples; `--json` or a comparison export retains
every value. Read limits apply to the two plan documents and their combined
retained identities. Comparison and export preserve parent artifacts.

## Interpret additional generation

Four additional designs do not rewrite the parent's target or mark it complete.
The child may also stop short if its effort budget or available designs are
exhausted. CLI generation then exits `3`, retaining its accepted prefix. Inspect
the child attempts and diagnostics to distinguish duplicate effort, rejected
candidates and exhausted enumeration.

## Combine saved runs

Pass runs in the order you want them to appear. The same design appears once,
even if its source is repeated. Equal DNA from different design histories retains
each full design reference. Repeated references with conflicting full records
fail, including when exporting only the sequence projection.

Continue the Python example above:

```python
# Read saved final sequences with design identities.
combined = da.inspect([parent, child, parent], view="sequences", all=True)
with combined.records() as records:
    sequences = list(records)
assert len(sequences) == 12
assert len({row.sequence_id for row in sequences}) == 12

# Publish the selected final sequences to a new destination.
receipt = da.export(
    [parent, child],
    view="sequences",
    all=True,
    format="fasta",
    out="combined.fasta",
)
assert receipt.records == 12
# Publish the selected coordinate annotations to a new destination.
da.export(
    [parent, child],
    view="placements",
    all=True,
    format="tsv",
    out="combined-placements.tsv",
)
```

The matching CLI accepts more than one source path:

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/parent runs/child --view sequences --all --json
# Write the declared selection or document to a new destination.
dense-arrays export runs/parent runs/child --view sequences --all \
  --format fasta --out cli-combined.fasta
```

Combined record views support `DesignFilter` and the existing CLI filter flags.
Use a full `run/cell/design` reference when a local design name is ambiguous,
`run/cell` for cells, and `collection_id/part_id` for parts. A part label shared
by multiple runs from the same collection is unambiguous. Groups are literal
supplied labels across collections. Unknown or ambiguous selectors fail before
the first result is emitted.

Each `LibraryView` binds all source revisions. Its `sources` property and JSON
output report those individual revisions; the union descriptor itself uses
revision zero. Continuation tokens preserve the bound revisions, source order,
view and filter even if an active run advances. Source paths can move provided
the same evidence remains available in the same order.

Pagination scans the earlier prefix again to reconstruct the exact identity set.
`ReadLimits` applies across all sources, including repeated sources, and export
receipts share the iterator's identity budget. `cost` reports this scan before
iteration. Use `all=True` / `--all` for a complete export. Very large source lists
can exceed the continuation-token limit; use a complete query or fewer sources.
Native record checks run as data is read, so a page does not certify later rows.
File exports publish only after every selected source has been read successfully;
stdout can retain a partial prefix on failure.

Combined sources support designs, sequences, placements and quality. Summary,
attempts and verification still operate on individual runs. Combining
accepted designs does not change either run's target or completion status.

## Assess a combined library

Quality reports use the same ordered sources and design filters as record
inspection. They include combined composition and per-cell metrics, with each
cell identified by `run/cell`. A shared part collection is counted once in the
eligible supply; equal local part names from different collections stay separate.
Group labels are combined literally, without inferring biological equivalence.

```python
# Use the typed requests and operations needed by this example.
from dense_arrays import reporting

# Read saved composition and search metrics.
quality = da.inspect([parent, child], view="quality")
metrics = quality.to_dict()
assert metrics["selection"]["designs"] == 12
assert metrics["selection"]["distinct_sequences"] == 12
assert metrics["attainment"] is None
assert [run["state"] for run in metrics["source_runs"]] == ["stopped", "completed"]
assert [run["attainment"]["shortfall"] for run in metrics["source_runs"]] == [4, 0]

# Choose designs by recorded identities and final-sequence properties.
child_only = reporting.DesignFilter(cells=(f"{child.run_id}/default",))
# Read saved composition and search metrics.
selected_quality = da.inspect([parent, child], view="quality", select=child_only)
assert selected_quality.to_dict()["selection"]["designs"] == 4
# Draw the selected saved evidence without generating new sequences.
da.render(selected_quality, view="library-quality", out="child-quality.png")
```

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/parent runs/child --view quality --json
# Render a figure from the selected saved records.
dense-arrays render runs/parent runs/child --view library-quality \
  --out combined-quality.png
```

The combined report has no merged target or completion state. `source_runs`
retains every distinct run's original attainment, status and search outcomes;
`selection` reports the selected design and sequence counts. Its aggregate search
context counts each source run once, including duplicate and rejected attempts.
Copies with conflicting design, plan or search evidence fail explicitly.

Filters change composition and usage denominators. Search context and eligible
source-plan supply remain visible even when a filter selects no designs from a
source. Combined requirement results retain their full cell references, so
identically named requirements do not silently merge. Empty selections return
null metrics with an explicit reason.

Usage pagination preserves all bound revisions and the filter. Each page retains
the complete aggregates and limits the top-level and per-cell usage tables.
Quality requires one consistent revision for each run within a query; compare
different snapshots in separate reports. Plots display the selected population
and retain the complete rendered report in PNG metadata. Text and JSON reports
need no graphics dependencies; PNG rendering uses the playback extra.

## Understand exclusions and effort

Extension retains the parent's resolved parts, requirements, length, assembly,
strand eligibility and policy versions. It does not accept overrides to those
rules. The new `limits` and `seed` are required, and `additional` is a positive
single-cell count. Changing the design rules requires a new design request.
Use [request revision](handoffs.md#save-and-revise-design-inputs) to edit the rules
and retain explicit exclusions from a resolved extension.

Planning requires a terminal, verifiable parent with no active writer. It pins
the committed revision and accepted-record digest, verifies its evidence, and
copies the parent and ancestor sequence exclusions into the child plan. It does
not rerun the parent's solver or require original source tables. A resolved child
plan can execute after the parent paths become unavailable.

Every candidate is checked against exact final-sequence identities. A parent or
ancestor duplicate consumes an attempt and records `parent_duplicate`, the
candidate sequence digest, and the matching full design reference. It does not
count toward the additional target. A second extension inherits those ancestor
exclusions, so it cannot regenerate earlier sequences. Enumeration restarts and
may revisit excluded candidates; allow enough attempts for that work.

Parent resolution accepts `reporting.ReadLimits` through `plan(..., read_limits=...)`.
CLI equivalents are `plan --max-read-records` and `--max-identity-entries`.
These bound evidence reading and retained lookup entries, not the child's
generation budget. [Recovery](recovery.md) continues eligible interrupted runs.
[Portable bundles](bundles.md)
carry the frozen exclusions needed to verify included designs.


## Extend a design matrix

For a matrix parent, provide an additional count for every cell, including
zeros. At least one count must be positive. An integer additional count is
accepted only for a single-cell parent. Each child retains the parent's cell
identities, parts, requirements and assembly policies.

```python
# Combine named choices and assign an explicit target to each combination.
matrix_request = planning.MatrixSpec(
    base=planning.DesignSpec(
        parts=[parts.Part(name, seq) for name, seq in list(sites.items())[:3]],
        length=planning.Length(maximum=16),
        strands="single",
    ),
    axes={"condition": {"a": planning.Variant(), "b": planning.Variant()}},
    allocation=planning.Allocation(per_cell=1),
    max_cells=2,
)
# Generate under the declared bounds into a new output directory.
matrix_parent = da.run(matrix_request, out="runs/matrix-parent")
# Request additional designs while excluding the saved parent library.
matrix_extension = planning.ExtensionSpec(
    parent=planning.ParentRun(matrix_parent.path),
    additional={"condition=a": 2, "condition=b": 0},
    limits=planning.Limits(attempts=30),
    seed=23,
)
# Resolve the request and bind its input records before execution.
matrix_plan = da.plan(matrix_extension)
# Save the resolved plan with its input bindings; keep the destination new.
matrix_plan.write("matrix-extension.plan.json")
# Generate under the declared bounds into a new output directory.
matrix_child = da.run(matrix_plan, out="runs/matrix-child")
# Read saved run state and attainment. Recount stored evidence before returning.
matrix_summary = da.inspect(matrix_child, verify=True)
assert matrix_summary.accepted == matrix_summary.target == 2
assert matrix_summary.cells["condition=a"].accepted == 2
assert matrix_summary.cells["condition=b"].state == "inactive"
```

Run the saved plan through the CLI with a new destination:

```bash
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run matrix-extension.plan.json --out runs/matrix-child-cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/matrix-child-cli --verify --json
```

A zero additional target keeps a cell visible and preserves its ancestor
exclusions for a later extension. It does not erase earlier designs. Equal DNA
in different cells remains independent: accepting a sequence in one cell never
excludes it from another. Duplicates within a cell consume the child's effort
and identify the matching parent or ancestor design.

The resolved plan freezes per-cell accepted libraries and source-part evidence.
It remains executable after the original source files or parent run are moved.
`export --view request` preserves those exclusions in the matrix's `exclude`
mapping, so editing and replanning cannot silently discard the accepted history.
Each frozen library records its cell's target, attainment, revision and original
design references. Exclusion mappings currently preserve cell identities; cell
plans containing named exclusions execute within their matrix.

Use [run recovery](recovery.md) to continue a measured interruption of the child.
Use combined inspection or export to select both parent and child results;
lineage alone does not add ancestors to a selected library.
