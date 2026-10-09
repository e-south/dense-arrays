---
title: Prepare curated binding sites
description: Import binding-site tables, require groups and reuse selected parts.
author: Eric J. South
---

# Prepare curated binding sites

The four synthetic 16-base sites below illustrate group requirements and
occurrence identity. Run commands in a new directory with Dense Arrays installed.

```python
from pathlib import Path  # Write the example input table.
import dense_arrays as da  # Prepare, generate and verify records.
from dense_arrays import parts, planning, reporting  # Construct typed requests.
```

## Bind curated parts and requirements

Create this synthetic `parts.csv`:

```csv
part_id,sequence,group
a,ACGTTGCAAGTCCTGA,A
b,AGTCCTGATCGTACCG,A
c,TCGTACCGATGCTTAG,B
d,ACGTTGCAAGTCCTGA,A
```

Save `design.yaml` beside it:

```yaml
schema: dense_arrays.design.v1  # Request type and wire-format version.
parts: {table: parts.csv, format: csv}
length: {maximum: 40}  # Length constraint for the assembled design.
requirements:  # Rules that every accepted design must satisfy.
  - {id: two-A, kind: occurrences, select: {groups: [A]}, min: 2, max: 2}
  - {id: both-groups, kind: group_coverage, groups: [A, B], min: 2}
target: {count: 1}  # Requested number of accepted final designs.
strands: single
limits: {attempts: 1000, active_seconds: 300, solver_seconds: 30}  # Global and per-solver work allowances.
seed: 7  # Seed for versioned candidate streams.
```

```bash
# Resolve inputs and save an executable plan.
dense-arrays plan design.yaml --out design.plan.json
# Generate into a new directory with explicit bounds.
dense-arrays run design.plan.json --out runs/curated
# Read or verify saved evidence without generating again.
dense-arrays inspect runs/curated --verify --json
```

The `a` and `d` rows have different occurrence identities despite equal DNA.
Counts apply to selected occurrences. Group coverage counts distinct supplied
group labels. Upper and lower bounds are inclusive and enforced in CBC, then
independently recounted from persisted placements.

CSV and TSV formats require a header and explicit format. Map other column names
through `columns: {part_id: site, sequence: bases, group: family}`. DNA is strict
uppercase A/C/G/T by default. Opt into
`normalization: {uppercase: true, trim_outer_whitespace: true}` when needed;
internal whitespace and ambiguous bases still fail. Use `id_policy: row` only
when no supplied ID column exists. Metadata columns require an explicit mapping.
For typed Parquet columns or a named Excel worksheet, use the optional
[table readers](../tables.md) with the same request and mappings.

Planning reads inputs, checks supported requirements, and binds their digests;
it does not solve or prove feasibility. Running a saved plan rejects changed
input bytes. Inspection uses the persisted snapshot and does not reopen source
tables. Unknown request fields, duplicate YAML/JSON keys and unsupported
requirements fail before generation.


## Prepare a reusable pool

Prepare once when several designs use the same curated inputs. The pool preserves
the supplied occurrence IDs, normalized sequences, import transformations, ignored
columns and retention policy. Identical sequences with different IDs remain distinct.

Save this `prepare.yaml` beside the `parts.csv` above:

```yaml
schema: dense_arrays.prepare.v1  # Request type and wire-format version.
source: {kind: table, table: parts.csv, format: csv}  # Input model or part collection.
retain: {select: {groups: [A]}}  # Number and policy for keeping eligible candidates.
```

```bash
# Resolve inputs and save an executable plan.
dense-arrays plan prepare.yaml --out prepare.plan.json
# Publish the prepared part pool.
dense-arrays prepare prepare.plan.json --out pools/curated
# Read or verify saved evidence without generating again.
dense-arrays inspect pools/curated --view parts --group A --all --json
# Read or verify saved evidence without generating again.
dense-arrays inspect pools/curated --verify
```

Planning reports exact curated retention counts without scoring, sampling, solving
or creating a pool. Execution rejects changed source bytes. `prepare` accepts
preparation requests/plans, or a resolved generation plan with an explicit
[batch sampling policy](../batches.md). `run` accepts generation
requests/plans.

The matching Python path uses typed values:

```python
# Use the typed requests and operations needed by this example.
from dense_arrays import reporting

table = Path("python-parts.csv")
# Write the example input or request so it can also be used from the CLI.
table.write_text(
    "part_id,sequence,group\na,ACGTTGCAAGTCCTGA,A\nb,AGTCCTGATCGTACCG,A\nc,TCGTACCGATGCTTAG,B\nd,ACGTTGCAAGTCCTGA,A\n"
)
# Keep group A occurrences, including equal DNA with distinct part IDs.
preparation = parts.PreparationSpec(
    source=parts.PartTable(table, "csv"),
    retain=parts.Retention(select=parts.PartFilter(groups=("A",))),
)
prepared_plan = da.plan(preparation)  # Preview retention without creating a pool.
assert prepared_plan.preview["retained_parts"] == 3
pool = da.prepare(prepared_plan, out="pools/python-curated")  # Publish reusable parts.
assert da.inspect(pool, verify=True).retained_parts == 3

# Read saved prepared part records.
view = da.inspect(
    pool, view="parts", all=True, read_limits=reporting.ReadLimits(records=10)
)
assert view.cost.records_estimate == 3
with view.records() as rows:
    assert [row.part_id for row in rows] == ["a", "b", "d"]
    assert rows.examined == rows.returned == 3

# Bind the immutable pool to a new generation request.
reused = da.run(
    planning.DesignSpec(
        parts=parts.PoolSource(pool),
        length=planning.Length(maximum=40),  # Permit shorter accepted sequences.
        strands="single",
    ),
    out="runs/from-pool",
)
assert da.inspect(reused, verify=True).accepted == 1
```

A design file references the same pool with
`parts: {pool: pools/curated, select: {groups: [A]}}`; relative paths resolve from
that file. Pool inspection and reuse do not need the original table. Move the
complete pool directory and reference its new location explicitly.

Part filters combine values within a field with OR and different fields with AND.
Length ranges are inclusive; unknown IDs/groups or unavailable metrics fail
without rescoring. For a saved predicate, pass `--selection FILE` using the
`dense_arrays.part-filter.v1` schema. It is exclusive with `--part-id`/`--group`.
