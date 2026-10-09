---
title: Library export contracts
description: Publish native documents, sequence projections and portable selected libraries with explicit receipts.
author: Eric J. South
---

# Library export contracts

Exports publish an explicitly selected population and return a receipt naming
its sources and outputs. [Selection](selection.md) defines membership;
[artifact schemas](artifacts.md#schema-compatibility) define source-manifest
bindings. These examples assume saved `runs/initial` and `runs/additional`
from the [extension recipe](generation.md#revise-or-extend-a-library), with
new output destinations. User-facing format examples are in the
[export guide](../../library-workflow/results/export.md).

## Tabular and sequence output

Export supports JSON, CSV/TSV, and FASTA where meaningful. Sequence tables contain
full design reference, native sequence digest, sequence, length, GC, and plan
identity. Placement tables join on the full design reference and add placement
ID, source collection/part ID, group, orientation, start/end, core geometry when
known, and lineage references. Missing values stay null/empty under a documented
schema; no stringified nested objects masquerade as scalar columns. Optional
caller metadata needs explicit column selection and a namespace. Every export
identifies its table schema and selection in the API receipt or CLI diagnostic.

Report/request/plan exports use JSON: reports retain their versioned envelope,
while editable requests and plans use their native input schema. Record exports
support JSON or the documented scalar CSV/TSV projection; FASTA requires the
sequences view. Bundle and selection formats require the designs view and a
complete selection declaration. Reject unsupported view/format combinations
before reading the population. Exported data and its receipt have separate
schemas; an editable request must never be wrapped as an inspection response.

```text
dense-arrays export runs/initial runs/additional --view sequences --all \
  --format csv --out library.csv
dense-arrays export runs/initial runs/additional --view placements --all \
  --format tsv --out placements.tsv
dense-arrays export runs/initial runs/additional --view sequences --all \
  --format fasta --out library.fasta
```

```python
import dense_arrays as da

sources = ["runs/initial", "runs/additional"]
receipt = da.export(
    sources, view="sequences", all=True, format="csv", out="python-library.csv"
)
placement_rows = da.inspect(sources, view="placements", all=True)
with placement_rows.records() as records:
    for placement in records:
        print(placement.design_ref, placement.start, placement.end)
```

Use full encoded design references in FASTA headers so repeated local IDs do
not collide; retain a separate sequence digest. Pure data stdout never contains diagnostic text or JSON error envelopes.
Use `--out -` for a stream; its receipt and diagnostics remain on stderr even
with `--json`. For file exports, the default receipt goes to stderr and `--json`
requests a versioned receipt on stdout. Data format and receipt format are
independent. Atomic
create-only publication is required for explicit file destinations; existing
paths fail. Shell redirection may leave partial bytes on failure and returns a
nonzero exit. Parquet output requires a separately specified and qualified
projection; optional Parquet input support does not imply export support.
## Portable bundle

```text
dense-arrays export runs/initial runs/additional --view designs --all \
  --format bundle --out library-bundle
dense-arrays inspect library-bundle --verify
```

```python
import dense_arrays as da

sources = ["runs/initial", "runs/additional"]
# Bundle the selected records and the evidence needed to read them after a move.
bundle_receipt = da.export(
    sources, view="designs", all=True, format="bundle", out="python-library-bundle"
)
verification = da.inspect("python-library-bundle", verify=True)
```

`bundle` requires an explicit new directory and a complete selection declaration
(`--all`, optionally filtered, or a bounded take/selection snapshot); a default page cannot accidentally become a whole
library. The versioned bundle manifest records selection, origin snapshots,
parent/exclusion identities, metric versions, and content digests. It includes
selected design/placement records, their requirement evaluations, and the bound
part/plan/pool evidence needed to inspect and verify those records without access
to the original paths. Required referenced bytes are copied, not silently left
as remote locators. Unavailable evidence fails the self-contained export.

The bundle represents a selected collection, not a new completed run. Unselected
designs, debug logs, scorer binaries, and environments are not implicitly shipped.
Origin-wide summaries, if included, retain their original scope; they are not
recomputed completion claims for the subset. Verification names the included
record/evidence boundary and does not claim to reproduce the original execution.
Exclusion-set identities can remain provenance references when their records
are unnecessary to verify selected placements; that limitation is reported.

Exported plans use relative locators and preserve semantic plan identity while
receiving new export-byte digests. Reserved machine-location fields are remapped;
raw configuration/log text is not copied as an unreviewed path leak. Caller
metadata follows the explicitly selected columns. Publication uses owned staging
and a manifest commit; interrupted staging is not a valid bundle. Moving the
completed bundle must preserve verification and record access using only its
contained evidence. Verification does not invoke the solver.
