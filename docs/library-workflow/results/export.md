---
title: Export sequences and placements
description: Save final DNA and annotated coordinates with shared design identities.
author: Eric J. South
---

# Export sequences and placements

Start with the saved run in [Generate a saved library](../../library-workflow.md).
The CLI examples read `runs/first`; the Python examples read `runs/python-first`.

```python
import dense_arrays as da  # Read and export the saved run.
from dense_arrays import reporting  # Build a typed filter.

result = "runs/python-first"  # Output of the first-library Python example.
with da.inspect(result, view="designs", limit=1).records() as records:
    design = next(records)  # Select one saved identity without rerunning generation.
```

## Select and export sequences

Use `sequences` for final DNA and `placements` for its annotated coordinates.
Both carry the same full `design_ref`, sequence digest and plan identity. A
placement includes its original part ID, group, orientation, source collection,
and zero-based, half-open coordinates in the final sequence. Known core
coordinates and orientations are transformed into that same frame; unknown
core annotations stay null.

```bash
# Read or verify saved evidence without generating again.
dense-arrays inspect runs/first --view sequences --limit 10 --json
# Export every matching record to a new file.
dense-arrays export runs/first --view sequences --all \
  --format fasta --out library.fasta
# Export every matching record to a new file.
dense-arrays export runs/first --view sequences --all \
  --format csv --out sequences.csv
# Export every matching record to a new file.
dense-arrays export runs/first --view placements --all \
  --format tsv --out placements.tsv
```

Python uses the same query and output formats:

```python
# Use the typed requests and operations needed by this example.
from dense_arrays import reporting

# Filter full design identities and final-sequence composition together.
selected = reporting.DesignFilter(
    design_ids=(design.reference,),
    metrics={"gc_fraction": reporting.Range(min=0.0, max=1.0)},
)
# Read saved final sequences with design identities.
sequences = da.inspect(result, view="sequences", select=selected, all=True)
print(sequences.cost.to_dict())  # Inspect the read budget before opening records.
with sequences.records() as records:
    for row in records:
        print(row.design_ref, row.sequence, row.length, row.gc_fraction)

# Write all matching sequences; the destination must be new.
sequence_receipt = da.export(
    result,
    view="sequences",
    select=selected,
    all=True,
    format="csv",
    out="python-sequences.csv",
)
# Export annotations with the same selected design identities.
placement_receipt = da.export(
    result,
    view="placements",
    select=selected,
    all=True,
    format="tsv",
    out="python-placements.tsv",
)
assert sequence_receipt.records == 1
assert sequence_receipt.design_refs == placement_receipt.design_refs
```

`DesignFilter` accepts bare or full design IDs, cell references, part IDs, groups,
and inclusive ranges over `length`, `gc_fraction`, `placement_count` and
`packing_density`. Part/group filters match recorded placements; an incidental
sequence match does not count. Values within a field mean OR and fields combine
with AND. Unknown identities or metrics fail; a valid filter matching no designs
returns an empty collection.

Use `--design-id`, `--cell`, `--part-id` and `--group` for equivalent CLI
predicates, or save `selected.to_dict()` as a declared
`dense_arrays.design-filter.v1` JSON/YAML file and pass `--selection FILE`.
Filter files and flags are exclusive. Filters apply to designs, sequences,
placements and quality reports. Placement pagination
can continue inside a design using the returned cursor, preserving row order.

| Format | Supported views | Content |
| --- | --- | --- |
| JSON | designs, sequences, placements, attempts, pool parts | Versioned record export with source revision and predicate |
| JSON documents | request, plan, summary, quality, diagnostics, plan/quality comparison | Native schema with the document's declared scope |
| CSV / TSV | sequences, placements | Scalar columns; null annotations become empty cells |
| FASTA | sequences | Encoded full design reference in each header, plus sequence digest |
| Bundle | designs | Selected records, resolved plan evidence and original source context |

Ordinary record exports require `all=True` / `--all`, including exports narrowed
by a `DesignFilter`. This exports every matching record rather than an inspection
page. A bounded `LibrarySelection` or saved `SelectionSnapshot` already declares
the exported membership and rejects `--all`; see
[saved panels](../selection.md).
Empty CSV/TSV exports retain column headers. Every file destination is
create-only and published after its complete write succeeds.
Receipts identify the source revision, filter, output count and design references.
Record export receipts and JSON headers also include `manifest_digest`, the
SHA-256 digest of each source's canonical committed manifest. The exporter
rechecks that manifest before publishing the completed file. File receipts
include a separate byte count and checksum of the exported file. A manifest
digest identifies committed metadata; it is not a checksum of every source
record or a substitute for `inspect --verify`.

Use `--out -` to stream pure data to stdout. Diagnostics and receipts then go to
stderr, including with `--json`; for a file export, `--json` puts the receipt on
stdout instead. A stream can retain partial bytes on failure and exits nonzero.
Python also accepts a writable text stream and leaves it open. Export reads
stored records without solving or scoring, and needs no plotting extra.

Design, sequence and placement queries also accept a list of runs in Python,
or multiple positional paths in the CLI. See
[combine saved runs](../extension.md#combine-saved-runs) for ordering,
identity checks and pagination. Attempt and pool queries still use one source.
Use [saved panels](../selection.md) for total or per-cell sampling
and reusable selection snapshots. Use
[portable bundles](../bundles.md) to inspect selected designs after
their original sources become unavailable. See
[save requests, plans and reports](../handoffs.md) for editable inputs,
preparation settings and scoped document exports.

For pool-backed designs, `collection_id` is the original pool ID. Inline/table
collections use a digest of ordered part records, import evidence and source byte
digests. Target, seed, source paths and extension lineage do not change this
identity.

To generate more designs with the same rules, follow
[add designs to a saved library](../extension.md). Extension creates
a child target and excludes parent/ancestor sequences while preserving earlier
results. The guide includes a stopped-parent example and plan comparison.
