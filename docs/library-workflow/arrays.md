---
title: Share supplied arrays
description: Save existing sequences, parts and placements as a verified, portable collection.
---

# Share supplied arrays

Save existing DNA arrays with their part identities, coordinates and source
metadata. A collection can be inspected, filtered, exported and rendered after
moving it to another directory or computer.

Use an `ArrayCollection` when the inputs are sequences and recorded placements.
Use a [generated-design bundle](bundles.md) when sharing a saved Dense Arrays
run with its plan and execution evidence.

## Save sequences and placements

Run this synthetic binding-site example in a new directory. Both parts are
16 bases long; their shared eight-base sequence joins them into a 24-base array.

```python
import dense_arrays as da  # Use the shared workflow operations.
from dense_arrays.arrays import ArrayCollection, ArrayFilter  # Bind and select arrays.
from dense_arrays.parts import Part  # Preserve original part IDs and annotations.
from dense_arrays.realized import (
    Orientation,
    Placement,
    RealizedArray,
)  # Saved geometry.

# Define the source-oriented parts, including an optional motif-core window.
parts = (
    Part(
        "site-a",
        "ACGTTGCAAGTCCTGA",
        group="regulator-a",
        core_start=2,
        core_end=13,
        core_orientation="forward",
    ),
    Part("site-b", "AGTCCTGATCGTACCG", group="regulator-b"),
)
# Positions are zero-based starts in the final sequence; end coordinates are exclusive.
array = RealizedArray(
    "promoter-001",
    "ACGTTGCAAGTCCTGATCGTACCG",
    (
        Placement(
            "occurrence-1", "site-a", "tfbs", parts[0].sequence, 0, Orientation.FORWARD
        ),
        Placement(
            "occurrence-2", "site-b", "tfbs", parts[1].sequence, 8, Orientation.FORWARD
        ),
    ),
    provenance={"family": "two-regulator"},
)
# Keep dataset metadata separate from individual array annotations.
collection = ArrayCollection(
    parts, (array,), provenance={"dataset": "binding-site-panel"}
)
# Publish into a new directory; failed validation leaves no completed collection.
receipt = da.export(collection, view="arrays", all=True, format="bundle", out="arrays")
# Verify every part, placement, record digest and database checksum.
summary = da.inspect("arrays", verify=True)
assert summary.arrays == 1 and summary.placements == 2
# Read matching occurrences with their final-sequence core coordinates.
selected = ArrayFilter(array_ids=("promoter-001",))
with da.inspect("arrays", view="placements", select=selected).records() as rows:
    assert [(row.start, row.end) for row in rows] == [(0, 16), (8, 24)]
# Export every matching placement as a scalar table for analysis.
da.export(
    "arrays",
    view="placements",
    select=selected,
    all=True,
    format="tsv",
    out="array-placements.tsv",
)
```

Each placement must agree with the final DNA and its source part in the declared
orientation. For a reverse placement, `Placement.sequence` contains the reverse
complement of `Part.sequence`. Core coordinates and strands are transformed from
the source part into the final sequence. Unspecified core annotations stay null.

Array IDs are unique within a collection. Equal DNA sequences can retain
different array IDs and annotations. Repeated occurrences of a part need distinct
placement IDs within their array. The part catalog retains unused parts too.

## Inspect, export and render

```bash
# Verify a saved collection and show its counts and provenance.
dense-arrays inspect arrays --verify --json
# Read a page of placements; --after continues from the returned cursor.
dense-arrays inspect arrays --view placements --array-id promoter-001 --json
# Export DNA with preserved array IDs in FASTA headers.
dense-arrays export arrays --view sequences --all --format fasta --out arrays.fasta
# Render exactly one array; this requires the playback extra.
dense-arrays render arrays --view array --array-id promoter-001 --out promoter.png
```

| View | Contents | Export formats |
| --- | --- | --- |
| `summary` | Counts, provenance, exporter and collection digest | JSON |
| `arrays` | Full realized arrays, including placement metadata | JSON, JSONL, bundle |
| `parts` | Complete source catalog, including unused parts | JSON |
| `sequences` | Array ID, sequence ID, DNA, length and GC fraction | JSON, CSV, TSV, FASTA |
| `placements` | Occurrences, source part IDs, groups, spans and core strands | JSON, CSV, TSV |

`ArrayFilter` accepts array IDs, part IDs and groups. Values within one field
mean OR; fields combine with AND. Part and group filters select whole arrays
containing a matching occurrence. Their placement view includes every occurrence
in those arrays. The equivalent CLI flags are `--array-id`, `--part-id` and
`--group`. For `parts`, filters apply directly to the complete catalog.

Record export requires `all=True` or `--all`. Inspection defaults to a page of
100 records. Full verification and export use explicit [read limits](resources.md);
large collections may need a higher `--max-read-records` allowance.

## Carry provenance with the data

Collection `provenance` stores supplied dataset metadata. Use array `provenance`
and part or placement `metadata` for annotations at those levels. Values must
be finite JSON data. Record original software, solver settings, seeds and source
checksums when those are known; keep missing values null or omit them. Represent
integers beyond the target reader's exact range as documented decimal strings.

The manifest's `exporter` records the package and version that wrote the
collection. Source metadata survives re-export unchanged. Verification checks
the supplied geometry and file integrity; it does not establish a solver outcome
or reproduce optimization. Execution views such as attempts, plans and generation
quality require a saved run.

To exchange supplied inputs through the CLI, first export `arrays` as `jsonl`.
The stream contains a versioned catalog header, one realized array per line and
a completion record with a count and checksum. Import it with
`dense-arrays export input.jsonl --view arrays --all --format bundle --out arrays`.
Truncated inputs fail before publication. Use JSONL or a bundle to retain full
metadata; scalar CSV/TSV and FASTA are analysis projections.
