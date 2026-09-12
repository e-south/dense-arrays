---
title: DenseArray results
description: Interpret sequence length, strand offsets, motif counts, and terminal output.
---

# DenseArray results

`Optimizer` returns an immutable `DenseArray`. Read `sequence` for the realized
DNA and the two offset lists for the supplied motif entries. The library and
offset properties return independent lists; editing them does not change the result.

| Field or property | Meaning |
| --- | --- |
| `library` | Original motif entries, in input order |
| `sequence_length` | Requested length limit |
| `sequence` | Realized sequence; it can be shorter than the limit |
| `offsets_fwd` | Zero-based starts of forward placements, or `None` |
| `offsets_rev` | Zero-based starts of reverse-complement placements, or `None` |
| `nb_motifs` | Number of placed entries across both strands |
| `compression_ratio` | Sum of lengths of selected library entries divided by the requested length limit |

`str(result)` displays both strands and pads unused trailing space with `-`.
Padding is not part of `sequence`. The CLI shows the realized length followed
by the requested limit: `length 3 / 6 limit`. In Python, read `len(result.sequence)`
and `result.sequence_length`, respectively.

When constructing a result directly, the offset lists must match the library
length. Motifs must be nonempty uppercase `A/C/G/T` strings, and offsets must
be non-negative integers or `None`. Booleans and fractional offsets are
rejected. Each library entry can use at most one orientation. Placed motifs must agree where they overlap, stay within the length
limit, and cover a contiguous sequence from coordinate zero. The constructor
rejects missing placements, conflicting bases, gaps, and out-of-bounds offsets.
Direct construction can describe compatible contained placements; exact and
greedy packing select entries under the [path-entry rules](../method.md).

The result has no automatic conversion to playback JSON. Callers that need
playback translate their selected placements into a
[RealizedArray](realized.md).

## Signatures

Direct construction uses `DenseArray(library, sequence_length, offsets_fwd,
offsets_rev)`. Most callers should read the result returned by `Optimizer`
rather than construct one.

::: dense_arrays.solution.DenseArray

Return to the [API index](../api.md).
