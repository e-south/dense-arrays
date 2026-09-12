---
title: DenseArray results
description: Interpret sequence length, strand offsets, motif counts, and terminal output.
---

# DenseArray results

`Optimizer` returns a `DenseArray`. Read `sequence` for the realized DNA and
the two offset lists for the supplied motif entries.

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
Padding is not part of `sequence`. The CLI's `length` field is the requested
limit; use `len(result.sequence)` for the realized length.

When constructing a result directly, the offset lists must match the library
length. Placed motifs must agree where they overlap, stay within the length
limit, and cover a contiguous sequence from coordinate zero. The constructor
rejects missing placements, conflicting bases, gaps, and out-of-bounds offsets.

The result has no automatic conversion to playback JSON. Callers that need
playback translate their selected placements into a
[RealizedArray](realized.md).

## Signatures

Direct construction uses `DenseArray(library, sequence_length, offsets_fwd,
offsets_rev)`. Most callers should read the result returned by `Optimizer`
rather than construct one.

::: dense_arrays.solution

Return to the [API index](../api.md).
