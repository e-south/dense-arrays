---
title: Input and record compatibility
description: Integrate Dense Arrays with explicit input validation, solver outcomes and saved-record contracts.
author: Eric J. South
---

# Input and record compatibility

Pin the Dense Arrays version used by your analysis. Check the inputs, exceptions
and record schemas your code consumes before changing that version. Package
versions and persisted schema versions are separate identifiers.

## Optimization inputs and outcomes

Pass actual integer lengths, counts, indices and interval bounds. Booleans and
fractional values are rejected. Positional intervals are an integer, `None`, or
an ordered two-item tuple. Position bounds are non-negative; negative spacers
request intentional overlap. Motifs contain nonempty uppercase `A/C/G/T` DNA.

The optimizer snapshots its input library. Returned configuration views are
independent copies and results are immutable. Create a fresh optimizer to change
the packing problem. Rejected bias, weight and forbid operations leave state
unchanged.

Catch `InfeasibleError` when the offered model has no feasible result.
`OptimizationError` subclasses distinguish backend failures, unproven optimality
and invalid solver output. Both names are exported from `dense_arrays`.
Enumeration propagates execution failures, including failures after an earlier
result. An exception does not mean normal enumeration exhaustion; see the
[solver outcome table](reference/optimizer.md#solver-outcomes).

`approximate()` rejects constraints, side biases and modified models. Use exact
methods for those requirements. Counts refer to selected entries; incidental
contained substrings add no entries. Each supplied entry can be selected in
only one orientation.

## Persisted records

Use nonblank string IDs, supported enum values and integer coordinates.
`RealizedArray` validates sequence alignment and references during construction.
Handle malformed input at that boundary, before reconstructing playback.

Metadata and provenance contain JSON values with string keys and finite numbers.
Nested values are immutable snapshots. Use public `dumps_*()` helpers for JSON
text or `*_to_dict()` for JSON-compatible dictionaries. A shallow
`dict(record.provenance)` does not recursively convert nested snapshots.

Saved playback plans must match coordinate ordering, predecessor references,
newly covered spans and constraint evaluations. Their authority is
`placement_reconstructed`; they do not claim solver-recorded chronology.
A legitimate `passed=False` requirement remains valid evidence. Keep that failure
when saving or reconstructing a record.

Coordinate conversion belongs to the producer of the record. If a conversion
requires qualification, supply explicit `PlaybackNotice` records with `notices=`.
Metadata labels alone do not establish how a record was recovered. Producers
also verify source bytes and digests before making provenance claims.

For generated-library compatibility, see the
[artifact schema contracts](architecture/library-workflow/artifacts.md#schema-compatibility).
Unknown schemas and policies fail explicitly rather than being reinterpreted.

## Presentation and export

Import `PlaybackDocument` from `dense_arrays.playback` or
`dense_arrays.playback.presentation`. Use public render functions; private
renderer helpers are not compatibility interfaces.

Supply placement-ID label and color maps and explicit legend entries. Available
profiles are `categorical`, `uniform` and `constraints`.
`graph_detail="none"` requires `graph_fraction=0`; `"reduced"` draws the traversal
layout. See [presentation settings](reference/playback-presentation.md).
Callback frames must satisfy the requested-step, shape and dtype contracts every
time they are evaluated, including during lazy rendering.

The playback CLI requires at least one of `--poster`, `--mp4` or `--gif`.
Overwriting an existing media file requires `--replace`. Outputs are rendered
before publication and published atomically per file. Check exit status and
stderr for any partial publication caused by a filesystem failure. An existing
file or earlier terminal output alone does not establish success.

Workflow exports use create-only destinations. Their receipts bind source
identity and selected membership; the [export guide](library-workflow/results/export.md)
describes checksums, coordinate joins and streamed-output failures.

## Check an integration

1. Record the package version and the schemas your application accepts.
2. Check public imports and representative saved records against that version.
3. Exercise a successful case and relevant invalid-input, solver and export failures.
4. Inspect the final sequence, placements or media together with their recorded evidence.

Keep interpretation of Dense Arrays records with the receiving application.
