---
title: Update an existing caller
description: Adapt integrations to strict inputs, distinct solver outcomes, and explicit playback presentation.
author: Eric J. South
---

# Update an existing caller

Before upgrading, check the interfaces your integration uses: optimization,
saved placements, or media export. Stricter validation can reject inputs and
settings that older versions accepted. The sections below identify what to
change and which failures to handle.

## Optimization callers

Pass actual integer lengths, counts, indices, and interval bounds; booleans
and fractional values are rejected. Positional intervals are an integer,
`None`, or an ordered two-item tuple. Position bounds are non-negative;
negative spacers remain available for intentional overlap. Use nonempty
uppercase `A/C/G/T` motifs and valid regulator labels.

The optimizer snapshots its input library. Its returned library and model
configuration views are independent copies, and results are immutable.
Create a fresh optimizer to change the packing problem. Rejected bias, weight,
and forbid operations leave state unchanged.

Catch `InfeasibleError` for absence of a feasible result. Catch
`OptimizationError` for backend execution, unproven optimality, or malformed
solver-output failures; its subclasses distinguish these cases. Both names are
exported from `dense_arrays`. Iteration now propagates execution failures,
including a failure after an earlier result. Do not treat every failure as an
empty iterator. See the [outcome table](reference/optimizer.md#solver-outcomes).

`approximate()` rejects constraints, side biases, and model weight/forbid
changes. Use exact methods for those requirements. Counts follow selected
entries and occurrences; incidental contained substrings do not add entries.
A library entry may be selected in only one orientation. Review code that
compared heuristic substring counts with exact counts.

## Persisted records and adapters

Use nonblank string IDs and supplied labels, supported enum values, and integer
coordinates. Numeric IDs and fractional or boolean coordinates are rejected.
`RealizedArray` validates sequence alignment and references during construction;
move error handling to that boundary instead of waiting for reconstruction.

Metadata and provenance must contain JSON values with string keys and finite
numbers. Nested values are immutable snapshots. Serialize through the public
`dumps_*()` helpers for JSON text or `*_to_dict()` for ordinary JSON-compatible
dictionaries. Do not assume `dict(record.provenance)` recursively thaws the snapshot.

Saved v1 plans must match coordinate ordering, predecessor references, exact
newly covered spans, and actual constraint evaluations. `solver_selected` is
reserved and rejected. A valid `passed=False` constraint remains acceptable;
do not replace it with an inaccurate success flag to make a plan load.

Adapters own coordinate conversion and evidence for recovery procedures.
Metadata such as `offset_raw` does not cause the compiler to assert that
recovery occurred. Pass explicit `PlaybackNotice` records with `notices=` when
that qualification is justified. Producers also retain responsibility for
verifying source bytes and digests.

## Presentation

Import `PlaybackDocument` from `dense_arrays.playback` or
`dense_arrays.playback.presentation`. Use public render functions; private
renderer drawing/export helpers are not compatibility interfaces.

Supply placement-ID label/color maps and caller-authored legend entries.
Generic profiles are `categorical`, `uniform`, and `constraints`; study-specific
profiles such as `secg` are rejected. `graph_detail="none"` requires
`graph_fraction=0`. Replace `graph_detail="inset"` with `"reduced"` and keep
the same `graph_fraction` to preserve the traversal-only layout. Supply
`PlaybackPresentation.legend_entries` directly; there is no profile-derived
legend helper. Review
[media presentation settings](reference/playback-presentation.md) and validate
producer raster frames against their requested-step, shape, and dtype rules.
Titles, subtitles, and full evidence are stored in
[native media metadata](reference/playback-presentation.md#read-the-evidence).
A producer callback can be evaluated lazily, so it must
return valid frames throughout the requested render.

## Export handling

The CLI requires at least one of `--poster`, `--mp4`, or `--gif`, and
`--replace` to overwrite existing files. It rejects aliases
and collisions, renders all requested outputs before publication, and publishes
atomically per file. Check its exit code and stderr, including reported partial
publication on a filesystem failure. Do not infer successful completion from
an existing artifact or earlier terminal output.

## External consumer checklist

Verify the migration in each consumer's environment:

1. Pin the Dense Arrays version or commit being adopted.
2. Check public imports, adapters, saved fixtures, palettes, and callback frames
   against that version. If an old plan fails validation, reconstruct it from
   the pinned realized placements, preserving legitimate failed requirements.
3. Run one successful case and relevant invalid-input and export-failure cases.
4. Inspect the exported media and its evidence before publishing recipe outputs.

Consumer owners perform these checks in their repositories. Updating Dense
Arrays alone does not migrate external adapters, notebooks, recipes, or artifacts.
