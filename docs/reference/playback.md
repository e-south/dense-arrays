---
title: Playback API
description: Reconstruct, serialize, and render validated saved feature placements.
---

# Playback API

The [playback guide](../playback.md) builds a complete example.
`reconstruct_playback()` derives an order from persisted coordinates and
returns a `PlaybackPlan`. The NetworkX/Matplotlib media pipeline consumes that
plan without running or importing a solver. Rendering requires the playback
extra; contracts, reconstruction, and serialization remain importable without
the solver or raster libraries.

## Interpretation and validation

Playback v1 supports `placement_reconstructed` authority and
`coordinate_precedence` relations. `solver_selected` is reserved and rejected;
there is no exact-trace contract in v1. Equal starts or containment produce
`ambiguous` order; internal gaps produce `layout_only` order.

Python construction and JSON loading share validation. They enforce supported
field types and enum values, sequence alignment, bounds, unique placement and
constraint IDs, known references, and contiguous step indices. Steps must use
the deterministic order `(start, length, placement_id)`, with each predecessor
identifying the preceding step. Reveal spans must exactly describe maximal
runs of newly covered bases, without repeated or unrelated positions.

Constraint results must match the actual placement coordinates and declared
range. A valid failed requirement stays representable as `passed=False`;
a contradictory flag or distance is rejected. JSON loaders also reject missing,
unknown, and duplicate object keys at every record level.

`source_digest` and `realization_digest` are checked as SHA-256 identifiers.
A loaded plan is not cryptographic proof of its source: v1 does not contain
all source fields required to recompute the realization digest. Producers
remain responsible for preserving source bytes and verifying their digests.

## Notices and presentation

`reconstruct_playback(realized, notices=())` accepts explicit caller-authored
`PlaybackNotice` records. The package does not infer how coordinates were
recovered from metadata names. Reserved notice codes cannot contradict the
plan's authority or order; arbitrary caller prose remains caller evidence.
Notice levels are `info` or `warning`.

Renderers derive the essential authority, order, and failure text from validated
plan fields. Optional detailed notices, labels, colors, graph choices, and
media inputs are described in the [presentation reference](playback-presentation.md).
`PlaybackDocument` is defined in `dense_arrays.playback.presentation` and
exported from `dense_arrays.playback`.

## Signatures

::: dense_arrays.playback
    options:
      show_root_heading: false
      members:
        - PlaybackPlan
        - PlaybackStep
        - CoordinateSpan
        - PlaybackAuthority
        - OrderingStatus
        - ConstraintResult
        - PlaybackNotice
        - NoticeLevel
        - reconstruct_playback
        - dumps_realized_array
        - loads_realized_array
        - realized_array_to_dict
        - realized_array_from_dict
        - dumps_playback_plan
        - loads_playback_plan
        - playback_plan_to_dict
        - playback_plan_from_dict

Return to the [API index](../api.md).
