---
title: Playback API
description: Reconstruct, serialize, and render saved feature placements.
---

# Playback API

The [playback guide](../playback.md) builds a complete example.
`reconstruct_playback()` derives an order from persisted coordinates and
returns a `PlaybackPlan`. Renderers consume that plan without solving again.

## Interpretation and validation

Reconstruction reports `placement_reconstructed` authority. Check
`ordering_status`, `constraint_results`, and `notices` when interpreting the
plan. Equal starts or containment can make the order ambiguous; internal gaps
produce `layout_only` status. Neither case supplies a recorded solver path.

JSON loaders reject unsupported schema versions and check keys on the main
records, but nested reveal spans currently accept extra fields. The loaders
currently coerce some field values and do not enforce every cross-field
invariant of a loaded plan. Prefer reconstructing from validated realized
placements when producing new plans. Loading a saved plan is not equivalent
to checking the underlying realization again.

The current HTML and default raster views do not display all plan qualifications.
Read the serialized status and constraint results, and supply a caption that
states any ambiguity or failed requirement. The
[hardening plan](../development/improvement-plan.md) tracks validation and
visible qualification work.

## Signatures

::: dense_arrays.playback
    options:
      members:
        - PlaybackPlan
        - PlaybackStep
        - PlaybackAuthority
        - OrderingStatus
        - ConstraintResult
        - PlaybackDocument
        - reconstruct_playback
        - dumps_realized_array
        - loads_realized_array
        - dumps_playback_plan
        - loads_playback_plan
        - render_playback_html
        - render_playback_collection_html

Return to the [API index](../api.md).
