# Dense-array solution playback

## Status

This document defines the first public playback contract. It supports persisted
realized placements now and leaves exact solver-selected traces as a compatible
future source.

## Ownership

`dense-arrays` owns renderer-independent realized-array and playback-plan
contracts, deterministic reconstruction, validation, and reference playback.
Producer packages own translation from their persisted schemas. A study owns
only selected record identities, domain labels, captions, and review evidence.

BaseRender may provide sequence-frame or video-publication integration. It does
not own graph semantics, solver claims, or the playback clock.

## Truth levels

`placement_reconstructed` means the order was derived from persisted placement
coordinates. Relations are coordinate precedence, not recorded solver-selected
edges. The ordering status further qualifies the result:

- `unique`: the persisted intervals imply one strict left-to-right order.
- `ambiguous`: equal starts or containment require a deterministic tie-break.
- `layout_only`: an internal uncovered span prevents a complete placement chain.

`solver_selected` is reserved for a future exact trace captured from the solver
result. A reconstructed plan must never use that authority value.

## Contracts

`RealizedArray` contains a realized sequence, stable feature and placement
identities, oriented sequences, zero-based half-open coordinates, declared
constraints, and source provenance.

`PlaybackPlan` contains immutable semantic steps, newly revealed sequence spans,
constraint evaluations, ordering status, authority, and notices. Renderers
consume only this plan and must not import optimizer or OR-Tools modules.

For placement reconstruction, the deterministic order is:

1. Start coordinate.
2. Shorter placement first for equal starts.
3. Stable placement ID.

The compiler validates sequence agreement, bounds, unique identities, constraint
references, and declared distance ranges before producing a plan.

## DenseGen integration

DenseGen translates `densegen__used_tfbs_detail` into `RealizedArray`. Display
coordinates use the persisted `offset`, while `offset_raw` and padding remain in
placement metadata. Fixed upstream/downstream pairs become declared distance
constraints. DenseGen labels remain producer metadata; endpoint recipes decide
whether those elements are neutral anchors or biological -35/-10 elements.

Endpoint recipes pin source-table digests and record IDs. They must not select
records randomly during publication. Generated publication bundles contain:

```text
manifest.json
playback.html
playback.mp4
poster.png
```

The manifest records source and realization digests, authority, ordering status,
label profile, renderer version, and output digests.

## Compatibility path

Existing DenseGen corpora compile through:

```text
DenseGen record -> RealizedArray -> reconstructed PlaybackPlan -> renderers
```

Future exact solves compile through:

```text
SolveResult -> ExactSolutionTrace -> solver-selected PlaybackPlan -> renderers
```

Both pathways share the same renderer contract. Historical arrays do not need to
be regenerated, and exact traces do not require a second animation system.
