---
title: Playback contract
description: Ownership, coordinate interpretation, implemented checks, and remaining enforcement gaps.
---

# Playback contract

Playback reconstructs an explanation from persisted placements. This page
defines the ownership and interpretation rules; the
[API reference](../reference/playback.md) describes callable interfaces.

## Implemented behavior and remaining gaps

Reconstruction checks placement bounds, sequence agreement, and constraint
references, then derives an order and evaluates declared distances. Invalid
placements fail; a valid layout that violates a distance requirement returns
a failed constraint result. Ordering and provenance qualifications are recorded
separately in `notices`.

Saved plans do not yet receive equivalent semantic validation. Renderers also
do not display every authority, ordering, or constraint qualification. The
[audit](../development/audit.md) records these gaps. The rules below state the
required interpretation; they are not evidence that every entrypoint enforces
it. Exact solver trace capture is not implemented.

## Ownership

`dense-arrays` owns renderer-independent realized-array and playback-plan
contracts, deterministic reconstruction, validation, and reference playback.
Producer packages own translation from their persisted schemas. A study owns
only selected record identities, domain labels, captions, and review evidence.

BaseRender may provide sequence-frame or video-publication integration. It does
not own graph semantics, solver claims, or the playback clock.

## What the order means

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

`PlaybackPlan` contains frozen semantic records, newly revealed sequence spans,
constraint evaluations, ordering status, authority, and notices. Renderers
consume this plan. They must not depend on optimizer or OR-Tools state.
The package root currently imports the optimizer eagerly; isolating those
imports is part of the hardening plan.

For placement reconstruction, the deterministic order is:

1. Start coordinate.
2. Shorter placement first for equal starts.
3. Stable placement ID.

Record construction checks local field invariants and unique placement and
constraint IDs. Reconstruction validates sequence agreement, bounds, and
constraint references, and evaluates distance ranges before producing a plan.

## Producer handoffs

Producer adapters and publication recipes remain outside this package. DenseGen
is one caller; installing Dense Arrays does not require DenseGen or a study
repository. The following describes that adapter's ownership, not a command
implemented by Dense Arrays.

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

## Data flow and future trace work

Existing DenseGen corpora compile through:

```text
DenseGen record -> RealizedArray -> reconstructed PlaybackPlan -> renderers
```

A future exact-trace design could use:

```text
SolveResult -> ExactSolutionTrace -> solver-selected PlaybackPlan -> renderers
```

The intended reuse point is the renderer's plan input. Any exact-trace extension
needs its own versioned validation and migration decision; it is outside the
current hardening plan.
