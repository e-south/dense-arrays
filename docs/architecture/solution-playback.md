---
title: Playback contract
description: Ownership, coordinate interpretation, validated v1 semantics, and producer handoffs.
---

# Playback contract

Playback reconstructs an explanation from persisted placements. Dense Arrays
validates the layout, derives its coordinate order, and evaluates declared
distances. A renderer presents those facts together with caller-owned labels
and evidence. Exact solver trace capture is not implemented.

## Ownership

`dense-arrays` owns renderer-independent realized-array and playback-plan
contracts, deterministic reconstruction, validation, and reference renderers.
Callers supply native `RealizedArray` records. A study owns
selected record identities, domain labels, captions, and review evidence.

Optional frame callbacks provide sequence illustrations. Dense Arrays retains
graph semantics, solver evidence and the playback clock.

## What the order means

Version 1 accepts only `placement_reconstructed` authority and
`coordinate_precedence` relations. The order is derived from placement
coordinates, with these qualifications:

- `unique`: the persisted intervals imply one strict left-to-right order.
- `ambiguous`: equal starts or containment require a deterministic tie-break.
- `layout_only`: an internal uncovered span prevents a complete placement chain.

The deterministic order is start coordinate, shorter placement first for an
equal start, then stable placement ID. Each step names the preceding step as
its predecessor. This reference describes coordinate ordering; it does not
claim a recorded solver edge. A `layout_only` view suppresses active traversal.

The reserved `solver_selected` enum value is unsupported in v1. Constructing
or loading a v1 plan with that authority is rejected.

## Record and plan validation

`RealizedArray` contains an IUPAC sequence, stable feature and placement
identities, oriented feature sequences, zero-based half-open coordinates,
declared constraints, and immutable JSON provenance. Construction validates
sequence agreement, bounds, unique IDs, and constraint references.
See the [realized-array reference](../reference/realized.md) for field domains.

`PlaybackPlan` contains frozen steps, newly revealed spans, constraint results,
ordering status, authority, and notices. Python constructors and JSON loaders
share validation. Plans must match actual layout order, exact newly covered
bases, predecessor references, and distance evaluations. JSON record keys are
checked at every level; numeric and textual fields are not coerced into other
types. See [plan validation](../reference/playback.md#interpretation-and-validation).

The distance between declared placements is
`downstream.start - upstream.end`. A layout that violates its declared range
remains valid data and returns `passed=False`. A result whose stored distance
or `passed` flag contradicts that layout is malformed and rejected.

Digests are SHA-256 identifiers with validated syntax. The source holder must
verify source bytes separately: a playback plan cannot authenticate external
source content or recompute every realization field from v1 alone.

## Evidence and presentation

Renderers consume validated plans without importing optimizer or OR-Tools
state. A `PlaybackDocument` resolves artifact metadata, labels, colors, and
compact visible evidence. NetworkX supplies the established graph layout;
Matplotlib draws the scene and writes PNG, MP4, or GIF through the same media
pipeline. A compact summary preserves ambiguous/gapped ordering qualifications
and failed requirements. Routine reconstructed authority and unique coordinate
order are metadata only. Full evidence is stored in native media metadata;
notices are an optional addition. Presentation settings and evidence retrieval
belong in the [presentation reference](../reference/playback-presentation.md).

Every animated scene begins with its complete configured graph and duplex in
neutral gray. Placement order controls progressive color and active emphasis,
not whether known sequence coordinates exist on the canvas. Uncovered bases
remain gray, and `layout_only` still suppresses active traversal. Geometry and
font sizes remain fixed across resting, intermediate, and completed states.
Scene changes preserve visible context instead of fading through white.

A producer frame callback receives `step_index: int | None`. `None` explicitly
means the complete pre-placement scene; integers identify the state after a
placement. Both use the same fixed image shape and two-frame validation/cache
contract. The producer owns its neutral annotations and fixed crop; Dense
Arrays blends successive images in place without moving nucleotide glyphs.

Callers can supply `PlaybackNotice` records through
`reconstruct_playback(realized, notices=...)`. Dense Arrays preserves explicit
caller evidence and rejects conflicting reserved authority/order codes. It
does not infer biological identity or coordinate-recovery methods from labels,
IDs, or metadata keys. Caller-authored prose remains the caller's evidence.

## Producer handoffs

Callers provide `RealizedArray` values with feature coordinates and fixed-element
relationships in the declared realized-sequence coordinate system. A recipe
selects labels, colors, captions and outputs. Other coordinate fields remain
metadata; labels and metadata do not change the coordinates Dense Arrays reads.

Endpoint recipes pin source-table digests and record IDs. They own source
verification, record selection, biological interpretation, and publication
bundles. A recipe can publish normalized input and plan JSON, MP4, a
poster, and a manifest of input/output digests and versions. The package CLI
renders requested files; it does not create that manifest or publish a site.

Use the [realized-array reference](../reference/realized.md) for input fields
and the [playback reference](../reference/playback.md) for construction and
validation failures.

## Trace boundary

Playback consumes this data flow:

```text
Producer record -> RealizedArray -> reconstructed PlaybackPlan -> renderer
```

An exact solver trace would require a separate versioned validation contract:

```text
Solver result -> exact solution trace -> solver-selected plan -> renderer
```

That extension is outside v1. The reusable renderer input remains an explicit,
validated plan.

## Future trace work

See the [trace boundary](#trace-boundary).
