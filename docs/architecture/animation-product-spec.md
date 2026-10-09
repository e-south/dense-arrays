---
title: Playback product brief
description: Visual, accessibility and publication requirements for figures from persisted placements.
---

# Playback product brief

Playback figures preserve placement evidence, ambiguous or gapped ordering,
and failed requirements. Review typography, contrast and publication context
alongside those data guarantees. The [playback contract](solution-playback.md)
defines validated semantics; the [development gate](../development.md) defines
repository review checks.

## Purpose

Dense-array playback explains how an ordered set of overlapping sequence
features realizes a compact DNA sequence. Its graph and duplex explain
persisted placement geometry.

The animation has two synchronized representations:

- a compact explanation graph showing the ordered feature relations;
- a linear duplex emphasizing each feature at its fixed realized coordinates.

Resting and completed frames must also work as legible stills in a presentation.

## Authority

The serialized contracts and authority language in
[`solution-playback.md`](solution-playback.md) are normative.

Version 1 plans use `authority=placement_reconstructed`. Their order is
derived from persisted placements. Playback must not describe this as the solver search, candidate
graph, or exact solver-selected path.

Exact solver-trace capture is unsupported. Realized arrays provide the sequence
and placement evidence needed for reconstructed playback.

## Product hierarchy

The public package owns neutral contracts, validation, reconstruction,
serialization, layout and reference renderers. Callers supply native realized
records; recipes own record selection and publication.

Three presentation tiers exercise the same public surface:

1. **Generic packing** teaches that an overlap-efficient order incrementally
   creates a compact sequence.
2. **Generic constraints** teaches that fixed anchors and their required span
   remain invariant while other sites pack around them.
3. **Study application** uses study-owned identities and labels to explain a
   specific promoter architecture.

Generic examples use supplied sequences and placements. Study applications add
selected records, biological labels and interpretation through caller-owned
recipes.

## Visual contract

- Keep the duplex visually dominant; the graph is an explanation scaffold.
- Begin with the complete graph, duplex, and placement context in neutral gray;
  progressively color the represented placements without moving nucleotide glyphs.
- Reveal one causal event at a time: traverse, place, settle, hold.
- Use one canonical curve for the visible route, progressive stroke, and point.
- Freeze layout for the full scene and keep Start and End as compact horizontal
  anchors.
- Do not add edge costs, categorical colors, legends, or annotations unless
  they serve the tier's premise.
- Use equal output-space nucleotide typography in the graph and duplex.
- Preserve fixed-element tracks and annotation lanes from the first frame so
  later annotations never cause existing features to move.
- Ship MP4 as the canonical motion artifact and a poster PNG as the still and
  reduced-motion alternative.

## Accessibility contract

- Essential text must retain at least 4.5:1 contrast against its fill.
- Essential distinctions must survive grayscale and cannot depend on hue
  alone.
- At intended slide placement, nucleotide cap height should be at least 18 px,
  Start/End labels at least 22 px, the active point at least 10 px, and the
  active route at least 3 px.
- Scene transitions must be restrained and the completed state must hold long
  enough to read.
- Presentation authors should attach concise alt text describing the ordered
  graph-to-duplex correspondence.

## Publication contract

Each endpoint publishes a digest-addressed bundle containing:

- normalized input and playback plan JSON;
- a manifest with schema, authority, ordering status, theme, renderer, and
  source provenance;
- MP4 playback;
- poster PNG.

Transient producer runs and dogfood media remain outside tracked source. The
maintained teaching assets `docs/assets/playback-example.mp4`,
`docs/assets/playback-opening.png`, and `docs/assets/playback-poster.png` are
reviewed documentation outputs and may be committed with their reproduction
instructions. Renderers consume validated playback plans without importing solver
or OR-Tools internals. Recipe owners provide the manifest and publication step;
the package CLI renders the requested files.
