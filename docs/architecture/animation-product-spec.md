# Dense-array playback product brief

## Purpose

Dense-array playback explains how an ordered set of overlapping sequence
features realizes a compact DNA sequence. It is a publication surface, not a
dashboard and not a visualization of optimizer search.

The animation has two synchronized representations:

- a compact explanation graph showing the ordered feature relations;
- a linear duplex showing each feature settle into its realized coordinates.

Every completed frame must also work as a legible still in a presentation.

## Authority

The serialized contracts and authority language in
[`solution-playback.md`](solution-playback.md) are normative.

Current historical DenseGen records produce
`authority=placement_reconstructed`. Their order is recovered from persisted
placements. Playback must not describe this as the solver search, candidate
graph, or exact solver-selected path.

A future solver-authoritative trace may use the same renderer once dense-arrays
captures the selected oriented path at solve time. Exact trace capture is not a
prerequisite for using existing realized arrays.

## Product hierarchy

The public package owns neutral contracts, validation, reconstruction,
serialization, layout, and reference renderers. Producer adapters and recipes
remain outside this repository.

Three presentation tiers exercise the same public surface:

1. **Generic packing** teaches that an overlap-efficient order incrementally
   creates a compact sequence.
2. **Generic constraints** teaches that fixed anchors and their required span
   remain invariant while other sites pack around them.
3. **Study application** uses study-owned identities and labels to explain a
   specific promoter architecture.

The first two recipes belong to DenseGen. Study selection, biological labels,
and interpretation belong to the owning study repository.

## Visual contract

- Keep the duplex visually dominant; the graph is an explanation scaffold.
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

Generated media remains producer-owned output and is not committed to this
repository. Renderers consume only validated playback plans and must not import
solver or OR-Tools internals.
