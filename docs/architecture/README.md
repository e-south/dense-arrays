# Playback architecture

The public playback system has two complementary documents:

- [Animation product specification](animation-product-spec.md) preserves the
  visual and publication intent that initiated the work.
- [Solution playback architecture](solution-playback.md) is the authoritative
  local contract for ownership, truth levels, schemas, and producer handoffs.

Implementation lives under `dense_arrays.realized` and
`dense_arrays.playback`. Renderers consume `PlaybackPlan`; they do not inspect
optimizer or OR-Tools state.
