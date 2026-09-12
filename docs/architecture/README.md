# Architecture

Dense Arrays separates optimization, persisted placement contracts, and
presentation. These surfaces have different inputs and responsibilities:

| Surface | Owns |
| --- | --- |
| `dense_arrays.optimizer` and `constraints` | Motif packing, positional requirements, regulator coverage, and solver interaction |
| `dense_arrays.solution` | The `DenseArray` result: sequence, motif offsets, and terminal representation |
| `dense_arrays.sequence` | Sequence and overlap utilities |
| `dense_arrays.realized` | Producer-neutral `RealizedArray` and placement contracts |
| `dense_arrays.playback.reconstruction` | Validation and deterministic reconstruction from saved placements |
| Playback renderers | HTML, still, and video views of a `PlaybackPlan` |

Optimization returns a `DenseArray`; it does not automatically serialize a
`RealizedArray`. Producer adapters own translation into the persisted-placement
contract. Study identities, biological labels, selected records, and publication
interpretation remain caller-owned.

## Playback contracts

The public playback system has two complementary documents:

- [Animation product specification](animation-product-spec.md) preserves the
  visual and publication intent that initiated the work.
- [Solution playback architecture](solution-playback.md) is the authoritative
  local contract for ownership, truth levels, schemas, and producer handoffs.

Implementation lives under `dense_arrays.realized` and
`dense_arrays.playback`. Renderers consume `PlaybackPlan`; they do not inspect
optimizer or OR-Tools state.

For runnable entrypoints, use [the playback guide](../playback.md). Interface
details belong in [the API reference](../api.md); contributor checks belong in
[development](../development.md).
