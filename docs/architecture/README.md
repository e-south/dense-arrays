---
title: Code map
description: Find the module, contract, and tests that own a change.
---

# Code map

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

## Find the files for a change

Paths below are relative to `src/dense_arrays/` and `tests/`, respectively.
The table identifies current owners; the
[improvement plan](../development/improvement-plan.md) describes proposed
splits that have not been implemented.

| Task | Implementation under `src/dense_arrays/` | Tests under `tests/` |
| --- | --- | --- |
| Motif validation or exact solving | `optimizer.py`, `constraints.py` | `test_optimize.py` |
| Sequence overlap calculations | `sequence.py` | `test_optimize.py` |
| Result construction or terminal layout | `solution.py` | `test_optimize.py`, `test_cli.py` |
| Optimizer CLI parsing or messages | `cli.py` | `test_cli.py` |
| Reject playback inputs before writing output | `playback/cli.py`, `playback/serialization.py` | `test_playback.py` covers loading; playback CLI coverage still needs to be added |
| Placement or plan JSON | `realized.py`, `playback/models.py`, `playback/serialization.py` | `test_playback.py`: round trips and invalid inputs |
| Coordinate reconstruction | `playback/reconstruction.py` | `test_playback.py`: sequence, ambiguity, gaps, constraints |
| HTML controls or labels | `playback/html.py` | `test_playback.py`, plus a browser check |
| Graph layout or media export | `playback/graph/`, `playback/matplotlib_renderer.py` | `test_playback.py`, plus rendered stills |
| Dependency imports | `__init__.py`, `playback/__init__.py`, `playback/graph/__init__.py` | `test_optional_playback_imports.py` |

Use `rg -n 'def |class ' <file>` to locate symbols. Run a focused test file with
`uv run pytest -q tests/<file>`, then finish with the
[full development gate](../development.md#local-verification).

## Change boundaries

Keep input validation independent of solver execution. Keep placements and
plan semantics independent of drawing. A useful module owns one decision that
can change without editing unrelated code; file length alone does not establish
that boundary.

The optimizer currently combines model construction, solution extraction,
enumeration, and a greedy heuristic. The raster renderer combines scene drawing
and export loops. Both need bounded refactoring after the
[documented contract gaps](../development/audit.md) have regression tests.

## Playback contracts

Use the contract for implementation and the brief for presentation intent:

- [Animation product specification](animation-product-spec.md) preserves the
  visual and publication goals, including requirements not yet fully enforced.
- [Solution playback architecture](solution-playback.md) is the authoritative
  local contract for ownership, interpretation, schemas, and producer handoffs;
  it distinguishes current checks from intended invariants.

Implementation lives under `dense_arrays.realized` and
`dense_arrays.playback`. Renderers consume `PlaybackPlan`; they do not inspect
optimizer or OR-Tools state. Package initialization currently imports the
optimizer eagerly; import isolation remains a hardening task.

For runnable entrypoints, use [the playback guide](../playback.md). Interface
details belong in [the API reference](../api.md); contributor checks belong in
[development](../development.md).
