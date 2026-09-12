---
title: Code map
description: Find the module, contract, and tests that own a change.
---

# Code map

Dense Arrays separates optimization, persisted placement contracts, and
presentation. Find the owner of the behavior before editing:

| Surface | Owns |
| --- | --- |
| `problem`, `constraints` | Immutable packing inputs and requirement validation |
| `model`, `optimizer`, `errors` | Model construction, solving, enumeration, and distinct solver outcomes |
| `greedy` | Greedy selection with explicit entry occurrences and orientations |
| `solution` | Immutable `DenseArray` results and terminal representation |
| `realized`, `_record_validation` | Producer-neutral placements, JSON snapshots, and sequence alignment |
| `playback.models`, `playback.validation` | Supported plan fields and cross-record geometry/evidence checks |
| `playback.reconstruction`, `playback.serialization` | Deterministic compilation and strict JSON boundaries |
| `playback.presentation`, renderers | Document choices, visible evidence, NetworkX graph layout, stills, and video |

Optimization returns a `DenseArray`; it does not automatically serialize a
`RealizedArray`. Producer adapters own that translation. Study identities,
biological labels, selected records, and publication interpretation remain
caller-owned.

## Find the files for a change

Paths below are relative to `src/dense_arrays/` and `tests/`, respectively.
Public entrypoints stay small; helper modules own the listed decisions.

| Task | Implementation under `src/dense_arrays/` | Tests under `tests/` |
| --- | --- | --- |
| Packing inputs, counts, ranges, or mutation checks | `problem.py`, `constraints.py`, `optimizer.py` | `test_core_contracts.py`, `test_optimize.py` |
| Exact model construction | `model.py` | `test_optimize.py`, `test_core_contracts.py` |
| Solver status or enumeration failure | `optimizer.py`, `errors.py` | `test_solver_outcomes.py`, `test_cli.py` |
| Greedy selected entries and occurrences | `greedy.py` | `test_greedy.py` |
| Sequence overlap calculations | `sequence.py` | `test_optimize.py` |
| Result construction or terminal layout | `solution.py` | `test_core_contracts.py`, `test_optimize.py`, `test_cli.py` |
| Optimizer command inputs or errors | `cli.py` | `test_cli.py` |
| Playback CLI preflight and publication | `playback/cli.py`, `playback/output.py` | `test_playback_cli.py` |
| Realized fields, nested provenance, or placement alignment | `realized.py`, `_record_validation.py` | `test_playback_contracts.py` |
| Plan JSON and evidence/geometry validation | `playback/models.py`, `playback/validation.py`, `playback/serialization.py` | `test_playback_contracts.py` |
| Coordinate reconstruction and caller notices | `playback/reconstruction.py` | `test_playback_contracts.py`, `test_playback.py` |
| Document labels, colors, and visible evidence | `playback/presentation.py`, `playback/theme.py` | `test_playback_presentation.py` |
| Graph projection, selected relations, layout, routing | `playback/graph/`, `playback/graph_drawing.py` | `test_playback_graph.py`, plus rendered stills |
| Raster scene and sequence frames | `playback/scene_drawing.py`, `playback/duplex_drawing.py`, `playback/duplex_frames.py` | `test_playback.py`, `test_playback_presentation.py` |
| Frame timing, writers, and figure cleanup | `playback/frame_schedule.py`, `playback/export.py`, `playback/matplotlib_renderer.py` | `test_playback_exports.py` |
| Dependency import boundaries | `__init__.py`, `playback/__init__.py`, `playback/graph/__init__.py` | `test_optional_playback_imports.py` |
| Runnable documentation and built links | `docs/`, `mkdocs.yml` | `test_documentation.py` |

Use `rg -n 'def |class ' <file>` to locate symbols. Run focused checks from the
repository root with
`uv run --no-sync pytest -q tests/<file>`, then finish with the
[full development gate](../development.md#local-verification).

## Change boundaries

Input validation is independent of solver execution. The optimizer holds an
immutable problem and delegates model construction and greedy realization.
A `DenseArray` can describe a valid contiguous layout more broadly than an
exact path; `forbid()` checks the narrower exact-path requirement separately.

Plan validation and evidence projection are independent of drawing. Graph
projection is independent of raster geometry. Default and injected graph
layout engines receive the same selected topology. Raster scene drawing,
producer-frame adaptation, timing, and writer lifecycle have separate owners;
public export functions compose them.

Keep these boundaries when adding behavior. A useful module owns a decision
that can change without editing unrelated solver or renderer code.

## Continue by task

- [Playback contract](solution-playback.md): coordinates, authority, and producer handoffs.
- [Product brief](animation-product-spec.md): visual and publication goals that need output review.
- [Playback guide](../playback.md): runnable entrypoints.
- [API reference](../api.md): inputs, results, and failures.
- [Caller migrations](../migration.md): deliberate input and integration changes.
- [Development](../development.md): verification; the [audit](../development/audit.md)
  and [improvement plan](../development/improvement-plan.md) retain implementation history and status.
