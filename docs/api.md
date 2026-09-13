---
title: Python API
description: Choose the interface for solving, reading a result, or rendering saved placements.
---

# Python API

Choose the interface needed for your task. For a runnable example, start with
the [first-array tutorial](quickstart.md) or [playback guide](playback.md).

## Optimization

[Optimizer](reference/optimizer.md): construct a problem, add requirements,
solve with CBC, enumerate arrangements, and handle distinct failure outcomes.

## Dense-array results

[DenseArray](reference/results.md): read the sequence, offsets, motif count,
and compression ratio. `Optimizer` and `DenseArray` are exported by
`dense_arrays`.

## Sequence utilities

[Sequence utilities](reference/sequence.md): complements, pairwise overlaps,
and sequence-display helpers.

## Realized arrays and playback

- [Realized arrays](reference/realized.md): describe an existing sequence and
  its persisted placements.
- [Playback](reference/playback.md): reconstruct, serialize, and render those
  placements through `dense_arrays.playback` and its rendering module.

For command-line use, see [CLI options and failures](reference/cli.md).
Existing integrators should review [caller migrations](migration.md).
