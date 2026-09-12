---
title: Python API
description: Choose the interface for solving, reading a result, or rendering saved placements.
---

# Python API

Start with the interface needed for your task. Each reference combines current
behavior and generated signatures; the [first-array tutorial](quickstart.md)
and [playback guide](playback.md) provide runnable examples.

## Optimization

[Optimizer](reference/optimizer.md): construct a problem, add requirements,
solve with CBC, or enumerate arrangements.

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
  placements. These interfaces live in explicit submodules.

For command-line use, see [CLI options and failures](reference/cli.md).
