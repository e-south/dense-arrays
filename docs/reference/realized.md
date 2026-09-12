---
title: Realized arrays
description: Persisted sequence placements, coordinates, and reconstruction checks.
---

# Realized arrays

Use `RealizedArray` to describe a sequence that already exists. The producing
package supplies placement identities, oriented feature sequences, and start
coordinates. Follow the [runnable playback example](../playback.md) first.

Coordinates are zero-based and half-open. A placement's end is its start plus
its sequence length. The sequence is already oriented to the realized array;
`orientation` records that choice and does not reverse-complement the input.
Unlike optimizer motifs, these contracts accept IUPAC DNA and normalize case.

Construction checks local field invariants, including non-negative integer
placement starts, unique placement IDs, and non-empty placements.
`reconstruct_playback()` checks agreement against the realized sequence,
placement bounds and constraint references. Keep `coordinate_space` set to
`realized_sequence`; arbitrary coordinate-space labels are not validated.
Call reconstruction before treating a constructed record as a validated layout.

Declared distance constraints are evaluated during reconstruction. A failed
distance requirement is retained as a failed result; it does not
prevent a plan from being returned. Inspect `plan.constraint_results` before
using the layout as evidence that a requirement passed.

## Signatures

::: dense_arrays.realized

Continue to the [playback reference](playback.md) or
[ownership and interpretation rules](../architecture/solution-playback.md).
