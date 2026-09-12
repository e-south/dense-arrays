---
title: Realized arrays
description: Persisted sequence placements, validated coordinates, and immutable provenance.
---

# Realized arrays

Use `RealizedArray` to describe a sequence that already exists. The producing
package supplies placement identities, oriented feature sequences, and start
coordinates. Follow the [runnable playback example](../playback.md) first.

Coordinates are zero-based and half-open. A placement's end is its start plus
its sequence length. The sequence is already oriented to the realized array;
`orientation` records that choice and does not reverse-complement the input.
Unlike optimizer motifs, these contracts accept IUPAC DNA and normalize case.

| Field | Accepted values |
| --- | --- |
| `kind` | `tfbs`, `fixed_element`, `other` |
| `orientation` | `fwd`, `rev`, `unspecified` |

## Construction and loading

Python construction and JSON loading enforce the same record invariants:

- IDs and supplied labels are nonblank strings; numeric identities are not converted
  to text. Kinds and orientations must be supported enum values or their strings.
- Coordinates and distance bounds are integers, excluding booleans. Placement
  starts and declared distance bounds are non-negative; ranges are ordered.
- Placements fit within and agree exactly with the realized sequence.
  Placement and constraint IDs are unique; constraints reference two known,
  different placements.
- Metadata and provenance are detached immutable JSON snapshots, including
  nested objects and arrays. Keys are strings, numbers are finite, and cyclic
  or non-JSON values are rejected.

`RealizedArray` requires at least one placement. It can include uncovered
sequence and overlapping placements. `coordinate_space` is a nonblank caller
label; use `realized_sequence` for sequence-relative positions. The package
does not translate coordinate systems.

Declared distance constraints are evaluated during reconstruction as
`downstream.start - upstream.end`. A valid layout can violate a requirement:
reconstruction retains that result with `passed=False`. The actual distance
may be negative when placements overlap, even though the declared allowed
range is non-negative.

`source_digest`, when provided, is a SHA-256 hexadecimal identifier. Parsing
checks its form; it does not authenticate external source bytes.

## Signatures

::: dense_arrays.realized
    options:
      show_root_heading: false
      members:
        - RealizedArray
        - Placement
        - PlacementKind
        - Orientation
        - DeclaredConstraint

Continue to the [playback reference](playback.md) or
[ownership and interpretation rules](../architecture/solution-playback.md).
