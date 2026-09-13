---
title: Constraints
description: Set positional requirements, regulator coverage, and side preferences before solving.
---

# Constrain an array

Require motifs at particular positions, require coverage of motif groups, or
favor one side of the array. Positional and group constraints are hard
requirements; side biases are preferences.

Run each independent example with `uv run python` from the
[installed checkout](quickstart.md#install-from-source). Configure a fresh
`Optimizer` before solving: adding constraints after the model is built raises
`RuntimeError`. Use `optimal()`, `solutions()`, or `solutions_diverse()`;
`approximate()` rejects configured constraints and side biases.

## Position two motifs

Use `add_promoter_constraints()` to require one motif upstream of another.
With the four motifs from the first-array tutorial, require the first motif
at start 0, 1, or 2 and the fourth motif seven to nine bases after its end:

```python
from dense_arrays import Optimizer

motifs = [
    "ACGTTGCAAGTCCTGA",
    "AGTCCTGATCGTACCG",
    "TCGTACCGATGCTTAG",
    "ATGCTTAGGACGTTCA",
]
optimizer = Optimizer(motifs, sequence_length=40, strands="single")
optimizer.add_promoter_constraints(
    upstream=motifs[0],
    downstream=motifs[3],
    upstream_pos=(0, 2),
    spacer_length=(7, 9),
)
best = optimizer.optimal()
print(best)
upstream_start = best.offsets_fwd[0]
downstream_start = best.offsets_fwd[3]
assert upstream_start in range(0, 3)
assert downstream_start is not None
assert 7 <= downstream_start - (upstream_start + len(motifs[0])) <= 9
```

Positions are zero-based start coordinates. A two-element range includes its
endpoints; a single integer fixes the value. `spacer_length` measures bases
between the end of the upstream motif and the start of the downstream motif.
`downstream_pos` can constrain the downstream start separately.
Use an integer, `None`, or an ordered two-item tuple whose bounds are integers
or `None`. An omitted bound is unrestricted. Position bounds must be
non-negative; negative spacers are allowed to request overlap. Booleans,
fractional bounds, lists, and reversed ranges are rejected before model allocation.

Both motifs must occur in the supplied library. Reusing a motif in another
pair requires another copy of that motif in the library. Positional
requirements are enforced by the solver and can make a request infeasible.

## Require regulator coverage

Map each motif entry to a regulator label to specify which groups must appear.
This example requires both entries labeled `R1` and at least one entry from
another group. Reducing the length limit to 35 bases leaves room for three of
the four motifs:

```python
from dense_arrays import Optimizer

motifs = [
    "ACGTTGCAAGTCCTGA",
    "AGTCCTGATCGTACCG",
    "TCGTACCGATGCTTAG",
    "ATGCTTAGGACGTTCA",
]
optimizer = Optimizer(motifs, sequence_length=35, strands="single")
optimizer.add_regulator_constraints(
    ["R1", "R1", "R2", "R3"],
    required={"R1"},
    min_count_by_regulator={"R1": 2},
    min_required_regulators=2,
)
best = optimizer.optimal()
print(best)
assert best.offsets_fwd[0] is not None and best.offsets_fwd[1] is not None
assert any(offset is not None for offset in best.offsets_fwd[2:])
assert best.nb_motifs == 3
```

Use `required` for named labels, `min_required_regulators` for a minimum number
of different labels, and `min_count_by_regulator` for minimum counts of entries
assigned to a label. Declare the regulator requirements together in one call.
Supply positive integers for minimum counts and nonempty regulator labels
without surrounding whitespace. Fractional and boolean counts are rejected;
the library must contain enough entries to meet every declared minimum.

The example selects the first three entries, representing `R1` and `R2`.
Minimum counts refer to entries in the supplied mapping; the minimum number
of regulators refers to distinct labels. The length limit determines how many
additional entries can fit.

## Prefer a side

Side biases favor left or right positions when other scoring terms are equal.
Use them for preferences; use positional constraints when a coordinate or span
is required.

```python
from dense_arrays import Optimizer

left_motif = "ACGTTGCAAGTCCTGA"
right_motif = "ATGCTTAGGACGTTCA"
optimizer = Optimizer([left_motif, right_motif], sequence_length=32, strands="single")
optimizer.add_side_biases(left=[left_motif], right=[right_motif])
best = optimizer.optimal()
print(best.sequence)
assert best.sequence.startswith(left_motif)
assert best.sequence.endswith(right_motif)
assert best.offsets_fwd == [0, len(best.sequence) - len(right_motif)]
```

All preferred motifs must belong to the original library. Side biases can be
combined with positional and regulator constraints before solving. A rejected
bias update leaves both side preferences unchanged.
The [optimizer reference](api.md#optimization) documents the accepted arguments
and validation failures.
