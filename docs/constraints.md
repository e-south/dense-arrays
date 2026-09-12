# Constrain an array

Configure constraints on a fresh `Optimizer`, then solve. Each example below
is independent and can be run with `uv run python`. Once a model has been
built, adding constraints raises `RuntimeError`; create a new optimizer to
try another specification.

## Position two motifs

`add_promoter_constraints()` requires two library entries in a declared
upstream/downstream relationship. Despite the method name, the positional
contract can be illustrated with synthetic motifs:

```python
from dense_arrays import Optimizer

optimizer = Optimizer(
    ["GCA", "CCC", "ATGC", "CATT"], sequence_length=10, strands="single"
)
optimizer.add_promoter_constraints(
    upstream="ATGC",
    downstream="CCC",
    upstream_pos=(0, 2),
    spacer_length=(0, 3),
)
best = optimizer.optimal()
print(best)
```

Positions are zero-based start coordinates. A two-element range includes its
endpoints; a single integer fixes the value. `spacer_length` measures bases
between the end of the upstream motif and the start of the downstream motif.
`downstream_pos` can constrain the downstream start separately.

Both motifs must occur in the supplied library. Reusing a motif in another
pair requires another copy of that motif in the library. Positional
requirements are enforced by the solver and can make a request infeasible.
Using biological motif sequences does not establish that the resulting
arrangement functions as a promoter.

## Require regulator coverage

Map every motif entry to a regulator label. Use `required` for named labels,
`min_required_regulators` for a minimum number of different labels, and
`min_count_by_regulator` for minimum counts of motif entries assigned to a label.

```python
from dense_arrays import Optimizer

optimizer = Optimizer(["AAA", "CCC", "GGG", "TTT"], sequence_length=9, strands="single")
optimizer.add_regulator_constraints(
    ["R1", "R1", "R2", "R3"],
    required={"R1"},
    min_count_by_regulator={"R1": 2},
    min_required_regulators=2,
)
best = optimizer.optimal()
print(best)
```

Here both motif entries labeled `R1` must appear, together with at least one
other regulator label. These are hard coverage requirements over the supplied
mapping. They do not measure binding or simultaneous occupancy.
Declare the regulator requirements together in one call.

## Prefer a side

Side biases favor left or right positions when other scoring terms are equal.
Use them for preferences; use positional constraints when a coordinate or span
is required.

```python
from dense_arrays import Optimizer

optimizer = Optimizer(["AAA", "CCC"], sequence_length=6, strands="single")
optimizer.add_side_biases(left=["AAA"], right=["CCC"])
best = optimizer.optimal()
print(best.sequence)  # AAACCC
```

All preferred motifs must belong to the original library. Side biases can be
combined with positional and regulator constraints before solving.
The [optimizer reference](api.md#optimization) documents the accepted arguments
and validation failures.
