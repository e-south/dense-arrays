---
title: First array
description: Pack four binding-site-sized sequences with CBC and read their positions and overlaps.
---

# Create your first array

Pack four 16-base motifs into a 40-base sequence, then read where each motif
starts. The compatible overlaps let the motifs share bases. For your own
library, the length limit may allow only a subset of the motifs to be placed.

## Set up an environment

Follow [installation](installation.md) to install the package in a virtual
environment. The base install includes CBC; FIMO is not needed for this example.
The sequences below are synthetic 16-base sites, chosen to illustrate overlaps.

## Solve from the terminal

Supply non-empty uppercase `A/C/G/T` motifs and a positive integer length limit:

```bash
# Search both strands for an arrangement within 40 bases.
dense-arrays optimize \
  --motif ACGTTGCAAGTCCTGA \
  --motif AGTCCTGATCGTACCG \
  --motif TCGTACCGATGCTTAG \
  --motif ATGCTTAGGACGTTCA \
  --length 40 --strands double
```

The terminal displays the sequence, its complement, and the placed motifs.
One optimal arrangement is shown below; its reverse complement is equally valid:

```text
ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGTTCA
```

In this orientation, the four input motifs start at 0, 8, 16, and 24. Each
consecutive pair shares eight bases, so the first motif contributes 16 bases and each later motif adds
eight: `16 + 8 + 8 + 8 = 40`. Read the [packing method](method.md) for the
overlap calculation and its graph interpretation.

`--strands double` searches both supplied motifs and their reverse complements
and is the default. Use `--strands single` to restrict placements to the supplied
orientation.

For a larger input, replace the repeated `--motif` options with
`--motifs-file motifs.txt`. Write one motif per line; blank lines and lines
starting with `#` are ignored. The two input forms cannot be combined.

Use `dense-arrays optimize --help` to see all options. CBC is the
default OR-Tools backend. `--solver` selects another backend only if the local
OR-Tools installation can create it. An unavailable backend produces an error
and a nonzero exit. See [CLI failures](reference/cli.md).

## Solve from Python

Run this code with `python` in your activated environment, or `uv run python`
in a uv project:

```python
from dense_arrays import Optimizer  # Configure and solve the packing problem.

# Four synthetic 16-base sites; adjacent sites share eight bases.
motifs = [
    "ACGTTGCAAGTCCTGA",
    "AGTCCTGATCGTACCG",
    "TCGTACCGATGCTTAG",
    "ATGCTTAGGACGTTCA",
]
optimizer = Optimizer(
    motifs, sequence_length=40, strands="double"
)  # Allow both strands.
best = optimizer.optimal()  # Require a proven optimum for this packing model.
print(best.sequence)  # The final DNA, without terminal display padding.
print(best.nb_motifs)  # Four selected input entries in this example.
print(best.offsets_fwd)  # Forward-strand starts in input order; None means absent.
print(best.offsets_rev)  # Reverse-complement starts on the same final sequence.
assert len(best.sequence) == 40  # This arrangement fills the length limit.
assert best.nb_motifs == 4  # Every supplied site was selected.
```

Both offset lists use zero-based starts in input order. For the forward
arrangement above, `offsets_fwd` is `[0, 8, 16, 24]` and each occupied span is
16 bases long. For its reverse complement, `offsets_rev` is `[24, 16, 8, 0]`.
Reverse offsets locate the reverse-complement motif on the returned sequence.
The right endpoint of each span is excluded.

`optimal()` returns a `DenseArray`. `sequence_length` is the requested limit;
`len(best.sequence)` is the realized length and may be shorter. Terminal
padding uses `-` for unused trailing space; those characters are not part of
`best.sequence`. An offset of `None` means that motif was not placed on that
strand. See the [result reference](api.md#dense-array-results) for the full
interface.

Configure every [constraint](constraints.md) before solving. Solving builds
the model; to change constraints afterward, create a new `Optimizer`.
If no motif can fit or the declared constraints are infeasible, `optimal()`
raises `InfeasibleError`, a `ValueError` subclass. Backend failure, an unproven
feasible result, and invalid solver output have distinct exceptions; see
[solver outcomes](reference/optimizer.md#solver-outcomes).

## Request further solutions

Limit enumeration to the number of results you intend to inspect. This bounds
the result count, not the time needed to solve each result:

```bash
# Print at most three arrangements, favoring sites used less often so far.
dense-arrays solutions \
  --motif ACGTTGCAAGTCCTGA \
  --motif AGTCCTGATCGTACCG \
  --motif TCGTACCGATGCTTAG \
  --motif ATGCTTAGGACGTTCA \
  --length 40 --strands double \
  --max-solutions 3 --diverse
```

Check the command's exit status: a solver failure exits nonzero even if an
earlier result was printed.

In Python, reuse the `motifs` library defined above with a fresh optimizer:

```python
from itertools import islice  # Bound how many results are consumed.

# Use the typed requests and operations needed by this example.
from dense_arrays import Optimizer

optimizer = Optimizer(
    motifs, sequence_length=40, strands="double"
)  # Allow both strands.
for solution in islice(optimizer.solutions_diverse(), 3):  # Stop after three results.
    print(solution.sequence, solution.nb_motifs)  # Report DNA and selected-entry count.
```

`solutions()` enumerates arrangements in decreasing score order.
`solutions_diverse()` adjusts motif weights to favor underrepresented motif
entries across returned arrangements. This is a motif-representation bias;
different arrangements can produce the same DNA sequence. The sequence pool
is not guaranteed to be unique or balanced across regulators.

Next, [add constraints](constraints.md) or [look up the optimizer API](api.md#optimization).
