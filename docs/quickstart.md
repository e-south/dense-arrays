---
title: First array
description: Install from source, solve a small CBC example, and read its sequence and offsets.
---

# Create your first array

Pack four 16-base motifs into a 40-base sequence, then read where each motif
starts. The compatible overlaps let the motifs share bases. For your own
library, the length limit may allow only a subset of the motifs to be placed.

## Install from source

Use Python 3.12 or later. With [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/e-south/dense-arrays.git
cd dense-arrays
uv sync --frozen
```

Run the commands below from that checkout. To install into an existing Python
environment instead, run `python -m pip install .` from the repository root.
The core installation includes the optimizer and persisted-placement
contracts. Rendering requires the extra dependencies described in
[playback](playback.md); contributor tools are listed in [development](development.md).

## Solve from the terminal

Supply non-empty uppercase `A/C/G/T` motifs and a positive integer length limit:

```bash
uv run dense-arrays optimize \
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

Use `uv run dense-arrays optimize --help` to see all options. CBC is the
default OR-Tools backend. `--solver` selects another backend only if the local
OR-Tools installation can create it. An unavailable backend produces an error
and a nonzero exit. See [CLI failures](reference/cli.md).

## Solve from Python

Run this code with `uv run python` in the checkout environment:

```python
from dense_arrays import Optimizer

motifs = [
    "ACGTTGCAAGTCCTGA",
    "AGTCCTGATCGTACCG",
    "TCGTACCGATGCTTAG",
    "ATGCTTAGGACGTTCA",
]
optimizer = Optimizer(motifs, sequence_length=40, strands="double")
best = optimizer.optimal()
print(best.sequence)
print(best.nb_motifs)  # 4
print(best.offsets_fwd)
print(best.offsets_rev)
assert len(best.sequence) == 40
assert best.nb_motifs == 4
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
uv run dense-arrays solutions \
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
from itertools import islice

from dense_arrays import Optimizer

optimizer = Optimizer(motifs, sequence_length=40, strands="double")
for solution in islice(optimizer.solutions_diverse(), 3):
    print(solution.sequence, solution.nb_motifs)
```

`solutions()` enumerates arrangements in decreasing score order.
`solutions_diverse()` adjusts motif weights to favor underrepresented motif
entries across returned arrangements. This is a motif-representation bias;
different arrangements can produce the same DNA sequence. The sequence pool
is not guaranteed to be unique or balanced across regulators.

Next, [add constraints](constraints.md) or [look up the optimizer API](api.md#optimization).
