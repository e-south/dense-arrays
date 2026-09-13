---
title: First array
description: Install from source, solve a small CBC example, and read its sequence and offsets.
---

# Create your first array

Pack four 16-base motifs into a 37-base sequence, then read where each motif
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
  --motif AAGTCCTGATCGTACC \
  --motif GATCGTACCGATGCTT \
  --motif CCGATGCTTAGGACGT \
  --length 37 --strands single
```

The terminal displays the sequence, its complement, and the placed motifs:

```text
ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGT
```

The four input motifs start at 0, 7, 14, and 21. Each consecutive pair shares
nine bases, so the first motif contributes 16 bases and each later motif adds
seven: `16 + 7 + 7 + 7 = 37`. Read the [packing method](method.md) for the
overlap calculation and its graph interpretation.

`--strands single` permits motifs only in their supplied orientation;
`--strands double` also permits reverse complements and is the default.

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
    "AAGTCCTGATCGTACC",
    "GATCGTACCGATGCTT",
    "CCGATGCTTAGGACGT",
]
optimizer = Optimizer(motifs, sequence_length=37, strands="single")
best = optimizer.optimal()
print(best.sequence)  # ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGT
print(best.nb_motifs)  # 4
print(best.offsets_fwd)  # [0, 7, 14, 21]
assert best.sequence == "ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGT"
assert best.nb_motifs == 4
assert best.offsets_fwd == [0, 7, 14, 21]
```

Offsets are zero-based starts in input order: the first entry occupies
`[0, 16)`, the second `[7, 23)`, the third `[14, 30)`, and the fourth `[21, 37)`.
The right endpoint is excluded. Use these positions to recover each motif
from the sequence.

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
  --motif AAGTCCTGATCGTACC \
  --motif GATCGTACCGATGCTT \
  --motif CCGATGCTTAGGACGT \
  --length 37 --strands single \
  --max-solutions 3 --diverse
```

Check the command's exit status: a solver failure exits nonzero even if an
earlier result was printed.

In Python, reuse the `motifs` library defined above with a fresh optimizer:

```python
from itertools import islice

from dense_arrays import Optimizer

optimizer = Optimizer(motifs, sequence_length=37, strands="single")
for solution in islice(optimizer.solutions_diverse(), 3):
    print(solution.sequence, solution.nb_motifs)
```

`solutions()` enumerates arrangements in decreasing score order.
`solutions_diverse()` adjusts motif weights to favor underrepresented motif
entries across returned arrangements. This is a motif-representation bias;
different arrangements can produce the same DNA sequence. The sequence pool
is not guaranteed to be unique or balanced across regulators.

Next, [add constraints](constraints.md) or [look up the optimizer API](api.md#optimization).
