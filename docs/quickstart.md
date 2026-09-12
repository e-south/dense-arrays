# Create your first array

Dense Arrays takes non-empty uppercase `A/C/G/T` motifs and a positive length
limit. It searches for an arrangement that packs as many supplied motifs as
possible into that space. The result may include only a subset of the library.

## Install from source

Use Python 3.12 or later. With [uv](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/e-south/dense-arrays.git
cd dense-arrays
uv sync --frozen
```

Run the commands below from that checkout. To install into an existing Python
environment instead, run `python -m pip install .` from the repository root.
The core installation includes the optimizer and HTML playback. Extra
dependencies are described in [playback](playback.md) and
[development](development.md).

## Solve from the terminal

```bash
uv run dense-arrays optimize \
  --motif CAG --motif AGC --motif CGT --length 6 --strands single
```

The synthetic motifs `CAG`, `AGC`, and `CGT` fit into `CAGCGT`, with starts at
0, 1, and 3. Its overlapping bases count once toward the six-base limit.
The terminal displays the sequence, its complement, and the placed motifs.
`--strands single` permits motifs only in their supplied orientation;
`--strands double` also permits reverse complements and is the default.

For a larger input, replace the repeated `--motif` options with
`--motifs-file motifs.txt`. Write one motif per line; blank lines and lines
starting with `#` are ignored. The two input forms cannot be combined.

Use `uv run dense-arrays optimize --help` to see all options. CBC is the
default OR-Tools backend. `--solver` selects another backend only if the local
OR-Tools installation can create it; an unsupported backend fails explicitly.

## Solve from Python

Run this code with `uv run python` in the checkout environment:

```python
from dense_arrays import Optimizer

optimizer = Optimizer(["CAG", "AGC", "CGT"], sequence_length=6, strands="single")
best = optimizer.optimal()
print(best.sequence)       # CAGCGT
print(best.nb_motifs)      # 3
print(best.offsets_fwd)    # [0, 1, 3]
```

`optimal()` returns a `DenseArray`. `sequence_length` is the requested limit;
`len(best.sequence)` is the realized length and may be shorter. Terminal
padding uses `-` for unused trailing space; those characters are not part of
`best.sequence`. An offset of `None` means that motif was not placed on that
strand. See the [result reference](api.md#dense-array-results) for the full
interface.

Configure every [constraint](constraints.md) before solving. Solving builds
the model; to change constraints afterward, create a new `Optimizer`.
If no motif can fit or the declared constraints are infeasible, `optimal()`
raises `ValueError`. The CLI reports an unsuccessful solve with a nonzero exit.

## Request further solutions

Limit enumeration to the number of results you intend to inspect:

```bash
uv run dense-arrays solutions \
  --motif CAG --motif AGC --motif CGT --length 6 --strands single \
  --max-solutions 3 --diverse
```

The equivalent bounded Python iteration is:

```python
from itertools import islice

from dense_arrays import Optimizer

optimizer = Optimizer(["CAG", "AGC", "CGT"], sequence_length=6, strands="single")
for solution in islice(optimizer.solutions_diverse(), 3):
    print(solution.sequence, solution.nb_motifs)
```

`solutions()` enumerates arrangements in decreasing score order.
`solutions_diverse()` adjusts motif weights to favor underrepresented motif
entries across returned arrangements. This is a motif-representation bias;
different arrangements can produce the same DNA sequence. The sequence pool
is not guaranteed to be unique or balanced across regulators.

Next, [add constraints](constraints.md) or [look up the optimizer API](api.md#optimization).
