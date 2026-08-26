# dense-arrays

[![CI](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml/badge.svg)](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-gitlab_pages-blue)](https://dunloplab.gitlab.io/dense-arrays)

**dense-arrays** finds short double-stranded DNA sequences that densely pack a
requested set of protein-binding motifs. It also emits explicit realized-array
records and optional playback views that explain how the motifs overlap.

For more detailed documentation, please visit our [documentation site](https://dunloplab.gitlab.io/dense-arrays).

<img src="images/SPP_overview.png" alt="Dense Arrays Diagram" width="600"/>

*Formulation of the nucleotide String Packing Problem (SPP) as an Orienteering Problem (OP). For more details, see the [associated paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012276).*

---

### Installation

From the repo root:

```bash
uv sync --extra dev --extra playback --extra docs
```

This creates a local `.venv` and installs dev tools (pytest/ruff).

If you prefer pip:

```bash
pip install .
```

---

### Usage modes

Use **dense-arrays** either via the CLI for quick runs or the Python API for scripting.

CLI:

```bash
dense-arrays optimize --motifs-file motifs.txt --length 30 --strands double
dense-arrays solutions --motifs-file motifs.txt --length 30 --max-solutions 5 --diverse
```

Motifs file format: one motif per line (blank lines and `#` comments ignored). If no motif can fit within the requested length, the CLI exits non-zero with an explicit error; in Python, `Optimizer.optimal()` raises a `ValueError`.

Python API:

```python
import dense_arrays as da
import dense_arrays.sequence as seq

motifs = [
    "ATAATATTCTGAATT",
    "TCCCTATAAGAAAATTA",
    "TAATTGATTGATT",
    "GCTTAAAAAATGAAC",
    "TGCACTAAAATGGTGCAA",
]

opt = da.Optimizer(motifs, sequence_length=30)
best = opt.optimal()
# Best (highest score) solution.
print(best)

# Enumerate all solutions in decreasing score order.
for solution in opt.solutions():
    print(solution)

print("Shift metric:", seq.shift_metric("ATGCATTA", "CATTATG"))
```

Constraints (promoter/regulator/side bias) must be configured before calling `optimal()` or `solutions()`. To change constraints, create a new `Optimizer`.

Regulator constraints (optional, solver-level):

```python
regulators = ["R1", "R1", "R2", "R3", "R4"]
opt.add_regulator_constraints(regulators, min_required_regulators=2)
```

---

### Solver Backends

The methods `Optimizer.optimal` and `Optimizer.solutions` allow you to specify a solver backend. They accept any solver supported by `ortools`. The available options include:

- `"CBC"` (default)
- `"SCIP"`
- `"GUROBI"`
- `"CPLEX"`
- `"XPRESS"`
- `"GLPK"`

---

### Development

Install dev tools and enable pre-commit hooks:

```bash
uv sync --extra dev
uv run pre-commit install
```

Run the full hook suite locally:

```bash
uv run pre-commit run --all-files
```

### Solution playback

Persisted placements can be compiled into a renderer-independent playback plan
without claiming that their order is the solver-recorded path:

```python
from dense_arrays.playback import reconstruct_playback, render_playback_html
from dense_arrays.realized import RealizedArray

realized: RealizedArray = load_realized_array()
plan = reconstruct_playback(realized)
html = render_playback_html(plan, title="Dense-array packing")
```

The standalone renderer accepts strict `RealizedArray` or `PlaybackPlan` JSON:

```bash
dense-arrays-playback realized.json \
  --html playback.html \
  --poster poster.png \
  --mp4 playback.mp4
```

HTML output is self-contained. Poster and MP4 export require the optional
`playback` extra and a local FFmpeg installation. Reconstructed plans always
declare `placement_reconstructed` authority; `solver_selected` is reserved for
future exact traces captured by the optimizer.
