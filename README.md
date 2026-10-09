# ![Dense Arrays — overlapping motifs within a sequence-length limit](https://raw.githubusercontent.com/e-south/dense-arrays/v0.2.1/docs/assets/dense-arrays-banner.png)

[![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue)](https://github.com/e-south/dense-arrays/blob/main/pyproject.toml)
[![CI](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-read-blue)](https://dunloplab.gitlab.io/dense-arrays)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](https://github.com/e-south/dense-arrays/blob/main/LICENSE)

Dense Arrays designs nucleotide sequences with densely packed DNA-protein binding
sites. It uses optimal string packing to arrange overlapping sites on both DNA
strands, with controls for sequence length, site placement, and library diversity.

## First array

Use Python 3.12 or later. This example packs four synthetic 16-base motifs into a
40-base array.

```bash
python -m venv .venv                       # Create an isolated Python environment.
source .venv/bin/activate                  # Activate it on macOS or Linux.
python -m pip install dense-arrays         # Install the released package and CBC solver.

# Fit four sites into at most 40 bases, considering both DNA strands.
dense-arrays optimize \
  --motif ACGTTGCAAGTCCTGA \
  --motif AGTCCTGATCGTACCG \
  --motif TCGTACCGATGCTTAG \
  --motif ATGCTTAGGACGTTCA \
  --length 40 --strands double
```

One optimum contains all four sites in 40 bases. The
[first-array guide](https://github.com/e-south/dense-arrays/blob/main/docs/quickstart.md)
explains the overlaps and shows the matching Python code.
[Installation](https://github.com/e-south/dense-arrays/blob/main/docs/installation.md)
covers optional features, version availability, Windows, uv projects, and Pixi.

## Choose a task

The [documentation index](https://github.com/e-south/dense-arrays/blob/main/docs/index.md)
routes by use case: prepare binding-site pools, constrain and generate libraries,
inspect shortfalls, compare results, export sequences, or render figures and playback.
It also links the API, method, and contributor guides.

## Cite and contribute

Andreani V, South EJ, Dunlop MJ (2024). Generating information-dense promoter
sequences with optimal string packing. *PLOS Computational Biology* 20(7):
e1012276. [DOI: 10.1371/journal.pcbi.1012276](https://doi.org/10.1371/journal.pcbi.1012276).

[Report an issue](https://github.com/e-south/dense-arrays/issues) ·
[Contribute](https://github.com/e-south/dense-arrays/blob/main/docs/development.md) ·
[Security](https://github.com/e-south/dense-arrays/blob/main/SECURITY.md) ·
[MIT license](https://github.com/e-south/dense-arrays/blob/main/LICENSE)
