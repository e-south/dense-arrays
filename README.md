# ![Dense Arrays — overlapping motifs within a sequence-length limit](https://raw.githubusercontent.com/e-south/dense-arrays/v0.2.0/docs/assets/dense-arrays-banner.png)

[![CI](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml/badge.svg)](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-gitlab_pages-blue)](https://dunloplab.gitlab.io/dense-arrays)

Pack overlapping DNA motifs into a sequence-length limit. Dense Arrays selects
an arrangement and returns the sequence and each motif's position, so you can
inspect the overlaps or generate further solutions. Choose single- or
double-strand placement, require motif groups, or constrain positions.

The package also renders saved placements as images and video. [Watch the
worked playback](https://dunloplab.gitlab.io/dense-arrays/playback/#watch-four-overlapping-motifs)
to see how four motifs share sequence space.

## First array

With Python 3.12 or later and [uv](https://docs.astral.sh/uv/), run from a
[source checkout](https://github.com/e-south/dense-arrays/blob/main/docs/quickstart.md#install-from-source):

```bash
uv sync --frozen
uv run dense-arrays optimize \
  --motif ACGTTGCAAGTCCTGA \
  --motif AGTCCTGATCGTACCG \
  --motif TCGTACCGATGCTTAG \
  --motif ATGCTTAGGACGTTCA \
  --length 40 --strands double
```

The search considers both strands. One optimal arrangement packs all four
16-base motifs into 40 bases:

```text
ACGTTGCAAGTCCTGATCGTACCGATGCTTAGGACGTTCA
```

In this orientation, the motifs start at 0, 8, 16, and 24. The solver may return
the reverse-complement arrangement. The
[quickstart](https://github.com/e-south/dense-arrays/blob/main/docs/quickstart.md) explains these offsets and shows Python use
and bounded enumeration of further solutions.

## Documentation

- [Create and read your first array](https://github.com/e-south/dense-arrays/blob/main/docs/quickstart.md).
- [Set positional and regulator constraints](https://github.com/e-south/dense-arrays/blob/main/docs/constraints.md).
- [Render saved placements as images or video](https://github.com/e-south/dense-arrays/blob/main/docs/playback.md).
- [Understand the packing method](https://github.com/e-south/dense-arrays/blob/main/docs/method.md).
- [Look up Python interfaces](https://github.com/e-south/dense-arrays/blob/main/docs/api.md).
- [Update an existing caller](https://github.com/e-south/dense-arrays/blob/main/docs/migration.md).

The [documentation index](https://github.com/e-south/dense-arrays/blob/main/docs/index.md) also routes integrators and
contributors to the relevant interfaces and checks.

## Citation

If you use Dense Arrays in your research, please cite:

Andreani V, South EJ, Dunlop MJ (2024). Generating information-dense promoter
sequences with optimal string packing. *PLOS Computational Biology* 20(7):
e1012276. [doi:10.1371/journal.pcbi.1012276](https://doi.org/10.1371/journal.pcbi.1012276).

Record the package version or commit alongside your results.

## Contribute

[Development](https://github.com/e-south/dense-arrays/blob/main/docs/development.md) covers local checks and documentation builds.
[Report bugs](https://github.com/e-south/dense-arrays/issues) and follow the
[security policy](https://github.com/e-south/dense-arrays/blob/main/SECURITY.md) for vulnerabilities.
Dense Arrays is available under the [MIT license](https://github.com/e-south/dense-arrays/blob/main/LICENSE).
