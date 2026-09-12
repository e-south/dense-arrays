# ![Dense Arrays — overlapping motifs within a sequence-length limit](https://raw.githubusercontent.com/e-south/dense-arrays/main/docs/assets/dense-arrays-banner.svg)

[![CI](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml/badge.svg)](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-gitlab_pages-blue)](https://dunloplab.gitlab.io/dense-arrays)

When many DNA motifs must fit into a short sequence, compatible overlaps can
save space. Dense Arrays searches for an arrangement of supplied motifs within
a length limit and returns the sequence and motif offsets. Choose single- or
double-strand placement, require motif groups, or constrain positions.

The package also accepts persisted feature placements for playback. These
views explain a realized array; their ordering is reconstructed from its
coordinates. Motif packing and playback describe sequence arrangements, not
binding, expression, or laboratory performance.

## First array

With Python 3.12 or later and [uv](https://docs.astral.sh/uv/), run from the
repository root:

```bash
uv sync --frozen
uv run dense-arrays optimize \
  --motif CAG --motif AGC --motif CGT --length 6 --strands single
```

This synthetic example packs three motifs into `CAGCGT`. The
[quickstart](docs/quickstart.md) covers installation, Python use, and bounded
enumeration of further solutions.

## Choose a task

- [Create and read your first array](docs/quickstart.md).
- [Set positional and regulator constraints](docs/constraints.md).
- [Render saved feature placements](docs/playback.md).
- [Look up Python interfaces](docs/api.md).

For the formulation and associated paper, see [the method](docs/method.md).
The [documentation index](docs/index.md) also routes integrators and
contributors to the relevant contracts.

## Contribute

[Development](docs/development.md) covers local checks and documentation builds.
[Report bugs](https://github.com/e-south/dense-arrays/issues), follow the
[security policy](SECURITY.md) for vulnerabilities, and read
[AGENTS.md](AGENTS.md) when working with a coding agent.
Dense Arrays is available under the [MIT license](LICENSE).
