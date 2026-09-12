# ![Dense Arrays — overlapping motifs within a sequence-length limit](https://raw.githubusercontent.com/e-south/dense-arrays/main/docs/assets/dense-arrays-banner.svg)

[![CI](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml/badge.svg)](https://github.com/e-south/dense-arrays/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-gitlab_pages-blue)](https://dunloplab.gitlab.io/dense-arrays)

Dense Arrays packs supplied DNA motifs into a sequence within a requested
length. Overlapping motifs share compatible bases. Choose single- or
double-strand placement, require particular motif groups or positional
relationships, and inspect the resulting sequence and motif offsets.

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
[quickstart](https://github.com/e-south/dense-arrays/blob/main/docs/quickstart.md) covers installation, Python use, and bounded
enumeration of further solutions.

## Choose a task

- [Create and read your first array](https://github.com/e-south/dense-arrays/blob/main/docs/quickstart.md).
- [Set positional and regulator constraints](https://github.com/e-south/dense-arrays/blob/main/docs/constraints.md).
- [Render saved feature placements](https://github.com/e-south/dense-arrays/blob/main/docs/playback.md).
- [Look up Python interfaces](https://github.com/e-south/dense-arrays/blob/main/docs/api.md).

For the formulation and associated paper, see [the method](https://github.com/e-south/dense-arrays/blob/main/docs/method.md).
The [documentation index](https://github.com/e-south/dense-arrays/blob/main/docs/index.md) also routes integrators and
contributors to the relevant contracts.

## Contribute

[Development](https://github.com/e-south/dense-arrays/blob/main/docs/development.md) covers local checks and documentation builds.
[Report bugs](https://github.com/e-south/dense-arrays/issues), follow the
[security policy](https://github.com/e-south/dense-arrays/blob/main/SECURITY.md) for vulnerabilities, and read
[AGENTS.md](https://github.com/e-south/dense-arrays/blob/main/AGENTS.md) when working with a coding agent.
Dense Arrays is available under the [MIT license](https://github.com/e-south/dense-arrays/blob/main/LICENSE).
