---
title: Development checks
description: Find focused checks, run the full local gate, and preview documentation.
---

# Develop Dense Arrays

Start from the [source checkout](quickstart.md#install-from-source). Read the
[architecture map](architecture/README.md) before changing module boundaries,
and the [playback contract](architecture/solution-playback.md) before changing
serialized placements, reconstruction, or rendering.

Choose a focused test from the [task-to-file map](architecture/README.md#find-the-files-for-a-change).
For documentation, use the [writing and routing guide](development/documentation.md).
The [audit](development/audit.md) records the original findings; the
[improvement plan](development/improvement-plan.md) tracks implementation and
verification status.

## Local verification

Run the full gate from the repository root before handoff:

```bash
uv sync --frozen --extra dev --extra playback --extra docs
uv run pre-commit run --all-files
uv run ruff check .
uv run ruff format --check .
uv run pytest -q
uv run mkdocs build --strict
uv export --frozen --all-extras --no-emit-project | uv run pip-audit -r /dev/stdin --require-hashes --disable-pip --progress-spinner off
uv build
```

To run checks when committing, install the hooks with
`uv run pre-commit install`. The optional `dev`, `playback`, and `docs` extras
supply the tools needed by the full gate.

The test suite executes the guides' Python examples and checks internal links
and anchors after a strict documentation build. The dependency audit uses the
complete locked export, including all extras and hashes, without resolving a
different environment.

## Documentation changes

Keep [the README](https://github.com/e-south/dense-arrays/blob/main/README.md)
short enough to choose a task. Use [the documentation index](index.md) to route
readers to a guide, and keep signature details in [the API reference](api.md).
Run changed examples in the locked environment. Examples that add constraints
must construct a fresh optimizer before solving.

`tests/test_documentation.py` executes the Python examples in `quickstart.md`,
`constraints.md`, and `playback.md`, then builds the site strictly and checks
built links, assets, and anchors. Keep examples on each page runnable in order;
temporary playback outputs are isolated by the test. The linked presentation
example continues from the playback guide's `plan` and output directory and
also needs a manual run when changed.

Preview documentation with `uv run mkdocs serve`; `uv run mkdocs build --strict`
writes the static site to `public/`. Shared documentation images live in
`docs/assets/`, so both GitHub and the built site use one source asset. Inspect
SVGs at their intended display size and retain accessible descriptions.

## Hosted checks and publication

The [GitHub workflow](https://github.com/e-south/dense-arrays/blob/main/.github/workflows/ci.yml)
checks the lock, source quality, tests, documentation, dependencies, and package
build on pushes to `main` and pull requests.
The [GitLab workflow](https://github.com/e-south/dense-arrays/blob/main/.gitlab-ci.yml)
also defines the GitLab Pages documentation build on its default branch.
Building documentation locally does not publish that site.
