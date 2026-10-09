---
title: Development checks
description: Find focused checks, run the full local gate, and preview documentation.
---

# Develop Dense Arrays

Start from the [source checkout](installation.md#install-from-source). Read the
[architecture map](architecture/README.md) before changing module boundaries,
and the [playback contract](architecture/solution-playback.md) before changing
serialized placements, reconstruction, or rendering.

Choose a focused test from the [task-to-file map](architecture/README.md#find-the-files-for-a-change).
For documentation, use the [writing and routing guide](development/documentation.md).

For optimizer changes, use the [model-building measurements](development/optimizer-performance.md)
and their constraint-parity checks before comparing end-to-end execution time.

## Local verification

Run these checks from the repository root. `uv sync` manages the checkout’s
`.venv` from `pyproject.toml` and `uv.lock`; `uv run` selects that environment:

```bash
uv sync --frozen --extra dev --extra playback --extra docs  # Reproduce the locked tool environment.
uv run pre-commit run --all-files  # Check tracked source, secrets and formatting.
uv run ruff check .  # Check Python correctness and style rules.
uv run ruff format --check .  # Verify formatting without rewriting files.
uv run pytest -q  # Exercise behavior, failure paths and guide examples.
uv run mkdocs build --strict  # Reject broken documentation configuration.

# Audit the exact locked versions and hashes across optional features.
uv export --frozen --all-extras --no-emit-project | uv run pip-audit -r /dev/stdin --require-hashes --disable-pip --progress-spinner off
uv build  # Build the wheel and source distribution.
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

`tests/test_documentation.py` executes the Python examples in the first-array,
constraints, playback and saved-library guides, including their linked export,
assembly and quality examples. Its page list names the additional workflow
guides covered by the test. It also builds the site strictly and checks links,
assets and anchors. Keep examples on each page runnable in order; test outputs
use temporary directories. Run changed preparation examples separately, with
FIMO when scoring is requested. The linked playback presentation example also
needs a manual run when changed.

Preview documentation with `uv run mkdocs serve`; `uv run mkdocs build --strict`
writes the static site to `public/`. Shared documentation images live in
`docs/assets/`, so both GitHub and the built site use one source asset. Inspect
SVGs at their intended display size and retain accessible descriptions.

## Hosted checks and publication

Use the [release procedure](development/releases.md) for immutable software tags,
qualified distributions and PyPI publishing. Documentation hosting is separate.

The [GitHub workflow](https://github.com/e-south/dense-arrays/blob/main/.github/workflows/ci.yml)
checks the lock, source quality, tests, documentation, dependencies, and package
build on pushes to `main` and pull requests.
The [GitLab workflow](https://github.com/e-south/dense-arrays/blob/main/.gitlab-ci.yml)
publishes GitLab Pages from its default branch after the test jobs pass. Its
slim Python image omits Git, so the pre-commit job installs Git before checking
the repository.

After merging on GitHub, fast-forward the GitLab default branch to the same
commit through the existing repository remote. Verify the test and Pages jobs,
then open a changed guide on the hosted site and check its content and assets.
Building documentation locally or passing GitHub checks does not publish that
site.
