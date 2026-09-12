# Develop Dense Arrays

Start from the [source checkout](quickstart.md#install-from-source). Read the
[architecture map](architecture/README.md) before changing module boundaries,
and the [playback contract](architecture/solution-playback.md) before changing
serialized placements, reconstruction, or rendering.

## Local verification

Run the full gate from the repository root before handoff:

```bash
uv sync --frozen --extra dev --extra playback --extra docs
uv run pre-commit run --all-files
uv run ruff check .
uv run ruff format --check .
uv run pytest -q
uv run mkdocs build --strict
uv export --frozen --all-extras --no-hashes --no-emit-project | uv run pip-audit -r /dev/stdin --progress-spinner off
uv build
```

To run checks when committing, install the hooks with
`uv run pre-commit install`. The optional `dev`, `playback`, and `docs` extras
supply the tools needed by the full gate.

## Documentation changes

Keep [the README](https://github.com/e-south/dense-arrays/blob/main/README.md)
short enough to choose a task. Use [the documentation index](index.md) to route
readers to a guide, and keep signature details in [the API reference](api.md).
Run changed examples in the locked environment. Examples that add constraints
must construct a fresh optimizer before solving.

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
