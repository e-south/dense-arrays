# Dense Arrays agent router

Dense Arrays is a public Python package for optimizing motif-dense DNA arrays
and rendering explicit realized-array playback records.

Start with:

- `README.md` for user-facing scope and commands;
- `docs/architecture/README.md` for module boundaries;
- `docs/architecture/solution-playback.md` for playback authority;
- `pyproject.toml` for supported Python, extras, and quality commands.

Keep optimizer semantics, realized-array contracts, and playback presentation
separate. Playback may explain persisted placements; it must not invent a
solver-recorded order. Run the full local gate before handoff:

```bash
uv sync --frozen --extra dev --extra playback --extra docs
uv run pre-commit run --all-files
uv run pytest -q
uv run mkdocs build --strict
uv export --frozen --all-extras --no-hashes --no-emit-project | uv run pip-audit -r /dev/stdin --progress-spinner off
uv build
```
