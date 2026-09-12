# Dense Arrays agent router

Dense Arrays is a public Python package for optimizing motif-dense DNA arrays
and rendering explicit realized-array playback records.

Start with:

- `README.md` for user-facing scope and commands;
- `docs/architecture/README.md` for module boundaries;
- `docs/architecture/solution-playback.md` for playback authority;
- `docs/development.md` for the required local verification gate;
- `pyproject.toml` for supported Python and extras.

Keep optimizer semantics, realized-array contracts, and playback presentation
separate. Playback may explain persisted placements; it must not invent a
solver-recorded order. Run the full local gate in `docs/development.md` before
handoff.
