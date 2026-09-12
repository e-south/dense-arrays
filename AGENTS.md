---
id: dense-arrays-router
intent: Route tasks to the relevant documentation and code owners.
---

# Dense Arrays agent router

Dense Arrays is a public Python package for optimizing motif-dense DNA arrays
and rendering explicit realized-array playback records.

Choose the route for the task; do not load every reference:

| Task | Start here |
| --- | --- |
| Install or create an array | `README.md` → `docs/quickstart.md` |
| Add positional or motif-group requirements | `docs/constraints.md` |
| Render saved placements | `docs/playback.md` |
| Change implementation or tests | `docs/architecture/README.md` task-to-file map |
| Change placement JSON or interpretation | `docs/architecture/solution-playback.md` → serialization row in the code map |
| Revise documentation | `docs/development/documentation.md` |
| Plan hardening work | `docs/development/improvement-plan.md` and its linked findings |

`pyproject.toml` owns supported Python and extras. Nested source and test
directories inherit this router unless a closer instruction file applies.

Keep optimizer semantics, realized-array contracts, and playback presentation
separate. Playback may explain persisted placements; it must not invent a
solver-recorded order. Run the full local gate in `docs/development.md` before
handoff. Preserve all existing author credits, including Virgile Andreani's.
Attribute new work to Eric J. South; do not replace joint authorship with a
single-author header. Keep generated dogfood media outside tracked source.
