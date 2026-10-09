---
id: dense-arrays-router
intent: Route tasks to the relevant documentation and code owners.
audience: runtime-agent
load: always
navigation:
  start:
    - docs/index.md
  implementation:
    - docs/architecture/README.md
  documentation:
    - docs/development/documentation.md
---

# Dense Arrays agent router

Dense Arrays is a public Python package for optimizing motif-dense DNA arrays
and preparing, generating, inspecting and exporting saved DNA libraries.

Choose the route for the task; do not load every reference:

| Task | Start here |
| --- | --- |
| Install or create an array | `README.md` → `docs/installation.md` → `docs/quickstart.md` |
| Generate or inspect a saved library | `docs/library-workflow.md` and its task guides |
| Add positional or motif-group requirements | `docs/constraints.md` |
| Render saved placements | `docs/playback.md` |
| Change library-workflow contracts | `docs/architecture/library-workflow/index.md` → domain, operations, artifacts or delivery; verification requirements are in delivery |
| Change implementation or tests | `docs/architecture/README.md` task-to-file map |
| Change placement JSON or interpretation | `docs/architecture/solution-playback.md` → “Plan JSON and evidence/geometry validation” in `docs/architecture/README.md` |
| Revise documentation | `docs/development/documentation.md` |
| Review historical engineering evidence | `maintenance/history/improvement-plan.md` and its linked audit |

`pyproject.toml` owns supported Python and extras. Nested source and test
directories inherit this router unless a closer instruction file applies.

Keep optimizer semantics, realized-array contracts, and playback presentation
separate. Playback may explain persisted placements; it must not invent a
solver-recorded order. Run the full local gate in `docs/development.md` before
handoff. Preserve all existing author credits, including Virgile Andreani's.
Attribute new work to Eric J. South; do not replace joint authorship with a
single-author header. Keep generated dogfood media outside tracked source;
the reviewed teaching example in `docs/assets/` is maintained with its guide.
