---
title: Dense Arrays hardening plan
description: Completed contract, usability, and module-boundary work with verification evidence.
---

# Dense Arrays hardening plan

Author: Eric J. South. Status: implemented and locally verified on
12 September 2026, following the [audit](audit.md).

The six slices below address every finding C1–C5 and P1–P4. Users can distinguish
infeasibility from execution failure, validate saved placements before export,
and interpret media without mistaking coordinate reconstruction for a recorded
solver path. The [architecture map](../architecture/README.md) routes each task
to its implementation and focused tests.

## Scope and decisions

Dense Arrays remains responsible for motif packing, explicit realized placements,
and playback. All seven original author-credit lines remain unchanged, including
Virgile Andreani's. Extracted model construction retains joint authorship; new
implementation and documentation are attributed to Eric J. South.

Playback uses the established NetworkX/Matplotlib publication pipeline and
producer-supplied duplex frames. The separate HTML renderer, unused SVG-frame
path, duplicate export loops, ignored presentation aliases, and label-based
biological inference were removed. Existing packing and constraint showcase
artifacts supplied the visual baseline.

Exact and greedy packing count selected library entries and occurrences.
Compatible contained placements remain representable, but they do not imply an
exact path; `forbid()` validates that path before changing the model. Constraints
require genuine integers and preserve intentional negative spacers. Inputs,
results, and nested provenance are immutable snapshots.

Python and JSON records share semantic validation. Supported v1 plans describe
coordinate reconstruction; reserved solver-recorded authority is rejected.
Producers supply evidence for recovered coordinates explicitly. The
[migration guide](../migration.md) records caller changes.

The DenseGen adapter migration belongs to its producer repository. Dense Arrays
has no dependency on DenseGen, Research Studies, HOP, their recipes, or their
data. Scientific claim expansion, exact solver-trace formats, framework
replacement, and new study dependencies were outside this work.

## Delivered slices

| Slice | Result | Regression coverage |
| --- | --- | --- |
| 1. Solver outcomes and approximation | Only proven infeasibility ends enumeration normally. Execution, unproven-optimality, and invalid-result failures remain distinct. Unsupported approximation requirements fail before work. | `test_solver_outcomes.py`, `test_greedy.py`, `test_cli.py` |
| 2. Configuration and result identity | Validated immutable problem inputs; integer counts, bounds, and indices; consistent entry/orientation semantics; rejected bias, weight, and forbid operations leave state unchanged. | `test_core_contracts.py` and existing exact/constraint tests |
| 3. Saved-plan validation | Python and JSON share nested shape, reference, sequence/span, reveal, order, and constraint checks. Valid failed constraints remain representable; contradictory results fail before rendering. | `test_playback_contracts.py`, including all nine audit mutations |
| 4. Truthful publication playback | Authority, order, and failed requirements stay visible. Layout-only scenes have no active placement chain. Explicit caller labels, colors, and legends replace study inference. Producer sizing and annotations are preserved. | `test_playback_presentation.py`, `test_playback_graph.py`, actual showcase media |
| 5. Module boundaries | Model construction and greedy packing have separate owners. Presentation, scheduling, drawing, graph routing, and writer lifecycle are separate. Default/injected graph engines agree; semantic imports avoid solver and raster dependencies. | `test_optional_playback_imports.py`, graph/curve tests, narrowed Ruff rules |
| 6. Commands and documentation | Both CLIs give concise errors. Export preflight rejects aliases/collisions, stages requested formats, and reports partial publication. Runnable guides, separate API references, and task routing have maintained checks. | `test_playback_cli.py`, `test_playback_output.py`, `test_playback_exports.py`, `test_documentation.py` |

## Publication contracts

Compact media preserves authority, order, and failed-requirement qualifications.
Long evidence and many native distance brackets use visible summaries with all
declared distance results in native metadata; optional notices remain
caller-selected. The [presentation reference](../reference/playback-presentation.md#read-the-evidence)
explains how to read that evidence. This policy keeps the established canvas
proportions while avoiding clipped text.

Export destinations are checked before rendering. Requested formats are staged
together, then published atomically per file. A publication failure reports any
files already written. A drawing failure preserves prior outputs, closes figures,
and retains the original error if writer cleanup also fails.

Producer callbacks supply a constant image shape within each scene, with at most
two frame snapshots cached for transitions. The DenseGen migration pins each
scene to its completed-frame crop, preventing a late terminus reveal from changing
frame dimensions. Both anchored showcase final frames remain pixel-identical to
the previous tight crop.

## Verification record

The [full local gate](../development.md#local-verification) passed:

- **323 tests**, including guide execution and built-link checks; three existing
  SWIG deprecation warnings remain.
- Pre-commit, Ruff, formatting, strict documentation build, and wheel/source builds.
  Distribution inspection confirms that the HTML player is absent.
- The complete frozen dependency export, including all extras and hashes, passed
  vulnerability auditing without resolving a different environment.
- All four fresh-reader tasks passed: local CBC packing, regulator coverage,
  saved-placement PNG interpretation, and nested contributor routing. The nested
  malformed-input check rejected invalid coordinates before creating output.
- Independent review matched **32 small CBC cases** against an exhaustive path
  oracle. Its two export/evidence findings were repaired and re-reviewed; no
  confirmed issue remains in the reviewed scope.
- Both existing producer showcase recipes exported MP4/poster bundles. Packing
  produced 602 frames at 2400×360; constraints produced 269 frames at 2400×450,
  both at 30 fps. Native video playback and metadata were inspected.
- Documentation was checked in plain Markdown and at wide and narrow browser
  widths. The shared banner, task routes, and readable reference pages were retained.

These checks establish bounded local behavior, not an exhaustive algorithm proof,
performance benchmark, or backend portability claim. Only local CBC was used;
no Gurobi or BU SCC job ran. A remote Gurobi check remains optional only if a
specific backend question arises. Local documentation builds do not publish the
hosted site; hosted CI and publication state belong to the PR/release record.
