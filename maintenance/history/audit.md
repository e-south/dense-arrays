---
title: September 2026 baseline audit
description: Historical findings from before the September 2026 contract and playback repairs.
---

# September 2026 baseline audit

Author: Eric J. South. Audited 12 September 2026 against
[`006b361e`](https://github.com/e-south/dense-arrays/tree/006b361e462c460d9fc398bd55ed7be174790f76).
This records the baseline findings, before runtime repairs. The
[completed hardening record](improvement-plan.md) describes the repairs and
verification, including removal of the separate HTML playback renderer.

The package has useful ownership boundaries and working small examples. Its
largest risks are incomplete input contracts and output that can hide a failure
or overstate what a saved layout proves. Passing existing tests does not close
these gaps. Rewriting modules before fixing those contracts would make the
behavior harder to assess.

## Findings that affect results

Priority 1 means a result or its interpretation can be wrong. Priority 2 means
an input, state, or presentation contract needs repair before further expansion.

| ID | Priority | Observed behavior | Code owner |
| --- | --- | --- | --- |
| C1 | 1 | Greedy approximation ignores configured hard requirements. Repeated entries can share one occurrence; some feasible inputs produce invalid gapped results. | [optimizer.py:756](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/optimizer.py#L756) |
| C2 | 1 | Iterators catch every solve `ValueError`, including abnormal backend status and invalid result construction. `optimal()` then calls the failure infeasibility. | [optimizer.py:542](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/optimizer.py#L542), [enumeration:667](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/optimizer.py#L667) |
| P1 | 1 | Saved plans bypass equivalent semantic validation. Invalid sequence, references, reveal spans, and contradictory constraint results load successfully. | [models.py:128](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/models.py#L128), [serialization.py:260](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/serialization.py#L260) |
| P2 | 1 | Visible playback omits ordering qualifications and failed constraints. A `layout_only` plan still draws an active edge. Metadata embedded in HTML is not visible evidence. | [html.py:151](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/html.py#L151), [matplotlib_renderer.py:147](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/matplotlib_renderer.py#L147) |
| C3 | 2 | Fractional quotas are truncated; malformed ranges and noninteger lengths fail late or change the requested problem. | [constraints.py:38](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/constraints.py#L38), [counts:87](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/constraints.py#L87) |
| C4 | 2 | Mutating the public library leaves cached overlaps stale. A rejected right-bias update can already have changed the left bias. Foreign results and invalid motif indices can mutate a model. | [optimizer.py:62](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/optimizer.py#L62), [biases:163](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/optimizer.py#L163), [forbid:593](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/optimizer.py#L593) |
| C5 | 2 | Direct results accept empty/invalid motifs, boolean offsets, and conflicting entry-orientation semantics. Realized records coerce identity values and accept enum states that later fail serialization. | [solution.py:37](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/solution.py#L37), [realized.py:30](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/realized.py#L30) |
| P3 | 2 | HTML ignores presentation settings; subtitles and label overrides have no renderer consumers. Generic theme code infers study-specific categories from labels. | [html.py:18](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/html.py#L18), [theme.py:63](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/theme.py#L63) |
| P4 | 2 | Playback imports the optimizer through package initialization; semantic graph projection imports Matplotlib. Default and injected layout engines also receive different context topology. | [package initialization](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/__init__.py#L16), [graph initialization](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/graph/__init__.py#L4), [layout.py:284](https://github.com/e-south/dense-arrays/blob/006b361e462c460d9fc398bd55ed7be174790f76/src/dense_arrays/playback/graph/layout.py#L284) |

The CLI adds operational gaps: unavailable backends and malformed playback
inputs can produce tracebacks, rejected solver options are ignored, and HTML
is written before optional exports are checked. Export paths can overwrite
inputs or one another. These belong to the existing CLI owners, not a new
orchestration layer.

## Reproduce the core failures

Run small examples with `uv run python` in the locked checkout; exact comparisons
below use CBC. These describe the audited behavior, not desired regression-test
assertions.

```python
from dense_arrays import Optimizer

optimizer = Optimizer(["AAA", "CCC"], 3, "single")
optimizer.add_regulator_constraints(["R1", "R2"], required={"R2"})
print(optimizer.approximate().sequence)  # AAA: misses required R2
print(optimizer.optimal().sequence)  # CCC: includes required R2

optimizer = Optimizer(["AAA", "AAA"], 3, "single")
print(optimizer.approximate().nb_motifs)  # 2, both offsets are 0
print(optimizer.optimal().nb_motifs)  # 1 under exact path-entry semantics

optimizer = Optimizer(["AAA", "CCC"], 3, "single")
optimizer.add_regulator_constraints(["R", "R"], min_count_by_regulator={"R": 1.9})
print(optimizer.optimal().nb_motifs)  # 1: the fractional quota was truncated
```

Additional probes established these cases:

| Input or action | Result |
| --- | --- |
| `CCC` upstream of `AAA`, fixed at zero with no spacer | Greedy returns `AAACCC`; exact returns `CCCAAA` |
| `['AAAAAC', 'ACG', 'CGT']`, length 5 | Greedy reconstruction reports a gap; exact returns valid `ACGT` |
| Synthetic backend returns `ABNORMAL` | Direct solve reports abnormality; enumeration returns empty; `optimal()` reports no feasible solution |
| Change library `['AAA', 'AAT']` to `['AAA', 'CCC']` after construction | Cached overlaps are reused and the failed reconstruction is reported as infeasibility |
| `set_motif_weight(-1, 2.0)` on a built model | Terminal arcs are modified instead of rejecting the index |
| Construct `DenseArray(['N'], 1, [0], [None])` | Construction succeeds; printing raises `KeyError` |

The abnormal-status probe used an isolated test double to exercise error
propagation. It is not evidence of an actual CBC backend failure.

## Playback negative controls

A valid reconstructed plan was serialized, then changed one field at a time.
All nine invalid cases loaded: sequence mismatch, duplicate placement IDs,
unknown predecessor, out-of-bounds reveal span, extra reveal-span field,
boolean coordinate, fractional coordinate, contradictory constraint `passed`,
and unknown constraint placement. The out-of-bounds span later raised
`IndexError` in raster rendering.

A real array `AAATTTCCC` with placements at `[0,3)` and `[6,9)` and a declared
zero-distance constraint correctly reconstructs as `layout_only` with
`passed=False`. An inspected Matplotlib graph still drew an active green edge.
Even enabling the optional authority notice displayed only “Realized order.”
HTML source inspection found no rendering of order, notices, or constraint
results. Browser-level interaction is a separate acceptance check in the plan.

## Documentation and module structure

The original API source was 43 lines but built to 482,209 bytes with 48
headings. GitHub readers saw mostly generator directives. The documentation
pass splits references by interface, adds source-readable behavior summaries,
and retains the API index's existing section anchors.

The largest responsibilities are concentrated in `optimizer.py` (819 lines),
`matplotlib_renderer.py` (952), and graph `routing.py` (695). The renderer's
graph-drawing function alone spans 214 lines; MP4 and GIF writers repeat their
frame loops. These are reasons to separate model construction, drawing, frame
scheduling, and encoding. They are not evidence that every long file needs
splitting or that performance is poor. Broad playback lint exceptions should
shrink as those owners become explicit.

The documentation pass also adds task-based agent routing, a file/test map,
consistent front matter, plainer openings, and current limitations. Existing
author credits and the shared banner remain intact. Presentation aspirations
are distinguished from implemented behavior.

## Evidence and limits

Three delegated readers worked independently: core contracts, playback, and
fresh task routing. All runtime code remained unchanged. Baseline focused
suites passed: 52 core/CLI tests and 26 playback/import tests, with three
existing SWIG deprecation warnings in each invocation.

The fresh reader completed the documented CBC array, regulator coverage, and
saved-JSON HTML examples. Root and nested routing inherited the repo router;
there were no deeper scoped instructions. Global and workspace-ancestor
instructions were identified separately. The missing contributor file/test
route and source-readable API were the main routing gaps.

This establishes local behavior for small cases. It is not an exhaustive
algorithm proof, performance benchmark, fresh installation test, remote-host
check, or Gurobi validation. No remote job or Gurobi solve was run. The full
documentation handoff gate is recorded in the accompanying PR; remaining
runtime repairs require the plan's own tests and review.
