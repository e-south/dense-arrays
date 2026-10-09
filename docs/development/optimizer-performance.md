---
title: Packing model construction
description: Measure model-building cost and preserve subtour semantics when changing OR-Tools constraint assembly.
author: Eric J. South
---

# Packing model construction

Build subtour constraints with direct OR-Tools coefficients to reduce Python
model-building work. Each selected edge requires its successor to have a higher
integer rank. For `N` oriented entries, the row is:

`u_i - u_j + N × X_ij ≤ N - 1`, with `1 ≤ u_i, u_j ≤ N`.

When `X_ij = 1`, the row requires `u_j ≥ u_i + 1`. When the edge is absent,
the row imposes no further restriction within these rank bounds. A directed
cycle therefore cannot satisfy the selected edges. Direct coefficients preserve
the existing inequality while avoiding temporary expression trees for each
of the `N × (N - 1)` possible edges.

## Measured construction cost

The [measurement record](measurements/optimizer-model-build.json) contains the
exact libraries, before/after methods, all timings, model identities and CBC
outputs. Measurements used Python 3.14.7, OR-Tools 9.15.6755 and macOS 26.5 arm64.
Each case used three paired, sequential runs after imports. Times include
optimizer construction, model building and deterministic model serialization.
They exclude solver search.

| Input parts, both strands | Expression assembly median | Direct coefficients median | Traced Python peak, before/after |
| --- | ---: | ---: | ---: |
| 30 | 65.17 ms | 41.75 ms | 3.44 / 3.44 MB |
| 70 | 364.39 ms | 236.13 ms | 18.44 / 18.43 MB |

These synthetic libraries contain 12-base motifs and use a maximum array length
of 100 bases. The solver has a one-second allowance, although the timed region
does not call it. Allocation measurements use a separate `tracemalloc` run and
exclude native solver memory. The results support this construction change;
they do not establish faster search or predict end-to-end throughput.

## Correctness checks

The complete serialized models are byte-identical before and after the change
for both timed cases. A separate six-part, single-strand example uses real CBC
with a one-second allowance per solve. Its first five solutions, occurrence
counts and subsequent path exclusions agree exactly.

`tests/packing/test_continuity.py` independently exhausts all 64 directed edge
subsets on three nodes. For each subset, it checks whether the generated rows
admit an integer rank assignment and compares that result with the existence
of a topological order. This verifies the absence of cycles without depending
on the expression used to construct the rows.

Run the focused optimizer checks from the repository root:

```bash
uv run pytest tests/packing tests/test_optimize.py tests/test_solver_outcomes.py -q
```

Retain [solver outcome distinctions](../library-workflow/search.md) and the
[full development gate](../development.md) when changing model construction.
