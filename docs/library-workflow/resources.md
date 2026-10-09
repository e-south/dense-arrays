---
title: Bound preparation and packing work
description: Set candidate-base and model-size limits before allocating work, alongside time and attempt budgets.
---

# Bound preparation and packing work

Set both size and effort limits when preparing parts or generating libraries.
Size limits reject oversized requests before candidate generation or packing
model allocation. Time limits bound ongoing work cooperatively: model building
and a complete preparation batch can exceed the remaining time allowance.

## Prepare candidate sequences

`parts.CandidateBudget` separates candidate count, batch size and base limits:

| Setting | Default | Meaning |
| --- | ---: | --- |
| `candidates` | Required | Maximum candidates processed by one recipe |
| `batch_size` | 1,000 | Maximum candidates proposed together |
| `batch_bases` | 1,000,000 | Maximum declared bases in one batch |
| `total_bases` | 100,000,000 | Maximum declared bases across all candidates in one recipe |
| `seconds` | Unset | Cooperative elapsed-time allowance for one recipe |

Admission uses the maximum sequence length, including flanks, even when shorter
lengths are possible. The two size checks are:

```text
min(candidates, batch_size) × maximum_length ≤ batch_bases
candidates × maximum_length ≤ total_bases
```

These are upper bounds on candidate bases, not predictions of retained yield.
Early stopping does not reduce the admitted upper bound. Every recipe in a
preparation set must pass before execution starts; the caps apply independently
to each recipe. Candidate evidence from earlier recipes remains in memory, so
the sum of their budgets matters when sizing a set.

Use the same names in Python and request files:

```python
# Use the typed requests and operations needed by this example.
from dense_arrays import parts

budget = parts.CandidateBudget(
    candidates=100_000,
    batch_size=500,
    seconds=120,
    batch_bases=50_000,
    total_bases=10_000_000,
)
```

```yaml
budget:  # Maximum sampling effort; separate from retained count.
  candidates: 100000
  batch_size: 500
  seconds: 120
  batch_bases: 50000
  total_bases: 10000000
```

Conditional background construction also has separate
[counting-table limits](background.md#bound-work-and-interpret-outcomes). FIMO has separate
[scoring limits](../reference/motif-scoring.md). All applicable limits must pass.

## Size a packing model

`planning.Limits(model_pairs=250_000)` limits the square of the number of
oriented input parts. A single-strand search has one node per offered part; a
double-strand search has two. The pair count includes diagonal entries because
the adjacency matrix stores them. It bounds the quadratic model dimension,
rather than claiming an exact count of solver variables or constraints.

For example, 250 parts on both strands have 500 nodes and 250,000 pairs. The
default admits this model. Offering 251 parts requires either a smaller
[candidate batch](batches.md) or an explicit increase to `limits.model_pairs`.
The limit applies to exact and greedy search because both use pairwise overlaps.

For a schedule, admission uses the largest offered batch. Runtime resampling
uses its declared batch size. Matrix admission checks each active cell's
offered input; the base recipe and inactive cells do not allocate models.
Planning exposes `oriented_nodes` and `path_variables` without building a model.

Planning a large collection remains available so you can prepare a smaller
candidate batch from it. Size admission happens when starting or resuming
generation, before creating a new output or allocating a packing model.
For example, `prepare(plan(request), sampling=planning.BatchSampling(size=80,
seed=7), out="batch.json")` selects 80 offered parts before search admission.

```yaml
limits:  # Global and per-solver work allowances.
  attempts: 1000
  active_seconds: 300
  solver_seconds: 30
  model_pairs: 400000
```

An exceeded size limit fails with the calculated bound and the setting to
change. No request is shortened automatically. Saved plans and results remain
readable; execution rechecks size admission. To change a saved plan's limits,
edit its request and create a new plan.

Base and pair limits are work measures, not RAM or wall-clock guarantees.
Python objects, sequence lengths, extra constraints and native solver state
also affect resource use. Use measured workloads to choose larger allowances;
see the [packing measurements](../development/optimizer-performance.md).
