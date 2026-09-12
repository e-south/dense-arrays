---
title: Optimizer
description: Inputs, lifecycle, solver behavior, and Python signatures for motif packing.
---

# Optimizer

Use `Optimizer` when you have exact motif strings and a sequence-length limit.
Start with [one array](../quickstart.md), then add
[positional or regulator requirements](../constraints.md).

## Inputs and lifecycle

- Supply a non-empty list of non-empty uppercase `A/C/G/T` strings, a positive
  integer length, and `single` or `double` strands. Each list entry is a
  selectable motif; repeated strings remain separate entries.
- Configure all requirements before building or solving. Adding requirements
  to a built model raises `RuntimeError`. Create a new optimizer to change the
  problem; do not mutate its library or cached model state directly.
- `optimal()`, `solutions()`, and `solutions_diverse()` build the model
  automatically. Results are [DenseArray objects](results.md).
- Bound enumeration with `itertools.islice`. The bound limits returned
  results, not the time spent solving each result.

## Current solver limitations

CBC is the default. Backend availability depends on the installed OR-Tools
build. A backend that cannot be created raises `RuntimeError`.
`solver_options` forwards backend-specific strings; their acceptance is not
checked by this package. Do not treat that argument as a verified portable
timeout control.

`solve()` accepts only an optimal solver status. It raises `ValueError` for
other statuses, including infeasibility, an unproven feasible result, and
solver failure. The enumeration methods currently catch these errors and stop;
`optimal()` then reports no feasible solution. An empty or shortened iterator
therefore does not distinguish exhaustion from failure.

`approximate()` is an unconstrained greedy heuristic. It does not apply promoter
requirements, regulator coverage, or side biases configured on the optimizer.
Use the exact solver methods for those requirements. Its substring-based result
extraction can also count repeated entries at the same location or leave gaps
that the result constructor rejects. See the [audit examples](../development/audit.md#reproduce-the-core-failures)
and [selected-entry meaning](../method.md) before comparing heuristic and exact
motif counts.

## Signatures

The main calls are readable here as well as in the generated reference:

```python
Optimizer(library, sequence_length, strands="double")
optimizer.optimal(solver="CBC", solver_options=None)
optimizer.solutions(solver="CBC", solver_options=None)
optimizer.solutions_diverse(solver="CBC", solver_options=None)
optimizer.approximate()
```

For requirement arguments, use the complete [constraint examples](../constraints.md).
The generated class reference below includes advanced model methods.

::: dense_arrays.optimizer

Return to the [API index](../api.md) or locate the implementation in the
[code map](../architecture/README.md).
