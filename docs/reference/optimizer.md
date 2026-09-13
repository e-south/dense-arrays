---
title: Optimizer
description: Inputs, lifecycle, solver outcomes, and Python signatures for motif packing.
---

# Optimizer

Use `Optimizer` when you have exact motif strings and a sequence-length limit.
Start with [one array](../quickstart.md), then add
[positional or regulator requirements](../constraints.md).

## Inputs and lifecycle

- Supply a non-empty sequence of non-empty uppercase `A/C/G/T` strings, a
  positive integer length, and `single` or `double` strands. Each library entry
  is selectable once, in at most one orientation; repeated strings retain
  separate entry identities.
- The optimizer keeps an immutable input snapshot. Its library, adjacency, and
  configured-requirement views return independent copies. Mutating a supplied
  list or returned view does not reconfigure the problem.
- Configure requirements before building or solving. Adding requirements to
  a built model raises `RuntimeError`; use a fresh optimizer for another problem.
- `optimal()`, `solutions()`, and `solutions_diverse()` build the model
  automatically and return [DenseArray objects](results.md).
- Bound enumeration with `itertools.islice`. This limits returned results,
  not the time spent solving each result.

## Solver outcomes

CBC is the default. A requested backend must be available in the installed
OR-Tools build. `solver_options` contains backend-specific strings; rejected
options raise `ValueError` before replacing an existing model. Their syntax
and support depend on that backend, so they are not a portable timeout API.

Only an optimal solver status returns a result. Exceptions are exported from
`dense_arrays` and `dense_arrays.errors`:

| Exception | Meaning |
| --- | --- |
| `InfeasibleError` | No feasible result exists; a `ValueError` subclass |
| `UnprovenSolutionError` | A feasible incumbent exists, but optimality was not proved |
| `SolverBackendError` | The backend is unavailable, fails to run, or returns an unsuccessful status |
| `InvalidSolverResultError` | Reported solver output fails path or result validation |
| `OptimizationError` | Base `RuntimeError` for the three execution failures above |

The enumeration methods stop normally only on `InfeasibleError`. Other failures
propagate even after a result has been yielded. An iterator that stops because
the caller reached an `islice` bound has not proved exhaustion.
`optimal()` raises `InfeasibleError` if enumeration has no first result.

Malformed problem inputs and rejected API operations use `ValueError` or
`TypeError`; they are separate from backend outcomes. Both optimizer commands
report failures on stderr and exit nonzero. See [CLI behavior](cli.md).

## Greedy approximation

`approximate()` builds a feasible path with a multi-start greedy heuristic;
it does not prove optimality. It rejects configured promoter constraints,
regulator requirements, side biases, and model changes made by `forbid()` or
`set_motif_weight()`. Use an exact solver method when those requirements matter.

The heuristic records selected entries and orientations as it builds the
path. Repeated motifs need separate placements, and an incidental contained
substring does not add a selected entry. See the [packing method](../method.md)
for this counting rule. If no entry fits, it raises `InfeasibleError`.

## Advanced model operations

`build_model()` creates the model; `solve()` uses it without rebuilding.
`forbid(result)` requires the same library, length limit, strand policy, and
exact path-placement semantics. `set_motif_weight(index, weight)` requires an
original library index and a finite real weight. Invalid indices, boolean
weights, and nonfinite values are rejected before model mutation.

## Signatures

Create an optimizer, then choose a solving method:

```text
Optimizer(library, sequence_length, strands="double")
optimizer.optimal(solver="CBC", solver_options=None)
optimizer.solutions(solver="CBC", solver_options=None)
optimizer.solutions_diverse(solver="CBC", solver_options=None)
optimizer.approximate()
```

For requirement arguments, use the complete [constraint examples](../constraints.md).

::: dense_arrays.optimizer.Optimizer

Return to the [API index](../api.md) or locate the implementation in the
[code map](../architecture/README.md).
