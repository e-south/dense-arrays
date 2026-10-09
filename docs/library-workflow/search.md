---
title: Choose packing search
description: Choose exact or greedy search, interpret their guarantees, and favor underused parts during exact enumeration.
author: Eric J. South
---

# Choose packing search

Generation defaults to exact CBC search within the declared oriented-path packing model.
The primary objective maximizes placed part occurrences. With no preference,
remaining paths are enumerated by that objective. Set
`packing_preference="underused_parts"` to favor parts used less often in earlier
proposed packings from the same offered batch.

This preference can reduce bias while enumerating alternatives. It does not
optimize nucleotide distance, measure binding function or guarantee equal part
usage in the accepted library. Screening and duplicate rejection can change the
composition of accepted designs independently of packing search.

| Request | Search behavior | Guarantee |
| --- | --- | --- |
| `search="exact"` (default) | Enumerate optimal remaining packing paths. | Proven optimality within each offered model. |
| `search="exact", packing_preference="underused_parts"` | Favor below-mean proposed usage while maximizing occurrences. | The occurrence objective remains primary. |
| `search="greedy"` | Choose one deterministic greedy packing per offered batch. | Validated placements and final screens; optimality is unproven. |

Exact enumeration searches the offered packing paths; it does not enumerate
arbitrary gaps or every possible DNA sequence. The seed controls sampling and
padding streams, while the solver determines the order of equal packing optima.
See [assembly](generation/assembly.md#interpret-padding-outcomes) for how each
packing proceeds through padding and final screening.

## Use the same request in Python and the CLI

Run this example in a new directory:

```python
import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

# Declare the part collection, sequence bounds and generation policy.
request = planning.DesignSpec(
    parts=[
        parts.Part("a", "ACGTTGCAAGTCCTGA"),
        parts.Part("b", "GATCAGTACCTAGGTC"),
        parts.Part("c", "TTGACCGATAGCTACG"),
    ],
    length=planning.Length(maximum=32),
    strands="single",
    target=planning.Target(4),
    packing_preference="underused_parts",
)
# Publish editable settings to a new destination.
da.export(request, view="request", out="balanced.json")
# Resolve the request and bind its input records before execution.
resolved = da.plan(request)
assert resolved.preview["packing_preference"] == "underused_parts"
# Save the resolved plan with its input bindings; keep the destination new.
resolved.write("balanced.plan.json")
# Generate under the declared bounds into a new output directory.
library = da.run(resolved, out="runs/python")
assert da.inspect(library, verify=True).accepted == 4

with da.inspect(library, view="attempts", all=True).records() as records:
    attempts = list(records)
first = attempts[0].evidence["packing_objective"]
assert first["weights"] == {"a": 1.0, "b": 1.0, "c": 1.0}
assert first["proposed_packings"] == 0
assert all(len(row.candidate.packed.placements) == 2 for row in attempts)
print(attempts[1].evidence["packing_objective"])
```

The preference lives in the request file; it adds no command or separate CLI
execution path:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan balanced.json --out balanced-cli.plan.json --json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run balanced-cli.plan.json --out runs/cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/cli --view attempts --all --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/cli --verify --json
```

## Meaning of the weights

For `N` offered parts, let `u_i` be the number of earlier proposed packings that
placed part `i`. Count a part once regardless of orientation. Its weight is:

```text
1 + 0.5/N, if N * u_i < sum(u)
1,         otherwise
```

The maximum total bonus is 0.5, below the contribution of one additional
occurrence. The preference therefore preserves occurrence count as the primary
objective. Ties remaining after weighting are resolved by the solver; a seed is
not a promise of identical solver tie choices across versions or environments.
The versioned preference is `underused_parts.v1`.

Every completed solve records its effective weights, per-part usage counts,
number of earlier proposed packings and population scope under
`attempt.evidence["packing_objective"]`. Accepted, rejected and duplicate
proposals all contribute to subsequent usage. A solve without a recorded packing
contributes nothing. Prior path exclusions use the same committed packing prefix.
Feasible-but-unproven and backend outcomes retain their distinct termination
semantics; the preference never permits an unproven candidate to be accepted.

## Cells, batches and recovery

Each matrix cell has independent history. A prepared batch uses its offered
identities and cardinality, rather than the full eligible collection. Scheduled
or resampled batches start new usage histories. Pool sampling and its
coverage/failure feedback remain separate controls over which parts are offered.

On supported resume, the engine restores both excluded paths and usage weights
from committed packing records. Unresolved interrupted attempts still consume
their declared effort; they do not invent a placement history. Existing batch
advancement and resource limits continue to apply.

`inspect(..., verify=True)` independently recounts usage from saved placements
and checks each objective before advancing the history. It does not rerun the
solver or independently prove optimality. Reader limits cover the objective's
usage/weight maps and the retained per-cell history. Saved plan comparisons
include preference changes, and portable libraries retain the declared policy.
An exported library without attempt history cannot establish the sequence of
search decisions.

## Generate a greedy proposal

Choose `search="greedy"` when a quick packing proposal is useful and an optimality
proof is unnecessary. The method extends a path greedily from each fitting
oriented entry, chooses the path with most placed occurrences, and breaks ties
by shorter realized length, then stable input order. It uses the same overlap
geometry as exact packing and selects each supplied part identity at most once.

Continue from the Python example above:

```python
greedy = request.with_changes(
    search="greedy", packing_preference=None, target=planning.Target(1)
)
# Resolve the request and bind its input records before execution.
greedy_plan = da.plan(greedy)
assert greedy_plan.preview["solver"] is None
assert greedy_plan.preview["proof_scope"] is None
# Save the resolved plan with its input bindings; keep the destination new.
greedy_plan.write("greedy.plan.json")
# Generate under the declared bounds into a new output directory.
greedy_library = da.run(greedy_plan, out="runs/greedy-python")
assert da.inspect(greedy_library, verify=True).accepted == 1

# One batch provides one proposal, even when a larger target is requested.
shortfall = da.run(
    greedy.with_changes(target=planning.Target(2)), out="runs/greedy-shortfall"
)
# Read saved run state and attainment. Recount stored evidence before returning.
shortfall_report = da.inspect(shortfall, verify=True)
assert shortfall_report.accepted == 1
assert shortfall_report.termination_reason == "heuristic_exhausted"
```

Use the same saved plan from the CLI:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan greedy.plan.json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run greedy.plan.json --out runs/greedy-cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/greedy-cli --view attempts --all --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/greedy-shortfall --view diagnostics --json
```

A proposal can be rejected by screening or duplicate checks. Greedy search does
not replace it with another path from that batch. Use [prepared or resampled
batches](batches.md) for further proposals. Cell quotas, batch caps and global
effort limits still apply; unused targets are not silently reallocated.

Greedy supports maximum length and single- or double-strand eligibility. Exact
final length requires explicit padding. Literal exclusions and GC screens run
on the assembled sequence. Counts, group coverage, fixed placements, spacing and
`packing_preference` fail planning with the unsupported rule or setting named.
Use exact search for these requirements. There is no automatic method fallback.

The saved `heuristic` evidence names `greedy_multistart.v1` and distinguishes a
candidate, exhausted offered search and a cooperative time limit. It carries no
solver status or proof. Exhaustion does not establish infeasibility or enumerate
all possible designs. `limits.solver_seconds` bounds each packing search in both
modes; greedy checks its deadline between starts and extensions. Adjacency
construction and one extension scan are not hard wall-clock bounds.

Resume restores whether the batch's proposal was consumed from its committed
packing. An interrupted attempt without a committed packing can retry that
proposal under the remaining budget. Verification binds evidence to the declared
method, checks the one-proposal rule per cell and batch, and recounts placements
and final requirements without rerunning the heuristic.
