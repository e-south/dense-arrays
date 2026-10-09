---
title: Bound preparation effort and explain shortfalls
description: Separate candidate effort, eligible supply, retained targets, and incomplete outcomes.
author: Eric J. South
---

# Bound preparation effort and explain shortfalls

The mining-target example below is standalone after installing Dense Arrays. To inspect an existing shortfall, begin with the saved pool produced by a [preparation recipe](../preparation.md).

## Interpret effort and shortfalls

`CandidateBudget.candidates` is independent of `Retention.count`; changing the
retained target never increases the effort cap. `batch_size` defaults to 1,000.
An optional positive `seconds` budget is checked between batches and bounds each
FIMO call by the remaining time. This is a cooperative budget, not a hard deadline
for the whole Python call. Scorer window and output limits are separate; planning
rejects a batch whose maximum window count exceeds the configured scorer limit.

Quality reports reconcile these stages:

```text
processed = eligibility_rejected + eligible + execution_error
eligible = duplicate_discarded + eligible_unique
eligible_unique = retained + not_selected
```

They also show requested retention, candidate budget, all rejection counts and
the stopping reason. A candidate that fails several rules contributes once to
`eligibility_rejected` and once to each corresponding rejection count.
Use [preparation-quality figures](../preparation-quality.md) to plot yield, recorded
MMR selection distances and declared score bands for each recipe.

The pool is incomplete when fewer than the requested parts are retained, a
declared mining target is unmet, or a scoring batch fails. Failed scoring candidates count as execution errors, never
as eligibility rejection. The CLI writes the pool result and exits **3** for a
retention or mining-target shortfall, or **4** for recorded scoring failure. A completed pool exits
**0**. Native output is committed together; interrupted pool creation has no
committed pool and cannot be resumed or overwritten.

Generation accepts completed pools. An incomplete pool remains available for
inspection and export, including when it retained every requested part but missed
its mining target. Review `inspect --view quality`, revise the preparation budget
or requirements, and prepare into a new destination before generating a library.

Conditional construction also saves an incomplete pool when it proves empty
sampling support or reaches a construction limit, with CLI exit **3**. Its
work and outcome stay separate from processed candidate counts.

Inspection and verification use recorded candidates and scores, without
rerunning the sampler or FIMO. Quality reports expose their scan cost and accept
`read_limits` in Python or the existing CLI read caps. Verification checks
selection, accounting, score geometry and retained-part joins. Background
candidates must use bases with positive declared probability and cannot carry
a primary motif core or score. Recorded scoring failures require a declared
scorer. These checks do not
independently recompute FIMO's statistical model. Exact or bounded-range lengths, `first_eligible`, `top_score` and `mmr`
are supported. Unsupported policies fail explicitly.

## Stop after enough eligible candidates

Add `MiningTarget` to stop once the observed eligible supply reaches a declared
size. `eligible_unique` names an absolute count. Alternatively,
`max_retained_fraction=f` requests at least `ceil(Retention.count / f)` unique
eligible candidates. For example, retaining eight parts with a maximum fraction
of 0.125 requires at least 64 unique eligible candidates.

```python
import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

# Eight retained parts at fraction 0.125 require at least 64 eligible unique candidates.
mining_recipe = parts.PreparationSpec(
    parts.Background(),
    parts.Retention(count=8, policy="first_eligible"),
    sampling=parts.Sampling(planning.Length(exact=20)),
    budget=parts.CandidateBudget(candidates=256, batch_size=16),
    mining_target=parts.MiningTarget(max_retained_fraction=0.125),
    seed=7,
)
# Resolve the request and bind its input records before execution.
mining_plan = da.plan(mining_recipe)
assert mining_plan.preview["mining_target"]["eligible_unique"] == 64
# Save the resolved plan with its input bindings; keep the destination new.
mining_plan.write("mining-target.plan.json")
# Prepare the declared parts or batch and save its identities for reuse.
mining_pool = da.prepare(mining_plan, out="pools/mining-python")
# Read saved run state and attainment. Recount stored evidence before returning.
mining_summary = da.inspect(mining_pool, verify=True)
assert mining_summary.state == "completed"
assert mining_summary.preparation.stop_reason == "mining_target"
assert mining_summary.preparation.counts["eligible_unique"] >= 64
# Publish the scoped quality report to a new destination.
da.export(mining_pool, view="quality", out="mining-quality.json")
```

Use the same saved plan through the CLI:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan mining-target.plan.json
# Prepare the declared pool or offered batch.
dense-arrays prepare mining-target.plan.json --out pools/mining-cli --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/mining-cli --view quality
```

The target is checked after each complete scored and screened batch. A batch may
overshoot the count; its candidates all participate in final representative
selection and retention. Changing `batch_size` can therefore change where mining
stops and which parts are retained. Candidate and time caps still bound work.
With no mining target, preparation uses the full candidate budget unless time or
execution failure stops it. `minimum_candidates` optionally postpones target
stopping until that much candidate effort has been processed.

Eligible supply uses the declared uniqueness rule after scoring and screening:
different flanks around the same oriented core count once under core uniqueness.
The fraction describes this observed eligible population. It does not establish
binding affinity. Retention remains a separate policy: top-score ranking and MMR
keep their existing selection rules, pool limits and score cutoffs.

Choose exactly one target form. The fraction must be in `(0, 1]`; the absolute
count must be positive. A fraction target with zero requested retention and no
minimum effort stops before generating candidates. Reaching the effort cap
before a target produces an inspectable incomplete pool, even when all requested
parts were retained. Quality reports preserve the resolved target and its
attainment; preparation sets report these independently for each recipe.

Verification recounts eligibility and uniqueness and checks that target stopping
occurred at the first eligible batch boundary. The policy is
`eligible_unique_batch_stop.v1`. Fraction targets resolve with decimal arithmetic;
they are planning rules for sample supply, rather than population estimates.


Return to [prepare a sampled pool](../preparation.md).
