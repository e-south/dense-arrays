---
title: Generate constrained background parts
description: Sample DNA under GC and literal exclusions with a defined distribution and bounded construction work.
author: Eric J. South
---

# Generate constrained background parts

Use `Sampling(strategy="conditional")` to draw background parts that already
satisfy sequence-scoped `GC` and `Avoid` rules. This is useful when independent
draws would rarely pass the constraints. Conditional sampling counts valid
suffixes once, then reuses those counts for every candidate in the recipe.

```python
import dense_arrays as da

from dense_arrays import parts, planning

recipe = parts.PreparationSpec(
    source=parts.Background(),
    sampling=parts.Sampling(
        planning.Length(exact=30),
        strategy="conditional",
        limits=parts.ConditionalLimits(seconds=10),
    ),
    screening=(
        planning.GC("gc", "sequence", 29 / 30, 1),
        planning.Avoid("no_triples", ("CCC", "GGG")),
    ),
    budget=parts.CandidateBudget(candidates=100, batch_size=20),
    retain=parts.Retention(count=20, policy="first_eligible"),
    seed=7,
)
# Check the GC and exclusion rules and bind the construction limits.
resolved = da.plan(recipe)
# Save the resolved plan with its input bindings; keep the destination new.
resolved.write("background.plan.json")
# Generate constrained candidates and retain 20 distinct parts.
pool = da.prepare(resolved, out="pools/python-background")
report = da.inspect(pool, view="quality").to_dict()
assert report["construction"]["status"] == "feasible"
assert report["counts"]["eligibility_rejected"] == 0
assert report["counts"]["retained"] == 20
```

Use the same saved plan through the CLI:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan background.plan.json
# Prepare a separate pool from the same saved plan.
dense-arrays prepare background.plan.json --out pools/cli-background
# Read construction status and candidate counts.
dense-arrays inspect pools/cli-background --view quality
# Verify saved decisions, constraints and retained-part records.
dense-arrays inspect pools/cli-background --verify
```

In a YAML or JSON request, the sampling declaration is:

```yaml
sampling:  # Proposal strategy and candidate-length bounds.
  strategy: conditional
  length: {exact: 30}
  limits: {seconds: 10}
```

`plan` records the distribution, compiled rule IDs and construction limits. It
does not build counting tables or predict retained yield. Each recipe builds its
own table during `prepare`; there is no global cache.

## Choose the distribution

`Background.base_probabilities` gives the ordered A/C/G/T probabilities before
constraints. Conditional sampling conditions that distribution on valid
sequences. Equal base probabilities at a fixed length give every valid sequence
equal probability; unequal probabilities preserve the sequences' relative
background weights. Bases with zero probability are never sampled.

`LengthRange(minimum, maximum)` supplies a uniform **prior** over lengths. Joint
conditioning on sequence validity can change that distribution: for lengths one
or two with forward `AA` forbidden, length probabilities become `16/31` and
`15/31`. The preview labels this as a uniform prior conditioned on constraints.
This strategy draws length and sequence together.

Draws are with replacement; uniqueness and retention remain separate stages.
Use the default `stochastic` strategy for independent candidate sampling followed
by screening. Conditional sampling requires a `Background` source; PWM proposals have their own
[construction strategies](preparation/motifs.md#choose-a-pwm-proposal-strategy).

## Bound work and interpret outcomes

Candidate proposals must pass the [batch and total base limits](resources.md)
before counting or sampling starts. Construction has additional limits below.

`ConditionalLimits` controls counting, independently of candidate effort:

| Field | Default | Bound |
| --- | --- | --- |
| `states` | 250,000 | Admitted counting states and length entries |
| `automaton_states` | 4,096 | Forbidden-pattern prefix states |
| `mass_bits` | 32,000,000 | Sum of retained integer-mass bit lengths; also bounds each temporary integer |
| `seconds` | 30 | Cooperative construction allowance |

The smaller remaining allowance from `CandidateBudget.seconds` also applies to
construction. Time is checked during table work and between complete candidate
batches; a batch may overrun the remaining allowance. These controls are not a
hard whole-operation deadline. The state caps bound table overhead separately
from integer mass storage. There is no silent switch to another generator.

| Construction status | Meaning | Preparation behavior |
| --- | --- | --- |
| `feasible` | Exact positive mass was counted | Draw candidates, apply remaining screens, deduplicate and retain |
| `infeasible` | Exact mass is zero | Save `construction_infeasible` with zero processed candidates |
| `limited` | A state, mass or time cap stopped counting | Save `construction_limited` with its specific cause; feasibility remains unknown |

The proof concerns the declared base support, length bounds and compiled rules.
Increase the named limit or revise the constraints in a new request when
counting is limited. Construction failures
produce an incomplete pool and CLI exit **3**. They do not become candidate
rejections. A zero mining target stops before construction begins.

`PWMExclusion` remains independent post-proposal FIMO screening and requires
explicit motif/scorer bindings. It can cause a retained-part shortfall even
after feasible construction. Duplicate draws can also limit retained yield.
Part validity does not establish that joins or padding in a later assembled
design satisfy final-sequence constraints.

## Inspect saved evidence

Exact integer masses preserve the normalized decimal probability ratios
recorded in the request. Randomness is versioned as
`conditional_background_shake256.v1` and bound to model, seed and candidate index.
Batch size, retention, effort limits and later candidate-budget increases
preserve the existing sequence prefix when counting completes. Changing the
compiled rules, length range or base probabilities changes the model and its stream.

Quality reports expose a separate `construction` record containing model ID,
policy, status, cause and work counts. Its hexadecimal `mass` is the exact
integer sampling mass; for weighted or ranged requests it is not a count of
distinct sequences. Candidate proposal evidence binds the same model.

`inspect(..., verify=True)` checks those bindings, independent sequence
constraints, positive base support, work allowances and reconciled candidate
decisions. It reads saved evidence without rebuilding the counting table or
running FIMO. The recorded zero-mass result is a construction outcome; inspection
does not independently repeat its proof. Pools and exported quality reports
remain portable through the existing [inspection and export operations](preparation/inspection.md#inspect-candidate-decisions).
