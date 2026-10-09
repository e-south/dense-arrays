---
title: Prepare a sampled pool
description: Choose a preparation task, generate reusable parts, and inspect their saved evidence.
author: Eric J. South
---

# Prepare a sampled pool

Prepare reusable DNA parts from a background distribution or a supplied motif.
Keep candidate effort, eligibility, uniqueness and retained count separate.
Planning binds inputs and limits; preparation saves candidates, decisions and
retained parts in an immutable pool. Use a new output directory for each run.

| Task | Guide | Inputs |
| --- | --- | --- |
| Generate background parts | [First recipe](preparation.md#generate-background-parts) | Declared base distribution and length |
| Generate motif-bearing candidates | [Motif recipes](preparation/motifs.md) | Motif artifact and FIMO |
| Vary candidate length or select a motif window | [Lengths and windows](preparation/windows.md) | Complete base recipe |
| Combine independently configured recipes | [Preparation sets](preparation/sets.md) | Named complete recipes |
| Retain high scores and core diversity | [Retention](preparation/retention.md) | Scored eligible candidates |
| Exclude qualifying motif hits | [Motif screens](preparation/screening.md) | Explicit exclusion motifs and scorer |
| Bound effort or diagnose a shortfall | [Effort and targets](preparation/effort.md) | Budget and retained target |
| Read decisions and reuse a pool | [Inspection](preparation/inspection.md) | Saved pool |
| Make a preparation figure | [Preparation quality](preparation-quality.md) | Saved pool or quality report |

For an existing binding-site table, start with [curated parts and pools](preparation/curated.md).
The examples use synthetic sequences and models to demonstrate the computation;
they do not establish binding or expression activity.

## Generate background parts

Save `background.yaml`:

```yaml
schema: dense_arrays.prepare.v1  # Request type and wire-format version.
source: {kind: background}  # Input model or part collection.
sampling: {strategy: conditional, length: {exact: 20}}  # Proposal strategy and candidate-length bounds.
budget: {candidates: 200}  # Maximum sampling effort; separate from retained count.
screening:
  - {id: gc, kind: gc, scope: sequence, min: 0.25, max: 0.75}
uniqueness: {key: sequence}  # Identity used to group duplicate candidates.
retain: {count: 8, policy: first_eligible}  # Number and policy for keeping eligible candidates.
seed: 7  # Seed for versioned candidate streams.
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan background.yaml --out background.plan.json
# Prepare the declared pool or offered batch.
dense-arrays prepare background.plan.json --out pools/background
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/background --view quality
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/background --view parts --all --json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/background --verify
```

The matching Python request is:

```python
import dense_arrays as da

# Use the typed requests and operations needed by this example.
from dense_arrays import parts, planning

recipe = parts.PreparationSpec(
    source=parts.Background(),
    sampling=parts.Sampling(length=planning.Length(exact=20), strategy="conditional"),
    # Examine at most 200 candidates, then retain up to eight unique parts.
    budget=parts.CandidateBudget(candidates=200),
    screening=(planning.GC("gc", "sequence", 0.25, 0.75),),
    uniqueness=parts.Uniqueness(key="sequence"),
    retain=parts.Retention(count=8, policy="first_eligible"),
    seed=7,
)
# Planning binds inputs and limits; preparation writes the immutable pool.
preview = da.plan(recipe)
assert preview.preview["retained_count_status"] == "unknown"
# Prepare the declared parts or batch and save its identities for reuse.
pool = da.prepare(preview, out="pools/python-background")
# Read saved composition and search metrics.
quality = da.inspect(pool, view="quality")
assert quality.to_dict()["counts"]["retained"] == 8
```

`Background` defaults to equal A/C/G/T probabilities and group `background`.
`base_probabilities=(0.1, 0.4, 0.4, 0.1)` supplies a different ordered distribution;
zero probabilities are allowed. Its group label is configurable. A background
pool declares how sequences were generated; it does not establish an experimental
control.

The `conditional` strategy generates candidates under the declared sequence GC
and literal exclusions. See [constrained background parts](background.md) for
distribution guarantees, construction limits and feasible/infeasible/limited
outcomes. The default `stochastic` strategy draws candidates independently and
screens them afterward.

Sequence-scoped `GC` and literal `Avoid` rules share generation's sequence
semantics. Every applicable rule is evaluated, so a candidate can carry several
rejection reasons. Placement exceptions and padding-scoped GC do not apply to
unassembled candidates.



<a id="sample-a-motif"></a>
<a id="import-a-meme-or-jaspar-motif"></a>
<a id="choose-a-pwm-proposal-strategy"></a>
<a id="vary-candidate-length"></a>
<a id="choose-a-motif-window"></a>
<a id="prepare-named-motif-windows"></a>
<a id="prepare-several-recipes-together"></a>
<a id="retain-score-and-core-diversity"></a>
<a id="describe-the-eligible-score-distribution"></a>
<a id="exclude-qualifying-motif-hits"></a>
<a id="interpret-effort-and-shortfalls"></a>
<a id="stop-after-enough-eligible-candidates"></a>
<a id="inspect-candidate-decisions"></a>
<a id="use-retained-parts"></a>
<a id="reopen-a-pool-report"></a>
