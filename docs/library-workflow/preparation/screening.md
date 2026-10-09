---
title: Screen prepared parts for motif hits
description: Exclude qualifying hits using explicit motif models, scoring thresholds, and resource limits.
author: Eric J. South
---

# Screen prepared parts for motif hits

Prerequisites: create the background `recipe` using the [first recipe](../preparation.md#generate-background-parts) and `motif.json` using the [motif guide](motifs.md). Continue in that directory and Python session with FIMO installed.

## Exclude qualifying motif hits

Add `PWMExclusion` to reject prepared parts with qualifying hits to specified
motif models. Each rule declares its scoring threshold and strand policy.
GC and literal rules continue to report their own failures.

These screens apply to individual parts before assembly. Generation checks
declared literal and GC rules on the final sequence; it does not rerun
`PWMExclusion`. Joins and padding can introduce additional motif hits. When
final-library motif exclusion is required, scan the exported sequences with
the declared motif models and scoring policy.

```python
# Screen complete prepared parts, including flanks, using an explicit model.
screened_recipe = recipe.with_changes(
    screening=(
        *recipe.screening,
        parts.PWMExclusion(
            id="motif_hit",
            motifs=(parts.PWMArtifact("motif.json"),),
            scoring=parts.FimoScoring(hit_pvalue_max=0.1, strands="double"),
            reject="any_hit",
        ),
    ),
)
# Bind the exclusion model and scorer, then check scoring limits.
screened_plan = da.plan(screened_recipe)
# Save the resolved plan with its input bindings; keep the destination new.
screened_plan.write("screened.plan.json")
# Generate candidates and exclude those with qualifying motif hits.
screened_pool = da.prepare(screened_plan, out="pools/python-screened")
print(da.inspect(screened_pool, view="quality").to_dict())
```

```bash
# Apply the same motif exclusion while preparing the CLI pool.
dense-arrays prepare screened.plan.json --out pools/cli-screened
# Read retention counts and the reasons candidates were rejected.
dense-arrays inspect pools/cli-screened --view quality
```

The recipe file field uses `kind: pwm_exclusion` inside `screening`. Motifs are
ordered objects with `path` and optional `motif_ids`; scoring uses the same
fields as a PWM source. Use distinct motif IDs within each rule. Separate named
rules expose separate rejection counts.

For a score cutoff, use `reject="score_above"`, declare
`score_field="raw"` or `score_field="fraction_of_max"`, and set `threshold`.
A score exactly equal to the threshold passes. The ratio uses each motif's
FIMO-calibrated theoretical maximum. A qualifying hit with a nonpositive maximum
makes a requested ratio undefined and records an execution error. It cannot be
counted as a passed or rejected candidate. `any_hit` needs no ratio.

Only hits admitted by `hit_pvalue_max` and the declared strand policy enter
these checks. Passing means no qualifying hit met the rule's rejection criterion.
Screens cover the complete candidate, including flanks, and preserve
forward/reverse hit coordinates separately from any core used to generate it.

Planning resolves all motif files and tools, checks per-batch scorer limits and
previews the total exclusion-window bound, including one calibration per motif
per batch. It does not score candidates. Execution gives each scoring call the
remaining preparation time budget. A failed call stops the batch with execution
errors; successful earlier observations remain recorded. An incomplete scan
never becomes a no-hit result.

Every completed candidate records a hit or an explicit no-hit observation for
every bound motif, even if another screen already rejects it. Verification
checks those bindings and hit geometry, then recounts all rejection reasons
from saved observations. Reading the pool does not require the original motif
files or scorer. The human report shows up to 20 named rejection counts;
`--json` includes all counts.


Return to [prepare a sampled pool](../preparation.md).
