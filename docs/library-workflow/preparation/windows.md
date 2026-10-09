---
title: Choose motif windows and candidate lengths
description: Keep motif-window selection separate from the length and background of complete candidate sequences.
author: Eric J. South
---

# Choose motif windows and candidate lengths

Candidate length controls the complete DNA sequence. Motif-window selection
controls which positions of a supplied model are used for sampling and scoring.
Choose a task below; each example links to its setup.

## Vary candidate length

Use `parts.LengthRange(minimum, maximum)` for inclusive, uniformly sampled
integer lengths under ordinary sampling. Conditional background sampling uses
that uniform length prior conditioned jointly on the sequence constraints;
see [length conditioning](../background.md#choose-the-distribution).
Use `planning.Length(exact=N)` for a fixed length. A packing
`Length(maximum=N)` does not declare a sampling distribution and is rejected
here. Candidate count and retained count remain independent of length.

Create `recipe` using the [background first recipe](../preparation.md#generate-background-parts),
then continue in that directory and Python session:

```python
# The maximum length times 200 candidates gives the 4,800-base admission bound.
ranged_recipe = recipe.with_changes(
    sampling=parts.Sampling(length=parts.LengthRange(16, 24), strategy="conditional"),
)
# Check the length range and calculate its candidate-base bound.
ranged_plan = da.plan(ranged_recipe)
assert ranged_plan.preview["candidate_bases_bound"] == 4800
# Save the resolved plan with its input bindings; keep the destination new.
ranged_plan.write("ranged.plan.json")
# Generate and retain parts with lengths from 16 through 24 bases.
ranged_pool = da.prepare(ranged_plan, out="pools/python-ranged")
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan ranged.plan.json
# Generate variable-length parts from the saved plan.
dense-arrays prepare ranged.plan.json --out pools/cli-ranged
# Inspect five candidates and their realized lengths.
dense-arrays inspect pools/cli-ranged --view candidates --limit 5
```

In a request file, the same field is
`sampling: {strategy: conditional, length: {minimum: 16, maximum: 24}}`. Both endpoints are required,
and `exact` cannot be combined with them. Equal endpoints are valid and produce
that one length. PWM sources support ranges with all three proposal strategies;
the minimum must accommodate the selected motif. Increase the candidate length
or [choose a motif window](#choose-a-motif-window) explicitly; requesting a shorter
sequence does not trim the motif.

Ordinary length draws have their own versioned, candidate-local random stream. Batch
size, retention settings and later candidate-budget increases do not change an
existing proposal prefix. Each saved sequence records the realized length;
verification checks it against the frozen interval without drawing again.
Retained length frequencies can differ from the uniform proposal distribution
because eligibility, uniqueness and retention change the population.

The preview bounds total candidate bases using the maximum length. Scorer
admission likewise uses maximum-length windows for every member of a batch,
including calibration. A batch that exceeds the declared window cap fails in
planning even if some realized lengths might be shorter. Execution accounts for
the actual windows scanned. An exclusion motif wider than a candidate has no
full-width window: that candidate records an explicit no-hit observation and
is omitted from that motif's scoring batch. Other screens still apply.

## Choose a motif window

Declare a fixed window when preparing parts from a longer motif. Planning
selects a contiguous interval with maximum relative entropy against a declared
background and reports the information discarded. Omitting the window keeps
the full motif. A request for the original width preserves its model.

Create `motif.json` and `pwm_recipe` using the [motif guide](motifs.md#sample-a-motif),
then continue in that directory and Python session. This example selects ten of
the model's twelve positions while keeping the complete candidate at 20 bases:

```python
# Selection changes the scoring model; the complete candidate stays 20 bases.
window_recipe = pwm_recipe.with_changes(
    source=parts.PWMArtifact(
        "motif.json", window=parts.MotifWindow(length=10, background="motif")
    ),
)
# Select the window and report retained and discarded information.
window_plan = da.plan(window_recipe)
print(window_plan.preview["motif_window"])
# Save the resolved plan with its input bindings; keep the destination new.
window_plan.write("window.plan.json")
# Sample and score candidates using the selected ten-position model.
window_pool = da.prepare(window_plan, out="pools/python-window")
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan window.plan.json
# Prepare parts from the selected ten-position motif window.
dense-arrays prepare window.plan.json --out pools/cli-window
# Verify candidate decisions, model bindings and retained parts.
dense-arrays inspect pools/cli-window --verify
```

Use `background="uniform"` to rank windows by base conservation (`2 − entropy`),
or supply four positive A/C/G/T probabilities for another reference.
The selection background does not replace sampling or scoring backgrounds.
The same `PWMArtifact.window` declaration works in an exclusion rule.
In a request file it is `source.window: {length: 10, background: motif}`.

Every candidate is scored against the selected model, with fresh calibration
for that model. Shorter sequences never trigger implicit trimming. Increasing
candidate length adds sampled flanks; it does not restore discarded motif
positions. Consult the [window reference](../../reference/motif-windows.md) for
coordinates, model identities and limits on biological interpretation.

## Prepare named motif windows

Use `PreparationSet.from_windows` to expand one PWM recipe into named windows
with separate scoring and retention. This example uses the twelve-position teaching motif
and requests one retained part per window:

```python
import json

window_base = pwm_recipe.with_changes(
    budget=parts.CandidateBudget(candidates=40),
    scoring=parts.FimoScoring(hit_pvalue_max=1.0),
    retain=parts.Retention(count=1, policy="top_score", rank_by="best_hit_score"),
)
# Compare a 10-position window with the unchanged full 12-position model.
named_windows = {
    "ten_base": parts.MotifWindow(10),
    "full_length": parts.MotifWindow(12),
}
window_set = parts.PreparationSet.from_windows(
    window_base,
    windows=named_windows,
    candidate_length="window",
    max_recipes=2,
)
# Expand the two windows and check their combined budgets.
windows_plan = da.plan(window_set)
assert windows_plan.preview["candidate_budget"] == 80
assert windows_plan.preview["requested_retention"] == 2
# Save the resolved plan with its input bindings; keep the destination new.
windows_plan.write("windows.plan.json")
# Prepare each window with its own scoring and retention.
windows_pool = da.prepare(windows_plan, out="pools/python-windows")

# Publish editable settings to a new destination.
da.export(window_base, view="request", out="window-base.json")
# Save the base recipe and named-window choices as a compact request.
Path("windows.json").write_text(
    json.dumps(
        {
            "schema": "dense_arrays.preparation_windows.v1",
            "base": json.loads(Path("window-base.json").read_text()),
            "windows": [
                {"id": name, "window": window.to_dict()}
                for name, window in named_windows.items()
            ],
            "candidate_length": "window",
            "max_recipes": 2,
        }
    )
)
```

The CLI accepts either the compact declaration or its resolved plan:

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan windows.json --out windows-cli.plan.json
# Prepare both named windows from the compact request.
dense-arrays prepare windows.json --out pools/cli-windows
# Compare candidate and retained counts for each window.
dense-arrays inspect pools/cli-windows --view quality
# Inspect three candidates from the ten-position window.
dense-arrays inspect pools/cli-windows --view candidates --recipe-id ten_base --limit 3
# Verify candidate decisions, model bindings and retained parts.
dense-arrays inspect pools/cli-windows --verify
```

`candidate_length="window"` makes each candidate exactly its selected window's
width. `"base"` preserves the base recipe's exact length or range, including
sampled flanks; its minimum must accommodate every declared window. Choose the
mode explicitly. The base must be an unwindowed PWM source, so all coordinates
refer to the original motif. Window selection uses each declared background.

`max_recipes` bounds expansion before source reads. **Candidate budgets and
retention targets apply per window**, and the preview shows their totals.
Seeds and other policies are copied from the base; names do not reseed sampling.
Use explicit recipe edits when windows need different budgets, targets or seeds.
The ordinary set's sequence-collision policy still applies.

Each window has its own scoring model, calibration, uniqueness decisions and
retention. MMR therefore compares cores of one width. Scores and score bands
remain recipe-local. Expansion does not draw random widths or rank different
models together, and the selected widths do not imply biological boundaries.

Saved plans and exported requests contain complete ordinary recipe sets. They
have the same semantic identities as explicitly written equivalent recipes;
execution requires no additional expansion step. Inspect their saved models
and coordinates without reopening the motif source or running FIMO.


Return to [prepare a sampled pool](../preparation.md).
