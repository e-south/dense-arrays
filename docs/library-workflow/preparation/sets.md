---
title: Prepare several recipes together
description: Combine independently configured motif and background recipes with explicit collision policies.
author: Eric J. South
---

# Prepare several recipes together

Prerequisites: create `pwm_recipe` and `example.meme` using the [motif guide](motifs.md), and the background `recipe` using the [first recipe](../preparation.md#generate-background-parts). Continue in that directory and Python session. Each output destination below must be new.

## Prepare several recipes together

Use `parts.PreparationSet` to prepare several motifs and background parts in one
pool. Give each recipe a stable name and a complete `PreparationSpec`. Budgets,
seeds, eligibility and retention remain independent; scores from different
models are never ranked together. Recipes can use different input formats,
candidate lengths and motif windows.

Continue from the linked motif and background setup examples:

```python
# A second synthetic motif supplies an independent scoring model.
second_consensus = "TGCAGTACCGAT"
# Write the example input or request so it can also be used from the CLI.
Path("second.jaspar").write_text(
    ">second\n"
    + "\n".join(
        f"{base} [ {' '.join('7' if base == expected else '1' for expected in second_consensus)} ]"
        for base in "ACGT"
    )
    + "\n"
)
# Each recipe keeps its own 200-candidate budget and eight-part target.
set_recipe = parts.PreparationSet(
    {
        "motif_a": pwm_recipe.with_changes(
            source=parts.PWMArtifact("example.meme", format="meme"),
        ),
        "motif_b": pwm_recipe.with_changes(
            source=parts.PWMArtifact("second.jaspar", format="jaspar"),
        ),
        "background": recipe,
    },
    sequence_collisions="preserve",
)
# Resolve the request and bind its input records before execution.
set_plan = da.plan(set_recipe)
assert set_plan.preview["candidate_budget"] == 600
assert set_plan.preview["requested_retention"] == 24
# Publish editable settings to a new destination.
da.export(set_recipe, view="request", out="set-request.json")
# Save the resolved plan with its input bindings; keep the destination new.
set_plan.write("set.plan.json")
# Prepare the declared parts or batch and save its identities for reuse.
set_pool = da.prepare(set_plan, out="pools/python-set")
print(da.inspect(set_pool, view="quality").to_dict()["recipes"])
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan set-request.json
# Prepare the declared pool or offered batch.
dense-arrays prepare set.plan.json --out pools/cli-set
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/cli-set --view quality
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/cli-set --view candidates --recipe-id motif_b --limit 5
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/cli-set --verify
```

The request schema is `dense_arrays.preparation_set.v1`. Its `recipes` is an
ordered array of `{id, request}` entries, each containing a complete
`dense_arrays.prepare.v1` request. Exporting a Python request writes this format
with relative input paths. Recipe names are distinct from motif IDs and part
groups. Mapping order controls execution and publication order. Reordering or
adding a recipe does not change another recipe's candidate-local random stream.

Each retained part has a qualified ID such as `motif_b/candidate_12`. Candidate
records have a pool-wide `index` and representative link, plus `recipe_id` and
`recipe_index`; ranks are local to the recipe. Select candidates with
`CandidateFilter(recipes=("motif_b",))` or repeatable CLI `--recipe-id`.

The default `sequence_collisions="error"` refuses publication if different
recipes retain the same complete sequence. The example explicitly preserves
such occurrences under distinct part IDs. There is no cross-recipe score-based
deduplication. Within each recipe, its existing uniqueness policy still applies.

Set `core_collisions="error"` to also reject equal observed core strings retained
by different recipes:

```python
distinct_cores = set_recipe.with_changes(core_collisions="error")
# Publish editable settings to a new destination.
da.export(distinct_cores, view="request", out="distinct-cores.json")
```

In JSON or YAML, set top-level `core_collisions: error` on a preparation-set or
named-window request. The default is `preserve`. Both collision choices appear
in the human plan preview; neither changes a recipe's candidate draws or ranking.

Cores are compared in their recorded motif orientation. In this minimal
two-base orientation example, a forward `AC` hit and
a reverse-strand `GT` slice both describe the oriented core `AC` and collide.
Forward `AC` and forward `GT` remain distinct. Parts without an observed core,
including background parts, are excluded from this comparison. Exact equality
does not establish biological redundancy or comparable scores across models.

A collision reports the two recipes and parts and prevents publication. Revise
the recipes or explicitly preserve the occurrences; preparation does not drop,
replace, rescore or resample them to satisfy the policy. Sequence and core
collision checks are independent and are repeated when verifying a saved pool.

All inputs and required tools are preflighted before output creation. Recipes
then execute sequentially under their own limits; a recorded scoring failure or
shortfall does not prevent another recipe from executing. The pool is complete
only if every recipe meets its target without execution errors. Reports retain
each stopping condition and stage count alongside totals. All records publish
in one transaction; interruption before commit leaves no readable pool. Source
changes and cross-recipe collision errors also prevent publication. Preparation
resume, nested sets and curated-table entries are not supported.


Return to [prepare a sampled pool](../preparation.md).
