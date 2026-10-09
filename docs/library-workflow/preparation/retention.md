---
title: Retain scores and core diversity
description: Select representatives within one scoring model and describe their eligible score population.
author: Eric J. South
---

# Retain scores and core diversity

Prerequisites: create `motif.json` and `pwm_recipe` using the [motif guide](motifs.md). Continue in that directory and Python session with a working FIMO installation.

## Retain score and core diversity

`policy="mmr"` greedily balances score relevance and similarity to the parts
already selected. Set an explicit `pool_size`, `relevance_weight` in `(0, 1]`,
and `score_scaling`. Larger relevance weights favor score more strongly.

```python
# MMR compares eligible cores within this one motif/scoring model.
mmr_recipe = pwm_recipe.with_changes(
    scoring=parts.FimoScoring(hit_pvalue_max=1.0),
    uniqueness=parts.Uniqueness(key="core"),
    retain=parts.Retention(
        count=4,
        policy="mmr",
        rank_by="best_hit_score",
        mmr=parts.MMR(
            pool_size=32,
            relevance_weight=0.5,
            score_scaling="fraction_of_max_clipped",
        ),
    ),
)
# Resolve the request and bind its input records before execution.
mmr_plan = da.plan(mmr_recipe)
# Save the resolved plan with its input bindings; keep the destination new.
mmr_plan.write("mmr.plan.json")
# Prepare the declared parts or batch and save its identities for reuse.
mmr_pool = da.prepare(mmr_plan, out="pools/python-mmr")
print(da.inspect(mmr_pool, view="quality").to_dict())
```

```bash
# Prepare the declared pool or offered batch.
dense-arrays prepare mmr.plan.json --out pools/cli-mmr
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/cli-mmr --view quality
```

This example uses a permissive hit threshold to admit additional cores. Choose
hit and score thresholds for the intended computation. MMR cannot create
variation absent from the eligible pool.

After eligibility and deduplication, representatives are ordered by descending
raw score, then core and full sequence lexically. The highest-ranked
`pool_size` representatives enter selection. Optional
`minimum_fraction_of_max` applies an inclusive score-ratio cutoff before that
cap. These are retention decisions: excluded representatives still count as
eligible. A cap smaller than the requested count produces a visible shortfall.

To scale the choice pool with the retained target, use an explicit multiplier and
maximum instead of a fixed integer:

```python
relative_mmr = mmr_recipe.with_changes(
    retain=parts.Retention(
        count=4,
        policy="mmr",
        rank_by="best_hit_score",
        mmr=parts.MMR(
            pool_size=parts.PoolSize(per_retained=10, maximum=24),
            relevance_weight=0.5,
            score_scaling="fraction_of_max_clipped",
        ),
    ),
)
# Resolve the request and bind its input records before execution.
relative_plan = da.plan(relative_mmr)
assert relative_plan.preview["retention"]["pool_limit"] == 24
# Save the resolved plan with its input bindings; keep the destination new.
relative_plan.write("relative-mmr.plan.json")
# Prepare the declared parts or batch and save its identities for reuse.
relative_pool = da.prepare(relative_plan, out="pools/python-relative-mmr")
print(da.inspect(relative_pool, view="quality").to_dict()["retention"])
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan relative-mmr.plan.json
# Prepare the declared pool or offered batch.
dense-arrays prepare relative-mmr.plan.json --out pools/cli-relative-mmr
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/cli-relative-mmr --view quality
# Write the declared selection or document to a new destination.
dense-arrays export pools/cli-relative-mmr --view quality --out relative-mmr-quality.json
```

This example requests a choice pool of `ceil(10 × 4) = 40` and caps it at 24.
The multiplier is at least one; the positive integer maximum always bounds
selection effort. Changing the retained count resolves a new limit during
planning. A zero retained target resolves to zero. Arithmetic uses the declared
decimal multiplier under `retained_multiplier.v1`.

The score cutoff remains fixed. The highest-ranked representatives above it
enter the pool up to the resolved limit; scarce supply can leave fewer candidates.
Reports distinguish the requested size, capped limit, available supply, pool
shortfall and `has_choice`. The latter is true only when the admitted pool is
larger than a positive retained target; it does not claim that distinct cores or
biological properties were achieved. A choice-pool shortfall can coexist with
successful retention. It does not make a completed retained pool incomplete.

The planning preview gives the resolved limit and a conservative core-distance
work bound. Candidate mining keeps its separate budget. Score bands do not control
admission, and no score threshold is relaxed automatically.

Two relevance scales are available:

- `fraction_of_max_clipped` uses the recorded raw-score/theoretical-maximum
  ratio, clipped to `[0, 1]`. It requires a positive maximum.
- `score_percentile` uses average ranks for equal raw scores within the admitted
  pool, scaled to `[0, 1]`. If all scores are equal, relevance is one. This scale
  works with a zero maximum unless a ratio cutoff is also requested.

The `pwm_tolerant_hamming` distance compares motif-oriented cores at matching
positions. For position distribution `P_i` and artifact background `B`, its
weight is `1 − clip(D_KL(P_i || B) / −log2(min(B)), 0, 1)`. A mismatch receives
less weight at a position that differs more strongly from that background.
With a uniform background this emphasizes more variable positions. With a
biased background, even an invariant common base can receive a large weight;
the weights do not measure binding tolerance. Weights revert to one when their
sum is at most `1e-6`. Core length must equal motif width.

For each choice, similarity is `1 / (1 + nearest_weighted_distance)` and utility
is `relevance_weight * relevance - (1 - relevance_weight) * similarity`.
The first choice has zero similarity penalty and no nearest-neighbor evidence.
Equal utilities use descending raw score, lexical core, then lexical full
sequence. The versioned `greedy_mmr.v1` policy keeps this order even when the
whole admitted pool is retained. Changing the requested count preserves its
selection prefix when the admitted pool is unchanged.

Saved candidates carry pool admission, relevance, selection rank, utility and
nearest-selected distance/similarity where computed. Quality reports show the
admitted count and exclusions from the score cutoff and pool cap. Plans expose
a conservative bound on core-position comparisons; verification can repeat
selection from saved evidence within that bound. Selection keeps one distance
per admitted candidate rather than an all-pairs matrix. Its work grows with
pool size, retained count and core width.

## Describe the eligible score distribution

Add `ScoreBands` to a PWM recipe to report where retained parts lie in the
observed score distribution. Declare cumulative upper fractions explicitly:
`(0.1, 0.5)` creates the upper 10%, the remainder through 50%, and the remainder
of the population. These fractions are illustrative analysis choices.

```python
import json

from dense_arrays.reporting import CandidateFilter

# Bands describe the eligible score population; they do not alter retention.
band_recipe = pwm_recipe.with_changes(score_bands=parts.ScoreBands((0.1, 0.5)))
# Resolve the request and bind its input records before execution.
band_plan = da.plan(band_recipe)
assert band_plan.preview["score_bands"]["counts"] is None
# Save the resolved plan with its input bindings; keep the destination new.
band_plan.write("bands.plan.json")
# Prepare the declared parts or batch and save its identities for reuse.
band_pool = da.prepare(band_plan, out="pools/python-bands")
band_quality = da.inspect(band_pool, view="quality").to_dict()
for band in band_quality["score_bands"]["bands"]:
    print(band["band"], band["count"], band["retained"], band["scores"])

upper_band = CandidateFilter(score_bands=(1,))
# Write the example input or request so it can also be used from the CLI.
Path("upper-band.json").write_text(json.dumps(upper_band.to_dict()))
# Publish the declared records to a new destination.
da.export(
    band_pool, view="candidates", select=upper_band, all=True, out="upper-parts.json"
)
```

```bash
# Prepare the declared pool or offered batch.
dense-arrays prepare bands.plan.json --out pools/cli-bands
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect pools/cli-bands --view quality
# Write the declared selection or document to a new destination.
dense-arrays export pools/cli-bands --view candidates --selection upper-band.json --all --out cli-upper-parts.json
# Write the declared selection or document to a new destination.
dense-arrays export pools/cli-bands --view quality --out bands-quality.json
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect bands-quality.json --view quality
```

The population is every eligible unique representative in this recipe, before
retention-pool limits. Scores use the best qualifying hit's raw FIMO log2-odds,
bound to the report's `scoring_id`. Each boundary takes the score at descending
rank `ceil(fraction × population)` and includes every tied score in the upper
band. Actual counts can exceed nominal fractions; repeated boundary scores can
leave an empty band. Reports show actual eligible and retained counts and each
band's minimum, median and maximum score. Empty populations remain explicit.

Bands annotate existing decisions; they do not change sampling, representatives
or retained membership. They describe empirical ranks within one scoring model,
not biological activity or affinity. Preparation sets report each recipe
separately. To query a set by score band, select exactly one recipe, for example
`CandidateFilter(recipes=("motif_b",), score_bands=(1,))`. A band number has
meaning only with that recipe's boundaries and scoring identity.

Planning records the policy `upper_rank_include_ties.v1` and leaves counts
unknown. Verification rederives membership and summaries from saved candidate
evidence. Exported quality reports retain the population and scoring identity
and check internal consistency without requiring the pool or invoking FIMO.


Return to [prepare a sampled pool](../preparation.md).
