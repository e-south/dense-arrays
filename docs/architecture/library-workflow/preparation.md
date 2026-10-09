---
title: Preparation recipes
description: Compare curated and sampled pool requests, effort bounds, retention and recorded candidate accounting.
author: Eric J. South
---

# Preparation recipes

Preparation creates reusable parts with explicit input, scoring and retention
evidence. The [domain contract](domain.md#part-inputs) defines part identity;
[operation rules](operations.md) define invocation and publication.

## Curated pool

Use the 16-base-site `parts.csv` from the
[curated binding-site example](../../library-workflow/preparation/curated.md#bind-curated-parts-and-requirements).
For a reusable curated pool, `curate.yaml` is:

```yaml
schema: dense_arrays.prepare.v1
source:
  kind: table
  table: parts.csv
  format: csv
```

```text
dense-arrays prepare curate.yaml --out pools/curated
dense-arrays inspect pools/curated --view parts --group A --all
```

```python
import dense_arrays as da
from dense_arrays import parts

curated_request = parts.PreparationSpec(
    source=parts.PartTable(
        table="parts.csv",
        format="csv",
    )
)
curated = da.prepare(curated_request, out="pools/python-curated")
selected_parts = da.inspect(
    curated, view="parts", select=parts.PartFilter(groups=("A",)), all=True
)
```

Reusing a prepared pool binds its snapshot and optional `PartFilter`; it never
remines or rescores. [Input rules](domain.md#part-inputs) define mappings, identity,
normalization and source-row diagnostics.

## PWM and background preparation

Preparation supports three named source kinds under
`dense_arrays.prepare.v1`: `table`, `pwm_artifact`, and `background`. Source
acquisition and motif inference remain external. Each recipe has independent
mining effort, eligibility, uniqueness, and retention settings; increasing a
retained target does not silently increase effort.

This parameterized recipe requires a native `motif.json` containing motif ID
`example`, a FIMO background file at `background.txt`, and the external FIMO
executable. [Motif inputs](../../reference/motif-scoring.md) defines those file
contracts; [sampled-pool recipes](../../library-workflow/preparation.md) supplies
the user workflow. Save the request as `prepare.yaml`:

```yaml
schema: dense_arrays.prepare.v1
source: {kind: pwm_artifact, path: motif.json, motif_ids: [example]}
sampling: {strategy: stochastic, length: {exact: 20}}
budget: {candidates: 10000}
scoring: {backend: fimo, hit_pvalue_max: 0.0001, background: background.txt}
eligibility: {best_hit_score_min_exclusive: 0.0}
uniqueness: {key: core}
retain: {count: 100, policy: top_score, rank_by: best_hit_score}
seed: 7
```

```python
import dense_arrays as da
from dense_arrays import parts, planning

preparation = parts.PreparationSpec(
    source=parts.PWMArtifact(path="motif.json", motif_ids=("example",)),
    sampling=parts.Sampling(strategy="stochastic", length=planning.Length(exact=20)),
    budget=parts.CandidateBudget(
        candidates=10000
    ),  # Cap proposals independently of retention.
    scoring=parts.FimoScoring(hit_pvalue_max=0.0001, background="background.txt"),
    eligibility=parts.Eligibility(best_hit_score_min_exclusive=0.0),
    uniqueness=parts.Uniqueness(key="core"),
    # Retain up to 100 eligible unique parts; this does not guarantee that yield.
    retain=parts.Retention(count=100, policy="top_score", rank_by="best_hit_score"),
    seed=7,
)
preparation_plan = da.plan(preparation)
pool = da.prepare(preparation_plan, out="pools/python-example")
pool_report = da.inspect(pool, view="quality")
```

```text
dense-arrays plan prepare.yaml --out prepare.plan.json
dense-arrays prepare prepare.plan.json --out pools/example
dense-arrays inspect pools/example --view quality
```

Run from the request directory with declared motif/background files present and
the optional scorer available. Preflight validates artifact IDs, score-policy
versions, tool availability and output ownership before mining. The report may
show fewer than 100 retained parts; it must identify incomplete status and the
stage responsible. It never labels that pool complete merely because effort ended.

Use `top_score` and MMR retention as explicit choices.
MMR requests must name score scaling, sequence/core distance, relevance weight,
eligible-pool bound, tie rule, and algorithm version. Preserve clipping,
trimming, strand, flank/core geometry, background model, and selection lineage
where supported. Unsupported combinations fail with a named field and a
supported next action.
Scores carry scorer/version, units, orientation and background identity.
Raw log-odds score, per-length normalization, fraction of theoretical maximum,
and p-value are distinct fields; no unlabeled `score_norm` in the native contract.
These settings describe computation, not binding affinity or biological activity.

Uniqueness applies within each preparation recipe. A
[`PreparationSet`](../../library-workflow/preparation/sets.md#prepare-several-recipes-together)
combines named recipes with independent collision policies for retained full
sequences and observed motif-oriented cores. Sequence collisions default to
`error`; core collisions default to `preserve`. Either policy may reject equal
strings retained by different recipes. Parts without an observed core are
excluded from core comparison. These checks prevent publication on a collision;
they do not alter selection or establish biological redundancy.

Background recipes declare length, GC range, fixed candidate budget, sequence
uniqueness, and optional literal-pattern/PWM-hit exclusions. They use the same
named screening definitions as final sequence checks where meanings match.
Scoring remains optional unless a request explicitly requires it. A background
pool is an input choice, not an inferred experimental control.

A background request uses `source: {kind: background}`, equal A/C/G/T base
probabilities unless an explicit distribution is supplied, a length contract
and fixed candidate budget. `strategy: stochastic` draws independently and
screens afterward. `strategy: conditional` samples under declared GC and literal
exclusions; its [distribution and resource contract](../../library-workflow/background.md)
defines zero-mass and resource-limited outcomes. Its `screening` list can reuse `avoid` and sequence-scoped `gc`
requirements; retention uses `policy: first_eligible` in stable candidate order
or a separately versioned supported selection policy. For PWM-hit exclusion,
the preparation-only `pwm_exclusion` screen names motif artifacts, scorer,
background, hit p-value threshold, strands and either zero-hit-only or a named
score cutoff. These fields must be explicit; it is not a literal `avoid` rule.
Absent optional scoring means no claim that PWM hits were screened.

Every pool quality report includes supplied/mined, eligibility rejection,
eligible, duplicate-discarded, eligible-unique, retained, and not-selected counts,
with requested retention, effort, and completion. The stage equations reconcile:
processed candidates = eligibility-rejected + eligible + execution-error;
eligible = duplicate-discarded + eligible-unique; eligible-unique = retained +
not-selected. Scoring/execution failures are separate from biological filters.
Per-candidate decisions contain all rejection reasons, representative mapping,
and rank/selection evidence where computed. Bounded diagnostics summarize even
when full rejected sequence retention is disabled; the retention policy and any
missing detail are explicit. Never rerun preparation to manufacture an explanation.
