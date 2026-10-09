---
title: Inspect and reuse prepared parts
description: Read saved candidates, inspect selection decisions, and reuse retained parts without rescoring.
author: Eric J. South
---

# Inspect and reuse prepared parts

Continue from the [background first recipe](../preparation.md#generate-background-parts),
where `pool` is the completed prepared collection.

## Inspect candidate decisions

Inspect a retained candidate to see its original sequence, candidate index and
selection outcome. Candidate indices follow sampling order; retained-part
ordinals follow retention order. Outcome and reason filters help investigate
shortfalls after you have seen which decisions the pool contains.

```python
from dense_arrays.reporting import CandidateFilter

# Inspect saved decisions; this does not sample candidates or run FIMO.
retained_page = da.inspect(
    pool,
    view="candidates",
    select=CandidateFilter(outcomes=("retained",)),
    limit=5,
)
with retained_page.records() as decisions:
    retained_examples = [row.to_dict() for row in decisions]
    continuation = decisions.next_cursor

# Export all candidate decisions, including those not retained.
da.export(pool, view="candidates", all=True, out="background-candidates.json")
```

```bash
# Show the first five retained candidates.
dense-arrays inspect pools/background --view candidates --outcome retained --limit 5
# Read the original candidate at sampling index one.
dense-arrays inspect pools/background --view candidates --candidate-index 1 --json
# Export every candidate decision as JSON.
dense-arrays export pools/background --view candidates --all --out cli-background-candidates.json
```

Each row includes its pool identity, complete sequence, rejection reasons,
representative index, rank and any recorded scoring or MMR evidence. Candidate
indices and representative links are local to that pool. To examine rejected
candidates, choose `outcomes=("eligibility_rejected",)` or CLI
`--outcome eligibility_rejected`; use the quality report's observed reason codes
to narrow the query further.

`CandidateFilter` accepts `indices`, `outcomes`, `reasons`, optional set
`recipes` and recipe-local `score_bands`. Pass a saved predicate with
`--selection`; the corresponding convenience CLI
options are repeatable `--candidate-index`, `--outcome` and `--reason`.
Values within a field are alternatives; different fields must all match.
Unknown indices and rejection reasons fail explicitly. Rejection reasons are
resolved against the pool's observed rejection counts. Outcomes are
`eligibility_rejected`, `execution_error`, `duplicate_discarded`, `retained`
and `not_selected`. The last outcome describes an eligible unique representative
that was not retained. Every rejected candidate preserves all applicable reasons.

A candidate page defaults to 100 rows. Pass its continuation token as Python
`after=continuation` or CLI `--after TOKEN` with the same filter. Read limits bound
examined records, including rows excluded by the filter. Close an iterator after
early termination, or use `with` as above. Unfiltered pages read indexed rows;
filtered queries scan candidate evidence. Notebook and terminal summaries stay
compact; `to_dict()` and `--json` expose complete records.

Pages validate the rows they read. Use `inspect(pool, verify=True)` or
`inspect POOL --verify` to check complete accounting, representatives and
retained-part joins. JSON candidate export streams the complete selected
population and publishes a new file only after a successful read. Curated table
pools expose supplied `parts`; they do not invent mined candidate evidence.

## Use retained parts

```python
# Reuse saved 20-base parts in a separate 40-base packing request.
design = planning.DesignSpec(
    parts=parts.PoolSource(pool),
    length=planning.Length(maximum=40),
    strands="single",
)
# Pack the retained parts into an array at most 40 bases long.
run = da.run(design, out="runs/prepared")
assert da.inspect(run, verify=True).accepted == 1
# Publish the scoped quality report to a new destination.
da.export(pool, view="quality", out="background-quality.json")
```

Pools retain normalized input and scoring evidence, so reading, verification
and generation reuse work after the original motif file or scorer is removed.
Use the pool's full identity when joining parts across collections; candidate
indices alone are local to their pool.

## Reopen a pool report

Export quality once to share its counts independently of the pool directory:

```python
# Reopen candidate-yield and retention counts from the exported report.
recorded_report = da.inspect("background-quality.json", view="quality")
assert recorded_report.to_dict()["pool_id"] == pool.pool_id
# Save a copy that preserves the report values and identities.
da.export(recorded_report, out="background-quality-copy.json")
```

```bash
# Read the quality report independently of the source pool.
dense-arrays inspect background-quality.json --view quality
# Save another copy of the recorded report.
dense-arrays export background-quality.json --view quality --out cli-background-quality-copy.json
```

The returned `reporting.PoolQualitySnapshot` checks the saved schema, source
identities, stage equations and completion state. It can be read after the pool,
motif files and scorer are unavailable. The source pool's candidate evidence is
not included in this report; `--verify` requires that native pool. Recorded
reports cannot be filtered, paginated or compared as library-quality metrics.
Re-export preserves the recorded values and their pool/plan identities.

Return to [prepare a sampled pool](../preparation.md).
