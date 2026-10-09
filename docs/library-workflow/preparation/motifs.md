---
title: Prepare parts from a motif
description: Create a PWM recipe, import equivalent motif formats, and choose how candidate sequences are drawn.
author: Eric J. South
---

# Prepare parts from a motif

## Create a motif artifact

This example prepares parts from a synthetic 12-position model with consensus
`ACGTTGCAAGTC`. Run the code in a new directory with the
[library workflow](../../installation.md#use-the-library-workflow) and
[optional FIMO scoring](../../installation.md#configure-fimo-for-motif-scoring) installed.

```python
import json
import math
from pathlib import Path

import dense_arrays as da

from dense_arrays import parts, planning

consensus = "ACGTTGCAAGTC"
bases = "ACGT"
# Each position favors its consensus base without making the draw deterministic.
probabilities = [
    {base: 0.7 if base == expected else 0.1 for base in bases} for expected in consensus
]
background = {base: 0.25 for base in bases}
motif_document = {
    "schema_version": "1.0",
    "producer": "dense-arrays-documentation",
    "motif_id": "example",
    "alphabet": bases,
    "matrix_semantics": "probabilities",
    "background": background,
    "probabilities": probabilities,
    # This supplied log2 matrix supports direct matrix-scoring examples.
    # FIMO independently derives scores from probabilities and its own settings.
    "log_odds": [
        {base: math.log2(row[base] / background[base]) for base in bases}
        for row in probabilities
    ],
}
# Save the model probabilities, background and score matrix.
Path("motif.json").write_text(json.dumps(motif_document, indent=2))
```

## Sample a motif

For the model above, each 20-base candidate contains twelve bases drawn from the
motif probabilities and eight from the background distribution. The motif's
position is sampled within the candidate. FIMO scans the complete sequence;
the example's permissive p-value threshold applies to each motif-width window.
You can substitute a sourced motif that fits within the declared candidate length.

Save `pwm.yaml`:

```yaml
schema: dense_arrays.prepare.v1  # Request type and wire-format version.
source: {kind: pwm_artifact, path: motif.json}  # Input model or part collection.
sampling: {strategy: stochastic, length: {exact: 20}}  # Proposal strategy and candidate-length bounds.
budget: {candidates: 200}  # Maximum sampling effort; separate from retained count.
scoring: {backend: fimo, hit_pvalue_max: 0.1, strands: double}  # Scorer settings applied to candidate windows.
eligibility: {best_hit_score_min_exclusive: 0}  # Requirements for a candidate to enter retention.
uniqueness: {key: sequence}  # Identity used to group duplicate candidates.
retain: {count: 8, policy: top_score, rank_by: best_hit_score}  # Number and policy for keeping eligible candidates.
seed: 7  # Seed for versioned candidate streams.
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan pwm.yaml --out pwm.plan.json
# Sample and score candidates, then retain the eight highest-scoring parts.
dense-arrays prepare pwm.plan.json --out pools/cli-pwm
# Read candidate yield and retention counts.
dense-arrays inspect pools/cli-pwm --view quality --json
# Export the pool quality report as JSON.
dense-arrays export pools/cli-pwm --view quality --out pwm-quality.json
```

The matching Python request is independent of the background example:

```python
import dense_arrays as da

from dense_arrays import parts, planning

pwm_recipe = parts.PreparationSpec(
    source=parts.PWMArtifact("motif.json"),
    # Twelve modeled positions fit within each 20-base candidate.
    sampling=parts.Sampling(length=planning.Length(exact=20)),
    # Candidate effort and retained count are separate limits.
    budget=parts.CandidateBudget(candidates=200),
    # The p-value threshold applies to each scanned motif-width window.
    scoring=parts.FimoScoring(hit_pvalue_max=0.1, strands="double"),
    eligibility=parts.Eligibility(best_hit_score_min_exclusive=0),
    uniqueness=parts.Uniqueness(key="sequence"),
    retain=parts.Retention(count=8, policy="top_score", rank_by="best_hit_score"),
    seed=7,
)
# Check the motif and FIMO, then resolve sampling and scoring settings.
pwm_plan = da.plan(pwm_recipe)
# Save the resolved plan with its input bindings; keep the destination new.
pwm_plan.write("pwm-python.plan.json")
# Sample, score and retain eight parts in the Python output pool.
pwm_pool = da.prepare(pwm_plan, out="pools/python-pwm")
print(da.inspect(pwm_pool, view="quality").to_dict())
```


The default `stochastic` strategy samples PWM positions from the artifact probabilities. Any additional
positions use its background distribution; the sampled motif's insertion offset
is drawn across the available positions. Each candidate has a separate seeded
stream, so changing batch size preserves proposals. Score and retention changes
do not redraw the same candidate stream. The sampling algorithm version is
recorded in the plan.

The best qualifying hit defines the retained part's motif-oriented core and its
zero-based, half-open coordinates. This hit can differ from the sampled motif's
insertion position. Raw score, score per core base, theoretical maximum,
fraction of maximum and p-value retain separate labels. See the
[scoring reference](../../reference/motif-scoring.md#score-candidates-with-fimo)
for background and numerical interpretation.

Planning reads the motif and checks FIMO's version before any sampling or
scoring. Execution rechecks the source and executable bytes. The
CLI reports scorer preflight failures with exit **4** and, with `--json`,
`code: scoring_error` and a separate `reason`, such as `unavailable`.
No output destination is created when preflight fails.

`Uniqueness(key="core")` groups equal oriented cores. `key="sequence"` groups
identical complete candidates. Representatives use the highest recorded score,
then earliest candidate index. `top_score` ranks those representatives by the
same order; `first_eligible` retains them by candidate index. Equal sequences in
different groups are an error. Each recipe contains one motif per source; use
a [preparation set](sets.md) to combine independently configured recipes.

## Import a MEME or JASPAR motif

Set the format on `PWMArtifact`; all sampling, window selection, scoring and
retention settings apply to the imported model. These synthetic files describe
the same 12-position probability matrix:

```python
from pathlib import Path

consensus = "ACGTTGCAAGTC"
bases = "ACGT"
# Synthetic counts sum to 10 at every position; normalization gives 0.7/0.1.
counts = {
    base: [7 if base == expected else 1 for expected in consensus] for base in bases
}
meme_rows = [
    " ".join(str(counts[base][position] / 10) for base in bases)
    for position in range(len(consensus))
]
# Write the probability matrix in MEME format.
Path("example.meme").write_text(
    "MEME version 5\nALPHABET= ACGT\n\nMOTIF example\n"
    "letter-probability matrix: alength= 4 w= 12 nsites= 10\n"
    + "\n".join(meme_rows)
    + "\n"
)
# Write the equivalent counts in JASPAR format.
Path("example.jaspar").write_text(
    ">example\n"
    + "\n".join(f"{base} [ {' '.join(map(str, counts[base]))} ]" for base in bases)
    + "\n"
)
for input_format in ("meme", "jaspar"):
    # Only the input representation changes; sampling and retention stay fixed.
    imported_recipe = pwm_recipe.with_changes(
        source=parts.PWMArtifact(
            f"example.{input_format}", format=input_format, motif_ids=("example",)
        ),
    )
    imported_plan = da.plan(imported_recipe)
    imported_plan.write(f"{input_format}.plan.json")
    da.prepare(imported_plan, out=f"pools/python-{input_format}")
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan meme.plan.json
# Prepare parts from the imported MEME model.
dense-arrays prepare meme.plan.json --out pools/cli-meme
# Prepare parts from the equivalent JASPAR model.
dense-arrays prepare jaspar.plan.json --out pools/cli-jaspar
# Verify saved candidates, scoring observations and retained parts.
dense-arrays inspect pools/cli-meme --verify
# Verify saved candidates, scoring observations and retained parts.
dense-arrays inspect pools/cli-jaspar --verify
```

A multi-record file requires one exact primary ID; alternate names do not select
records. Omit `motif_ids` for a single-record file. To prepare another motif,
create a separate request and destination. A list of IDs does not implicitly
create multiple pools. Request files use
`source: {kind: pwm_artifact, path: example.meme, format: meme, motif_ids: [example]}`.

Import preserves probabilities without rounding them to inferred counts.
JASPAR counts are normalized per position without pseudocounts and retained as
source evidence. Missing backgrounds explicitly use uniform frequencies.
The preview reports the input format and absence of a supplied score matrix;
FIMO derives scores using the recipe's configured calibration. Consult the
[input reference](../../reference/motif-scoring.md#import-meme-or-jaspar) for accepted
formats, saved evidence and source-statistic interpretation.

## Choose a PWM proposal strategy

`Sampling.strategy` separates how a candidate is constructed from how FIMO
scores it. Every PWM-source strategy uses the same hit eligibility, sequence
screens, uniqueness and retention rules.

| Strategy | Candidate construction |
| --- | --- |
| `stochastic` | Draw an embedded motif from its position probabilities, then draw flanking bases from the sampling background. |
| `consensus` | Embed the highest-probability base at each motif position, breaking exact ties in A/C/G/T order; sample the offset and flanks. |
| `background` | Draw the whole candidate from the sampling background, with no embedded motif; then score it against the declared motif. |

The sampling background defaults to the motif artifact's background. Set
`Sampling.base_probabilities` to four explicit A/C/G/T probabilities to override
it. This affects flanks or whole-background proposals; it does not change the
FIMO background model declared in `scoring`. Plain `Background` sources own
their distribution directly and support `stochastic` and `conditional` sampling.

Continue from `pwm_recipe` above:

```python
proposal_pools = {}
for strategy in ("consensus", "background"):
    proposal_recipe = pwm_recipe.with_changes(
        sampling=parts.Sampling(
            length=planning.Length(exact=20),
            strategy=strategy,
            base_probabilities=(0.4, 0.1, 0.1, 0.4),
        ),
    )
    proposal_plan = da.plan(proposal_recipe)
    proposal_plan.write(f"pwm-{strategy}.plan.json")
    proposal_pools[strategy] = da.prepare(
        proposal_plan, out=f"pools/python-{strategy}-proposals"
    )
```

```bash
# Validate inputs and inspect or save the resolved plan.
dense-arrays plan pwm-consensus.plan.json
# Prepare consensus cores with sampled flanks and offsets.
dense-arrays prepare pwm-consensus.plan.json --out pools/cli-consensus-proposals
# Draw complete background candidates and score them against the motif.
dense-arrays prepare pwm-background.plan.json --out pools/cli-background-proposals
```

Consensus can produce several full sequences through different flanks and
offsets. It does not increase the candidate budget to satisfy retention.
At motif width, its sequence is deterministic; repeated proposals are removed by
uniqueness. Core uniqueness can likewise leave a shortfall despite varied flanks.
Background proposals with no qualifying motif hit are recorded as
`no_qualifying_hit`, not silently retained.

These explicit proposals record their intended forward interval as zero-based,
half-open `part.metadata["proposal"]["start"]` and `end`. Background proposals
record both as `None`. The best scored hit may occur elsewhere or on the reverse
strand; that observed hit owns the part's core coordinates. Verification checks
proposal geometry, consensus content and permitted bases from saved evidence,
without redrawing random candidates.


Return to [prepare a sampled pool](../preparation.md).
