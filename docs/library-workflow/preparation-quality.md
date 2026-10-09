---
title: Plot preparation quality
description: Inspect candidate yield, recorded selection distances and score bands from a saved pool.
---

# Plot preparation quality

Use `preparation-quality` to see where candidates were lost and what was retained.
The figure reads a sampled pool or its exported quality report. Install the
`playback` extra to write PNG figures.

```bash
# Render a figure from the selected saved records.
dense-arrays render pool --view preparation-quality --out preparation.png
# Write the declared selection or document to a new destination.
dense-arrays export pool --view quality --out preparation-quality.json
# Render a figure from the selected saved records.
dense-arrays render preparation-quality.json --view preparation-quality --out shared.png
```

```python
import dense_arrays as da

# Read saved composition and search metrics.
report = da.inspect("pool", view="quality")
print(report.cost)
# Draw the selected saved evidence without generating new sequences.
da.render(report, view="preparation-quality", out="preparation-python.png")
# Publish the declared records to a new destination.
da.export(report, out="preparation-quality-python.json")
```

Each recipe has its own row:

- **Candidate yield:** processed, eligible, representative and retained counts.
  These are nested populations. The dashed line marks requested retention;
  completion and stopping reason stay visible even when no candidates were made.
- **Selection distance:** recorded MMR distance to earlier selected cores, in
  selection order. The first choice has no comparison. The weighted distance
  depends on that recipe's motif model; it is neither a base count nor an
  all-pairs diversity statistic.
- **Score bands:** eligible representative counts and their retained subsets.
  Labels show declared rank intervals. Ties can increase actual band sizes.
  Counts describe each recipe's own scoring model and population.

A report without MMR distances or declared score bands shows that evidence as
unavailable. It does not substitute zeros. Sequence uniqueness, core uniqueness
and selection distances are separate properties; representative counts alone
do not measure core diversity.

Live reports verify saved candidate decisions and retained-part joins. MMR
verification repeats selection comparisons under `ReadLimits.pairs`, without
sampling or invoking FIMO. Use `--max-pairs` in the CLI, or bind `read_limits`
when inspecting a pool in Python. Exported reports retain their recorded
metrics and can be plotted after source files and scoring tools are unavailable;
they check internal consistency without re-verifying the source pool.

PNG output includes the complete report and its digest as metadata. The receipt
identifies the pool, plan and number of retained parts. Publication creates a
new destination and never overwrites an existing file. A figure supports up to
12 recipes and 24 score bands per recipe; larger reports remain available as
quality JSON. Long recipe labels are shortened only in the figure, with full
identities retained in its metadata.

See [preparation inspection](preparation/inspection.md#inspect-candidate-decisions) for
individual candidates and [saved report exports](handoffs.md) for portable data.
