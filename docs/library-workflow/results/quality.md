---
title: Explain shortfalls and assess a library
description: Inspect rejection evidence, composition and source-run attainment.
author: Eric J. South
---

# Explain shortfalls and assess a library

Use the saved run from [Generate a saved library](../../library-workflow.md).
Install the [playback extra](../../installation.md#optional-features) for PNG output.

```python
import dense_arrays as da  # Inspect saved evidence and render reports.
from dense_arrays import reporting  # Select a consistent report population.

result = "runs/python-first"  # Output of the first-library Python example.
summary = da.inspect(result)  # Original run attainment remains separate from selection.
with da.inspect(result, view="designs", limit=1).records() as records:
    design = next(records)  # Bind one identity for the filtered-report example.
selected = reporting.DesignFilter(design_ids=(design.reference,))  # Reusable predicate.
```

## Explain shortfalls and assess a library

Inspect the run's recorded failures and the composition of accepted designs:

```bash
# Read or verify saved evidence without generating again.
dense-arrays inspect runs/first --view diagnostics --limit 20
# Read or verify saved evidence without generating again.
dense-arrays inspect runs/first --view attempts --outcome rejected --limit 10
# Read or verify saved evidence without generating again.
dense-arrays inspect runs/first --view quality --limit 20 --json
# Render saved evidence to a new image.
dense-arrays render runs/first --view library-quality --out quality.png
```

Diagnostics report stable codes, requirement references, observed/expected
values, proof scope, and a next action. A forbidden match includes its final
coordinates, strand, and intersecting part or padding intervals. An attempt
limit explains a shortfall; it does not prove global infeasibility. Reason totals
can overlap, while attempt outcomes are mutually exclusive.

```python
# Use the typed requests and operations needed by this example.
from dense_arrays import reporting

diagnostics = da.inspect(
    result, view="diagnostics", limit=20
)  # Bound displayed reasons.
print(diagnostics.cost.to_dict())  # descriptor available before the scan
print(diagnostics.to_dict())
# Inspect rejected attempts separately from accepted library members.
rejected = da.inspect(
    result,
    view="attempts",
    select=reporting.AttemptFilter(outcomes=("rejected",)),
    all=True,
)
with rejected.records() as records:
    for attempt in records:
        print(attempt.attempt_id, attempt.outcome)
        if attempt.candidate is not None:
            print(attempt.candidate.packed.sequence)
            if attempt.candidate.final is not None:
                print(attempt.candidate.final.sequence)

quality = da.inspect(result, view="quality", limit=20)  # Limit displayed usage rows.
print(quality.cost.to_dict())
metrics = quality.to_dict()  # Compute exact aggregates for the declared population.
assert metrics["attainment"]["accepted"] == summary.accepted
# Draw the selected saved evidence without generating new sequences.
da.render(quality, view="library-quality", out="python-quality.png")
```

`attempt.candidate` contains the original packing and the last evaluated final
sequence, with their placements. Accepted, rejected and duplicate candidates
retain this evidence. An active-time stop before assembly has a packing but no
final sequence; attempts without a candidate, and older records lacking this
evidence, return `None`. Rejected candidates are not accepted designs and do not
appear in library exports.

`inspect(..., verify=True)` recounts the saved final checks, validates the packing
and assembly coordinates, and joins accepted candidates to their design records.
It reads recorded bytes without solving or drawing padding. Earlier padding
trials and the solver's optimality proof are not replayed by this check.

`AttemptFilter` accepts attempt ordinals, bare/full cell references, and outcome
codes. Its CLI flags are `--attempt-id`, `--cell`, and `--outcome`. Multiple
values within a field mean OR; fields combine with AND. It also applies to
diagnostics, whose report names the filtered attempt population while retaining
the original run target. A declared `dense_arrays.attempt-filter.v1` file can be
passed through `--selection`; filter files and convenience flags are exclusive.

Quality reports cover the designs selected by their `DesignFilter`, or all
accepted designs when no filter is supplied. They name
eligible and unused parts/groups, occurrence and design denominators, GC and
length distributions, placement counts, padding, and recorded requirement evidence.
`placement_count` counts all supplied part occurrences, including background parts;
it does not count detected motif hits. Reports identify this definition as
`library_composition.v2`. Saved reports with an unsupported metric policy remain
readable but are not compared numerically.
Density is the union of placement intervals divided by final length; compression
is summed part lengths divided by the pre-padding span. Positional occupancy
uses the number of designs reaching each position as its denominator. Empty
denominators are null with a reason. A missing recorded check remains missing;
inspection does not repeat screening.

`selection.designs` and `selection.distinct_sequences` describe the selected
population in `dense_arrays.quality.v3` reports.
`attainment` describes the original single run; filtering does not
change its target, accepted count or shortfall. Search outcomes cover all attempts
in the supplied native snapshots, including attempts outside the design selection.
For bundles, `search.availability` identifies complete, partial or absent attempt
history; unavailable counters are null. Each origin reports included and selected
design counts separately from its original attainment. Eligible
supply remains the parts and groups in the source plans, so filtered reports can
show more unused supply. Reuse the same filter for a report and a record export:

```python
selected_quality = da.inspect(
    result, view="quality", select=selected
)  # Reuse the export filter.
assert selected_quality.to_dict()["selection"]["designs"] == 1
# Draw the selected saved evidence without generating new sequences.
da.render(selected_quality, view="library-quality", out="selected-quality.png")
```

The quality CLI accepts the same `--selection`, `--design-id`, `--cell`,
`--part-id` and `--group` flags as record inspection. The same flags are available
on `render --view library-quality`. A precomputed Python report already binds its
filter and limits; rendering it does not replace those settings.

`limit` bounds displayed diagnostics or rows in each ranked part/group usage
table, including per-cell tables, not the
aggregate population. Quality usage tables support `after` cursors; every page
retains the same exact aggregates and recomputes the declared scan. Reports cache
their first successful computation in memory. `repr` and `cost` do no record scan.
Aggregate reports reject `all=True` because their displayed tables are bounded.

Reports charge source-plan rows plus each examined design/attempt; part/group
filter resolution also reads the relevant plans. Repeated source arguments are
checked again without multiplying the reported design or attempt counts. Quality
also caps retained lookup entries for source identities, distinct sequences,
distribution values, and occupancy boundaries. These are work/state limits,
not hard byte or time guarantees. CLI prints cost before computation and returns
an explicit error if a cap prevents an exact report. The quality render accepts
the same read-limit flags; graphics dependencies are checked before scanning.

The optional quality PNG shows ranked part use, GC, density, and search outcomes
from the same report. Its metadata contains the report and snapshot fingerprint.
Install the playback extra for rendering; text/JSON reports use the base install.
Quality also accepts [combined runs](../extension.md#assess-a-combined-library)
and [portable bundles](../bundles.md#combine-artifacts-and-assess-included-designs).
Use a [saved panel](../selection.md) to report and render the same
sampled identities at pinned source revisions.
