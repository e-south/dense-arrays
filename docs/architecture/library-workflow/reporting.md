---
title: Library reporting contracts
description: Interpret diagnostics and quality metrics with declared populations, work limits and reader lifetimes.
author: Eric J. South
---

# Library reporting contracts

Reports explain saved evidence without repeating preparation or generation.
Use [quality tasks](../../library-workflow/results/quality.md) for runnable
examples and [outputs](../../reference/outputs.md) to choose a figure or summary.
The examples below read the run created by the
[constrained generation recipe](generation.md#constrained-design). Rendering
requires the `playback` extra.

## Diagnostics

Diagnostics are records, shared by human text and Python/JSON: stable `code`,
`severity`, `stage`, requirement/input reference, observed/expected values,
evidence references, and suggested next action. Final sequence violations include
zero-based half-open match coordinates, strand, and intersecting placement/pad
intervals. Bounded output states how many additional diagnostics were omitted.

Distinguish three levels of evidence: a statically proven contradiction, a
packing-model proof, and an observed search bottleneck. Only the first two can
support their respective infeasibility claim. Repeated failures involving a part
do not prove that part caused infeasibility. Recommendations never edit requirements
or trigger a solver run. Minimal conflicting sets and automatic relaxation are
not promised.

For example, requesting three occurrences from a group with only two eligible
part IDs is a static contradiction when each ID can occur once. Planning
identifies the requirement and the requested and available counts before solving.
A literal match across a placement/padding join belongs to the `avoid`
requirement. An exhausted attempt budget reports `attempt_limit`; it does not
blame a constraint without evidence.

```text
dense-arrays inspect runs/curated --view diagnostics --limit 20
dense-arrays inspect runs/curated --view attempts --outcome rejected --limit 10
```

```python
import dense_arrays as da
from dense_arrays import reporting

result = "runs/curated"  # A saved run from the generation recipe.
diagnostics = da.inspect(result, view="diagnostics", limit=20)
rejections = da.inspect(
    result,
    view="attempts",
    select=reporting.AttemptFilter(outcomes=("rejected",)),
    limit=10,
)
```
## Quality reports

Every run requires a quality report, including empty and partial libraries.
`inspect --view quality` and its typed Python equivalent expose the same metrics;
no plotting extra is needed to understand them. Native completion still means
target attainment, not a judgment of experimental merit. No arbitrary quality
threshold silently changes acceptance.

| Required metric | Definition / interpretation |
| --- | --- |
| Attainment | Requested and accepted designs per cell; distinct final sequences reported separately. |
| Supply and use | Eligible part IDs/groups, selected occurrence counts, designs containing each part/group, and unused eligible parts. Both occurrence and design denominators are named. |
| Concentration | Ranked part-use shares and highest single-part share of selected occurrences, separately per cell and combined. Report the eligible reference pool; do not call this a sequence-diversity guarantee. |
| Geometry/composition | Final length and GC distributions; placement count; positional occupancy from supplied placements; amount and side of padding. |
| Packing density | Fraction of final bases covered by the union of selected placement intervals. Overlaps count once; named fixed parts are included. Also report pre-padding span and an explicitly separate sum-of-part-lengths/span compression ratio. |
| Search loss | Mutually exclusive accepted/rejected/duplicate/no-candidate/error outcomes and effort, plus multi-label reason totals. Include duplicates against a referenced parent separately. |
| Requirement evidence | Accepted requirement checks, explicit relaxations, proof scope and missing evidence. Unknown values are not zeros. |

Metrics carry a version, snapshot/selection identity, population definition,
denominator, and exact/sampled/not-computed status. Empty denominators yield
null with a reason. Default aggregate metrics cover the declared population,
not just the displayed page. Detailed per-part tables remain paginated. Library reports do not compute
all-pairs final-sequence distances. Preparation MMR reports record sequential
core distances within their declared selection process. Part-use balance,
core diversity and final-sequence distance describe different populations.

The rendering extra supplies named `library-quality`, `preparation-quality` and
`design` views. Quality views visualize the same reports available to inspection;
designs show selected placements, orientations, coordinates, padding and supported
requirements.
Rendering adds no acceptance calculations. Palette/title changes cannot change
metrics. [Preparation figures](../../library-workflow/preparation-quality.md)
show yield, recorded sequential MMR distances and declared score bands, separately
for each recipe. Missing metrics stay unavailable. Live pool report verification
replays MMR selection under the declared pair-work allowance; detached reports
reuse their saved metrics without that computation.
Study-specific figures and cohort interpretation remain external.

```text
dense-arrays render runs/curated --view library-quality --out quality.png
```

```python
import dense_arrays as da

visual = da.render("runs/curated", view="library-quality", out="python-quality.png")
```

`library-quality` explicitly selects the aggregate report over its declared
snapshot. A `design` view of a multi-design run requires design selection; it
cannot silently choose its first record.

## Cost and reader lifetime

| Operation | Cost contract |
| --- | --- |
| Unfiltered summary | Read bounded manifest/committed summary data; no full design scan, solver or optional imports. |
| Record iteration | Stream requested fields in bounded chunks; filtering may scan beyond a displayed page and reports records examined. No eager dataframe. |
| Filtered quality | Declare an exact scan or compatible derived index. One page of records is not the denominator. |
| Verification | State the artifact/evidence boundary and bytes/records checked; cost grows with that boundary, not the displayed limit. |
| MMR selection verification | A pair-work cap bounds distance recomputation on live pool reports. Detached reports use recorded distances. |
| Selection/export | Disclose candidate scan and requested output scope; copy only selected records and required evidence. |

Reports and record views expose a cost descriptor before iteration: source
snapshot, indexed/scan mode, known record/byte estimates, selected projection
and any explicit work cap. Unknown estimates stay unknown. CLI emits the same
descriptor before expensive work; it does not require interactive confirmation.
A caller may set a record/pair work bound. Reaching it returns explicit
not-computed/partial evidence or an error as appropriate; never relabel a
sampled or truncated result exact. Operational limits do not prove infeasibility.

The shared `read_limits` argument accepts `reporting.ReadLimits` with positive
integer `records`, `pairs` and `identities` caps. CLI equivalents are
`--max-read-records`, `--max-pairs` and `--max-identity-entries`. These apply to
inspect/export/render work, not generation attempts or output page size.
Omitted caps resolve to documented defaults and appear in the
cost descriptor. Unsupported cap enforcement fails before work; exporting a
partial stream after a limit still returns a nonzero status and pure data stdout.

A `RecordView` holds descriptors, not an open file or eagerly loaded rows.
Each `records()` call opens an independent iterator over the same revision
and repeats the declared page or full query; it never advances shared state.
The iterator is a context manager; exhaustion closes resources, and early-stop
callers close it explicitly. It does not follow new commits. Cursors bind the
query, ordering, schema and revision; changing any of them invalidates reuse.
Native runs retain committed revisions for the life of the run and
does not offer automatic pruning. Missing externally removed evidence fails
visibly; no implicit cursor reset or repair is permitted.

```python
import dense_arrays as da

view = da.inspect("runs/curated", view="placements", all=True)
with view.records() as records:
    for placement in records:
        print(placement.design_ref, placement.start, placement.end)
```

Row buffering is bounded; exact multi-source deduplication and a materialized
selection may need additional identity state. Declare and cap that state,
or use a qualified committed index. Hitting the cap fails explicitly before
publication; do not claim total constant memory merely because rows stream.
Selection snapshots expose bounded summaries and lazy reference access, not an
unbounded notebook representation. Inspect creates no spill files.

Derived indexes declare source revision and schema. Stale indexes cannot serve
a query; an explicit scan can provide the same result with its cost disclosed.
Rebuilding an index is a writer-owned operation, never a side effect of inspect.
Native records remain authoritative. Benchmark model construction, generation,
inspection, verification, selection and export separately before choosing caches,
parallelism or another storage backend. Measurements determine documented
default caps; unsupported hard memory/time guarantees must not be advertised.

## Rendering boundary

Per-design rendering requires a design or explicit bounded selection; a
multi-design artifact never silently chooses its first row. The library-quality
view instead renders the declared aggregate scope. Rendering consumes persisted
records and reporting metrics, never adds acceptance calculations or reruns
generation. Playback evidence describes the order reconstructed from placements. It does
not record a solver trace. Unsupported geometry/requirements fail before publication.
Visual failure leaves native results valid and preserves existing per-file
publication qualifications.
