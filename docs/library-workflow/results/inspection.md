---
title: Inspect saved records
description: Read saved designs, page through results, and choose limits for larger queries.
author: Eric J. South
---

# Inspect saved records

Use the run from [Generate a saved library](../../library-workflow.md). Read a
page of accepted designs:

```bash
# Display up to ten designs while allowing at most 1,000 decoded records.
dense-arrays inspect runs/first --view designs \
  --limit 10 --max-read-records 1000 --json
```

## Bound record inspection

Record views expose `cost` before opening their iterators. CLI inspection prints
the same descriptor to stderr before streaming rows. It identifies the source
revision, projection, indexed/scan access, record estimate and resolved work caps.
Byte estimates remain unknown. Each iterator reports its own `examined` and
`returned` counts; rejected filter rows count toward examined work.

`--limit` controls displayed rows. `--max-read-records` controls data records
decoded, including filtered-out rows. The initial default is 100,000 records;
`--all` removes pagination, not the read bound. Set a larger bound explicitly for
a larger scope. Identity lookup state is capped at 100,000 entries by default
(`--max-identity-entries`). Record views do no pairwise work; the reserved
`--max-pairs` cap does not request a pairwise calculation.

These caps bound record decoding and identity state, not bytes in a single row
or database metadata/index I/O. A filtered page can require a scan of the full
pool. Exceeding a work cap is an explicit failure, never a complete smaller
result. A streamed response may contain a prefix before that failure; check its
exit status. Verification uses the same explicit read limits and checks the
complete declared evidence boundary. Its `verification_cost` describes the
additional work before the scan; exceeding a cap fails verification.

A completed record page includes `next_cursor` when another page may exist.
Pass it as `--after TOKEN` (Python `after=token`) with the same view and filter.
The token retains the source revision even if generation subsequently advances.
A full final page can return a token whose next page is empty; a short final
page returns null. Closing an iterator early preserves a cursor after its last
returned row. Changed queries, substituted sources, and missing revisions fail
explicitly.

## Read a stable revision

A record view opens its reader when you call `records()`. Each call creates an
independent iterator over the same committed revision. Use `with view.records()`
when stopping early; exhausting the iterator also closes it. A Python
`artifacts.RunHandle(path, run_id, revision=N)` pins inspection, export and
rendering to that revision while generation continues.

`inspect --verify` checks checksums, coordinates, part identities, requirements,
design/attempt joins and reconciled counts. Its `verification_cost` describes
the scan before execution. The returned `verification` identifies the checked
record families and counts decoded records and UTF-8 JSON bytes. These are
logical record costs rather than physical database I/O; verification checks
saved results rather than reproducing generation.

## Read software versions

Run and pool summaries retain the Dense Arrays, Python and OR-Tools versions,
operating-system family and machine architecture. A run also records the solver
name and version once its model is built. Use `summary.producer.to_dict()` in
Python or `inspect --json` to read the full producer record.

These values help compare execution environments. They identify reported
software versions, not unpublished source edits or deterministic solver tie order.
For search and random-stream behavior, see [packing search](../search.md).
