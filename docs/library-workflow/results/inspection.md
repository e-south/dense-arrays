---
title: Bound record inspection
description: Read saved records with explicit work limits and revision-bound continuation.
author: Eric J. South
---

# Bound record inspection

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

```bash
# Read or verify saved evidence without generating again.
dense-arrays inspect pools/curated --view parts --group A \
  --limit 10 --max-read-records 1000 --json
```

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
