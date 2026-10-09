---
title: Library workflow artifacts
description: Persist versioned evidence, reconcile run accounting and recover committed work.
author: Eric J. South
---

# Library workflow artifacts

A run owns committed execution evidence. Native schemas, manifests and
recovery rules preserve its identity, accounting and accepted designs.
[Reporting](reporting.md), [selection](selection.md) and [exports](exports.md)
read those snapshots without changing the source records.

## Native evidence

Native artifacts use versioned logical records for pools, preparation and
generation plans, realized designs, attempt outcomes, selection snapshots and
manifests. Runtime stores and JSON exports carry those contracts without
creating duplicate authoritative state. A manifest references a coherent committed generation
of records. Columnar exports and indexes are derived projections.

Plans, scientific design records and display configuration have separate identities.
A title or palette cannot change a sequence digest or search identity. Canonical
encodings declare semantic ordering, included fields, algorithm and namespace
version; timestamps and machine paths do not enter semantic identity.
One `RealizedArray` owns each design's sequence/placements.

Inspection pins a committed revision and includes it in reports and cursors.
Reject an expired/missing revision rather than mixing pages across snapshots.
An active-run reader may see a coherent older revision; it never acquires writer
ownership or repairs state. Reports include schema version, operation, artifact
identity, scope, completion/verification status, diagnostics and data.

`inspect --verify` checks available bytes, schemas, digests, references,
placement alignment and accounting at the declared artifact boundary, regardless
of displayed page size. It never reruns a solver, fetches missing evidence or
claims biological validity. A cached valid flag is not evidence of integrity.
Missing evidence and unsupported versions receive distinct diagnostics.

## Schema compatibility

Each artifact family declares its reader support and writer version independently;
do not infer compatibility from the package version or file extension. Reader
and writer support names exact schemas. Realized-array and playback contracts
retain their own declared versions and geometry capabilities.

Unknown fields outside declared metadata namespaces, duplicate mapping keys,
unknown policy versions and
unknown schema versions fail with the observed value and supported versions.
Caller metadata is permitted only in explicitly designated namespaces; accepting
an extra metadata key must not add a constraint or override a core field.
Readers never apply today's defaults to a previously resolved plan.

Saved library-quality reports preserve an unknown metric policy as recorded
values. Comparisons mark those metrics incomparable, and rendering rejects the
unsupported policy. Preserving a report does not assign current metric semantics
to its values.

`export` writes supported native records, projections and bundles in the chosen
format. Receipts bind the source revision, selected identities, format, output
counts and file digests. Record exports also bind each source to its canonical
committed manifest digest and recheck it before publication; later run commits
do not change a pinned revision. These bounded metadata reads are separate from
the data-record scan allowance. Unsupported view/format combinations fail before
publication. Export preserves its source artifact and the declared native schema.

Wire encodings, supported versions and conformance fixtures define each native
artifact family. Readers accept declared supported schemas and reject
unrecognized formats.

## Accounting and completion

A candidate attempt is one request for the next design from an offered batch,
followed by evaluation or a terminal no-candidate/error outcome. Solver calls,
batches, attempts, designs, and output sequences have separate identifiers and
counts. Multiple constraint violations may annotate one rejected attempt; they
do not multiply the rejected-design count.

For a committed snapshot:

```text
started attempts = accepted + rejected + duplicate + no_candidate + error
                   + interrupted_unresolved + in_progress
```

Categories are mutually exclusive current outcomes. Resume changes a persisted
`in_progress` outcome to `interrupted_unresolved`; retrying starts a new attempt
ID linked to that interruption. Do not rewrite the original as a success or
reuse its random-stream identity. Reserve each attempt durably before work,
including requests that yield no candidate. Reason-level totals are clearly
labeled because a rejected sequence can have several violations.

Run state is `created`, `running`, `completed`, `stopped`, or `failed`.
`completed` means every declared target is met. A valid partial library can be
`stopped` by a limit, operator interruption, or established search exhaustion.
Execution/integrity failures are `failed`. Report per-cell attainment, stop
reason, proof scope, and whether continuation is supported. An interrupted or
failed run can retain a valid committed prefix without claiming completion.

In Python, anticipated search outcomes return structured results; malformed
requests and execution failures raise typed errors. A failure after publication
carries the run reference so callers can inspect its committed prefix. Keyboard
interrupt remains an interrupt, with best-effort checkpointing. The low-level
optimizer has its own typed exception contract.

Retry admission is an explicit policy over typed reasons and remaining budgets.
Invalid solver results, corrupt committed records, and unsupported capabilities
are failures, not reasons to resample. Infeasible batches and screening rejection
can continue only in a mode whose declared search policy permits that next step.
Never catch all exceptions and convert them into a successful partial library.

| Observation | Attempt outcome | Exact-mode response |
| --- | --- | --- |
| Proven packing result passes all final checks | `accepted` | Commit; continue until target or limit. |
| Final checks fail, or final sequence already exists in the cell | `rejected` or `duplicate` | Record reasons; request another candidate within limits. |
| Initial batch infeasible, or enumeration exhausted | `no_candidate` | Stop below target with the corresponding batch-scoped proof. |
| Feasible but unproven, or backend explicitly reports its limit | `no_candidate` | Stop below target; preserve the distinct proof and termination fields. |
| Unknown backend failure, invalid result, or integrity failure | `error` | Fail with an inspectable prefix; do not automatically retry. |

Sampled modes continue to another batch only under their explicit
policy. They do not weaken required proof for an accepted design.

CLI exits: `0` for a completed requested operation, `2` for invalid
usage/request, `3` for a run stopped below target, `4` for dependency/execution/
integrity failure, and `130` for interruption. A successful `inspect` exits `0`
even when it reports a stopped run; `inspect --verify` fails on invalid integrity.
Code `3` also covers a preparation operation stopped short of its retention target.
With `--json`, domain failures emit a versioned error/result envelope to stdout
and diagnostics to stderr. Parser errors before output-mode selection may remain
stderr-only and must have a documented nonzero exit.
## Resume and randomness

Retain single-writer workspace locking and input/configuration digest checks.
Native runs support one host and qualified local filesystems. It does
not claim distributed locking or safe simultaneous writers through a sync
service. Use atomic destination creation and an OS-held lock with a stable file
identity; do not unlink a lock path while another process may hold or acquire
its inode. Stale metadata is diagnostic, not authority to break a live lock.
Resume verifies the complete semantic plan and bound inputs before writing.
Do not share virtual environments, run roots, or accepted-output locations
between comparison executions.

Resuming a completed run verifies it and returns the existing result without
generating more designs. Resuming a stopped/failed run is permitted only when
its recorded reason is recoverable and its budgets remain available. The
inspection report states this explicitly. Exhausted budgets, changed inputs,
and damaged committed data fail with a specific recovery route; they cannot
be bypassed by repeating `--resume`.

A process can die while the saved state still says `created` or `running`.
After acquiring exclusive ownership, resume classifies that abandoned execution
from its committed records and applies the same interruption/integrity/budget
rules. It must not require a final state that the dead process could not write,
or infer ownership merely from a recorded PID.

The minimum guarantee is preservation of every committed design, no duplicate
publication, coherent quotas/accounting, and explicit handling of unfinished
work. The persistence implementation must commit records and state through one
recoverable frontier; file existence alone is not a commit. Physical storage
layout is an implementation choice to validate with crash injection.

A commit couples accepted records, attempt outcomes, uniqueness/quotas, and
consumed-budget state. Derived indexes may be rebuilt only from that committed
source of truth. Recovery covers process interruption on qualified storage; power-loss
durability and distributed storage require separate evidence. On recovery, uncommitted staging never contributes accepted designs.
If integrity cannot be established, keep the artifact inspectable and fail
resume without silently repairing or deleting evidence.

Published native artifacts are the handoff boundary. Reading them in another
application does not change the run's completion or ownership. A downstream
consumer owns any delivery, retry or storage protocol it adds.

Random streams derive from the root seed, algorithm version, cell, batch, and
attempt identities, with a specified stable derivation algorithm. No built-in
process-randomized hash or worker-completion order participates. Preparation
tie-breaking and representative-selection policies are versioned separately.
Changing a stream definition requires a new policy version and explicit
reproducibility evidence.

Exact continuation of a live solver iterator is not promised. Interrupted
search may restart from its last committed batch boundary with duplicate
suppression. Sampled work can be deterministic while tied solver solutions
differ. Byte-identical future sequences require a separately demonstrated
environment/backend/ordering guarantee. Resource limits count retried work;
checkpoint state cannot reset consumed budgets. Interrupted time that cannot
be recovered is labeled unknown or conservatively bounded, never fabricated.
If remaining active-time allowance cannot be established after a crash, fail
resume and offer a new linked run. Preserving a committed prefix does not imply
that every crashed run is continuable under its original limits.
`limits.attempts` and `limits.active_seconds` apply across the whole run;
`limits.solver_seconds` applies to one solve call. Paused time is excluded from
active time. Count active runtime from input staging through finalization,
including non-solver work; report conservative charging separately from measured
duration. Admission checks enforce the attempt count before work. Active-time
checks are cooperative at documented operation boundaries; clamp the requested
backend timeout to the remaining allowance and stop admitting work once spent.
Cleanup may outlast that allowance. Report the observed limit(s), checks, and
overrun rather than asserting which limit happened first inside opaque code.

These time controls are not a hard process deadline. If a caller requires one,
reject that unsupported capability until supervised cancellation has been
qualified. Do not add a process supervisor merely to rename a cooperative limit.
Resource caps are operational limits, not feasibility evidence. Seed equality
does not promise identical output under timing-dependent stopping conditions.

## Diagnostics

See [diagnostic evidence](reporting.md#diagnostics).

## Quality reports

See [populations and metric definitions](reporting.md#quality-reports).

## Result selection

See [filters and source unions](selection.md#result-selection).

### Bounded selection

See [counts, quotas and saved panels](selection.md#bounded-selection).

### Tabular and sequence output

See [record formats and publication](exports.md#tabular-and-sequence-output).

### Portable bundle

See [self-contained selected evidence](exports.md#portable-bundle).

## Cost and reader lifetime

See [read limits, iterators and cursors](reporting.md#cost-and-reader-lifetime).

## Rendering boundary

See [visuals from persisted evidence](reporting.md#rendering-boundary).
