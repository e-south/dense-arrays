---
title: Continue an interrupted run
description: Resume a measured interruption with the same target, inputs and remaining budgets.
author: Eric J. South
---

# Continue an interrupted run

After a clean interruption, inspect the saved run and continue its original
request when `resumable` is true. Run these commands from the directory containing
`runs/first`, using the same Dense Arrays, Python and solver environment:

```bash
# Read saved evidence; --verify also checks its integrity.
dense-arrays inspect runs/first --verify --json
# Generate into a new directory, or continue the explicitly named run.
dense-arrays run --resume runs/first --json
```

The Python operation is equivalent:

```python
import dense_arrays as da

# Generate under the declared bounds into a new output directory.
result = da.run(resume="runs/first")
# Read saved run state and attainment. Recount stored evidence before returning.
summary = da.inspect(result, verify=True)
print(summary.state, summary.accepted, summary.target)
```

The run keeps its identity, target, plan, seed and committed designs. Interrupted
attempts remain counted. Continuation reconstructs packing exclusions from saved
accepted, duplicate and rejected candidates before searching again. A resumed
run can still stop below target when its original effort is spent or the offered
packing model is exhausted.

`resume` is exclusive with a request, output destination and design options.
To change the target, seed or limits, create a
[new request or extension](extension.md). Resuming a completed run verifies its
inputs and persisted records, then returns it without building a solver or
changing any run bytes.

## What can continue

Recovery supports single-cell and matrix runs whose recorded interruption is
resumable. Each matrix cell retains its original target, counters and accepted
designs. Completed, inactive, exhausted and failed cells stay closed; only cells
stopped by the interruption reenter generation. The rotation continues after the
last recorded attempt, including an unresolved interrupted attempt. Per-cell
attempt ordinals continue to advance, preserving their random-stream positions.

For an [ordered batch schedule](batches.md#search-an-ordered-schedule), exhaustion
closes the offered batch. Resume preserves that boundary and continues the saved
schedule. If its final batch is already exhausted, resume closes the cell without
another attempt. Backend failures and unresolved solver outcomes cannot be
converted into a retry by interrupting execution.

| Recorded situation | Result |
| --- | --- |
| Clean interruption, measured elapsed time, remaining attempts/time and complete packing evidence | Continue after verification and exclusive ownership. |
| Completed target | Verify and return the existing result. |
| Exhausted attempts or active-time allowance | Refuse; repeating resume cannot replenish either budget. |
| Terminal search exhaustion, backend failure or unavailable candidate evidence | Refuse; inspect the recorded outcome before making another request. |
| Abrupt process death with state still `created` or `running` | Refuse with `active_time_unknown`; the committed prefix remains available for inspection or a linked request. |
| Changed inputs or damaged records | Refuse; restore the original inputs or create a new request. |
| A different producer environment | Refuse further generation; use the original environment or a new request. Verification of completed runs remains available. |
| Another process holds the writer lock | Refuse with `writer_busy`; retry after that writer exits. |

Inspection's `resumable` flag describes the recorded stop. Resume also checks
the current input bytes, environment, integrity and writer ownership before
admitting work. Original package/runtime/backend version records stay unchanged.
They describe versions, not an attestation of unpublished source changes.

## Time, snapshots and ownership

The active-time allowance covers successful resume verification, model rebuilding,
exclusion replay and generation. Time between executions is excluded. Checks are
cooperative; backend or cleanup work can exceed the allowance. An abrupt crash
does not provide a trustworthy measurement of its last active interval, so it
cannot silently reset that interval to zero.

Saved design records and earlier revisions are immutable. Repeated resume of a
completed run is idempotent. A clean interruption during continuation creates
another measured stop, retaining all earlier effort. Future tied solver choices
and final bytes are not promised to match an uninterrupted execution.

Recovery verifies the complete saved prefix and streams its packing records;
its work grows with the saved attempts and designs. Matrix replay scans the
history once and rebuilds only models for cells that will continue. Scheduled
batch replay uses two streaming passes to recount positions and restore only
unfinished batches. Rebuilding
and replay consume the same remaining active-time allowance as generation. Inspection remains available
with explicit read limits when only a summary or selected evidence is needed.

The writer holds the existing lock file for its entire lifetime. Do not delete
or replace that file to bypass a live writer. Qualification covers process
interruption on local macOS storage. Distributed filesystems, simultaneous access
through a synchronization service and power-loss durability require separate
qualification. External deliveries remain outside the native run transaction.

CLI recovery refusals use exit `4` with a structured code under `--json`.
Conflicting resume options use exit `2`; a continued run stopped below target
uses exit `3`, and interruption uses exit `130`.
