---
title: Dense Arrays hardening plan
description: Ordered contract, usability, and module-boundary work with measurable completion criteria.
---

# Dense Arrays hardening plan

Author: Eric J. South. Status: proposed implementation sequence, following the
[September 2026 audit](audit.md). The documentation and routing pass is applied;
the runtime slices below remain open. This plan is ready to become a separate
execution goal when selected.

Success means a user can formulate a packing problem, distinguish a valid result
from a solver failure, and interpret saved placements without hidden caveats.
A contributor should find the owner and relevant tests from one task route and
change that behavior without editing unrelated solver or renderer code.

## Scope and constraints

Keep the package's public purpose: motif packing, explicit realized placements,
and playback. Preserve all Virgile Andreani credits and joint authorship;
attribute new modules and documentation to Eric J. South. Retain existing
public entrypoints where their behavior is sound. Changes to invalid-input
acceptance must be deliberate and documented.

Use local CBC for bounded positive and negative cases. Do not run Gurobi on
this workstation. A BU SCC Gurobi check is optional follow-up only if a
backend-specific question remains: authenticate interactively, use the host's
approved job route, bound runtime and output, and remove transient work after
preserving the result. It is not a prerequisite for docs or core-contract work.

Exclude scientific claim expansion, new study dependencies, exact solver-trace
formats, framework replacement, and bulk rewrites. Producer packages own study
labels and adapters. Inventory consumers before changing an interface; migrate
them explicitly rather than adding compatibility shims or private imports.

## Delivery sequence

Each slice is independently reviewable. Add failing tests for the intended
contract, implement the smallest change, then run focused checks and the
[full gate](../development.md#local-verification). Update current-behavior docs
when the runtime changes; do not mark later work complete in advance.

1. **Preserve solver failures and reject unsupported approximation.** Address
   audit C1/C2 and backend-option handling. Distinguish infeasibility, feasible
   but unproven results, backend failure, and invalid result reconstruction.
   Only infeasibility may end enumeration normally. Reject configured hard
   requirements at `approximate()` entry until that method can honor them.
   Test all outcomes, including failure after one yielded result, through the
   iterator and CLI. Keep CBC success and genuine infeasibility controls.

2. **Validate configuration and result identity before mutation.** Address
   C3–C5. Normalize integer counts/lengths, interval shape/order, motif alphabet,
   enums, and IDs without lossy coercion. Hold an immutable problem snapshot;
   rejected bias, weight, or forbid operations leave state unchanged. Settle
   entry/orientation semantics and share result validation. Test malformed
   inputs before model allocation, caller-list mutation, foreign results,
   boolean indices, nonfinite weights, and failed-update atomicity.

3. **Validate saved plans before rendering.** Address P1 and direct realized
   contracts. One semantic validator must serve Python and JSON boundaries.
   Check every nested shape, reference, sequence/span agreement, unique ID,
   reveal span, and constraint evaluation. A valid failed constraint remains
   representable; a contradictory `passed` flag is rejected. Turn the nine
   audit mutations into negative controls and assert zero output for invalid
   CLI input. Define supported v1 authority/relation combinations explicitly.

4. **Make playback qualifications visible.** Address P2/P3. Derive concise
   authority/order/constraint text and permitted relations from the plan once,
   then render it in HTML and raster output. A layout-only scene must not
   animate a complete placement chain. Support presentation fields consistently
   or reject unsupported settings; carry explicit caller-owned labels and
   colors instead of inferring biology from names. Check actual DOM text,
   artists, labels, and active geometry for unique, ambiguous, gapped, and failed
   cases. Include a still and narrow-screen keyboard/reduced-motion review.

5. **Separate responsibilities behind the public interfaces.** After the
   relevant contracts pass, extract solver model construction and greedy
   occurrence realization from `Optimizer`; extract document/presentation data,
   frame scheduling, drawing components, and writer lifecycle from the raster
   renderer. Split tests along those same owners. Correct graph selection parity
   across default and injected engines. Keep semantic projection independent of
   raster imports and playback independent of eager solver imports. Test those
   import boundaries in isolated subprocesses. Narrow lint exceptions as each
   changed module earns the checks; do not build a generic plugin system.

6. **Finish command and documentation ergonomics.** Add concise failures and
   output preflight to both CLIs. Reject input/output aliasing and colliding
   outputs; define whether a requested export set is published together or as
   explicitly reported partial results. Share media scheduling and clean up
   figures on failure. Verify timing domains and optional dependencies before
   writing. Add maintained smoke examples and built-link checks through existing
   test/CI tooling. Repeat fresh-reader tasks from root and nested directories,
   check plain Markdown and the rendered site, and review the final API surface.
   Normalize module-header and docstring formatting while retaining the exact
   author credits; document accepted inputs, returns, and errors where generated
   signatures alone do not explain them. Verify the publication route before
   claiming that a local documentation build is the live site.

Slices 1–2 and 3 can run in parallel after assigning file ownership. Slice 4
depends on 3. Slice 5 follows its subsystem's contract repairs; slice 6 closes
the user workflow after those interfaces settle. Keep at most three workers
and one integration owner, with no recursive delegation by default. Separate
worktrees or nonoverlapping edit ownership prevent shared-file races.

## First implementation slice

**Goal:** an unsuccessful solver call cannot masquerade as exhausted enumeration,
and a constrained problem cannot receive an unconstrained greedy result.

**Touched owners:** `optimizer.py`, `cli.py`, their focused tests, and the
optimizer/constraints/CLI references. Do not move graph or renderer code.

**Done criteria:** normal CBC results remain valid; genuine infeasibility ends
enumeration; abnormal/not-solved/unproven statuses retain distinct meaning;
invalid result reconstruction propagates; failure after a yielded result exits
nonzero at the CLI; unsupported approximation requirements fail before work;
rejected solver options fail explicitly. No silent backend substitution.

**Verification:** independent coordinate/coverage checks for small CBC cases,
isolated solver-status injection for error propagation, a real CLI negative
case, focused tests, then the full gate. Test doubles establish boundary behavior,
not backend reliability. Preserve exact author-credit text in touched files.

**Handoff:** record changed behavior, commands/results, consumer impact, and
remaining audit IDs in the PR. Retain existing successful examples. Replan if
the error contract requires a public return-type change or a new dependency.

## Decisions to settle within their slice

- Keep exact path-entry packing as the baseline. Define how greedy occurrence
  selection treats repeats and contained motifs before repairing it; do not
  quietly change the optimization objective to substring coverage.
- Define allowed negative spacers separately from interval type validation so
  an input cleanup does not accidentally forbid intentional overlap.
- Decide whether nested provenance is an immutable JSON snapshot and which
  renderer-neutral document fields each backend implements. Test that declared
  contract instead of adding broad fallback behavior.

These decisions do not block the first slice. Escalate if a consumer needs a
conflicting contract, two fix attempts fail, or the scope expands beyond the
named owner. Keep proposed exact traces and new biological interpretation out
of this effort.

## Completion evidence

- Every priority-1 finding has a regression test and demonstrated fix; remaining
  priority-2 items are fixed or explicitly scoped and accepted before completion.
- CLI failure and success paths are distinguishable without reading a traceback.
- Python/JSON inputs share the same invariants, and failed mutations are atomic.
- HTML and still/video views preserve visible qualifications and caller labels.
- The large owners have semantic boundaries with focused tests; no duplicated
  export loop, alternate contract source, or study-specific inference remains.
- A fresh reader completes all [four routing tasks](documentation.md#dogfood-the-route),
  and a contributor identifies the edit/test scope without loading all docs.
- Examples, strict docs build, internal links, lint, tests, dependency audit, and
  distributions pass. Visual/accessibility evidence and author-credit preservation
  accompany the PR. Record any remote test separately from local CBC evidence.
