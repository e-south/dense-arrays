---
title: Library workflow qualification
description: Qualify native DNA library capabilities through reproducible tasks, scientific checks and installed-package evidence.
author: Eric J. South
---

# Library workflow qualification

A usable library workflow takes declared parts and requirements through bounded
generation, explains its outcomes and publishes results with traceable sequence
and placement identities. These gates establish that behavior through both
Python and the CLI. The [overview](index.md) routes to the domain, operation and
artifact contracts; the [native guides](../../library-workflow.md) describe
available interfaces.

## Workflow coverage

The complete path is curated parts → constrained designs → diagnostics and
quality → revision or extension → reproducible selection and portable handoff.
Recovery, matrix generation, PWM preparation and background preparation each
require additional method-specific evidence. A successful packing run alone
does not qualify these other capabilities.

| Required capability | Contract |
| --- | --- |
| Bounded counts, named geometry, final checks and proof scope | [Domain](domain.md#requirements) |
| Mapped inputs and reusable curated pools | [Part inputs](domain.md#part-inputs), [curated recipe](preparation.md#curated-pool) |
| Preview and generation through either interface | [Operations](operations.md#shared-operations) |
| Diagnoses, quality metrics and visible shortfalls | [Reporting](reporting.md#diagnostics) |
| Immutable extension and editable requests | [Revision](generation.md#revise-or-extend-a-library) |
| Fixed-size selection and portable publication | [Result selection](selection.md#result-selection) |
| Typed results, bounded iteration and discoverable tasks | [Python](operations.md#python-results), [cost](reporting.md#cost-and-reader-lifetime), [help](operations.md#help-and-discovery) |
| Recovery, matrices and sampled preparation | [Recovery](../../library-workflow/recovery.md), [matrices](../../library-workflow/matrices.md), [preparation](../../library-workflow/preparation.md) |

## Change discipline

Begin each behavior change with a failing contract check, implement the smallest
complete path and verify it before expanding scope. Structural refactoring,
algorithm changes and external integration have different evidence requirements;
review them separately. A module move must not silently change selection,
random draws, quotas or identifiers.

One module owns each policy. Use typed solver outcomes and supported controls;
never infer termination from exception text or patch private solver internals.
Core uniqueness, complete-sequence uniqueness, part-use balancing and library
selection remain distinct policies. Version random streams and representative
selection independently so recorded provenance can explain changed results.

Use isolated worktrees, environments and output destinations. Recheck live state
before writing; a clean checkout does not prove that no execution is active.
Fixtures must be immutable synthetic or explicitly authorized data. Development
checks must not start study jobs, modify scientific results or change another
project's environment.

## Capability qualification

| Capability | Required evidence |
| --- | --- |
| Table inputs | Explicit mappings, identity preservation, typed metadata, normalization and row/field diagnostics for each supported format. Optional readers remain optional. |
| PWM preparation | Declared motif/window and scorer, score units and background, core/sequence uniqueness, score or MMR retention, independent effort/target bounds, versioned selection rules and reconciled reports. |
| Background preparation | Exact conditional sampling under declared GC and literal rules, separate FIMO exclusions, exhaustive small oracles, distribution checks and distinct zero-mass/resource-limited outcomes. |
| Candidate batches | Bound ordered part identities, executable replay without drawing again, explicit resampling and feedback policies, per-batch and run-wide limits. |
| Matrices | Bounded expansion, named cells, complete allocation, explicit zero targets, independent cell identities and shared execution effort. |
| Extension | Additional targets, ancestor exclusions, immutable parent records, explicit lineage and semantic plan comparison. |
| Reporting | Coverage, concentration, occupancy, yield and declared diversity metrics with correct populations, denominators and work limits. |
| External handoff | Documented native export schemas, exact identity/coordinate joins, create-only publication and portable reads with original sources unavailable. |

Each capability needs implementation and test pointers in the [code map](../README.md),
working public examples and the applicable acceptance evidence below. Command
names and matching field names do not establish semantic equivalence.

## Acceptance by capability

Each capability must produce a usable, inspectable result. Review its behavior
through the domain owner and the acceptance evidence below.

| Capability | Reviewable outcome | Acceptance gate |
| --- | --- | --- |
| Native contracts | Shared vocabulary, six operations, identities and complete first-use examples. | Public types, wire encodings and behavior requirements are explicit. |
| Curated generation | One-cell planning, maximum length, occurrence/group requirements, bounded exact enumeration, placements and attempts. | Real CBC, equivalent Python/CLI semantics, coherent commits and bounded termination. |
| Assembly and presentation | Named fixed occurrences, exact counts/length, padding, final screens and rendering. | Independent identity/strand/coordinate checks, junction/padding fixtures, persisted-evidence rendering. |
| Reporting and revision | Mapped pools, diagnostics, quality, typed readers, request edits, extension and portable exports. | D1–D3, D5–D6, stable joins, parent preservation and explicit costs. |
| Reproducible selections | Total/per-cell first or seeded-random policies, saved membership and shortfalls. | D4, pinned revisions and stable identities across inspection, export and rendering. |
| Recovery | Locks, committed-prefix preservation, interrupted attempts and consumed-budget reconciliation. | Competing processes and crash injection; unsupported recovery states fail explicitly. |
| Matrices and sampling | Allocation, replay, resampling, quotas and coverage/failure weighting. | Independent cell and batch fixtures, uniqueness scope and bounded work. |
| Sampled preparation | PWM and background recipes, scoring, retention and complete pool evidence. | Small independent oracles, real-tool scoring checks, missing-tool preflight and memory/time measurements. |
| Distribution | Independently installable base and optional features with discoverable examples. | Built-wheel tasks, platform checks and the full development gate. |

Exclusive output ownership and coherent publication apply to every capability. Background preparation
must meet the [conditional-generation method](../../development/background-generation.md)
and its exhaustive validation cases. Preparation resume, additional geometry
and new numerical methods require their own acceptance evidence before support
is advertised.

## Product acceptance demonstrations

These demonstrations qualify complete user tasks alongside the runtime tests.
Execute them against an installed distribution and retain the commands, inputs,
results and environment. A written demonstration is an acceptance requirement,
not evidence of a passing implementation.

### D1. Express a constrained design through either interface

From installed help, import the same mapped synthetic table in CLI and Python;
require exactly two occurrences from A, coverage of A/B, the named fixed pair, final
length and screening rules in the [constrained recipe](operations.md#constrained-design). Produce the same normalized plan
and diagnostics. Compare accepted constraints/objectives where solver ties
prevent sequence equality. Check every placement and requirement on final bytes.
Use an additional fixture with three eligible A identities to prove that an
upper bound changes acceptance. Verify minimum-only, maximum-only, equal bounds,
zero maximum, overlapping selectors, repeated sequences with distinct IDs, fixed
parts exceeding a bound, and reversed/boolean/fractional bounds. Independent
placement recounts must agree with CBC-enforced requirements; sampling caps are
not sufficient evidence. Unsupported heuristic combinations fail planning.
Then add an impossible count and a malformed source row: both must identify
their source/rule before solving. Use public native requests and imports throughout. Record the actual commands, wrong turns and prerequisites.

### D2. Explain a shortfall and assess the accepted library

Use a deterministic controlled fixture yielding eight accepted designs toward
twelve before its attempt bound, with at least one duplicate and one final
screen rejection. Inspect the same persisted outcome through CLI and Python.
Reconcile the attempt categories; locate the rejected match and join/pad interval;
show the missing target, observed limits and proof scope without claiming global
infeasibility. Verify usage denominators, unused eligible parts, density/GC
calculations, and the quality plot against independently calculated fixture
values. Neither inspection nor rendering may run a solver. Pair this controlled
case with a real CBC integration run; a fake producer alone does not qualify it.

### D3. Add four designs and export a combined selection

Starting from D2's eight accepted designs, resolve an extension for four more.
The saved comparison must identify new target/effort/seed/exclusions and unchanged
requirements/inputs. Produce four additional sequences distinct from the parent
and each other in the cell, or report a truthful shortfall in the negative test.
The positive fixture contains enough feasible distinct sequences to reach four.
Check parent bytes unchanged, child target four, combined count twelve, and
stable sequence/placement joins in CSV/TSV/FASTA. Repeat with colliding local
design IDs, repeated source arguments and a second-generation extension.
Move the bundle to a new location; inspect/verify it after making the original
run paths unavailable in the isolated test. No source mutation, implicit ancestor
selection, original-path dependency, or fabricated completed-parent status passes.

D1–D6 qualify the complete workflow. Each capability needs its own applicable
recipe and reporting evidence before release. Document review and unrelated
passing tests do not establish task usability or runtime guarantees.

### D4. Select and reuse twenty-four designs

Use four qualified native run fixtures, each with one cell and ten accepted
designs, then declare six per full cell reference and a seed. This isolates result
selection from matrix generation. Python and CLI must
produce the same 24 ordered full references for the same ordered snapshots and
policy version. A total-count request must not imply those per-cell quotas.

Persist a selection snapshot, advance a source with a later valid commit, then
export/render the snapshot: it must use the original revision and chosen
identities. Repeating a source cannot multiply designs. Reusing a snapshot must
not draw random values or invoke the generator. Verify first-policy ordering,
zero allocations, ambiguous cell IDs, mutually exclusive count/per-cell fields,
and exact sequence duplicates with different design provenance.

Reduce one fixture to five candidates. Default shortfall fails before publication
with the affected cell and counts; explicit allow_partial produces 23 selected
designs, a qualified receipt and export exit 3, without borrowing from another
cell. A missing/changed source revision fails rather than substituting new data.
Move a complete selected bundle and verify it with original paths unavailable.

### D5. Preview preparation without doing it

For curated, PWM and background requests, Python and CLI
resolve equal preparation plans: source digests, effort, retention, filters,
scorer/background/algorithm versions and required capabilities. Assert zero
sampler/scorer/solver invocations and no workspace writes during planning.
Requested retention stays separate from unknown yield. Malformed inputs and
unsupported policies fail with field diagnostics.

Execute the matching plan through prepare. A generation plan requires an
explicit batch policy when passed to prepare; run rejects preparation plans. Replace an input between planning and execution
and verify failure before expensive work. Change defaults in a controlled fixture
and prove a saved normalized plan is unchanged. Missing tools have an explicit
capability error and installation route; no automatic installation or fallback.

### D6. Qualify cost and artifact contracts

Instrument the real artifact reader on small and larger fixtures. Unfiltered
summary must read committed summary/manifest data without iterating designs.
Verify streamed rows, early-close resource release, repeatable iteration,
cursor/query mismatch, coherent active snapshots and missing revisions.

For filtered quality and multi-source union, measure records examined and peak
memory, exercise explicit caps, and require declared scan/index state. Pairwise
work is bounded and opt-in; a stopped calculation must not claim exact metrics.
Verify the descriptor is available before iteration and CLI reports cost before
expensive work. Benchmark preparation, model construction, generation, verification,
selection and export independently; record fixture size, versions and resources.

For every native family, test a supported read/write round trip, unknown schema/
policy version, unexpected required field, malformed metadata namespace and
unchanged-source export. An unsupported view/format pair fails before publication.
For supported export projections, verify field changes and source/output
digests in the receipt; reject evidence loss. A valid checksum with invalid
placement/accounting data still fails verification.

## Verification matrix

| Risk | Evidence required |
| --- | --- |
| Python/CLI drift | Equivalent normalized plans, defaults, diagnostic codes/fields, and semantic records. Declare the comparison projection: exclude run IDs, destinations, timestamps, durations, and display text; retain part identity, constraints, lineage, proof and acceptance. Use deterministic fixtures for exact equality and validity/objective comparisons for tied solver results. |
| Product utility | D1–D6 and their linked domain/operation/artifact contracts: native inputs, requirement examples, explanations, independently checked quality metrics, extension without parent mutation, stable table joins, portable bundles and help-led task completion. |
| Invalid requests | Reject unknown fields, duplicate keys/IDs, boolean/fractional counts, invalid DNA, impossible bounds, stale inputs, unsupported mode/constraint combinations, and missing capabilities before heavy work. |
| Geometry and identity | Duplicate strings with distinct IDs, reverse complements, palindromes, repeated occurrences, overlaps, final coordinates after padding, core/part offsets, and narrow fixed-element exceptions. |
| Search interpretation | Real optimal and infeasible cases; controlled unproven/backend/invalid-result failures; exhausted offered batch versus unknown global feasibility; no false success after a printed prefix. |
| Matrix and uniqueness | Small totals, remainders, explicit inactive cells, expansion caps, per-cell uniqueness, and one sequence with multiple cell provenance. |
| Assembly/screening | Exact versus maximum length, discrete GC feasibility, junction-spanning forbidden motifs, pad-created motifs, and explicit relaxation records. |
| Accounting/recovery | Crash before and after the commit frontier, failure after output, repeated resume, competing writers through lock handoff, stale inputs, coherent active-run reads, interrupted attempts, unknown time and consumed budgets. |
| Randomness | Stable logical task streams across scheduling order where promised; preparation tie fixtures; documented solver tie limitations. |
| Algorithm changes | Fixed native fixtures, independently calculated results, versioned expected differences and readable prior artifact schemas. |
| Packaging/performance | Built-wheel base and optional installs, optional-import isolation, bounded inspection, small/medium preparation and generation benchmarks, and checks on each platform claimed by the release. |

### Match verification to the dependency

| Boundary | Verification approach |
| --- | --- |
| Normalization, allocation, identity and sequence checks | Exercise real deterministic functions with small fixtures and negative cases; compare invariants and semantic results, not private call order. |
| Packing | Real CBC solves for feasible/infeasible/tied cases; controlled backend statuses only for otherwise hard-to-reach failures. Never substitute a fake solver for all integration evidence. |
| Native runtime | Temporary directories with the actual persistence code; separate processes for locks; kill/fault injection around commit points. In-memory stores cannot qualify filesystem guarantees. |
| Time and randomness | A narrow injected monotonic clock for boundary/overrun tests, stable logical stream fixtures, real bounded smoke runs. Do not test exact sequence equality under wall-time stopping. |
| External scoring | Small command/parse/error fixtures plus a separately marked real-tool smoke test when installed. Missing tools must fail preflight before expensive work or partial publication. |
| Rendering | Reuse persisted realized records and existing renderer tests; prove rendering imports/execution do not invoke generation. |
| Artifact handoff | Export and reopen native records and bundles with matching identities, geometry and provenance; downstream protocols remain consumer-owned. |

Test a contract at the lowest boundary that proves it, then use a few complete
paths to check composition. Avoid duplicating every domain case through both
front ends or asserting mock call sequences that merely reproduce implementation.
Performance improvements require a measured bottleneck and before/after results;
measure algorithm changes separately from structural refactoring.

Before calling the ergonomics complete, run the
[product demonstrations and help recipes](delivery.md#product-acceptance-demonstrations)
from installed help, without giving the implementation file path. Record wrong
turns, undocumented prerequisites, commands required and unsupported requests.
Include recovery and preparation tasks wherever the release claims support.
The first array should need one execution command, and inspection/rendering one
command each; no workspace scaffolding or study knowledge is necessary.

Keep examples executable and record which acceptance cases have passed.
Qualify capability claims with runtime evidence.
Complete the [development gate](../../development.md) for changes
to this repository; consumer owners run their own checks before adoption.

## Release gates

Before publishing a schema, freeze its wire encoding, digest projection and
supported reader/writer versions. Before claiming recovery, qualify the commit
protocol and failure behavior on the declared storage. Benchmark default effort,
scan/state and model-size allowances with the workloads described in the
[resource guide](../../library-workflow/resources.md).

Release qualification uses built distributions in clean environments, both
base-only and with each supported extra. Record tested Python/platform/tool
versions and measured workloads. External systems own their deployment and
scientific acceptance; successful local generation does not establish either.

## Reader acceptance

A reader should find how to request exactly two parts, preview preparation
without scoring, explain eight accepted designs out of twelve, add four without
changing the parent, reuse a selected panel, and choose FASTA or a portable
bundle. A maintainer should find the owner and tests for a scorer or metric
through the code map. Verify task routes and documentation links; installed
command execution supplies separate usability evidence.

## Implementation routes

| Native owner | Evidence and interface |
| --- | --- |
| [Code map](../README.md) | Module ownership and corresponding tests. |
| [CLI reference](../../reference/cli.md) | Commands, request inputs, outputs and failure behavior. |
| [Optimizer](../../reference/optimizer.md) | Typed outcomes, solver controls and proof boundaries. |
| [Playback contract](../solution-playback.md) | Coordinate evidence and presentation requirements. |
| [Motif scoring](../../reference/motif-scoring.md) | Models, calibration, hit geometry and score interpretation. |
| [Background generation](../../development/background-generation.md) | Conditional distribution, counting method and resource outcomes. |
| [Packing measurements](../../development/optimizer-performance.md) | Measured model construction and solve behavior. |
| [Development gate](../../development.md) | Required repository validation. |

## Product milestone

See [workflow coverage](#workflow-coverage).

## Delivery slices and acceptance

See [acceptance by capability](#acceptance-by-capability).
