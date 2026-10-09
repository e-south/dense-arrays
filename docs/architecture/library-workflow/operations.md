---
title: Library workflow operations
description: Invoke the six Python and CLI operations with explicit inputs, output ownership and typed results.
author: Eric J. South
---

# Library workflow operations

Use the six operations to prepare parts, generate designs and read or publish
saved evidence. [Domain](domain.md) defines field meanings;
[artifacts](artifacts.md) defines persistence and recovery.
[Generation](generation.md) and [preparation](preparation.md) provide paired
contract examples. Runnable user workflows start in the
[saved-library guide](../../library-workflow.md).

## Shared operations

| Operation | Input and result | Effects |
| --- | --- | --- |
| `prepare` | Preparation request/set/plan → `PoolHandle`; generation/matrix plan plus batch policy → saved plan | Create a reusable pool or freeze candidate batches. |
| `plan` | Preparation request/set → `PreparationPlan`; design/extension request → `GenerationPlan`; matrix request → `MatrixPlan` | Read declared inputs and resolve policies. No sampling, scoring, solving or workspace creation. |
| `run` | Design/extension/matrix request or generation/matrix plan → `RunHandle` | Create one run, or resume its unchanged plan through the exclusive resume form. |
| `inspect` | Artifact(s) → typed report, record view or selection snapshot | Read only; no exports, repairs, scoring or solving. |
| `export` | Artifact(s) plus selection/representation → `ExportReceipt` | Publish data to an explicit destination through shared readers and serializers. |
| `render` | Persisted evidence plus selection/view → `ExportReceipt` | Publish a visual without regeneration or acceptance calculations. |

```text
prepare(request_or_plan, *, out, sampling=None, batch_count=1,
        attempts_per_batch=None, accepted_per_batch=None)
plan(request)
run(request_or_generation_plan, *, out)
run(*, resume)
inspect(artifact_or_artifacts, *, view="summary", select=None,
        compare=None, verify=False, limit=None, after=None, all=False, read_limits=None)
export(artifact_or_artifacts, *, out, view=None, select=None,
       format="json", all=False, read_limits=None, compare=None)
render(artifact_or_artifacts, *, out, view="design", select=None, read_limits=None)
```

```text
dense-arrays prepare SOURCE --out POOL
dense-arrays plan SOURCE [--out PLAN]
dense-arrays run SOURCE --out RUN
dense-arrays run --resume RUN
dense-arrays inspect ARTIFACT [ARTIFACT ...] [--view VIEW] [--verify]
dense-arrays export ARTIFACT [ARTIFACT ...] --out OUTPUT [--view VIEW] [--format FORMAT]
dense-arrays render ARTIFACT [ARTIFACT ...] --out OUTPUT [--view VIEW]
```

Python accepts typed values without YAML, shell calls or dataframes. File
dispatch uses a declared schema, never a guessed filename or directory.
Front ends converge before defaults, validation and policy resolution.
Specialized types have one documented import path in their owning package;
equal capabilities do not require Python users to emulate CLI flags.

The resume overload accepts no source, destination, changed seed, target or
limits. Extension is a new request and run, not a resume override.
`run` rejects preparation plans. `prepare` accepts a resolved generation or
matrix plan only with an explicit [batch sampling policy](../../library-workflow/batches.md);
it then freezes offered parts rather than generating designs.

## Invocation and output rules

- Sources and destinations are explicit. Relative input paths resolve from
  the request file; in-memory Python calls supply paths or handles explicitly.
  Never discover an active study or environment-selected request.
- `--out` names a destination to create. Collision checks precede expensive work.
  Native workflow outputs are create-only.
- `--json` requests a versioned report or receipt on stdout for any operation.
  Default progress/receipts use stderr; inspect/plan display human reports on
  stdout. Python always returns the corresponding typed object.
- `export --format` selects the data representation, not the receipt.
  `--out -` explicitly streams JSON/CSV/TSV/FASTA data to stdout; it rejects
  bundles and other directory formats. Diagnostics and receipts remain on stderr,
  including when `--json` is present.
  Explicit files use atomic create-only publication; redirected stdout can
  retain partial bytes on failure.
- Python export accepts a new path or an explicit caller-owned writable text
  stream for supported stream formats. CLI `--out -` lowers to stdout at this
  boundary; Python does not print data implicitly. Export never closes a
  caller-owned stream, and a stream failure reports partial publication.
- `inspect` has no `out` or file-format switch. Saving a report/request uses
  `export`; supporting Python serializers remain available.
- Use the six names in help and examples. No generate/execute/build aliases,
  workspace init/reset family, or separate status/doctor execution path.
- Inline motifs are a quick `run` form. Complex inputs use a typed/file request.
  Both resolve through the same planner; explicit `plan` is optional.
- Progress is bounded, non-interactive and off on non-TTY streams. Color carries
  no unique meaning. Errors name the field/artifact, observed and expected
  values, code and recovery action; normal errors omit a traceback.

Before publication a failure leaves no completed artifact. After a native commit,
errors expose its run reference and intact committed prefix. Preparation shortfalls
return explicit partial status; generation planning rejects incomplete pools.
See [accounting](artifacts.md#accounting-and-completion) for statuses and exits.

The human run receipt leads with accepted/target designs, length and cell scope,
termination reason, destination and one actionable next command. For example:

```text
Stopped: 8 / 12 designs accepted; attempt limit reached.
Results: runs/example
Details: dense-arrays inspect runs/example --view diagnostics
```

The report links to the extension recipe when additional effort requires a new
run. It never suggests that repeating resume restores an exhausted budget.

## Preview a request

Curated inputs may be parsed and normalized directly during planning; prepare
is optional for them. PWM sampling, external scoring and retention optimization
require an explicit preparation operation and reusable pool. Planning cannot
invoke those operations as a side effect.

A generation plan binds input digests, normalized requirements, resolved cells
and targets, algorithm versions, limits, and seed policy. It also records
required capabilities. Planning checks syntax, identity, geometry bounds, and
supported combinations; it does not prove feasibility or certify that a solver
license will be available later. Run preflight checks the actual dependencies,
input digests, and output ownership before starting work.

Defaults apply when resolving a request, not when loading an already normalized
plan. Unknown schema/policy versions fail; a new release cannot reinterpret an
old plan using new defaults. Preflight and execution must consume the same
verified input bytes, through an immutable snapshot or owned staging. Hashing a
mutable path and reopening it later is insufficient. Staged input capture is a
runtime responsibility; it must not make `plan` create a workspace.

The plan preview must show: input and retained-part counts; required occurrences
and groups; length interpretation; expanded cells and quotas; uniqueness scope;
search/proof/diversity policy; bounded effort; required external tools; and any
explicit relaxation. Separate definite contradictions from unknown feasibility.

For matrices, choose exactly one target declaration: per-cell counts, explicit
cell allocations, or a total with a named allocation rule. A total smaller than
the number of active cells fails unless zero-target cells are explicitly
declared. Allocation order, remainder handling, and inactive cells appear in
the saved plan. Expansion is bounded before execution.

A `PreparationPlan` binds source/scoring/background identities, normalization,
algorithm versions, effort and retention policies. Its preview shows source
part/motif counts, candidate budget, requested retention, screening stages,
required tools and retained-count uncertainty. It reads and validates inputs
and detects tool availability without invoking a scorer or sampling even once.
Unknown eligibility/yield remains unknown; requested retention is not a forecast.
Execution rechecks tools, bytes and destination ownership. Editing a request
invalidates the old plan; loading a normalized plan never reapplies defaults.

The generation preview also reports candidate-batch bounds, oriented part count
and formulation-specific size indicators. For the current path model, `n`
oriented nodes imply `n(n-1)` internal directed pair variables plus boundary
variables. Report an upper bound if a batch has not yet been sampled. Do not
build the solver model merely to preview it, infer runtime from that count, or
silently shrink a batch to meet a hidden cap. [Cost contracts](reporting.md#cost-and-reader-lifetime)
cover reporting; [delivery](delivery.md#verification-matrix) requires measurements.

## Python results

`prepare`, `plan`, `run`, `inspect`, `export`, and `render` remain the application operations.
Constructors, iteration, and serialization are supporting APIs, not separate
execution engines. Typed input does not require YAML or pandas. `run` retains an
explicit output destination; low-level packing remains available for transient
in-memory use. No hidden temporary workspace is introduced for notebooks.

| Return value | Required public behavior |
| --- | --- |
| `PoolHandle`, `RunHandle` | Stable artifact reference and schema/identity; accepted by other operations. No eager dataframe or dependence on private filenames. |
| `PreparationPlan`, `GenerationPlan` | Immutable resolved fields, capability requirements and public `write(path)` serializer; accepted by their matching executor and `inspect`. |
| Summary/quality reports | Typed counts, metrics, diagnostics, scope, and completion fields. JSON serialization uses the same schema as CLI. |
| `RequestReport` | Full editable typed request with parent lineage; native request-file serialization preserves explicit policies and remains accepted by `plan`. |
| `RecordView`, `LibraryView`, `BundleView`, `SelectionView` | Lazy `records()` iteration, declared selection/snapshot, page size and optional next cursor. Repeated access remains tied to the same snapshot. |
| `SelectionSnapshot` | Ordered full design references, pinned input revisions, filter/take policy and version, requested/selected counts and shortfalls; reusable by inspect/export/render. |
| `ExportReceipt` | Paths, format/schema, selection identity, counts, source snapshots and publication status. Errors report any published outputs. |
| `PlanComparison`, `QualityComparison` | Typed semantic field changes or comparable metric differences, with incompatible fields explicitly labeled. |

Default notebook/text representations show artifact identity, bounded summary,
completion and follow-up operations. Merely displaying an object never solves,
scans all records, renders media, or imports plotting/dataframe packages. Optional
dataframe conversion is a projection over explicit record selection. Exceptions
expose code/field/evidence and, after publication, a run reference; users should
not parse exception strings to recover a partial run.

`inspect` uses view-specific overloads: summary/quality return their report
types, record queries return the matching record view, request returns
`RequestReport`, comparisons return `PlanComparison` or `QualityComparison`,
and selection returns `SelectionSnapshot`.
There is no flag-dependent transition to file publication. `export` always
returns a receipt; filesystem paths and text formatting never replace typed data.
Views reject filters for a different record kind before reading records.

## Help and discovery

Use `dense-arrays --help` to choose an operation, then its `--help` for supported
inputs, views, formats and selection flags. The [CLI reference](../../reference/cli.md)
and [task guides](../../library-workflow.md) supply examples and failure routes.
Plan reports expose resolved defaults, required tools, static contradictions
and unknown feasibility. Tool checks do not run a solver or scorer.

Missing optional tools report an installation route without installing them.
Non-TTY output remains bounded and scriptable; live progress reports counts
and stage without predicting completion time.

## First array

See [maximum-length generation](generation.md#first-array).

## Constrained design

See [paired constrained requests](generation.md#constrained-design).

## Curated pool

See [curated preparation](preparation.md#curated-pool).

## PWM and background preparation

See [sampled preparation](preparation.md#pwm-and-background-preparation).

## Revise or extend a library

See [new requirements and additional designs](generation.md#revise-or-extend-a-library).
