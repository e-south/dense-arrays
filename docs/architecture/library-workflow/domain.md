---
title: Library workflow domain
description: Define parts, designs, count bounds, geometry, preparation policies and proof scope.
author: Eric J. South
---

# Library workflow domain

Use [operations](operations.md) for invocation, [artifacts](artifacts.md) for
persistence and [selection](selection.md) for accepted-record queries. The
[native guides](../../library-workflow.md) provide examples; [delivery](delivery.md)
defines the evidence required for capability claims.

The product is a bounded search over realizable sequence designs. A candidate
must satisfy hard requirements; optional objectives rank otherwise valid
candidates. Sampling determines which candidates are offered to search. Storage
and rendering report the results of that process. These are separate decisions.

## Vocabulary

| Term | Meaning |
| --- | --- |
| Part | A supplied sequence with stable identity and optional source/core annotations. |
| Pool | A reusable, provenance-bound collection of eligible parts. |
| Design specification | The user's desired sequences, constraints, selection policy, targets, and limits. |
| Generation plan | An immutable resolved specification with bound inputs and explicit cells/targets. A `PlaybackPlan` separately describes presentation of realized placements. |
| Matrix | Named axes and choices expanded into design combinations with explicit allocations. |
| Matrix cell | One resolved combination of axis choices, with eligible parts, requirements and a target. Its readable ID is scoped to the matrix; the cell plan fingerprint binds its content. |
| Design | A final sequence, selected placements, requirements evaluation, and construction provenance. |
| Selected library | An explicitly ordered collection of design references from declared snapshots, potentially spanning runs. It is a view or portable artifact, not another execution. |
| Run | One execution of a generation plan, with its original targets, work history, accepted designs and termination reason. |

A selected library carries selection attainment; it never inherits run completion.
Origin targets and attempt counts keep their run scope. A design refers to one
`RealizedArray`; sequence equivalence never erases placement or run provenance.

Preparation retention chooses eligible parts; generation sampling offers a
candidate batch to packing; result selection chooses accepted designs. These
are distinct policies. Share identifiers and range predicates where meanings
match, not a universal policy object. `PartFilter`, `DesignFilter` and
`AttemptFilter` query their own record types; `LibrarySelection` adds bounded
selection of designs. See [selection](selection.md#result-selection).

A **candidate batch** is the subset of parts offered to one packing search.
A pool or selected library requires caller-supplied criteria to define an
experimental cohort.

## Primitive contracts

| Decision | Required behavior |
| --- | --- |
| Identity | Separate part identity, sequence identity, placement identity, design identity, and run identity. Equal DNA strings do not establish equal construction provenance. |
| Repeated parts | Count selected occurrences. A repeated string may have distinct part IDs; incidental substring matches do not add supplied occurrences. |
| Length | Distinguish `length.maximum` from `length.exact`. Exact length requires a declared assembly policy, which may prohibit padding. Never interpret a packing bound as a final-length guarantee. |
| Geometry | Publish zero-based, half-open coordinates on the final sequence. Keep original orientation/core coordinates and any padding transform as provenance. |
| Constraints | Separate hard requirements from preferences and from sampling policies. Use typed selectors and explicit units. Unknown or unsupported constraints fail before search. |
| Fixed elements | Reference named part occurrences, not a first string match. General pair-spacing requirements can represent promoter anchors without hard-coding a biological role. |
| Group requirements | Distinguish number of selected occurrences from number of distinct represented groups. Labels are caller-supplied metadata, not inferred biology. |
| Diversity | Name the object and metric: core diversity during preparation, part-use balancing during generation, or sequence distance in a library. A single `diverse` switch cannot promise all three. |
| Feasibility | Name the offered parts, packing formulation, exclusions, and proof scope. A failed sampled batch cannot prove the entire request infeasible. |
| Screening | Check the final assembled sequence, both configured orientations, junctions, and padding. Exceptions attach to intended placement intervals. |
| Randomness | Bind randomness to logical work identities and algorithm versions. Record tie-breaking and representative-selection rules. |
| Targets | Make matrix allocation, unique-design policy, and incomplete completion visible before and after execution. |

Spacing uses an explicit relation: downstream start minus upstream end in final
sequence coordinates. Negative values mean intentional overlap. The
playback v1 distance contract has narrower capabilities. Generation requirements
are not playback distance records: reuse its placement representation, and keep
the complete requirement evaluations in the design record. Rendering that cannot
represent a declared requirement fails before publication; it must not drop that
requirement or convert a negative distance to zero. Any additional rendering geometry requires an explicit contract and
acceptance evidence.

Fixed requirements identify supplied part occurrences and their orientations.
Planning rejects unsupported combinations. Lower final-coordinate requirements
into packing coordinates explicitly, then recheck them after any padding
transform. Passing a final-coordinate value unchanged to a pre-padding solver
is not a valid translation.

A part is eligible for at most the declared number of occurrences. The
contract supports one occurrence per part ID; multiple copies use explicit
distinct IDs. Single/double-strand eligibility and a named fixed-element strand
are separate fields. Palindromes and equal strings do not remove identity.

Identity distinctions do not require five new global identity services. Part
IDs are scoped to their input collection; placement IDs to a design; native
design IDs to a run/cell; and run IDs to an execution. Keep operational identity
separate from content equivalence. Before a wire schema ships, document which
fields form each digest, its namespace/version, and which ordering is semantic.
Constructed design records contain one `RealizedArray` plus requirement and
lineage references; do not duplicate authoritative sequence/placement fields.

Uniqueness means exact final-sequence equality within a plan cell.
The same sequence may satisfy multiple cells, with distinct design provenance.
Report both design count and distinct sequence count. Reverse-complement equivalence and global uniqueness do not apply
implicitly. External sequence stores own their identifier conventions;
a native sequence digest must not be substituted for a store identifier.

## Requirements

Use a small typed requirement catalog. Every requirement has a unique `id`
within the specification; errors, reports, and placement annotations refer to
that ID. Counts refer to selected supplied occurrences, never incidental string
matches. Selectors contain either explicit `part_ids` or caller-supplied
`groups`, not executable predicates. A part has at most one group;
multiple independent grouping dimensions are unsupported.

| Requirement kind | Fields and meaning |
| --- | --- |
| `occurrences` | `select`, optional inclusive `min`/`max`: bounds on selected occurrences among the named parts/groups. Equal bounds express an exact count; `max: 0` excludes the selector. At least one bound is required. Bounds are enforced during packing and independently recounted from placements. |
| `group_coverage` | `groups`, `min`: at least this many distinct named groups represented. Distinct from occurrence counts. |
| `fixed` | `part_id`, `orientation`, optional inclusive `start: {min, max}` window in final coordinates. Selection refers to the named occurrence and its declared orientation. |
| `spacing` | `upstream`, `downstream`, inclusive `min`, `max` for downstream start minus upstream end. One declared fixed pair is supported; negative values mean overlap, subject to formulation and rendering support. |
| `avoid` | Literal A/C/G/T `patterns`, `strands: forward|both`, optional `except_placements` naming fixed part IDs. Final-sequence screening exempts a match only when wholly inside the specified realized interval. |
| `gc` | `scope: sequence|padding`, inclusive fraction `min`, `max`. Final checks; padding scope requires a padding policy. Zero added bases produce explicit `not_applicable`, not a fabricated GC value. |

Occurrence bounds are nonnegative integers, never booleans or fractional values;
an omitted minimum is zero and an omitted maximum is unbounded by this rule.
Reject reversed bounds and missing both bounds. Group-coverage minimums remain
positive and cannot exceed the number of distinct eligible groups. A minimum
above available eligible identities and a fixed occurrence forbidden by a zero
maximum are static contradictions. Overlapping selectors impose independent
bounds; do not merge them into one count or double-count one placement within
a selector. Sampling caps cannot stand in for accepted-design count bounds.
The packing model enforces these bounds, and final validation recomputes them
from selected placement identities. Fixed references must exist and declare one
orientation. All starts are zero-based; placement intervals are half-open.
Spacing is measured in forward final-sequence coordinates regardless of motif
orientation. `upstream` and `downstream` identify that coordinate relation, not
an inferred biological role. Fixed parts count toward a count rule only when
its selector includes them.

This catalog does not promise arbitrary pairwise ordering,
multiple interacting anchor pairs, selectable gaps,
or IUPAC pattern interpretation. These need named capability decisions
and tests; they must not be approximated by a minimum count or a literal match.
The planner reports the missing capability and an available supported route.
Length is the top-level `maximum`/`exact` choice. Objectives and
preferences remain separate from this list of hard requirements.

## Part inputs

Inputs include explicit CSV/TSV tables and typed Python `Part` collections. Native
pool artifacts are reusable input handles. CSV/TSV require a header; table
format is declared. Map `part_id`, `sequence`, optional `group`, `source`,
`core_start`, `core_end`, and `core_orientation` through `columns`. Import errors identify source, logical
one-based data-row number, original column, and normalized field. Blank optional
groups mean no group; blank sequence/required ID values are errors.

Default column names are the canonical names. Unmapped columns are ignored
with an inventory in the import report; `metadata_columns` explicitly preserves
chosen columns as namespaced caller metadata. They cannot override validated
core fields. Python accepts the same metadata mapping. Core coordinates require
both bounds and orientation and must align to the supplied sequence; absent
core annotations mean unknown, not an inferred PWM hit.

Default import policy preserves distinct rows/part IDs and rejects malformed
DNA. With `id_policy: provided`, IDs must be nonempty and unique in the input
collection. When IDs are absent, explicit `id_policy: row` creates IDs from
logical row order in the bound collection, preserving repeated sequences as
distinct occurrences. These IDs are stable within that immutable collection;
reordering or replacing a source is not a cross-file identity guarantee.

Normalization defaults to strict uppercase A/C/G/T with no whitespace. Optional
`normalization: {uppercase: true, trim_outer_whitespace: true}` records every
transformation; internal whitespace and ambiguous bases still fail. Neither
case conversion nor duplicate removal is implicit. `duplicates: retain` is the
default; preparation can explicitly deduplicate by sequence/core with a named
representative rule and a discarded-to-retained lineage map. Distinct source
identities must remain recoverable even when one sequence is retained.

Direct planning and curated `prepare` share this importer. Curated preparation
uses `schema: dense_arrays.prepare.v1`, `source.kind: table`, the same table
fields, and an optional `retain: {select: ...}` using a typed `PartFilter`.
Without retention selection it keeps all valid rows. It does not require
a mining budget or scorer because it performs neither task. Reusing a prepared
pool binds its immutable snapshot and selection; it does not rescore or remine.

A design can replace its table input with `parts: {pool: pools/curated}` or
`parts.PoolSource(pool=curated)`; optional `select` uses `PartFilter`.
Planning resolves that selection and checks requirements against its
eligible population. Selecting only A cannot silently satisfy a requirement
for group B. Inspection filters do not mutate or rebuild the pool. Import reports retain
row/field diagnostics with a bounded displayed sample and total error counts;
invalid input is not partially accepted as a valid pool.

Parquet and XLSX inputs use optional table readers with the same identity,
normalization and metadata contracts as CSV/TSV. Qualify each supported format
with typed and malformed inputs; see [table readers](../../library-workflow/tables.md).
FASTA sequence export and annotated placement export have distinct purposes.
FASTA input requires an explicit identity and annotation policy before support
can be claimed.

## Packing and acceptance

The resolved plan names its packing objective, placement preferences, and
library diversity policy separately. The default maximizes selected
part occurrences within the length bound. This is not an objective to minimize
sequence length or maximize measured function. Each attempt records the
effective objective, weights, and exclusions when usage balancing changes
them. PWM scores retain their declared scoring meaning; they are not silently
renamed binding affinity or functional performance.

The current model searches ordered paths using one predetermined overlap shift
per pair of oriented parts. Exact search means proven results within that model;
it does not enumerate arbitrary spacers, every possible overlap length, or all
padding realizations. Record the formulation and path-exclusion policy alongside
the offered parts. Path identity and final-sequence uniqueness remain different
concepts. Extending the feasible space to selectable gaps/overlaps is a separate
mathematical feature with its own fixtures and cost measurements.

The engine consumes a resolved candidate batch and typed constraints. Time,
thread, and backend-specific options are declared capabilities, validated once,
and passed through supported APIs. Unsupported thread control fails rather than
being silently ignored. Exact, enumerated, usage-balanced, and heuristic modes
retain their distinct mathematical guarantees.

Solver controls and typed solve reports belong at the packing boundary.
The workflow consumes typed outcomes rather than inferring statuses from
exception messages.
Report proof status separately from termination reason. The current optimizer
does not expose a reliable timeout reason: `NOT_SOLVED`, for example, becomes
`SolverBackendError`. Record an unknown termination reason when a backend cannot
establish it. A configured timeout or elapsed duration alone does not prove why
a solve ended. Likewise, distinguish an infeasible initial batch from the end
of its constrained enumeration without claiming global library exhaustion.

Exact mode requires proof of optimality for each packing result.
This is local optimality for the offered batch and declared packing objective;
it is not a claim that the final library is a globally optimal diverse set.
Feasible-but-unproven incumbents are reported separately and not silently
accepted. Exact mode does not accept unproven incumbents. Heuristics must reject
requirements they cannot enforce; they cannot discard constraints to proceed.

Successful packing passes through assembly/padding, final sequence checks,
uniqueness, and publication. Final hard constraints always apply after padding.
GC limits name their scope (pad or final sequence). Relaxation requires a named
policy and a recorded reason; strict behavior is the default. Joining valid
parts does not imply that the final sequence is valid.
Bound inner searches too: padding proposals and preparation draws need explicit
finite limits. Multiple assembly trials within one packing attempt do not count
as multiple solver attempts; record their effort separately. An exhausted
padding policy does not prove that no possible padding could pass.

An optimizer proves facts about its offered batch, not every possible batch or
the entire requested library. Exhausted enumeration, infeasible batch, global
search completion, budget exhaustion, and unavailable execution are distinct.
Do not infer global exhaustion from a resampling cap or repeated duplicates.
