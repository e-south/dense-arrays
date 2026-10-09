---
title: Library workflow architecture
description: Find the contracts and module owners for preparing parts and generating inspectable DNA libraries.
author: Eric J. South
---

# Library workflow architecture

Dense Arrays prepares reusable parts, generates constrained DNA libraries and
preserves the evidence needed to inspect, select and share designs. Python and
the CLI use the same application operations and request contracts.

Use the [saved-library guide](../../library-workflow.md) for working examples and
the [API reference](../../api.md) for supported interfaces. The contracts below
define operation behavior and evidence requirements. [Delivery](delivery.md)
specifies the checks required for capability and release claims.

| Task | Contract owner |
| --- | --- |
| Define inputs, identities, constraints or proof scope | [Domain](domain.md) |
| Invoke the six operations and read typed results | [Operations](operations.md) |
| Compare paired generation, revision and extension requests | [Generation recipes](generation.md) |
| Compare curated and sampled preparation requests | [Preparation recipes](preparation.md) |
| Store committed evidence and recover work | [Artifacts](artifacts.md) |
| Interpret diagnostics, metrics and read limits | [Reporting](reporting.md) |
| Filter records or save a bounded panel | [Selection](selection.md) |
| Publish projections or move selected evidence | [Exports](exports.md) |
| Implement a capability or qualify a release | [Delivery and acceptance evidence](delivery.md) |
| Locate the implementation and tests for a change | [Code map](../README.md) |

## Scope

Curated inputs, PWM and background preparation, plan expansion, candidate
sampling, packing, assembly, final screening, attempt accounting, recovery and
portable results belong to the library workflow. Packing and realized-placement
contracts supply its sequence geometry and proof boundaries.

Experimental design, biological interpretation, motif acquisition and inference,
sequence stores, campaign fitting, schedulers, notifications, study-specific
figures and publication decisions remain caller-owned. This architecture does not require a GUI, service,
distributed scheduler, plugin platform or arbitrary-gap packing formulation.

## Interface decision

Use six operations: **prepare, plan, run, inspect, export, render**. Preparation
and explicit plan persistence are optional. These names describe tasks, not a
mandatory six-stage pipeline. [Operation contracts](operations.md#shared-operations)
define their effects and return types.

`inspect` reads reports and records; `export` publishes data; `render` publishes
visuals. They share selection and artifact readers. Comparisons are inspection,
and resume is a run mode. Typed constructors, iteration and serializers support
these operations without introducing another execution system.

## Code ownership

Create a cohesive package when it owns behavior that can change independently.
Keep stable packing modules in place; avoid empty scaffolding and parallel
implementations of the same policy.

```text
src/dense_arrays/
├── parts/          # identities, importing, preparation, scoring, retention
├── planning/       # requests, requirements, resolution, capability checks
├── generation/     # batch sampling, assembly, screening, attempt decisions
├── artifacts/      # wire formats, commits, recovery, readers, publication
├── reporting/      # queries, quality, comparisons, diagnostic projections
├── workflow/       # public operations, lifecycle composition, CLI translation
├── playback/       # visual presentation of persisted placements
└── …               # packing inputs, models, optimizer and sequence geometry
```

The root package exposes the operation facade and common request/handle types.
Specialized public types live with their semantic owner. Front ends call the
same application functions. CLI parsing must not determine defaults, selection,
acceptance or recovery. Tests follow those owners, with focused complete-path
checks for their composition.

Reporting produces typed results; artifacts serializes and publishes them.
Artifacts does not calculate quality or select experimental cohorts. Workflow
composes planning and generation with input reads and native commits. Generation
emits attempt results without calling storage or orchestration. Domain records
do not import storage, presentation or study packages.

| Knowledge | Single owner |
| --- | --- |
| Part validity and preparation policies | Parts |
| Requirement values, defaults and plan resolution | Planning |
| Packing formulation and solver controls | Packing core |
| Final acceptance and attempt decisions | Generation |
| Committed evidence, schema versions and recovery | Artifacts |
| Query meaning, result selection and metric definitions | Reporting |
| Lifecycle and front-end translation | Workflow |
| Visual encoding | Playback |

A scorer change belongs in preparation and its scoring tests; a storage change
belongs in artifacts and its conformance tests; a metric change belongs in
reporting. CLI formatting cannot change acceptance. Optional tools must remain
removable without preventing base imports. Introduce interfaces at demonstrated
substitution boundaries, such as a solver, external scorer or monotonic clock.
One shared validator owns each input invariant; verification independently
recomputes evidence instead of trusting cached validity flags.

The base installation supports curated generation and inspection. Optional
scoring and rendering remain lazy: importing the workflow must not require
FIMO or plotting libraries. Packing retains its declared solver dependency.

## Documentation and agent routes

Each contract page owns its definitions; other pages link to it. Runnable recipes
belong in task guides. Keep `AGENTS.md` as a task-to-owner-and-checks router and
add scoped instructions only when a subtree needs different rules. Pages retain
title, description and authorship frontmatter. Capability limits belong beside
the affected operation, with a direct route to its supported alternative.

[Delivery](delivery.md) defines acceptance evidence and release gates.
[Development](../../development.md) defines the repository checks.
