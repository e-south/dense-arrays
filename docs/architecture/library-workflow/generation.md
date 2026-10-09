---
title: Generation and revision recipes
description: Compare Python and CLI requests for bounded generation, additional designs and requirement changes.
author: Eric J. South
---

# Generation and revision recipes

Use these paired examples to check operation semantics. For a complete user
workflow, start with the [saved-library guide](../../library-workflow.md). Run
from a new directory with the declared inputs; output paths must not exist.
[Operation rules](operations.md) define shared effects and failures.

## First array

Use an installed CBC backend for generation and the `playback` extra for the
PNG. The [installation route](../../installation.md) covers those prerequisites.
These synthetic sequences illustrate packing; they are not a TFBS model.
`--length` explicitly means a maximum:

```text
dense-arrays run --motif ACGTTGCAAGTCCTGA --motif AGTCCTGATCGTACCG \
  --motif TCGTACCGATGCTTAG --motif ATGCTTAGGACGTTCA \
  --length 40 --count 1 --seed 7 --out runs/first-array
dense-arrays inspect runs/first-array --verify
dense-arrays render runs/first-array --out first-array.png
```

Inline motifs receive stable occurrence IDs within that ordered request. The
quick path resolves to the same defaults as a one-cell specification. Defaults are a target of one design, exact search with required
optimality, double-strand eligibility, seed `0` when omitted, at most 1,000
candidate attempts across the run, 300 accumulated active seconds, and 30
seconds per solver call. [Resource bounds](../../library-workflow/resources.md)
explain additional size limits and measured work. Limits do not promise completion.

The result is a persistent run handle with inspectable accounting. Low-level
packing remains available for transient in-memory use. Workflow runs use the explicit output directory.

## Constrained design

Use the [curated binding-site example](../../library-workflow/preparation/curated.md#bind-curated-parts-and-requirements)
for paired Python/CLI requests with 16-base sites, occurrence counts and group
coverage. It creates `parts.csv` and a saved `runs/curated` run.

The [exact-length assembly example](../../library-workflow/generation/assembly.md#assemble-and-render-an-exact-length-design)
adds named 16-base anchors, declared spacing, right padding and final screens
for a 40-base sequence. It checks the persisted coordinates and renders the
selected design. The [domain contract](domain.md#requirements) defines those
requirements and their proof limits.

Planning binds normalized requirements and input evidence; bounded generation
may stop below its target. Read attainment and termination from the saved run.

## Revise or extend a library

Resume continues an unchanged plan. Extension creates a new plan/run with a
frozen parent selection, new effort budget, explicit seed, and an **additional**
target. These are different user intentions under the same `plan`/`run`
operations; no overloaded meaning of `--resume` is introduced.

The examples below require an integrity-verifiable terminal single-cell run at
`runs/initial`. Use that destination when creating the parent, or replace the
path consistently with an existing run. Save this `extend.yaml` beside `runs/`:

```yaml
schema: dense_arrays.extension.v1
parent: {run: runs/initial}
additional: 4
limits: {attempts: 1000, active_seconds: 300, solver_seconds: 30}
seed: 19
```

```text
dense-arrays plan extend.yaml --out extension.plan.json
dense-arrays inspect runs/initial --view plan --compare extension.plan.json
dense-arrays run extension.plan.json --out runs/additional
dense-arrays inspect runs/initial runs/additional --view quality
```

```python
import dense_arrays as da
from dense_arrays import planning

extension = planning.ExtensionSpec(
    parent=planning.ParentRun(run="runs/initial"),
    additional=4,  # Four new designs; the parent count is not part of this target.
    limits=planning.Limits(attempts=1000, active_seconds=300, solver_seconds=30),
    seed=19,
)
extension_plan = da.plan(extension)  # Pin the parent revision and exclusion set.
difference = da.inspect("runs/initial", view="plan", compare=extension_plan)
additional = da.run(extension_plan, out="runs/python-additional")
combined = da.inspect(["runs/initial", additional], view="quality")
```

Plan resolution pins the parent's committed revision and accepted-design digest.
Only terminal, integrity-verifiable parent runs are eligible for extension;
an actively written run must first stop. The child reuses immutable input pools,
requirements, cell identities, assembly, and selection policies. Limits and seed
are new explicit values. Parent counters or consumed budgets never transfer as
fresh allowance within the old run. The parent remains unchanged.

Every child candidate is checked for exact final-sequence duplicates against
the frozen parent library and the child accepted prefix within the corresponding
cell. A duplicate consumes effort and records its matching source reference;
it does not satisfy the additional target. The exclusion index must be usable
without replaying the parent's solver. Restarted enumeration may spend effort
on rejected duplicates; efficiency is measured, not promised by the contract.

If the parent accepted eight of twelve and the child accepts four, the child
is complete at four of four. The combined selection contains twelve unique
sequences for that cell; the original run remains stopped at eight of twelve.
The report must display those distinct facts. A short child is another explicit
shortfall, not an implied completion of the parent. For successive extensions,
the parent exclusion set includes the verified ancestor exclusions bound into
the parent plan, so a later child cannot regenerate an earlier ancestor's sequence.
This does not automatically add ancestors to a displayed/exported selection.

Integer `additional` is single-cell only. Matrices use explicit
per-cell additional allocations over unchanged cell identities; no implicit
quota redistribution is permitted. Changing requirements uses a new full
`DesignSpec`, not an extension request. Its optional `lineage.parent` references
the original run and may identify an explicit accepted-library exclusion input.
Such exclusion requires a declared cell mapping and uniqueness scope; it cannot
be inferred across changed matrices.

For ordinary revision, `export --view request --format json --out` serializes a
full editable `DesignSpec` reconstructed from resolved fields, not a checkpoint
or a new implicit inheritance chain. It preserves explicit policies and bound
input references and adds parent lineage. Relative locators resolve from the
new file location. The file uses the native design schema and can be passed to
`plan` after editing; it is not a JSON report envelope disguised as an input.
The Python `RequestReport.request` supports immutable `with_changes(...)`, which
returns a validated request value and does not plan or execute anything.

```text
dense-arrays export runs/initial --view request --format json --out revised.design.json
# Edit the desired requirement or length in revised.design.json.
dense-arrays plan revised.design.json --out revised.plan.json
dense-arrays inspect runs/initial --view plan --compare revised.plan.json
```

```python
import dense_arrays as da
from dense_arrays import planning

prior_request = da.inspect("runs/initial", view="request").request
# Editing creates a request value; planning and execution remain explicit.
revised_request = prior_request.with_changes(length=planning.Length(exact=48))
revised_plan = da.plan(revised_request)
revision_report = da.inspect("runs/initial", view="plan", compare=revised_plan)
```

This revision does not imply sequence exclusion. To prevent reuse of prior
sequences, explicitly retain/add the accepted-library exclusion described above;
the plan preview must show its scope. Extending unchanged requirements uses the
safer dedicated extension request and includes that exclusion automatically.

The plan comparison reports added/removed/changed inputs, requirements, targets,
policies, limits, seed, and exclusions. Paths, timestamps and display settings
appear separately from semantic changes. Run quality comparisons use the same
metric versions and expose population/denominator changes; incompatible metrics
are marked incomparable. They do not imply a causal effect or statistical
significance. Comparisons read artifacts only and never solve.
