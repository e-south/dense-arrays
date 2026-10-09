---
title: Constrained background generation
description: Choose a sampling method for GC bounds and forbidden DNA patterns, with explicit distribution and resource guarantees.
author: Eric J. South
---

# Constrained background generation

Use direct conditional sampling for background parts subject to sequence
GC bounds and literal exclusions. Choose each next base according to its declared
probability and the probability mass of valid suffixes. Keep OR-Tools available
for feasibility, enumeration and objectives requiring more general coupled
constraints.

Choose
`Sampling(strategy="conditional")`, with shared Python/CLI preparation, saved
construction outcomes and independent sequence verification. Use the stochastic
strategy when candidate sampling followed by screening meets the required yield.
See the [user guide](../library-workflow/background.md).

## Method choice

Forbidden patterns define a deterministic finite-state language. Combining its
state with sequence position and GC count gives a finite dynamic program:

`P(next base | prefix, valid) ∝ P(base) × valid suffix mass`.

For uniform DNA at a fixed length, exact integer counts and unbiased integer
rank selection produce uniform valid sequences. Determinism matters: counting
paths through an ambiguous representation can overcount strings. The connection
between deterministic automata and dynamic-programming string counts is
described by [Antonopoulos et al.](https://www.cs.yale.edu/homes/antonopoulos-timos/ICDT-2011.pdf).
The GC-count extension and Dense Arrays sampling contract are local design
decisions, not claims made by that paper.

CP-SAT supports [automaton constraints](https://or-tools.github.io/docs/python/classortools_1_1sat_1_1python_1_1cp__model_1_1CpModel.html)
and integer constraints. The comparison tested an automaton plus GC counts and
one Boolean per base/position with forbidden-word clauses plus GC counts. The
latter was markedly better for distinct-sequence enumeration in these fixtures.
A prefix of feasible solver solutions does not establish uniform or
background-weighted sampling.

| Need | Chosen approach |
| --- | --- |
| Sample under GC and literal exclusions | Completion counting and conditional draws |
| Find or enumerate feasible sequences | CP-SAT with explicit limits and sequence-level uniqueness |
| Optimize additional coupled requirements | Evaluate a CP-SAT model with a declared objective and proof scope |
| Apply FIMO qualification or other external screens | Preserve independent screening after proposals |
| Insufficient counting or solver resources | Return a limited/unknown outcome; no silent method substitution |

## Observed comparison

The [measurement record](measurements/background-generation.json) contains exact
synthetic inputs, outcomes and environment. These are one-run measurements on
macOS arm64, Python 3.14.7 and OR-Tools 9.15.6755, not general throughput claims.
Times below include model/table construction and generation of 100 sequences.
CP-SAT used one worker and a two-second solve allowance; its rows are feasible
enumeration results, not random samples.

| Fixture | Rejection draws | Counted sampling | CP-SAT clauses |
| --- | --- | --- | --- |
| 20 bases, 6–14 GC bases, four forbidden 6-mers | 100 accepted in 100 draws; 0.24 ms | 100 accepted; 2.78 ms | 100 distinct; 15.97 ms |
| 30 bases, 29–30 GC bases, no `CCC` or `GGG` | 0 accepted in 10,000 draws; 24.95 ms | 100 accepted; 1.32 ms | 100 distinct; 39.19 ms |
| 100 bases, 48–52 GC bases, six exclusions | 24 accepted in 10,000 draws; 74.30 ms | 100 accepted; 39.94 ms | 100 distinct; 207.91 ms |
| All four bases individually forbidden | No accepted draw; no infeasibility proof | Exact count zero | Proven infeasible |

The automaton CP-SAT encoding emitted repeated sequence projections from distinct
callback results in this environment. After deduplication it produced 14, 28 and
6 distinct sequences in the three two-second feasible cases. The clause encoding
produced 100 distinct sequences in each. This concerns these encodings and
enumeration settings; it does not establish a general CP-SAT limitation.

Three independent exhaustive oracles contained 514, 44 and 2 valid strings.
Both solver encodings matched each complete set. Counting matched each count,
and unranking every integer produced the lexically ordered oracle exactly.
This establishes the uniform fixed-length prototype's mapping. Native tests
also exhaust weighted integer masses and the ranged-length `16:15` case.

## Native implementation evidence

`parts/background/` separates resource/outcome contracts, forbidden-pattern
transitions, completion counting and compilation. Its iterative graph evaluation
does not depend on Python recursion depth. Existing native constraints remain
the source of GC and strand semantics. No solver or sampling framework was added.

The [native measurement record](measurements/conditional-background.json)
separates uninstrumented timings from a second run measuring traced Python
allocations. It uses the same synthetic inputs as the method comparison and
adds nonuniform probabilities and a case that reaches a resource cap.

| Case | Construction | 100 draws | Peak traced allocation |
| --- | --- | --- | --- |
| Ordinary 20 bases | 8.8 ms | 3.7 ms | 0.39 MB |
| Rare GC, 30 bases | 1.3 ms | 5.2 ms | 0.034 MB |
| Constrained 100 bases | 193 ms | 18.6 ms | 8.2 MB |
| Weighted 100 bases | 104 ms | 18.1 ms | 6.6 MB |
| 1,000 bases, mass cap reached | 737 ms | No draws | 30.8 MB |

These are single-run macOS arm64 measurements, not throughput guarantees or a
whole-process memory cap. They support separate default limits for state
overhead and stored integer bits. The limited case retained no candidates and
reported unknown feasibility.

Tests cover weighted multiplicity, conditional lengths, zero base support,
overlapping patterns on both strands, native GC boundary agreement, 1,500-base
sequences, resource outcomes, reusable tables, candidate-stream stability and
source-free saved evidence. A bounded CP-SAT clause model agrees with the
complete native support on a small fixture; 150 additional deterministic
exhaustive cases matched independent sequence predicates. FIMO remains a
separate screen and can still cause retained-part shortfalls.

### Counting memory

Completion counting keeps completed masses and an explicit traversal stack.
Each child has fewer remaining bases; depth-first traversal completes a state
before another queued copy reaches it. The mass table therefore also identifies
completed states, without a second set retaining every admitted state.

Three paired runs on Python 3.14.7, macOS arm64 measured the effect of removing
that set. The [measurement record](measurements/conditional-admission.json)
includes individual timings, traced allocation peaks, workloads and correctness
checks. Memory measurements ran separately after clearing Python freelists.

| Workload | Previous peak | Current peak |
| --- | ---: | ---: |
| Constrained 100-base compilation and 100 draws | 8.248 MB | 6.169 MB |
| Preparation, verification and report reads | 8.385 MB | 6.444 MB |
| 1,500-base single-support sequence | 0.502 MB | 0.371 MB |
| 10,000-base single-support sequence | 2.549 MB | 2.024 MB |

The reduction is 21–26% in these fixtures. Timings improved slightly; this does
not establish a general throughput gain or whole-process memory bound. Masses,
work counters, seeded sequences and native records matched. Independent weighted
and ranged-length oracles, 1,500 resource-cap comparisons and six deterministic
interruption points also matched.

## Distribution and outcome contracts

1. **Condition the declared proposal distribution.** Preserve base probabilities
   under the constraints. Uniform valid strings apply only to a uniform base
   distribution at one length. Zero-probability bases stay outside sampling
   support. Use exact integer/rational mass or a separately qualified numerical
   method; underflow must not become an infeasibility proof.
2. **Make length conditioning visible.** A uniform prior over lengths is generally
   nonuniform after conditioning on validity. For lengths one or two with `AA`
   forbidden, conditional length probabilities are `16/31` and `15/31`, rather
   than `1/2` each. Preserve and report that distinction. Explicitly uniform
   sampling over feasible lengths would be a different policy.
3. **Keep one constraint owner.** Compile existing sequence-scope `GC` and `Avoid`
   rules, including reverse complements. Reuse those rules for independent
   verification. FIMO screens remain separately bound and can still cause a
   shortfall. Part validity does not establish that later joins or padding
   satisfy final-sequence constraints.
4. **Bound construction and draws.** Limit automaton/counting states, integer
   mass size and elapsed work. Account for table construction separately from
   drawing. Reuse tables within one execution; avoid global caches. A time/state
   limit has no infeasibility meaning. A zero exact count proves no supported
   sequence satisfies the compiled constraints, not biological impossibility.
5. **Preserve solver outcomes.** If CP-SAT is used, distinguish feasible, proven
   infeasible, invalid model, unknown/limited and backend error. Model building,
   external screening and publication need separate accounting; a solver
   allowance is not a whole-operation hard deadline. Follow the
   [documented statuses](https://developers.google.com/optimization/cp/cp_solver)
   and [time-limit controls](https://developers.google.com/optimization/cp/cp_tasks).
6. **Version randomness and artifacts.** Bind draws to model, seed and candidate
   identity. Keep candidate effort, construction work, duplicates and retention
   separate. Persist method/outcome evidence through existing pool artifacts and
   the shared six operations. Preserve earlier requests and random streams.

## Screening cost

Preparation eligibility and conditional verification need a pass/fail result
from literal screens. They stop at the first forbidden match through a shared
iterator. Detailed final-design screening consumes the same iterator fully and
keeps every sorted match, including overlaps and placement intersections.

The [screening measurements](measurements/sequence-screening.json) compare three
paired native preparation, verification and record-reading runs after warmup on
macOS arm64 with Python 3.14.7. Both fixtures contain 64 candidates, use batches
of 16 and declare a retained target of one.

| Synthetic input | Full-match collection median | Predicate median | Traced Python peak, before/after |
| --- | ---: | ---: | ---: |
| Uniform background, 100 bases, six literal exclusions | 28.49 ms | 27.21 ms | 0.332 / 0.332 MB |
| All-A negative control, 1,000 bases, overlapping `AAA`/`AAAA`/`AAAAA` exclusions | 392.61 ms | 91.93 ms | 1.44 / 0.39 MB |

The ordinary case shows a small, noisy difference. The negative control isolates
the cost of collecting thousands of matches that eligibility discards. Allocation
measurements run separately and cover traced Python allocations, not total process
memory. Neither case measures solver search or FIMO throughput.

Plan and pool identities, candidate records and accounting are identical before
and after the change. Exhaustive checks compare boolean and detailed results
for 8,190 sequence/rule pairs. Independent overlap, strand, palindrome, exception
and GC-boundary checks preserve the detailed evidence contract.

## Stochastic proposal cost

Stochastic proposals prepare each required categorical base distribution once
and reuse it for the candidate's draws. Zero-support filtering, cumulative
probability arithmetic, boundary selection and random bytes remain unchanged.
The tables are local to one proposal; they create no persistent cache.

[Three paired measurements](measurements/stochastic-proposals.json) include
planning, preparation, saved-pool verification and quality-report serialization
for 1,000 candidates. Each workload retains at most 100 parts and screens GC
fractions between 0.25 and 0.75.

| Candidate length | Per-base preparation median (range) | Per-proposal preparation median (range) |
| --- | ---: | ---: |
| 20 bases | 148.78 ms (148.35–151.91) | 129.63 ms (128.07–130.90) |
| 100 bases | 278.69 ms (276.28–280.24) | 176.28 ms (171.64–179.61) |

These local macOS/Python 3.14 measurements show reductions of 13% and 37%.
Traced Python allocation peaks remained approximately 1.35 MB and 1.44 MB.
All native identities, accounting and quality records matched, as did 9,000
seeded proposals across motif, background and zero-support inputs. The result
does not measure conditional sampling, solver search or FIMO throughput.

## Distribution and boundary checks

A uniformly chosen next base need not produce uniformly distributed valid
strings. For two-base DNA excluding `AA`, there are fifteen valid strings.
Their first-base probabilities are `3/15`, `4/15`, `4/15` and `4/15` for
`A`, `C`, `G` and `T`. Completion masses supply these probabilities directly.

Resource exhaustion cannot establish infeasibility. Report a limited outcome
when counting reaches its bound without establishing total valid mass.
Similarly, converting GC fractions into integer counts must agree with the
sequence predicate. At length 50 and GC fraction `0.14`, seven GC bases satisfy
the declared fraction. Floating-point multiplication followed by rounding can
produce contradictory bounds; the implementation checks the integer boundary
against the predicate without adding a tolerance.

## Verification and module ownership

`parts/background/` owns constraint compilation, completion masses and
deterministic draws. Preparation contracts own request parsing; execution and
publication retain their existing owners. The method uses the same GC and
literal constraints as final sequence screening.

Independent exhaustive fixtures cover overlapping patterns, both strands,
exact/ranged GC counts, impossible support, weighted masses and ranged-length
conditioning. Regression checks cover integer boundaries, resource limits,
seed/prefix stability and zero-probability support. Bounded CP-SAT checks small
feasible sets independently; it does not serve as the sampling-distribution
oracle.

Qualify saved plans, matching Python/CLI results, self-contained pool
verification, optional FIMO screens, uniqueness, shortfalls and installed-package
execution whenever the method changes. Record resource measurements separately
from correctness evidence.
