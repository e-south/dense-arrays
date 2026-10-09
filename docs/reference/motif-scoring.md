---
title: Motif inputs and scoring
description: Import JSON, MEME or JASPAR motifs and score candidates with explicit units and bounded effort.
---

# Read motifs and score candidates

Use `dense_arrays.parts.motifs` to read motif models from JSON artifacts,
minimal MEME files or JASPAR frequency matrices. Use
`dense_arrays.parts.scoring` to evaluate candidates with
FIMO. These Python interfaces operate on explicit candidates. To generate and
retain candidates through Python or CLI, use a
[sampled preparation recipe](../library-workflow/preparation.md).

## Read a motif artifact

`parts.PWMArtifact(path, motif_ids=())` names one JSON file containing one motif.
Its optional [`window`](motif-windows.md) is applied by preparation planning;
`read_artifact` always returns the complete original model.
The optional selector must match that motif's ID. `read_artifact` returns a
`MotifInput` with an immutable `motif`, an absolute source path and its SHA-256
fingerprint. Call `.verify()` to check the source bytes before subsequent work.

The required fields are `schema_version: "1.0"`, `producer`, `motif_id`,
`alphabet: "ACGT"`, `matrix_semantics: "probabilities"`, `background`,
`probabilities` and `log_odds`. Matrix rows and background objects have exactly
`A`, `C`, `G`, `T` keys. Both matrices have the same nonzero width. Probability
rows must sum to one within 0.001; near-unit rows are normalized. Background
frequencies must be positive. Score entries must be finite numbers. An optional
`length` must match the width; additional producer annotations remain metadata.
Duplicate JSON keys, numeric strings and boolean matrix entries fail validation.

Use the 12-position `motif.json` from [the motif recipe](../library-workflow/preparation/motifs.md#create-a-motif-artifact).

```python
# Use the typed requests and operations needed by this example.
from dense_arrays import parts
from dense_arrays.parts.motifs import read_artifact, score_core, best_hit

source = read_artifact(
    parts.PWMArtifact("motif.json")
)  # Load the complete declared model.
source.verify()  # Confirm the bound input bytes are unchanged.
motif = source.motif  # Keep the immutable model for the two scoring examples.
score = score_core(motif, "ACGTTGCAAGTC")  # Score exactly one motif-width core.
hit = best_hit(motif, "GATTACGTTGCAAGTCTCGA", strands="double")  # Scan both strands.
```

Run in the directory containing `motif.json`. `score_core` requires exactly one
motif-width uppercase A/C/G/T core. `best_hit` checks every full-width window;
`single` checks only the supplied strand and `double` checks both orientations.
Its coordinates are zero-based, half-open intervals in the candidate. The core
is oriented to the motif. Ties choose the earliest interval, then forward strand.

These scores use the supplied `log_odds` entries exactly. Their units are
`declared_log_odds`; no logarithm base, pseudocount or p-value is inferred.
`raw`, `per_base`, `theoretical_max` and `fraction_of_max` have separate fields.
The fraction is `None` when the maximum is nonpositive. `Motif.model_id` binds
both matrices and the background, independently of labels or producer metadata.

## Import MEME or JASPAR

Choose the input format explicitly with
`parts.PWMArtifact(path, format="meme", motif_ids=("MA0001.1",))` or
`format="jaspar"`. File extensions do not select a parser. Both readers support
multiple records with unique IDs; preparation selects exactly one. Omit
`motif_ids` only for a single-record file. Selection is case-sensitive and uses
the primary ID, never an alternate name. Every record is validated before
selection, including records not selected. See the
[preparation example](../library-workflow/preparation/motifs.md#import-a-meme-or-jaspar-motif)
for Python and CLI use.

| Input | Supported content | Preserved evidence |
| --- | --- | --- |
| Minimal MEME | ACGT letter-probability matrices, optional width, `nsites`, `E`, `S`, alternate name and URL | Source probabilities, statistics, declared background and strands |
| JASPAR | Four bracketed, labeled A/C/G/T count rows or four unlabelled rows in A/C/G/T order | Source counts and position totals; probabilities normalized per position |

MEME probabilities are never converted to rounded counts. JASPAR normalization
adds no pseudocounts. Counts and `nsites` are source annotations, not evidence
that observations were independent. Missing backgrounds use uniform A/C/G/T
frequencies with `background_origin="uniform_default"` in metadata. An omitted
MEME alphabet uses ACGT and records that default. The supplied strand annotation
does not replace the explicit scoring setting. These choices follow the
[minimal MEME format](https://meme-suite.org/meme/doc/meme-format.html) and
[JASPAR matrix conventions](https://biopython.org/docs/latest/Tutorial/chapter_motifs.html#jaspar).

Both imports produce a probability-only `Motif`: `log_odds` and `score_units`
are `None`. `score_core` and `best_hit` fail clearly without a supplied matrix;
use the configured FIMO backend below. Probability-only models use motif
schema v2; models with supplied scores retain schema v1 and their identities.
Window selection preserves the absence of a supplied score matrix.

The readers accept UTF-8 text, blank separators and whole-line `#` comments.
Non-DNA alphabets, duplicate IDs, invalid probabilities, inconsistent matrix
widths and zero-total count positions fail. Full MEME discovery reports,
log-odds-only files and unlabeled matrices without a motif header are not
supported. Export a minimal probability or JASPAR count file before importing.

## Score candidates with FIMO

Use the optional [FIMO installation route](../installation.md#configure-fimo-for-motif-scoring).
Importing Dense Arrays does not discover or invoke FIMO. From the same directory:

```python
from dense_arrays.parts.scoring import bind_fimo, scan_fimo

settings = parts.FimoScoring(
    hit_pvalue_max=0.1,
    strands="double",
    limits=parts.ScoringLimits(seconds=30, windows=1000, output_bytes=1048576),
)
binding = bind_fimo(motif, settings)  # Bind the executable and scorer settings.
# Score two complete 20-base candidates against the bound 12-position motif.
result = scan_fimo(binding, ("GATTACGTTGCAAGTCTCGA", "TCGAGATCCGTAAGCTGTCA"))
for hit in result.hits:
    if hit is not None:
        print(hit.start, hit.end, hit.strand, hit.raw, hit.pvalue)
```

`executable=path` selects a binary explicitly; otherwise preflight resolves
`fimo` on `PATH`. Preflight queries only the version and fingerprints the binary.
No candidate is scored. Each scan verifies the executable and optional background
file before and after execution. Saved bindings can be loaded with
`FimoBinding.from_dict(...)` without accessing either file or invoking the tool.

The scorer uses motif probabilities, a declared background, pseudocount 0.1 and
site count 20. This is an explicit scoring convention; source `nsites` and
JASPAR position totals do not override it. The pseudocount is configurable.
Probability rows are written
with 17 significant digits. The artifact's supplied `log_odds` matrix is not
used by FIMO. `background=path` accepts a zero-order A/C/G/T frequency file;
omitting it uses the motif background. Double-strand scans explicitly average
complementary frequencies; the binding records original and effective values.

Each candidate receives one best qualifying hit or `None`. Thresholding occurs
inside FIMO; the adapter never reclassifies a borderline hit using the rounded
TSV p-value. Text mode provides no q-value. Geometry and tie rules match the
supplied-matrix scanner. Unknown IDs, duplicate hits, invalid coordinates,
nonfinite scores and mismatched cores fail the entire call.

`hit_pvalue_max` applies to each oriented motif-width window under the declared
zero-order background. The reported p-value belongs to that hit; selecting the
best hit does not turn it into a whole-candidate p-value or a library false
discovery rate. Longer candidates and additional strands create more windows
that can pass the threshold. This scope follows
[FIMO's occurrence-level scoring](https://pmc.ncbi.nlm.nih.gov/articles/PMC3065696/).
Generated or selected candidates need not follow the scoring background, so
their observed hit frequency can differ from its null expectation.

FIMO hits label their log2 likelihood-ratio scores `fimo_log2_odds`. A separate
maximizing core is scored using the same model and backend to obtain its reported
theoretical maximum. The result separates candidate windows from this calibration work.
A flat motif can have maximum zero and an unavailable fraction.
`fraction_of_max` divides two log scores; it is not a fraction of binding
probability or affinity. Scores and p-values describe the computation, not
biological activity.

## Work limits and failures

`ScoringLimits` defaults to 60 seconds, 1,000,000 oriented windows and 64 MiB of
combined standard output/error. Window admission includes the calibration core
and happens before scoring. Both subprocesses share the time and output budget.
The child is killed and reaped on timeout or excess output; parsing and input
verification check elapsed time between operations. This is not a hard deadline
for the entire Python call. Output buffering is bounded by the byte cap; total
Python and backend memory also depend on the admitted inputs and motif width.

`ScoringError.reason` distinguishes `unavailable`, `timeout`, `output_limit`,
`backend` and `malformed`. A missing qualifying hit is a successful scan with
`None` for that candidate. Changed bound files raise `ValueError`. Failed calls
return no partial result. `FimoHit` and `FimoResult` support `to_dict`/`from_dict`
for saved evidence; decoding checks derived fields and work accounting.
