---
title: Motif windows
description: Choose a bounded motif window with explicit information, background and coordinate semantics.
author: Eric J. South
---

# Select a motif window

`parts.PWMArtifact(path, window=parts.MotifWindow(length=N))` selects a
contiguous, fixed-width window during preparation planning. Omitting `window`
preserves the complete motif. Length must be positive and cannot exceed the
source width. Sampling length independently controls the complete candidate,
including any flanks. See the [preparation guide](../library-workflow/preparation/windows.md#choose-a-motif-window)
for matching Python and CLI execution.

For several widths, [prepare named windows](../library-workflow/preparation/windows.md#prepare-named-motif-windows).
Each expanded recipe selects directly from the original model and records its
own window, calibration and retained population. Budgets and targets are per
recipe; scores from different windows are not pooled into one ranking.

## Selection objective

The selector maximizes the sum of per-column relative entropies in bits:

```text
I_i = sum over A,C,G,T of P_i(base) * log2(P_i(base) / B(base))
selected start = earliest argmax_start sum(I_i over the requested window)
```

Zero-probability terms contribute zero. The calculation uses the supplied
probability model with no inferred observation counts or additional smoothing.
Exact score ties choose the earliest source position. Supplied log-odds rows
are sliced at the same coordinates; they do not determine the window.

| `MotifWindow.background` | Reference distribution for selection |
| --- | --- |
| `"motif"` (default) | The artifact's declared background. |
| `"uniform"` | Equal A/C/G/T probabilities; the metric reduces to `2 − entropy`. |
| Four A/C/G/T probabilities | An explicit positive, normalized distribution. |

This distribution belongs to window selection. Sampling probabilities and
FIMO's scoring background retain their separate declarations; changing either
does not silently alter the selected window.

Relative entropy measures expected discrimination between a known probability
model and its background. Under the independent-position model, maximizing its
fixed-width sum retains the most of that expected discrimination among the
allowed contiguous windows. This is an application of the log-odds interpretation
in [Yu et al., Methods §2.2](https://pmc.ncbi.nlm.nih.gov/articles/PMC4318935/).
It does not optimize binding affinity, sensitivity at a chosen threshold, or
an experimentally established functional core.

## Saved evidence and scoring

`plan.preview["motif_window"]` reports source and selected model identities,
zero-based half-open `start`/`end`, the effective background, selected and total
information, discarded bits and retained fraction. The fraction is `None` when
the source has zero information against the reference. A full-width request
reports its actual information and leaves the computational model unchanged.
Exclusion windows appear in `screening_windows`, identified by rule and motif.

The native plan retains the original source model and fingerprint. Proposals,
FIMO calibration and core-diversity selection use the selected model. Scored-hit
coordinates refer to the candidate sequence, independently of source-window
coordinates. Source-free inspection verifies this transformation without
running FIMO. Its log-odds scores and p-values are calculated for the selected
motif, as required by [FIMO's model-dependent scoring](https://meme-suite.org/meme/doc/fimo.html).
They are not scores against the discarded full-width model.

## Limits of the interpretation

Automatic biological boundary inference needs additional evidence. Probability
matrices alone do not provide the independent observation counts needed for
finite-sample corrections or the count-based boundary methods in
[Yu et al., Methods §§2.3–2.7](https://pmc.ncbi.nlm.nih.gov/articles/PMC4318935/).
The fixed site count used to serialize a FIMO model is a scoring convention,
not an estimate of the motif's independent observations.
Unpenalized relative entropy is nonnegative: allowing arbitrary width and
maximizing its sum would simply retain the whole motif.

Flanking sequence can affect specificity: experiments on two yeast bHLH factors
demonstrated context-dependent binding around shared E-box motifs.
[Gordân et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3640701/)
Information retained therefore describes the selected probability model, not
retained biological function. Background-sampled flanks do not reconstruct
discarded motif positions. Assess binding-preservation claims with independent
binding data and an explicitly defined evaluation population.
