# Dense Arrays

This library facilitates the design of double-stranded nucleotid sequences made of many motifs.

## Installation

Clone the repository with

```
$ git clone https://gitlab.com/dunloplab/dense-arrays
$ cd dense-arrays
```

Then create conda or pip environment (optional but recommended), and install the package with

```
$ pip install .
```

## Quick usage

In the simplest use case, the user can provide a list of motifs (e.g., transcription factor binding sites), along with a maximum length of the output sequence. The solver will return a sequence that accommodates a maximal number of the provided binding sites. If the specified length is too short to accommodate all provided sequences, then the returned sequence will only contain a subset of the provided binding sites.

``` python
import dense_arrays as da

motifs = [
    "ATAATATTCTGAATT",
    "TCCCTATAAGAAAATTA",
    "TAATTGATTGATT",
    "GCTTAAAAAATGAAC",
    "TGCACTAAAATGGTGCAA",
    "AATATGTAACCAAAAGTAA",
    "ACTGAATTTTTATGCAAA",
    "CGGGGATGAG",
]

opt = da.Optimizer(library=motifs, sequence_length=75)

best = opt.optimal()
print(f"Optimal solution, score {best.nb_motifs}")
print(best)

print("List of all solutions")
for solution in opt.solutions():
    print(solution)
```

If you want to generate multiple sequences composed of binding sites from the same library, and not just find the "densest solution," you can list all solutions found.

## Generating a Library of Diverse Solutions

Sometimes, depending on the motifs provided, the creation of dense arrays inherently tends to favor the inclusion of specific motifs, which is due to the unique interplay among sites within each specific library. To mitigate this, the user can stipulate either the default solver order of solutions returned or a modified diversity-driven order, such that the sequence of returned solutions tends to contain dense arrays composed of new combinations of motifs.

``` python
print("Return solutions in an order that transiently improves inclusivity in binding site representation.")
for solution in opt.solutions_diverse():
    print(solution)
```

## Incorporating Positional Constraints into Solutions

You can designate certain binding sites within the provided library as special, and thereby stipulate added constraints to the solver, such as these sites appearing at specific or rough positions within the final dense array sequence. While overly stringent constraints risk yielding no feasible solution, you can generate sequences with pairs of fixed binding sites, such as sigma factor recognition elements appearing at the -35 and -10 relative to the transcription start site (i.e., the terminus of L). Sigma factor recognition sites are crucial elements in bacterial promoters.

``` python
motifs = [
    "ATAATATTCTGAATT",
    "TCCCTATAAGAAAATTA",
    "TAATTGATTGATT",
    "GCTTAAAAAATGAAC",
    "TGCACTAAAATGGTGCAA",
    "AATATGTAACCAAAAGTAA",
    "ACTGAATTTTTATGCAAA",
    "CGGGGATGAG",
    "TTGACA",
    "TATAAT",
]

opt = da.Optimizer(
    library=motifs,
    sequence_length=75,
    strands="double",
)

opt.add_promoter_constraints(
    upstream="TTGACA",
    downstream="TATAAT",
    upstream_pos=(5, 10),
    spacer_length=(3, 13),
)

best = opt.optimal()
```

In some cases, motifs may come with labels, such as being associated with transcription factor activators or repressors. In such cases, these labels typically suggest a motif's natural placement within specific subregions of a gene's cis-regulatory region. For instance, activator binding sites are typically found upstream of the -35 sigma factor recognition element in bacterial promoters. Their effectiveness in gene regulation tends to decrease when positioned downstream of this element. The solver can be adapted to reflect these biological realities. By adjusting the initial edge weights in the SPP graph, the solver can be configured to preferentially select activator motifs for upstream positions. This adaptive weighting approach enhances the ability to create sequences that more accurately represent the complex organization of natural, multi-factor promoters, with site-specific enrichment of different types of binding motifs.

``` python
# Example code for adaptive weighting to favor activator binding sites
# [placeholder for usage code]
```

## Citation

If you use this package in your research, please cite the following preprint:

> *Generating information-dense nucleotide sequences with optimal string packing* Virgile Andreani, Eric J. South, Mary J. Dunlop, 2023, https://doi.org/10.1101/2023.11.01.565124
