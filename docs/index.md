# Dense Arrays

This library facilitates the design of double-stranded nucleotide sequences made of many motifs.

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

library = [
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

opt = da.Optimizer(library=library, sequence_length=100)

best = opt.optimal()
print(f"Optimal solution, containing {best.nb_motifs} motifs")
print(best)

# Optimal solution, score 8
# -->          GCTTAAAAAATGAAC     ACTGAATTTTTATGCAAA            AATATGTAACCAAAAGTAA
# --> CGGGGATGAG              TTGACA                          TATAAT
# --> CGGGGATGAGCTTAAAAAATGAACTTGACACTGAATTTTTATGCAAATCAATCAATTATAATATGTAACCAAAAGTAATTTTCTTATAGGGA--------
# <-- GCCCCTACTCGAATTTTTTACTTGAACTGTGACTTAAAAATACGTTTAGTTAGTTAATATTATACATTGGTTTTCATTAAAAGAATATCCCT--------
# <--                                              TTAGTTAGTTAAT                 ATTAAAAGAATATCCCT

print("List of all solutions")
for solution in opt.solutions():
    print(f"Solution containing {solution.nb_motifs} motifs")
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

You can designate certain binding sites within the provided library as special, and thereby stipulate added constraints to the solver, such as having these sites appear at specific or rough positions within the final dense array sequence. While overly stringent constraints risk yielding no feasible solution, you can generate sequences with pairs of fixed binding sites, such as sigma factor recognition elements appearing at the -35 and -10 sites relative to the transcription start site (e.g., the terminus of L). Sigma factor recognition sites are crucial elements in bacterial promoters.

``` python
# Adding another constraint to the previous optimizer
opt.add_promoter_constraints(
    upstream="TTGACA",
    downstream="TATAAT",
    upstream_pos=(40, 80),
    spacer_length=(16, 18),
)

best = opt.optimal()
print(best)

#                                                                position 60  spacer length: 17
#                                                                v    |<--------------->|
# -->                                                            TTGACA                 TATAAT
# --> TCCCTATAAGAAAATTA                               TAATTGATTGATT   ACTGAATTTTTATGCAAA
# --> TCCCTATAAGAAAATTACTTTTGGTTACATATTGTTCATTTTTTAAGCTAATTGATTGATTGACACTGAATTTTTATGCAAATATAATTCAGAATATTAT
# <-- AGGGATATTCTTTTAATGAAAACCAATGTATAACAAGTAAAAAATTCGATTAACTAACTAACTGTGACTTAAAAATACGTTTATATTAAGTCTTATAATA
# <--               AATGAAAACCAATGTATAA                                                    TTAAGTCTTATAATA
# <--                                  CAAGTAAAAAATTCG

```

In some cases, motifs may be associated with labels indicating their role as binding sites for transcription factor activators or repressors. Studies analyzing natural sequences have shown that activators or repressors often occupy specific areas within cis-regulatory regions. For example, activator binding sites are typically located upstream of the -35 sigma factor recognition element in bacterial promoters. Notably, the regulatory effectiveness of these sites diminishes when they are positioned downstream of this element. The solver can be adapted to mirror these biological patterns by allowing for the specification of motifs with upstream or downstream preferences. As the solver traverses the graph to assemble the output sequence, a position variable of each motif encountered is integrated into the scoring function. Here, position[i] signifies the starting point of motif i in the sequence. By dividing these position variables by a large constant K, the solver is configured to find solutions that favor positioning activator motifs upstream and repressor motifs downstream. This method, which can be combined with other promoter constraints, enables the generation of sequences that more accurately mimic the intricate arrangement of motifs in natural, multi-factor promoters.

``` python
motifs = [
    "GAAATAACATAATTGA",
    "TGTTAATAATAAGTAAT",
    "TTATATTTTACCCATTT",
    "AGGTTAATCCTAAAA",
    "ATTGAAACGATTCAGC",
    "CTCTGTCATAAAACTGTCATAT",
    "TTACGCATTTTTAC",
    "ATTTGTACACAA",
    "AAGGCATAACCTATCACTGT",
    "ACGCAAACGTTTTCTT",
    "TACATTTAGTTACA",
    "TTAATAAAACCTTAAGGTT",
    "CCTTTTAGGTGCTT",
    "TACTGTATATAAAAACAGTA",
    "TAAAATTCATGGTAATTAT",
    "AATGAGAATGATTATTAT",
    "TGTTTATATTTTGTTTA",
    "CATAAGAAAAA",
    "CATTCATTTG",
    "TTGACA",
    "TATAAT",
    "TATACT",
    "TGGCAGG",
    "TTGCA"
]

# Subset of motifs which you would prefer to be placed upstream
upstream = [
        "GAAATAACATAATTGA",
        "TGTTAATAATAAGTAAT",
        "TTATATTTTACCCATTT",
        "AGGTTAATCCTAAAA",
        "ATTGAAACGATTCAGC",
        "CTCTGTCATAAAACTGTCATAT",
        "TTACGCATTTTTAC",
        "ATTTGTACACAA",
        "AAGGCATAACCTATCACTGT",
        "ACGCAAACGTTTTCTT",
]

# Subset of motifs which you would prefer to be placed downstream
downstream = [
        "TTAATAAAACCTTAAGGTT",
        "CCTTTTAGGTGCTT",
        "TACTGTATATAAAAACAGTA",
        "TAAAATTCATGGTAATTAT",
        "AATGAGAATGATTATTAT",
        "TGTTTATATTTTGTTTA",
        "CATAAGAAAAA",
        "CATTCATTTG",
]

# Adding a promoter constraint: sigma D -35 and -10 motifs
opt.add_promoter_constraints(
    upstream="TTGACA",
    downstream="TATAAT",
    upstream_pos=(40, 80),
    spacer_length=(16, 18),
)

# Adding a promoter constraint: sigma E -24 and -12 motifs
opt.add_promoter_constraints(
    upstream="TGGCAGG",
    downstream="TTGCA",
    upstream_pos=(40, 80),
    spacer_length=(3, 5),
)

# Declare upstream and downstream preferences
opt.add_side_biases(
    left=upstream,
    right=downstream
)

best = opt.optimal()
print(best)

# Upstream preference: *                       position 42  spacer length: 18
# Downstream preference: +                     v    |<---------------->|     position 72  spacer length: 5
#                                                   |                  |     v     |<--->|
#                                                   |                  |           |     |
#                                                   |                  |           |     |
# -->                       TACATTTAGTTACA     TTGACA                  |     TGGCAGG     TTGCA
# -->           *TTACGCATTTTTAC        +CATTCATTTG+CATAAGAAAAA         TATAAT      TATACT
# --> TTGTGTACAAATTACGCATTTTTACATTTAGTTACATTCATTTGACATAAGAAAAATGGGTAAAATATAATGGCAGGTATACTTGCAAGCACCTAAAAGG
# <-- AACACATGTTTAATGCGTAAAAATGTAAATCAATGTAAGTAAACTGTATTCTTTTTACCCATTTTATATTACCGTCCATATGAACGTTCGTGGATTTTCC
# <-- AACACATGTTTA*                                        TTTACCCATTTTATATT*               TTCGTGGATTTTCC+
```

## Citation

If you use this package in your research, please cite the following preprint:

> *Generating information-dense nucleotide sequences with optimal string packing* Virgile Andreani, Eric J. South, Mary J. Dunlop, 2023, https://doi.org/10.1101/2023.11.01.565124
