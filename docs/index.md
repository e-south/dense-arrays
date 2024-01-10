# Dense Arrays

This library facilitates the design of double-stranded nucleotid sequences made of many overlapping motifs.

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

``` python
import dense_arrays as da

opt = da.Optimizer(["ATGC", "CGT", "ATTA", "TTATTA"], sequence_length=8)

best = opt.optimal()
print(f"Optimal solution, score {best.nb_motifs}")
print(best)

print("List of all solutions")
for solution in opt.solutions():
    print(solution)
```

## Citation

If you use this package in your research, please cite the following preprint:

> *Generating information-dense nucleotide sequences with optimal string packing* Virgile Andreani, Eric J. South, Mary J. Dunlop, 2023, https://doi.org/10.1101/2023.11.01.565124
