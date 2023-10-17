# dense-arrays

This library helps with the design of double-stranded nucleotide sequences made of many overlapping motifs.

Documentation available on [https://dunloplab.gitlab.io/dense-arrays](https://dunloplab.gitlab.io/dense-arrays)

## Installation

1. (optional but recommended) Create an empty conda environment: `conda create -n dense-arrays python` and activate it `conda activate dense-arrays`
2. Install the package with `pip install .`

## Paper figures

To regenerate the paper figures:

1. Recreate the data by executing the `benchmarks/main.py` script (will take a long time, so the result files are already present for convenience).
2. Plot the figures by executing the appropriate plotting script `figures/plot_*.py` as described in `figures/README.md`.

## Simple usage

``` python
import dense_arrays as da

opt = da.Optimizer(["ATGC", "CGT", "ATTA", "TTATTA"], sequence_length=8)

best = opt.optimal()
print(f"Optimal solution, score {best.nb_motifs}")
print(best)

print("List of all solutions")
for solution in opt.solutions():
    print(f"Solution with score {solution.nb_motifs}:")
    print(solution)
```

## Different solver backends

The methods `Optimizer.optimal` and `Optimizer.solutions` take an optional argument for the solver backend. It takes what `ortools` accepts:

* `"CBC"` (default)
* `"SCIP"`
* `"GUROBI"`
* `"CPLEX"`
* `"XPRESS"`
* `"GLPK"`
