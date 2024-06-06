# Data files

All the csv files in this directory contain results and benchmarks of the
algorithm on synthetic, randomly-generated promoter-like DNA sequences.

The association of these files to paper figures is explained by the README.md
in the `figures` directory.

The randomness in the benchmarks is always reproducible due to seeding the
pseudo-random number generator.

The data files produced are the following:

## Function `benchmark`

* Creates the files `benchmarks_single.csv` and `benchmarks_double.csv`

Each line represents solving the SPP on a different randomly generated library,
both approximately and exactly.

The columns include characteristics of the initial library (individual motif
sizes, library size), of the desired sequence L, and of the approximate and
exact solution as well as their timing.

## Function `topsols`

* Creates the files `topsols10_{i}_solver.csv` and `topsols10_{i}_control.csv`

Each file indexed by `i` corresponds to generating all top-scoring solutions
on a different promoter library. The first line corresponds to the motifs of
the library, and each of the following lines corresponds to whether a given
motif is or not included in the solution.

This function has been called with three different solvers, to investigate
solver transient bias.  Gurobi was the most efficient solver, followed by SCIP
and CBC:

+-----------+----------+----------+-----------+
|           |  Gurobi  |   SCIP   |    CBC    |
+-----------+----------+----------+-----------+
| Library 0 |   11m20s | 2h46m46s | 16h02m48s |
+-----------+----------+----------+-----------+
| Library 1 |   22m04s | 6h19m44s |   > 20h   |
+-----------+----------+----------+-----------+
| Library 2 | 3h29m28s |  > 20h   |   > 20h   |
+-----------+----------+----------+-----------+

## Function `multiple_promoters`

* Creates the file `multiple_promoters.csv`

Each couple of lines corresponds to solving the SPP on a different library.
The first column is the index of the problem (some are missing because of
solver timeouts).  The odd-numbered lines describe the library promoters, and
the even-numbered lines describe their offsets from the start of the desired
sequence, if present in the best solution.

## Function `side_bias`

* Creates the file `side_bias_both.csv`

Each couple of lines corresponds to solving the SPP on a different library, with
promoter constraints. The last two motifs are always the upstream promoter and
the downstream promoter.  The first column is the index of the problem (some
are missing because of solver timeouts).  The second is whether the problem was
solved with or without activating the side bias.  The rest are the offsets of
the motifs from the start of the sequence in the best solution.

## Function `size_bias`

* Creates the files `size_bias_{i}_control.csv` and `size_bias_{i}_solver.csv`

Same format as for `topsols`.
