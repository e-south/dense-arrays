"""Script to run to generate benchmarks."""

import itertools as it
import time
from pathlib import Path

import dense_arrays as da
import numpy.random as npr

ATGC = list("ATGC")


def benchmarks(*, double: bool) -> None:
    """
    Create the "benchmarks_single.csv" and "benchmarks_double.csv" files.

    They are the data plotted in Figure 2 of the paper.
    """
    rng = npr.default_rng(seed=42)
    solver = "Gurobi"
    motif_min_size = 5
    motif_max_size = 15
    library_sizes = [10 * i for i in range(1, 11)]
    sequence_lengths = [20 * i for i in range(1, 16)]
    with Path("benchmarks_" + ["single", "double"][double] + ".csv").open("w") as out:
        print(
            "motif_min_size,motif_max_size,library_size,nb_motif_letters,sequence_length,approx_nb_motifs,approx_solve_time,opt_nb_motifs,opt_solve_time,double,solver",
            file=out,
            flush=True,
        )
        for library_size, sequence_length, _replicate in it.product(
            [library_sizes, sequence_lengths, range(10)]
        ):
            motifs = [
                "".join(
                    rng.choice(
                        ATGC,
                        rng.integers(motif_min_size, motif_max_size, endpoint=True),
                    )
                )
                for _ in range(library_size)
            ]
            nb_motif_letters = sum(len(motif) for motif in motifs)
            tic = time.time()
            optimizer = da.Optimizer(
                motifs, sequence_length, strands=["single", "double"][double]
            )
            opt = optimizer.optimal(solver, solver_options=["TimeLimit=600"])
            opt_solve_time = time.time() - tic
            tic = time.time()
            approx = optimizer.approximate()
            approx_solve_time = time.time() - tic
            opt_nb_motifs = opt.nb_motifs if opt else 0
            print(
                f"{motif_min_size},{motif_max_size},{len(motifs)},{nb_motif_letters},{sequence_length},{approx.nb_motifs},{approx_solve_time},{opt_nb_motifs},{opt_solve_time},{double},{solver}",
                file=out,
                flush=True,
            )


def topsols(*, control: bool) -> None:
    """
    Create the "topsols10_{i}_solver.csv" and "topsols10_{i}_control.csv" files.

    Each is a csv with motifs in columns and solutions in rows.
    They are the data plotted in Figures 3 and 4 of the paper.
    """
    rng = npr.default_rng(seed=42)

    solver = "Gurobi"
    library_size = 10
    sequence_length = 50

    for i in range(1):
        motifs = ["".join(rng.choice(ATGC, 10)) for _ in range(library_size)]
        with Path(
            f"topsols10_{i}_" + ["solver", "control"][int(control)] + ".csv"
        ).open("w") as out:
            nodes = motifs + [da.optimize.reverse_complement(motif) for motif in motifs]
            print(",".join(nodes), file=out)

            optimizer = da.Optimizer(motifs, sequence_length)
            opt = optimizer.optimal(solver)
            sols = (
                optimizer.solutions_diverse(solver)
                if control
                else optimizer.solutions(solver)
            )
            for n, sol in enumerate(sols):
                print(f"{n}th solution, score {sol.nb_motifs}")
                if sol.nb_motifs < opt.nb_motifs:
                    break
                print(
                    *(int(x is not None) for x in sol.offsets_fwd + sol.offsets_rev),
                    sep=",",
                    file=out,
                    flush=True,
                )


def size_bias() -> None:
    """
    Create the "size_bias_{i}_control.csv" files.

    Each is a csv with motifs in columns and solutions in rows.
    They are the data plotted in Figures S2 and S3 of the paper.
    """
    rng = npr.default_rng(seed=42)
    library_size = 10
    sequence_length = 30
    solver = "Gurobi"
    for replicate in range(10):
        motifs = ["".join(rng.choice(ATGC, 5 + i)) for i in range(library_size)]
        with Path(f"size_bias_{replicate}_control.csv").open("w") as out:
            nodes = motifs + [da.optimize.reverse_complement(motif) for motif in motifs]
            print(",".join(nodes), file=out, flush=True)
            optimizer = da.Optimizer(motifs, sequence_length, strands="double")
            opt = optimizer.optimal(solver)
            for n, sol in enumerate(optimizer.solutions_diverse(solver)):
                print(f"{n}th solution, score {sol.nb_motifs}")
                if sol.nb_motifs < opt.nb_motifs:
                    break
                print(
                    *(int(x is not None) for x in sol.offsets_fwd + sol.offsets_rev),
                    sep=",",
                    file=out,
                    flush=True,
                )


benchmarks(double=True)

topsols(control=False)
topsols(control=True)

size_bias()
