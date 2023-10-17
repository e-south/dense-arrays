"""
Plot "performance.svg" (Figure 2) and "approximation.svg" (Figure S1)
from the file "benchmarks_double.csv".
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data = pd.read_csv("../benchmarks/benchmarks_double.csv").rename(
    columns={
        "library_size": "Library size |R|",
        "opt_solve_time": "Solve time (s)",
        "sequence_length": "Sequence length L (bp)",
    }
)

print(data[data.opt_nb_motifs == 0])

data = data[data["opt_nb_motifs"] > 0]


def plot_performance() -> None:
    """Plot the performance of the exact algorithm."""
    fig, axs = plt.subplots(1, 2, sharey=True, figsize=(8, 5))

    sns.lineplot(
        data=data,
        x="Library size |R|",
        y="Solve time (s)",
        hue="Sequence length L (bp)",
        ax=axs[0],
    )

    sns.lineplot(
        data=data,
        x="Sequence length L (bp)",
        y="Solve time (s)",
        hue="Library size |R|",
        ax=axs[1],
    )

    axs[0].set_yscale("log")

    fig.savefig("performance.png", bbox_inches="tight")
    fig.savefig("performance.svg", bbox_inches="tight")


def plot_approximation() -> None:
    """Plot the performance of the approximation algorithm."""
    fig, axs = plt.subplots(1, 2, sharey=True, figsize=(8, 5))

    sns.lineplot(
        data=data,
        x="Library size |R|",
        y=data["approx_nb_motifs"] / data["opt_nb_motifs"],
        hue="Sequence length L (bp)",
        ax=axs[0],
    )

    sns.lineplot(
        data=data,
        x="Sequence length L (bp)",
        y=data["approx_nb_motifs"] / data["opt_nb_motifs"],
        hue="Library size |R|",
        ax=axs[1],
    )

    axs[0].set_ylabel(
        "Proportion of binding sites included by the greedy algorithm\n"
        "compared to the optimal solution"
    )

    fig.savefig("approximation.png", bbox_inches="tight")
    fig.savefig("approximation.svg", bbox_inches="tight")


plot_performance()
plot_approximation()
