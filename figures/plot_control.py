"""
Plot files "control.svg" (Figure 4) and "size_bias_control.svg" (Figure S3).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from scipy import stats as st


def homogeneous() -> None:
    """
    Plot file "control.svg" (Figure 4) from the files "topsols10_{i}_control.csv".
    """
    fig, axs = plt.subplots(3, 3, sharex=True, sharey="row", figsize=(10, 8))

    for i in range(3):
        data = pd.read_csv(f"../benchmarks/topsols10_{i}_control.csv").to_numpy()
        data = data[:, : data.shape[1] // 2] + data[:, data.shape[1] // 2 :]
        cumdist = np.cumsum(data, axis=0).astype(np.float64)
        cumdist /= cumdist.sum(axis=1, keepdims=True)

        ps = np.arange(cumdist.shape[0]) + 1
        max_entropy = [-np.log(1 / data.shape[1])] * len(ps)
        axs[2, i].plot(ps, max_entropy, "--", label="Maximum entropy", color="gray")

        axs[1, i].plot(ps, cumdist)
        entropy = st.entropy(cumdist, axis=1)
        axs[2, i].plot(ps, entropy, label="Controlled order", color="C1")

        data = pd.read_csv(f"../benchmarks/topsols10_{i}_solver.csv").to_numpy()
        data = data[:, : data.shape[1] // 2] + data[:, data.shape[1] // 2 :]
        cumdist = np.cumsum(data, axis=0).astype(np.float64)
        cumdist /= cumdist.sum(axis=1, keepdims=True)

        axs[0, i].plot(ps, cumdist)
        entropy = st.entropy(cumdist, axis=1)
        axs[2, i].plot(ps, entropy, label="Solver order", color="C0")

        axs[0, i].set_title(f"Collection {i + 1}")

    axs[0, 0].set_xscale("log")
    axs[0, 1].set_xscale("log")
    axs[0, 2].set_xscale("log")
    axs[0, 0].set_ylabel("Solver order\nBinding sites representation")
    axs[1, 0].set_ylabel("Controlled order\nBinding sites representation")
    axs[2, 0].set_ylabel("Distribution entropy")
    axs[2, 0].set_xlabel("Number of generated solutions")
    axs[2, 1].set_xlabel("Number of generated solutions")
    axs[2, 2].set_xlabel("Number of generated solutions")
    axs[2, 2].legend()

    fig.savefig("control.png", bbox_inches="tight")
    fig.savefig("control.svg", bbox_inches="tight")


def inhomogeneous() -> None:
    """
    Plot file "size_bias_control.svg" (Figure S3)
    from the files "size_bias_{i}_control.csv".
    """
    fig, axs = plt.subplots(3, 3, sharex=True, sharey="row", figsize=(10, 8))

    for i in range(3):
        data = pd.read_csv(f"../benchmarks/size_bias_{i}_control.csv").to_numpy()
        data = data[:, : data.shape[1] // 2] + data[:, data.shape[1] // 2 :]
        cumdist = np.cumsum(data, axis=0).astype(np.float64)
        cumdist /= cumdist.sum(axis=1, keepdims=True)

        ps = np.arange(cumdist.shape[0]) + 1
        max_entropy = [-np.log(1 / data.shape[1])] * len(ps)
        axs[2, i].plot(ps, max_entropy, "--", label="Maximum entropy", color="gray")

        axs[1, i].plot(ps, cumdist)
        entropy = st.entropy(cumdist, axis=1)
        axs[2, i].plot(ps, entropy, label="Controlled order", color="C1")

        data = pd.read_csv(f"../benchmarks/size_bias_{i}_solver.csv").to_numpy()
        data = data[:, : data.shape[1] // 2] + data[:, data.shape[1] // 2 :]
        cumdist = np.cumsum(data, axis=0).astype(np.float64)
        cumdist /= cumdist.sum(axis=1, keepdims=True)

        axs[0, i].plot(ps, cumdist)
        entropy = st.entropy(cumdist, axis=1)
        axs[2, i].plot(ps, entropy, label="Solver order", color="C0")

        axs[0, i].set_title(f"Collection {i + 4}")

    axs[0, 0].set_xscale("log")
    axs[0, 1].set_xscale("log")
    axs[0, 2].set_xscale("log")
    axs[0, 0].set_ylabel("Solver order\nBinding sites representation")
    axs[1, 0].set_ylabel("Controlled order\nBinding sites representation")
    axs[2, 0].set_ylabel("Distribution entropy")
    axs[2, 0].set_xlabel("Number of generated solutions")
    axs[2, 1].set_xlabel("Number of generated solutions")
    axs[2, 2].set_xlabel("Number of generated solutions")
    axs[2, 2].legend(loc="upper right")

    fig.savefig("size_bias_control.png", bbox_inches="tight")
    fig.savefig("size_bias_control.svg", bbox_inches="tight")


def _plot_cumdist(
    solver: str, library: int, ax_cumdist: Axes, ax_entropy: Axes
) -> None:
    if solver == "control":
        data = pd.read_csv(
            f"../benchmarks/topsols10_{library}_control_Gurobi.csv"
        ).to_numpy()
    else:
        data = pd.read_csv(
            f"../benchmarks/topsols10_{library}_solver_{solver}.csv"
        ).to_numpy()
    data = data[:, : data.shape[1] // 2] + data[:, data.shape[1] // 2 :]
    cumdist = np.cumsum(data, axis=0).astype(np.float64)
    cumdist /= cumdist.sum(axis=1, keepdims=True)
    ps = np.arange(cumdist.shape[0]) + 1
    ax_cumdist.plot(ps, cumdist)
    entropy = st.entropy(cumdist, axis=1)
    isolver = {"Gurobi": 0, "control": 1, "CBC": 2, "SCIP": 3}[solver]
    solver_label = "Diversity-driven" if solver == "control" else solver
    ax_entropy.plot(ps, entropy, label=f"{solver_label} order", color=f"C{isolver}")
    if solver == "control":
        max_entropy = [-np.log(1 / data.shape[1])] * len(ps)
        ax_entropy.plot(ps, max_entropy, "--", label="Maximum entropy", color="gray")


def solver_bias() -> None:
    """
    Plot file "solver_bias.svg" from the files "topsols10_{i}_solver_{solver}.csv".
    """
    fig, axs = plt.subplots(5, 3, sharex=True, sharey="row", figsize=(10, 13))

    for i in range(3):
        for isolver, solver in enumerate(["Gurobi", "SCIP", "CBC", "control"]):
            _plot_cumdist(solver, i, axs[isolver, i], axs[4, i])
            solver_label = "Diversity-driven" if solver == "control" else solver
            axs[isolver, 0].set_ylabel(
                f"{solver_label} order\nBinding sites representation"
            )

        axs[0, i].set_title(f"Collection {i + 1}")
        axs[0, i].set_xscale("log")
        axs[4, i].set_xlabel("Number of generated solutions")

    axs[4, 0].set_ylabel("Distribution entropy")
    axs[4, 2].legend()

    fig.savefig("solver_bias.png", bbox_inches="tight")
    fig.savefig("solver_bias.svg", bbox_inches="tight")


solver_bias()
# homogeneous()
# inhomogeneous()
