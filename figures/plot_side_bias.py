"""
Plot "side_bias.svg" from the file "side_bias_both.csv".
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

activators = np.zeros((2, 100), dtype=np.float64)
repressors = np.zeros((2, 100), dtype=np.float64)
upstream = np.zeros((2, 100), dtype=np.float64)
downstream = np.zeros((2, 100), dtype=np.float64)

n = 0
nrng = -1
rng = np.random.default_rng(seed=42)
library_size = 20
motif_min_size = 5
motif_max_size = 15
sequence_length = 100
ATGC = list("ATGC")

with Path("../benchmarks/side_bias_both.csv").open("r") as f:
    for line in f.readlines():
        s = line.rstrip().split(",")
        bias = int(s[1] == "bias")
        n += bias
        while nrng < int(s[0]):
            # Advancing the prng state to get the same as from the simulation
            motifs = [
                "".join(
                    rng.choice(
                        ATGC,
                        rng.integers(motif_min_size, motif_max_size, endpoint=True),
                    )
                )
                for _ in range(library_size)
            ]
            _ = "".join(rng.choice(ATGC, 6))
            _ = "".join(rng.choice(ATGC, 6))
            nrng += 1
        for i, x in enumerate(s[2:]):
            if x == "None":
                continue
            if i % 22 == 20:  # noqa: PLR2004
                upstream[bias, int(x) : int(x) + 6] += 1
            elif i % 22 == 21:  # noqa: PLR2004
                downstream[bias, int(x) : int(x) + 6] += 1
            elif i % 2 == 0:
                activators[bias, int(x) : int(x) + len(motifs[i % 22])] += 1
            else:
                repressors[bias, int(x) : int(x) + len(motifs[i % 22])] += 1

fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)

axs[0].set_title("No side bias")
axs[0].plot(range(100), activators[0], label="upstream pref.")
axs[0].plot(range(100), repressors[0], label="downstream pref.")
axs[0].plot(range(100), upstream[0], label="-35")
axs[0].plot(range(100), downstream[0], label="-10")
axs[1].set_title("With side bias")
axs[1].plot(range(100), activators[1], label="upstream pref.")
axs[1].plot(range(100), repressors[1], label="downstream pref.")
axs[1].plot(range(100), upstream[1], label="-35")
axs[1].plot(range(100), downstream[1], label="-10")
axs[1].set_xlabel("Sequence position (bp)")
axs[0].legend(loc="upper left")
for ax in axs:
    ax.set_ylabel("Number of motifs\ncovering this position")
    ax.fill_between(
        [89 - 18 - 6, 91 - 16 - 1],
        0,
        1,
        color="C2",
        alpha=0.3,
        transform=ax.get_xaxis_transform(),
    )
    ax.fill_between(
        [89, 91 + 6 - 1],
        0,
        1,
        color="C3",
        alpha=0.3,
        transform=ax.get_xaxis_transform(),
    )
    ax.grid()
plt.savefig("side_bias.svg", bbox_inches="tight")
plt.savefig("side_bias.png", bbox_inches="tight")
