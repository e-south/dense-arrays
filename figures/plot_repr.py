"""
Plot the file "representation.svg" (Figure 3) from files "topsols10_{i}_solver.csv".
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

rng = np.random.default_rng(1)

fig = plt.figure(figsize=(15, 10), layout="constrained")
axs = fig.subplot_mosaic(
    [[f"dist_{i}", f"rand_{i}", "ordered"] for i in range(3)],
    width_ratios=[1, 1, 3],
    sharey=True,
    sharex=True,
)

all_data = []

for i in range(3):
    data = pd.read_csv(f"../benchmarks/topsols10_{i}_solver.csv").to_numpy()
    data = data[:, : data.shape[1] // 2] + data[:, data.shape[1] // 2 :]
    shuffled_data = np.array(data)
    for j in range(len(data)):
        rng.shuffle(shuffled_data[j, :])
    dist = data.sum(axis=0).astype(np.float64)
    shuffled_dist = shuffled_data.sum(axis=0).astype(np.float64)
    axs[f"dist_{i}"].set_title(
        f"Library {i + 1}: {round(dist.sum() / 6)} seqs with {6} motifs"
    )
    axs[f"dist_{i}"].set_ylabel("Frequency of motifs in sequences")
    axs[f"dist_{i}"].set_xticks(range(1, 11))
    axs[f"rand_{i}"].set_xticks(range(1, 11))
    axs[f"rand_{i}"].set_title(f"{round(dist.sum() / 6)} random samples of {6} motifs")
    dist /= dist.sum()
    shuffled_dist /= shuffled_dist.sum()
    axs[f"dist_{i}"].bar(np.arange(len(dist)) + 1, dist)
    axs[f"rand_{i}"].bar(np.arange(len(shuffled_dist)) + 1, shuffled_dist, color="C1")
    dist[::-1].sort()
    shuffled_dist[::-1].sort()
    for j, (freq, freq_shuffled) in enumerate(zip(dist, shuffled_dist, strict=False)):
        all_data.extend(
            [
                {
                    "library": i + 1,
                    "ordered index": j + 1,
                    "frequency": freq,
                    "distribution": "dense array",
                },
                {
                    "library": i + 1,
                    "ordered index": j + 1,
                    "frequency": freq_shuffled,
                    "distribution": "uniform",
                },
            ]
        )
    # axs[f"rand_{i}"].sharey(axs[f"dist_{i}"])

axs["dist_2"].set_xlabel("Motif ID")
axs["rand_2"].set_xlabel("Motif ID")

df_data = pd.DataFrame(all_data)

sns.lineplot(
    df_data,
    x="ordered index",
    y="frequency",
    hue="distribution",
    ax=axs["ordered"],
    style="distribution",
    markers=True,
    dashes=False,
    units="library",
    estimator=None,
)

axs["ordered"].set_xlabel("Motifs sorted from most to least represented")
axs["ordered"].set_ylabel("Frequency of motifs in sequences")
axs["ordered"].set_xticks(range(1, 11))

fig.savefig("representation.png", bbox_inches="tight")
fig.savefig("representation.svg", bbox_inches="tight")
