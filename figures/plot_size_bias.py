"""
Plot the file "size_bias.png" (Figure S2) from the files "size_bias_{i}_solver.csv".
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

replicates = []
ids = []
frequencies = []

for replicate in range(10):
    data = pd.read_csv(f"../benchmarks/size_bias_{replicate}_solver.csv")
    print(len(data), data.sum().to_numpy().sum() / len(data))
    s = data.sum().to_numpy()[:10].astype(np.float64)
    replicates += [replicate] * 10
    ids += list(range(5, 15))
    frequencies += list(s / s.sum())

data = pd.DataFrame(
    {
        "Binding site length (bp)": ids,
        "library": replicates,
        "Frequency in top-scoring solutions": frequencies,
    }
)

sns.lineplot(
    data=data, x="Binding site length (bp)", y="Frequency in top-scoring solutions"
)

plt.savefig("size_bias.png", bbox_inches="tight")
plt.savefig("size_bias.svg", bbox_inches="tight")
