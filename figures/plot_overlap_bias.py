"""
Plot file "overlap_bias.svg" (Figure S4).
"""
import arviz as az
import bambi as bmb
import dense_arrays as da
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

# from scipy.special import logit

fig, ax = plt.subplots()

pss = []
avgdists = []

for ilib in range(3):
    # The solver or control order doesn't matter, we just want the full solution set
    sols = pd.read_csv(f"../benchmarks/topsols10_{ilib}_solver_Gurobi.csv")
    data = sols.to_numpy()
    data = data[:, : data.shape[1] // 2] + data[:, data.shape[1] // 2 :]
    proms = sols.columns.to_numpy()[: data.shape[1]]
    ps = data.sum(axis=0) / data.sum()
    matrix = np.array(
        [[da.optimize.shift_metric(pa, pb) for pb in proms] for pa in proms]
    )
    avgdist = (matrix.mean(axis=0) + matrix.mean(axis=1)) / 2
    ax.scatter(avgdist, ps, label=f"Collection {ilib + 1}", marker="oxP"[ilib])
    pss.append(ps)
    avgdists.append(avgdist)

data = pd.DataFrame(
    {
        "abundance": np.concatenate(pss),
        "distance": np.concatenate(avgdists),
        "collection": np.concatenate([[i] * len(ps) for i, ps in enumerate(pss)]),
    }
)
# model = bmb.Model("abundance ~ (distance|collection)", data)
model = bmb.Model("abundance ~ distance", data)
print(model)
idata = model.fit()
print(az.summary(idata))
posterior = az.extract(idata).thin({"sample": 20})

distances = xr.DataArray([9.3, 10])
intercept_common = posterior["Intercept"]
slope_common = posterior["distance"]

bmb.interpret.plot_predictions(model, idata, "distance", ax=ax)

ax.set_ylabel("Relative abundance in the top-scoring solutions")
ax.set_xlabel("Average $d_{ij}$ to and from other binding sites")
ax.legend(loc="lower left")

fig.savefig("overlap_bias.png", bbox_inches="tight")
fig.savefig("overlap_bias.svg", bbox_inches="tight")
