# dense-arrays

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
