

This folder contains 3 sub-experiments:

1. exhaustive (all-pairs)
This benchmarks all-genomes sketching and all-pairs exhaustive comparisons.
2. topk
This is top-k calculated, using an LSH index to avoid quadratic comparisons.
This contains both speed and accuracy results.
3. thresholded
This is filtered for J > threshold
We use J = {.9, .8, .7, .6, .5, .4, .3, .2, .1}
