# CAUSE 1.2.0 source patches

CAUSE 1.2.0 does not run unmodified under R 4.3.3 with loo 2.10.1. Two patches to
`R/in_sample_elpd_loo.R` were required. Both are faithful to the published method;
neither changes the estimator. The patched source is not redistributed here.

## Patch 1: list-slice null test

```r
# before (a list slice is never NULL, so a missing model falls through)
} else if (is.null(loos[i2])) {
# after
} else if (is.null(loos[[i2]])) {
```
Applied at two sites (lines 23 and 108). The parallel `loos[[i1]]` test above each
was already correct, so this restores the intended symmetry.

## Patch 2: elpd difference computed directly

`loo_compare()` output indexing fails under loo 2.10.1 with the error
`non-numeric argument to binary operator` at `-1 * comp[2,1]`. The call was replaced
with the pointwise computation that `loo_compare` itself performs, and that the same
function already uses in its other two branches:

```r
# before
comp <- loo_compare(loos[[i1]], loos[[i2]])
if (rownames(comp)[1] == "model2") { est <- comp[2,1] } else { est <- -1*comp[2,1] }
se <- comp[2,2]
# after
diff <- loos[[i1]]$pointwise[,1] - loos[[i2]]$pointwise[,1]
est  <- sum(diff)
se   <- sd(diff) * sqrt(length(diff))
```
Applied at two sites (lines 28 and 113).

## Reproducing

```bash
curl -sSL -o cause.tar.gz https://github.com/jean997/cause/archive/refs/tags/v1.2.0.tar.gz
tar -xzf cause.tar.gz && cd cause-1.2.0
sed -i 's/is\.null(loos\[i1\])/is.null(loos[[i1]])/g; s/is\.null(loos\[i2\])/is.null(loos[[i2]])/g' R/in_sample_elpd_loo.R
# then apply patch 2 as shown above
R CMD INSTALL .
```

Environment: R 4.3.3, loo 2.10.1, CAUSE 1.2.0, 3 GB RAM. The nuisance-parameter
subset was reduced to 100,000 variants (CAUSE's default suggestion is 1,000,000)
to fit available memory; this affects the precision of the nuisance parameters,
not the model comparison.
