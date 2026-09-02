# Posterior class membership of an `lca()` fit

One row per subject, one column per latent class, rows summing to one:
the probability that a subject belongs to each class given its observed
item responses and its gating covariates. This is poLCA's `$posterior`.

## Usage

``` r
lca_probs(fit)
```

## Arguments

- fit:

  A `frmtmb_fit` with an
  [`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md) family.

## Value

A `n x K` matrix of posterior class probabilities with an `entropy`
attribute.

## Details

The relative entropy of the classification is attached as the `entropy`
attribute: `1 - sum(-p log p) / (n log K)`, which is 1 for a partition
with no ambiguity and 0 when every subject is equally likely to be in
any class. Values above about 0.8 are the usual rule of thumb for
classes worth naming.

This is
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
under an LCA-specific name and check; the two return the same matrix for
an [`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md) fit.

## See also

[`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md),
[`lca_profiles()`](https://aforren1.github.io/frmtmb/reference/lca_profiles.md),
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)

## Examples

``` r
set.seed(3)
n <- 200
cl <- rbinom(n, 1, 0.5) + 1
pr <- rbind(c(0.9, 0.85, 0.8, 0.9), c(0.1, 0.15, 0.2, 0.1))
Y <- matrix(0L, n, 4)
for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(row = seq_len(n))
dd$Y <- Y
fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)

p <- lca_probs(fit)
head(p)
#>           class1       class2
#> [1,] 0.050260948 0.9497390516
#> [2,] 0.999644971 0.0003550287
#> [3,] 0.001387876 0.9986121235
#> [4,] 0.041479368 0.9585206322
#> [5,] 0.991678548 0.0083214516
#> [6,] 0.986656917 0.0133430834
attr(p, "entropy")            # classification quality, 0 to 1
#> [1] 0.8561927
table(max.col(p), truth = cl) # the modal assignment, up to relabeling
#>    truth
#>      1  2
#>   1 11 95
#>   2 93  1
```
