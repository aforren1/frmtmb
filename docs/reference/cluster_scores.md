# Per-cluster score matrix

The gradient of each cluster's contribution to the marginal
log-likelihood, evaluated at the fitted estimates. This is the quantity
[`sandwich::estfun()`](https://zeileis.codeberg.page/sandwich/reference/estfun.html)
would return if the marginalized objective had per-observation
contributions, which it does not - see
[`vcov_cluster()`](https://aforren1.github.io/frmtmb/reference/vcov_cluster.md)
for the argument and the guards.

## Usage

``` r
cluster_scores(object, cluster)
```

## Arguments

- object:

  A `frmtmb_fit`.

- cluster:

  The clustering factor: a one-sided formula such as `~ g` (looked up in
  the model frame, then in the data the model was fitted to), a variable
  name, or a vector with one entry per fitted row.

## Value

A matrix with one row per cluster level and one column per estimated
outer parameter, named as in `vcov(object, full = TRUE)`. Signs follow
[`sandwich::estfun()`](https://zeileis.codeberg.page/sandwich/reference/estfun.html):
these are scores of the log-likelihood, not of the objective.

## Details

Rows sum to the gradient of the full objective, which is (near) zero at
a converged optimum; that identity is the cheapest check that the
clustering factor really does split the likelihood.

## See also

[`vcov_cluster()`](https://aforren1.github.io/frmtmb/reference/vcov_cluster.md),
which sandwiches these between two copies of the inverse observed
information.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:25, 8)))
dd$y <- rnorm(200, 1 + 0.5 * dd$x + rnorm(25, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = FALSE)

S <- cluster_scores(fit, ~ g)
dim(S)
#> [1] 25  4
# the scores add up to the gradient at the optimum
max(abs(colSums(S)))
#> [1] 8.139203e-05
```
