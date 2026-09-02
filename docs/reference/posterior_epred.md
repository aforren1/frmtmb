# Expected-value and predictive draws from sampled parameters

`posterior_epred()` evaluates the response-scale expectation per draw;
`posterior_predict()` additionally simulates responses from the family,
giving the posterior predictive distribution. Both condition on each
draw's own random effects (`re.form = NA` drops them).

## Usage

``` r
posterior_epred(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_epred(
  object,
  newdata = NULL,
  resp = NULL,
  re.form = NULL,
  ndraws = NULL,
  ...
)

posterior_linpred(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_linpred(
  object,
  transform = FALSE,
  newdata = NULL,
  resp = NULL,
  re.form = NULL,
  dpar = NULL,
  ndraws = NULL,
  ...
)

posterior_predict(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_predict(
  object,
  newdata = NULL,
  resp = NULL,
  re.form = NULL,
  ndraws = NULL,
  ...
)
```

## Arguments

- object:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Unused.

- newdata, resp, re.form:

  As in
  [`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md).

- ndraws:

  Number of draws to use (default: all).

- transform:

  For `posterior_linpred()`: if `TRUE`, apply the inverse link (the
  value of the `mu` dpar on its natural scale, brms's convention; unlike
  `posterior_epred()` this is not the response mean for zero-inflated
  and similar families).

- dpar:

  For `posterior_linpred()`: which distributional parameter's linear
  predictor to evaluate.

## Value

A draws-by-observations matrix; for a categorical outcome
`posterior_epred()` returns a draws-by-observations-by-categories array
(see the section below).

## Categorical outcomes

An ordinal family predicts a DISTRIBUTION per observation, not one
number: each draw's `predict(type = "response")` is an `n x K` matrix of
category probabilities. Those stack into a 3-D
`draws x observations x categories` array. `dimnames` are
`list(NULL, <observation names or NULL>, <category levels>)`, so
`ep[, , "high"]` is the draws-by-observations matrix for one category
and `ep[k, , ]` is draw `k`'s own `n x K` prediction, the matrix
`predict(type = "response")` returns. Every `ep[k, i, ]` sums to 1 for
an ordinal family.

This is brms's convention:
[`?brms::posterior_epred.brmsfit`](https://paulbuerkner.com/brms/reference/posterior_epred.brmsfit.html)
documents "an S x N x C array" for categorical and ordinal models and an
S x N matrix otherwise, and frmtmb follows brms spelling for brms-origin
functions. Any family whose per-draw response-scale prediction is a
matrix takes the array shape; every family that predicts one number per
observation keeps the plain `draws x observations` matrix.

`posterior_predict()` is unaffected - it draws one category per
observation - and so is `posterior_linpred()`, which is a statement
about one distributional parameter and stays an `n`-column matrix of the
latent predictor.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rpois(80, exp(0.3 + 0.4 * dd$x + rnorm(8, 0, 0.5)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

nd <- data.frame(x = c(-1, 0, 1),
                 g = factor(1, levels = levels(dd$g)))

# the expected response per draw: uncertainty in the mean
ep <- posterior_epred(ds, newdata = nd)
apply(ep, 2, quantile, c(0.025, 0.5, 0.975))

# the predictive distribution adds the family's own noise, so its
# intervals are wider
pp <- posterior_predict(ds, newdata = nd)
apply(pp, 2, quantile, c(0.025, 0.5, 0.975))

# the linear predictor itself, on the link scale by default
head(posterior_linpred(ds, newdata = nd, ndraws = 5))
}
#> Warning: There were 1 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.07, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>             [,1]        [,2]      [,3]
#> [1,]  0.42521058  0.69582071 0.9664308
#> [2,] -0.33351239 -0.05276348 0.2279854
#> [3,] -0.05381016  0.40258886 0.8589879
#> [4,] -0.13356468  0.19153765 0.5166400
#> [5,] -0.06035127  0.34516496 0.7506812
# }
```
