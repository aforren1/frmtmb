# Pointwise log-likelihood of posterior draws

The `ndraws x nobs` matrix of per-observation log-densities, each row
evaluated at one draw's own parameter vector. Every leave-one-out and
information-criterion quantity in this package is a function of this
matrix, and
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md),
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md) and
[`loo::psis()`](https://mc-stan.org/loo/reference/psis.html) take it
directly.

## Usage

``` r
log_lik(object, ...)

# S3 method for class 'frmtmb_draws'
log_lik(object, ndraws = NULL, resp = NULL, ...)
```

## Arguments

- object:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Unused.

- ndraws:

  Number of draws to use, evenly spaced through the matrix (default: all
  of them).

- resp:

  For a multivariate model without `rescor`, the response whose
  contribution to report; the default sums over responses.

## Value

A numeric matrix with one row per draw and one column per observation
(the rows the model was fitted on).

## What is conditioned on

The density is CONDITIONAL on the draw's own group-level values:
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
samples `b` alongside everything else, so each row of the draws matrix
is a complete parameter vector and no integration is left to do. This is
exactly brms's convention, where the Stan model also samples the
group-level parameters. The consequence is worth stating plainly:
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md) on such a
matrix is leave-one-OBSERVATION-out with the groups held fixed, not
leave-one-group-out, and for a model with few observations per group the
two differ.

Addition terms enter exactly as they enter the fitted objective, because
they enter through the same code: `cens()` replaces a row's density by
the matching CDF difference,
[`trunc()`](https://rdrr.io/r/base/Round.html) divides by the window
mass, [`weights()`](https://rdrr.io/r/stats/weights.html) multiplies the
row's contribution, and `trials()` is the family's own argument. brms
composes them in the same order (`log_lik_censor()`,
`log_lik_truncate()`, `log_lik_weight()`).

## Multivariate models

With `set_rescor(TRUE)` a column is the joint density of the row's
response VECTOR, so the matrix keeps one column per observation. Without
`rescor`, the responses are independent given the predictors and the
default sums their log-densities per row; pass `resp` to get one
response's contribution alone.

## Likelihoods with no per-observation column

A model whose smallest independent unit is a group has no
per-observation column to leave out, and this refuses rather than
inventing one: R-side residual correlation
([frmtmb-autocor](https://aforren1.github.io/frmtmb/reference/frmtmb-autocor.md)),
a [`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md)
sequence, and a group-level mixture (`mixture(groups = )`). An
[`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md) subject is
one row, so its column is well defined and is not refused. In-model
imputation (`mi()`, `me()`) is refused for the same kind of reason: a
latent value is a parameter, not an observation. Use
[`AIC()`](https://rdrr.io/r/stats/AIC.html) on the maximum-likelihood
fits or
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
for those.

## See also

[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md),
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md),
[`bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.md)

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

  ll <- log_lik(ds)
  dim(ll)
  # the column means are the per-observation expected log-densities
  head(colMeans(ll))
}
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] -0.8789809 -1.9343176 -1.3297286 -2.7670318 -0.8878848 -3.2526304
# }
```
