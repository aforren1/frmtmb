# Pointwise log-likelihood

`log_lik()` is the log-likelihood of every observation at every
posterior draw, the matrix
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md) and
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md) are
computed from. A `frmtmb_fit` is one maximum-likelihood parameter vector
rather than a posterior, so this refuses and names the route to draws,
as [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md) and
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md) do; the
estimator itself is in the `frmtmb.sample` package, on this same
generic.

## Usage

``` r
log_lik(object, ...)

# S3 method for class 'frmtmb_fit'
log_lik(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`, or (with `frmtmb.sample` loaded) draws.

- ...:

  Passed to methods.

## Value

This method signals an error on a maximum-likelihood fit.

## Details

The generic lives here rather than in `frmtmb.sample` so that a script
ported from brms is told to sample before it has sampled. Without it
`log_lik(fit)` stopped at "could not find function", which names neither
the reason nor the route.

## See also

[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md) and
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.md), which
average this matrix over draws;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) for the
maximum-likelihood total.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)

# the maximum-likelihood total is available directly
logLik(fit)
#> 'log Lik.' -52.08647 (df=3)
# the pointwise posterior matrix needs draws, and says so
try(log_lik(fit))
#> Error : log_lik() is the pointwise likelihood at every posterior draw, and this is a maximum-likelihood fit with one parameter vector. Sample first, with frmtmb.sample::log_lik(frmtmb.sample::frm_sample(fit)) once that package is installed; logLik() is the maximum-likelihood total already on the fit
```
