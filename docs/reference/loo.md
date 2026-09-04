# Approximate leave-one-out cross-validation

`loo()` and `waic()` estimate the expected log predictive density of a
model from posterior draws: `loo()` by Pareto-smoothed importance
sampling, `waic()` by the widely applicable information criterion.
`loo_compare()` ranks several of them. `LOO()` and `WAIC()` are brms's
deprecated capitalized spellings.

## Usage

``` r
loo(x, ...)

# S3 method for class 'frmtmb_fit'
loo(x, ...)

waic(x, ...)

# S3 method for class 'frmtmb_fit'
waic(x, ...)

loo_compare(x, ...)

# Default S3 method
loo_compare(x, ...)

LOO(x, ...)

# S3 method for class 'frmtmb_fit'
LOO(x, ...)

WAIC(x, ...)

# S3 method for class 'frmtmb_fit'
WAIC(x, ...)
```

## Arguments

- x:

  A `frmtmb_fit`, or (with `frmtmb.sample` loaded) draws.

- ...:

  Passed to methods.

## Value

These methods signal an error on a maximum-likelihood fit.

## Details

All of these average a likelihood over posterior draws, so on a
`frmtmb_fit` - which is one maximum-likelihood parameter vector - they
refuse and name the two routes to an answer:
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html), which are the
maximum-likelihood analogues already available on the fit, or sampling
the model first.

## Sampling

The estimators for posterior draws are in the `frmtmb.sample` package,
which also provides `frm_sample()`:

    # install.packages("remotes")
    remotes::install_github("aforren1/frmtmb",
                            subdir = "extensions/frmtmb.sample")
    library(frmtmb.sample)
    loo(frm_sample(fit))

It registers `frmtmb_draws` methods on these generics, so the spellings
on this page are the ones that keep working once it is loaded.

## See also

[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) for the maximum-likelihood
comparison,
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
for a resampling one.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)

# the maximum-likelihood comparison is available directly
AIC(fit)
#> [1] 110.1729
# the predictive one needs draws, and says so
try(loo(fit))
#> Error : loo() is a posterior quantity and this is a maximum-likelihood fit: an elpd averages the likelihood over draws. Sample first, with frmtmb.sample::loo(frmtmb.sample::frm_sample(fit)) once that package is installed, or compare maximum-likelihood fits with AIC() or BIC()
```
