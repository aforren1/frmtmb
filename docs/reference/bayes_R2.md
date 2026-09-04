# Bayesian R-squared

The proportion of the outcome's variance a model explains, computed per
posterior draw. A maximum-likelihood fit has one parameter vector rather
than a posterior, so this refuses on a `frmtmb_fit` and names the route
to draws; the estimator itself is in the `frmtmb.sample` package, on the
same generic.

## Usage

``` r
bayes_R2(object, ...)

# S3 method for class 'frmtmb_fit'
bayes_R2(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`, or (with `frmtmb.sample` loaded) draws.

- ...:

  Passed to methods.

## Value

This method signals an error on a maximum-likelihood fit.

## See also

[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x) + gaussian(), data = dd)
try(bayes_R2(fit))
#> Error : bayes_R2() is computed per posterior draw and this is a maximum-likelihood fit. Install frmtmb.sample and sample first: bayes_R2(frm_sample(fit))
```
