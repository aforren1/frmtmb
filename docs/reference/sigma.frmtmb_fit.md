# Residual standard deviation

Returns the estimated `sigma` distributional parameter on the response
scale when it is constant across observations (intercept-only or fixed).
When `sigma` is modeled with covariates the scalar summary does not
exist; the method warns and returns `NA` (use `predict(dpar = "sigma")`
for the per-observation values). Families without a `sigma` parameter
return 1, following glmmTMB.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
sigma(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

## Value

A scalar, or a named vector for multivariate fits.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# the residual SD, on the response scale
sigma(fit)
#> [1] 0.9658235
# which is what the standardized residuals divide by
max(abs(residuals(fit, type = "pearson") -
          residuals(fit) / sigma(fit)))
#> [1] 0

# a poisson fit has no dispersion parameter, so sigma() is 1
dd$cnt <- rpois(100, exp(0.5 + 0.3 * dd$x))
sigma(frm(bf(cnt ~ x) + poisson(), data = dd))
#> [1] 1
```
