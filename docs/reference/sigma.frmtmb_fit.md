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
