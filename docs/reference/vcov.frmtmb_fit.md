# Covariance matrix of the fixed-effect estimates

Covers the estimated coefficients of every linear predictor; dpars fixed
to constants are excluded.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
vcov(object, full = FALSE, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- full:

  If `TRUE`, include covariance parameters (`theta`), named as in
  [`confint()`](https://rdrr.io/r/stats/confint.html) (the glmmTMB
  `vcov(full = TRUE)` convention).

- ...:

  Unused.

## Value

A covariance matrix.
