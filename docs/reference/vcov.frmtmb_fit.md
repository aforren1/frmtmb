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

  If `TRUE`, include covariance parameters (`theta`).

- ...:

  Unused.

## Value

A covariance matrix.
