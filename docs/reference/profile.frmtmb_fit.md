# Likelihood profiles

Wraps [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html)
per parameter. The returned objects have
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) methods (from TMB).

## Usage

``` r
# S3 method for class 'frmtmb_fit'
profile(fitted, parm, ...)
```

## Arguments

- fitted:

  A `frmtmb_fit`.

- parm:

  Parameter names (as in
  [`confint()`](https://rdrr.io/r/stats/confint.html) rownames) or
  indices. Required; profiling is not free, so there is no
  all-parameters default.

- ...:

  Passed to
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html).

## Value

A `tmbprofile` data frame, or a named list of them when `parm` has
length above one.
