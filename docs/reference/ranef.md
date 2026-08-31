# Extract random-effect modes

Extract random-effect modes

## Usage

``` r
ranef(object, ...)

# S3 method for class 'frmtmb_fit'
ranef(object, condVar = FALSE, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

- condVar:

  If `TRUE`, attach the conditional SDs of the modes (from the Laplace
  posterior) as a `"condSD"` attribute on each matrix, in matching
  layout.

## Value

A named list of levels-by-coefficients matrices, one per random-effect
term. [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
gives the long form (with a `condsd` column when `condVar = TRUE` was
used).
