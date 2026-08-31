# Factor with numeric-coded levels for coordinate covariance structures

`ou()` and the spatial structures
([`exp()`](https://rdrr.io/r/base/Log.html), `gau()`, `mat()`) need the
positions of the term levels. `num_factor(x)` (one dimension) or
`num_factor(x, y)` (planar coordinates) encodes them in the level labels
the same way
[`glmmTMB::numFactor()`](https://rdrr.io/pkg/glmmTMB/man/numFactor.html)
does, so factors created by either function work.

## Usage

``` r
num_factor(x, y = NULL)
```

## Arguments

- x:

  Numeric positions (times, coordinates).

- y:

  Optional second coordinate.

## Value

A factor whose levels encode the unique positions.
