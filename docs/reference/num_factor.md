# Factor with numeric-coded levels for coordinate covariance structures

`ou()` (and future spatial structures) need the positions of the term
levels. `num_factor(x)` encodes them in the level labels the same way
[`glmmTMB::numFactor()`](https://rdrr.io/pkg/glmmTMB/man/numFactor.html)
does, so factors created by either function work.

## Usage

``` r
num_factor(x)
```

## Arguments

- x:

  Numeric positions (times, coordinates).

## Value

A factor whose levels encode the sorted unique positions.
