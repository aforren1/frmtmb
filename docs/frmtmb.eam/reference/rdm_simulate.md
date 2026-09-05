# Simulate from a racing diffusion model

Draws from the generative process directly: a uniform start point per
accumulator, then that accumulator's inverse-Gaussian first-passage
time, and whichever arrives first. This is the joint draw of choice and
time, which is what the model produces and what a fit needs;
[`simulate()`](https://rdrr.io/r/stats/simulate.html) on a fitted object
instead redraws only the time, holding each row's observed choice.

## Usage

``` r
rdm_simulate(n, v, A = 0.5, k = 0.5, ndt = 0.2)
```

## Arguments

- n:

  Number of trials.

- v:

  Drift rates, one per accumulator, all positive. Its length sets the
  number of accumulators. May be a matrix with one row per trial, for a
  drift that varies with a covariate.

- A:

  Upper end of the start-point range.

- k:

  Distance from the top of the start-point range to the threshold, so
  the threshold is `A + k`.

- ndt:

  Non-decision time.

## Value

A data frame with `choice` (the accumulator that reached the threshold,
from 1) and `rt`.

## Examples

``` r
set.seed(1)
dat <- rdm_simulate(500, v = c(3.0, 2.0, 1.2), A = 0.5, k = 0.5,
                    ndt = 0.2)
table(dat$choice)
#> 
#>   1   2   3 
#> 247 164  89 
tapply(dat$rt, dat$choice, mean)
#>         1         2         3 
#> 0.3541706 0.3660322 0.3619880 
```
