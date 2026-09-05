# Simulate from a linear ballistic accumulator

Draws from the generative process directly: a uniform start point and a
normal drift rate per accumulator, and whichever reaches the threshold
first. This is the joint draw of choice and time, which is what the
model produces and what a fit needs;
[`simulate()`](https://rdrr.io/r/stats/simulate.html) on a fitted object
instead redraws only the time, holding each row's observed choice.

## Usage

``` r
lba_simulate(n, v, A = 0.5, k = 0.4, ndt = 0.2, sd_v = 1, posdrift = TRUE)
```

## Arguments

- n:

  Number of trials.

- v:

  Drift means, one per accumulator. Its length sets the number of
  accumulators. May be a matrix with one row per trial, for a drift that
  varies with a covariate.

- A:

  Upper end of the start-point range.

- k:

  Distance from the top of the start-point range to the threshold, so
  the threshold is `A + k`.

- ndt:

  Non-decision time.

- sd_v:

  Drift standard deviation, fixed. One number or one per accumulator.

- posdrift:

  Truncate the drift distribution at zero? Matches the
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/lba.md)
  argument of the same name.

## Value

A data frame with `choice` (the accumulator that responded, from 1) and
`rt`.

## Examples

``` r
set.seed(1)
dat <- lba_simulate(500, v = c(2.5, 1.5, 1.0), A = 0.5, k = 0.4,
                    ndt = 0.2)
table(dat$choice)
#> 
#>   1   2   3 
#> 321 126  53 
tapply(dat$rt, dat$choice, mean)
#>         1         2         3 
#> 0.4374168 0.4691844 0.4590923 
```
