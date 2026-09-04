# Simulate a drift-diffusion data set

Draws response times and boundary choices jointly from the Wiener
first-passage process, which is what an experiment produces: the
boundary is an outcome, not a design variable. The family itself treats
the boundary as data and simulates the response time conditional on it,
so this is the function to use when you want a whole data set rather
than a posterior predictive draw.

## Usage

``` r
ddm_simulate(n, mu, bs, ndt, bias = 0.5)
```

## Arguments

- n:

  Number of trials.

- mu, bs, ndt, bias:

  Drift rate, boundary separation, non-decision time and relative start
  point, on the response scale. Each may be a single value or a vector
  of length `n`, which is how a two-condition design is built.

## Value

A data frame with `rt`, the response time, and `upper`, 1 for a response
at the upper boundary and 0 for the lower one, in the coding `vint()`
expects.

## Examples

``` r
set.seed(1)
cond <- rep(c(0, 1), each = 100)
dat <- ddm_simulate(200, mu = 0.2 + 1.1 * cond, bs = 1.4,
                    ndt = 0.3, bias = 0.5)
dat$cond <- factor(cond)
str(dat)
#> 'data.frame':    200 obs. of  3 variables:
#>  $ rt   : num  0.52 1.397 0.448 0.575 0.681 ...
#>  $ upper: int  1 0 1 1 0 0 0 1 0 0 ...
#>  $ cond : Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
```
