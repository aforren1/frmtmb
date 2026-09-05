# How many rows the grid could not represent

The generalized likelihood is floored: where the solved density at a
trial's own response time underflows, the log density is a large finite
negative number rather than `NaN`, so the optimizer gets a value it can
use and a mixture's log-sum-exp is not poisoned by one component. That
is the right behavior and it is silent, which is why this exists: it
says, once and after the fit, how many rows were answered by the floor
rather than by the solver.

## Usage

``` r
gddm_floored(fit)
```

## Arguments

- fit:

  A fitted
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)
  model.

## Value

The number of rows at the floor, with the row indices in the `"rows"`
attribute and the number of observations in `"n_obs"`.

## Details

A count of zero is the ordinary case and means the grid represented
every observation. A small count means a few trials sit in the leading
edge, within a few time steps of the fitted non-decision time, where a
fixed grid cannot resolve a first-passage density that is climbing
through orders of magnitude; those rows contributed a constant instead
of information. A large count means the fit is not to be trusted: shrink
`dt` in
[`gddm_control()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_control.md),
or give the model a lapse component with `gddm(lapse = "uniform")`,
which floors the density in the model rather than in the arithmetic.

## See also

[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md),
[`gddm_control()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_control.md)

## Examples

``` r
# \donttest{
set.seed(3)
ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 101)
dat <- gddm_simulate(200, mu = 2, bs = 2.5, ndt = 0.25, control = ctl)
dat$cond <- 1L
fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
           family = gddm(control = ctl), data = dat)
gddm_floored(fit)
#> [1] 0
#> attr(,"rows")
#> integer(0)
#> attr(,"n_obs")
#> [1] 200
# }
```
