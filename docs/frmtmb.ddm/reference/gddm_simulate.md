# Simulate from a generalized drift-diffusion model

Draws choices and response times by solving the model's own
Fokker-Planck equation and sampling the resulting defective densities,
so the draws come from the density
[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)
fits rather than from a discretized forward simulation, whose first
passages are late by however much the step misses excursions between
monitoring times.

## Usage

``` r
gddm_simulate(
  n,
  ...,
  coh = 0,
  drift = gddm_drift_constant(),
  bound = gddm_bound_constant(),
  start = gddm_start_point(),
  lapse = c("none", "uniform"),
  control = gddm_control()
)
```

## Arguments

- n:

  Number of trials.

- ...:

  Parameter values by name: `mu`, `bs`, `ndt` and whatever else the
  chosen components need (`alpha`, `leak`, `tau`, `kappa`, `bias`, `sz`,
  `lapse`). Anything not given takes the component's own starting value.

- coh:

  Coherence covariate, recycled to length `n`. Only meaningful with
  [`gddm_drift_coherence()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md).

- drift, bound, start, lapse, control:

  As in
  [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md).

## Value

A data frame with `rt`, `upper`, `cond` and, when a coherence drift is
used, `coh`.

## See also

[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

## Examples

``` r
set.seed(2)
head(gddm_simulate(20, mu = 1.5, bs = 2, ndt = 0.2,
                   control = gddm_control(t_max = 2, dt = 0.02,
                                          ny = 101)))
#>          rt upper cond
#> 1 0.8319035     1    1
#> 2 0.5872643     1    1
#> 3 1.1134061     1    1
#> 4 0.4327229     1    1
#> 5 0.5594095     1    1
#> 6 0.6637092     1    1
```
