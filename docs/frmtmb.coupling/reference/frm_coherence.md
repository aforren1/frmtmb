# Coherence and phase from a fitted coupling model, with intervals

Reads the `coh` and `phase` distributional parameters off a
[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
fit and returns them with confidence intervals.

## Usage

``` r
frm_coherence(
  fit,
  newdata = NULL,
  re.form = NULL,
  level = 0.95,
  allow_new_levels = FALSE
)

frm_phase(
  fit,
  newdata = NULL,
  re.form = NULL,
  level = 0.95,
  allow_new_levels = FALSE
)
```

## Arguments

- fit:

  A
  [`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
  fit whose family is
  [`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md).

- newdata:

  Optional data frame of predictor values. Defaults to the data the
  model was fitted to.

- re.form:

  `NULL` keeps the random effects, so the answer is per group; `NA`
  drops them, so the answer is the population one. Passed through to
  [`stats::predict()`](https://rdrr.io/r/stats/predict.html).

- level:

  Confidence level.

- allow_new_levels:

  Passed through to
  [`stats::predict()`](https://rdrr.io/r/stats/predict.html).

## Value

A data frame with one row per row of `newdata` and columns `.estimate`,
`.se`, `.lower` and `.upper`. `frm_coherence()` adds `.eta`, the value
on the logit scale the interval was built on.

## Why the interval is built on the link scale

Coherence lives on `(0, 1)` and a Wald interval on that scale walks off
the end whenever the estimate is near a boundary, which is exactly where
few segments put it. The interval here is a Wald interval on the logit
scale pushed through
[`plogis()`](https://rdrr.io/r/stats/Logistic.html), so it cannot leave
`(0, 1)` and it is asymmetric in the direction the sampling distribution
actually is.

Measured coverage of that interval, 40 subjects with a coherence random
effect, 150 replicates: 0.912 at 4 segments per subject and 0.893 at 16,
when every distributional parameter carries a random effect. With the
two channel powers left as free per-subject parameters instead, coverage
falls to 0.622 and 0.820, because two free numbers per subject that gain
no information as subjects are added starve the variance component.
`dev/xspec-findings.md` in the repository has the table. **Put
`(1 | id)` on every dpar, not only on `coh`.**

## Phase is an angle

`frm_phase()` returns radians and does no wrapping: the linear predictor
is unbounded and a phase of `-3.1` and one of `3.18` are the same angle
reached from two sides. When a group's phase is near the wrap point its
interval is wide for a reason that is a parameterization artifact rather
than data, and the fit is worth re-centering with an offset. Nothing
here detects that for you.

## Examples

``` r
set.seed(5)
src <- rnorm(4096)
xs <- frm_cross_spectrum(src + rnorm(4096), 0.9 * src + rnorm(4096),
                         segments = 16)
xs <- xs[xs$freq < 0.1, ]
fit <- frmtmb::frm(
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1, coh ~ 1, phase ~ 1),
  family = cross_wishart(), data = xs)
frm_coherence(fit, newdata = xs[1, ])
#>   .estimate       .se    .lower    .upper      .eta
#> 1  0.267868 0.1366232 0.2187026 0.3235095 -1.005466
frm_phase(fit, newdata = xs[1, ])
#>      .estimate        .se     .lower    .upper
#> 1 -0.003440959 0.05845064 -0.1180021 0.1111202
```
