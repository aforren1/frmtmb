# frmtmb.survival

Survival families for [frmtmb](https://aforren1.github.io/frmtmb/).

- `cox()` is the proportional hazards model brms calls `cox()`. Its
  baseline hazard is an M-spline in time with simplex weights.
  `cox_baseline()` gives the fitted weights.
- `royston_parmar()` is the flexible parametric model of Royston and
  Parmar (2002). It writes the log cumulative hazard as a natural cubic
  spline in log time, parameterized as `flexsurv::flexsurvspline()`
  parameterizes it. `rp_floored()` reports where a fit's likelihood is
  a floor and not a density.

Censoring, truncation and random-effect frailties come from the frmtmb
formula grammar.

## Install

The package is in the frmtmb repository. Install frmtmb first, then
this package:

```r
remotes::install_github("aforren1/frmtmb")
remotes::install_github("aforren1/frmtmb",
                        subdir = "extensions/frmtmb.survival")
```

## Example

```r
library(frmtmb.survival)
set.seed(2)
dd <- data.frame(x = rnorm(200))
dd$time <- rexp(200, exp(-0.5 + 0.7 * dd$x))
dd$out <- rbinom(200, 1, 0.3)
fit <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
fixef_by_dpar(fit)$mu
cox_baseline(fit)
```

`library(frmtmb.survival)` attaches frmtmb too. Without this package,
`family = "cox"` and `brmsfamily("cox")` stop with a message that names
it.
