# Report the deep censored rows and the non-monotone rows of a fit

Two things in
[`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
used to be floors rather than answers. One of them is gone; the other
refuses.

## Usage

``` r
rp_floored(object, action = c("error", "report"), max_nlogS = 19.2)
```

## Arguments

- object:

  A `frmtmb_fit` with a
  [`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
  family.

- action:

  `"error"`, the default, refuses when the MONOTONICITY floor was used.
  The censored count never refuses under either value; since frmtmb
  0.52.0 it is a diagnostic. `"report"` returns the same numbers without
  refusing.

- max_nlogS:

  The `-log S` on a censored row above which the row is reported as
  barely constrained. The default 19.2 is where the OLD
  probability-scale arithmetic passed 1e-8 of error; it is kept as the
  threshold so that the two versions report the same rows.

## Value

A list with `n_censored_floored`, `max_nlogS`, `threshold`,
`n_nonmonotone`, `scale` and `n_obs`, returned invisibly when nothing
was floored. The offending row indices are the `"rows"` attribute, a
list with elements `censored` and `nonmonotone`.

## The censored rows

a report, no longer a floor: Up to frmtmb 0.51.0 core formed a
right-censored contribution as `log(1 - F(y))` on the probability scale,
so the scored `log S` carried absolute error about
`.Machine$double.eps / S`: past `-log S` of about 19.2 that error passed
1e-8, and past 30 the term was FLAT with a gradient of exactly zero.

frmtmb 0.52.0 added the `lccdf` slot and this family supplies it, in
closed form on all three scales, so a right-censored row is scored from
`log S` directly and no complement is formed. Measured on the hazard
scale, `-log S = 40` was scored as -35.127363 and is now scored as -40
exactly.

The count is therefore a DIAGNOSTIC and never refuses. It still says
something: a censored row whose fitted survival probability is
`exp(-40)` is one the data barely constrain, whatever the arithmetic
does. The quantity is `-log S` at the fitted parameters, and it is one
quantity for all three scales: `exp(eta)` on `"hazard"`,
`log1p(exp(eta))` on `"odds"` and `-log(Phi(-eta))` on `"normal"`. On
the hazard scale it is the cumulative hazard `H`.

## The monotonicity floor

The cumulative hazard has to increase, so the spline's derivative in log
time has to stay positive; nothing enforces it and flexsurv does not
enforce it either. Where it goes non-positive there is no hazard and the
true log density is `-Inf`, and this family replaces it with a large
finite number so that the optimizer has something to work with. That
keeps the fit alive and makes
[`logLik()`](https://rdrr.io/r/stats/logLik.html) a pseudo-likelihood: a
60 percent cure-fraction dataset has been measured converging with 6
such rows and a reported log likelihood 3952 units away from the
density's.

## What this cannot do

[`logLik()`](https://rdrr.io/r/stats/logLik.html) reads
`object$opt$objective` directly, so a check that runs after the fit
cannot make [`logLik()`](https://rdrr.io/r/stats/logLik.html) or
[`AIC()`](https://rdrr.io/r/stats/AIC.html) refuse on their own. What
frmtmb 0.52.0 does provide is a fit-end hook: this family declares
`post$fit_check`, so a fit with a non-monotone row warns as it is
returned rather than only when someone calls this function.
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
and its two companions still call it for you.

## See also

[`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)

## Examples

``` r
set.seed(1)
n <- 300
dd <- data.frame(trt = rep(0:1, each = n / 2))
dd$t <- rweibull(n, shape = 1.4, scale = exp(1 - 0.5 * dd$trt))
dd$censored <- as.integer(dd$t > 3)
dd$t <- pmin(dd$t, 3)
fit <- frmtmb::frm(frmtmb::bf(t | cens(censored) ~ trt),
                   family = royston_parmar(df = 2), data = dd)
str(rp_floored(fit, action = "report"))
#> List of 6
#>  $ n_censored_floored: int 0
#>  $ max_nlogS         : num 2.56
#>  $ threshold         : num 19.2
#>  $ n_nonmonotone     : int 0
#>  $ scale             : chr "hazard"
#>  $ n_obs             : int 300
#>  - attr(*, "rows")=List of 2
#>   ..$ censored   : int(0) 
#>   ..$ nonmonotone: int(0) 
```
