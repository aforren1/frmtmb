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

  `"error"`, the default, refuses when the MONOTONICITY floor was used
  on an event row, and warns when a censored row sits where the fitted
  survival function rises. The deep censored count never refuses under
  either value; since frmtmb 0.52.0 it is a diagnostic. `"report"`
  returns the same numbers silently.

- max_nlogS:

  The `-log S` on a censored row above which the row is reported as
  barely constrained. The default 19.2 is where the OLD
  probability-scale arithmetic passed 1e-8 of error; it is kept as the
  threshold so that the two versions report the same rows.

## Value

A list with `n_censored_deep`, `max_nlogS`, `threshold`,
`n_nonmonotone`, `n_nonmonotone_censored`, `max_survival_rise`, `scale`
and `n_obs`, returned invisibly when nothing refuses. The offending row
indices are the `"rows"` attribute, a list with elements `censored`,
`nonmonotone` and `nonmonotone_censored`.

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

## Censored rows where the survival function rises

The floor above is only ever used on an EVENT row, because only an event
row has a density. A censored row contributes `log S`, which the family
scores exactly whatever the sign of the derivative. So a censored row
with a non-positive `d(eta)/d(log t)` does not make
[`logLik()`](https://rdrr.io/r/stats/logLik.html) wrong. It makes the
MODEL wrong: the fitted survival function rises with time there, and a
survival function that rises is not one.

This can happen with no event row affected at all. A random effect or a
covariate on `gamma1` gives each group its own slope, the
`log(gamma1 + u)` barrier that holds a slope positive lives in the
density, and a group whose rows are ALL censored contributes no density.
Measured on 40 centres of 10 with five centres followed to a common
administrative time with no deaths, seeds 20260910 to 20260915: on 4 of
6 seeds those five centres converge at slopes of -0.18 to -0.31, and one
centre's fitted survival goes from 4.2e-50 at `t = 1e-12` to 0.774 at
`t = 2.8`.

`n_nonmonotone_censored` counts those rows. They are counted apart from
`n_nonmonotone` because the two answer different questions:
`n_nonmonotone` says the reported likelihood is not the model's, and
`n_nonmonotone_censored` says the fitted model is not a survival
distribution where the data are. An interval-censored row is tested at
both ends.

`action = "error"` REFUSES on `n_nonmonotone` and WARNS on
`n_nonmonotone_censored`, and
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
and its two companions do the same, so on a fit with censored rows only
they answer with the warning. The warning gives `max_survival_rise`, the
largest amount the fitted survival climbs above its running minimum on
any flagged row's own coefficients, because the size is what a user
needs to judge it. The line is drawn where a negative hazard contradicts
an observed event, as flexsurv draws it:
[`flexsurv::dsurvspline()`](http://chjackson.github.io/flexsurv-dev/reference/Survspline.md)
sets the density to 0 where the derivative is not positive, so an event
row there makes the likelihood `-Inf` and the fit avoids it, while a
censored row is never checked and a rise there passes silently. frmtmb
draws the same line and says so. brms's `cox()` family cannot produce a
rising survival at all, because its baseline hazard is an M-spline basis
times non-negative weights; rstpm2 penalizes a negative hazard during
the fit.

It also fires on an ordinary design, and a user should expect it there:
a cure fraction with a time-varying effect. With `gamma1 ~ arm` and part
of each arm never having the event, the spline past an arm's last event
is fitted to censored rows only, and it can turn over there. Measured on
two arms of 250 with 40 and 55 percent cured, `df` 3 to 6, 20 seeds
each: 27 of 80 fits fire, every flagged row lies past its own arm's last
event, and the fitted survival rises across them by 5.2e-04 to 0.13. The
count is right, since that fitted survival function does rise, and the
rise is where the data carry no event to prevent it. With no
time-varying effect this cannot happen past the last event: every row
then shares one spline, and past the boundary knot its derivative is the
one at the last event row. A fit whose flagged rows are ALL censored
warns and does not refuse. A fit with a flagged EVENT row refuses, as in
0.7.0, whatever else it flags, and a cure design fitted with
proportional hazards can: its shared spline can turn over at an event
row. Fewer knots or a single `gamma1` remove a censored-row rise.

The test is at each row's own time. For `df = 1` the derivative is the
same at every time, so a row test is exact. For `df >= 2` the derivative
can dip between two observed times of one group and recover at both;
this check does not search for that.

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
#> List of 8
#>  $ n_censored_deep       : int 0
#>  $ max_nlogS             : num 2.56
#>  $ threshold             : num 19.2
#>  $ n_nonmonotone         : int 0
#>  $ n_nonmonotone_censored: int 0
#>  $ max_survival_rise     : num 0
#>  $ scale                 : chr "hazard"
#>  $ n_obs                 : int 300
#>  - attr(*, "rows")=List of 3
#>   ..$ censored            : int(0) 
#>   ..$ nonmonotone         : int(0) 
#>   ..$ nonmonotone_censored: int(0) 
```
