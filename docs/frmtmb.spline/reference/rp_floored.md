# Rows this family answered with a floor rather than with a density

Two things in
[`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
are floors rather than answers, and both are silent in the fitted
object: [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`AIC()`](https://rdrr.io/r/stats/AIC.html) report the floored value
with nothing to say it is floored. This function is where that goes to
be read, and by default it REFUSES rather than reports, because a fit in
either region is one whose numbers are not the model's.

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

  `"error"`, the default, refuses when either floor was used. `"report"`
  returns the same numbers without refusing.

- max_nlogS:

  The largest `-log S` on a censored row that is still scored
  accurately. The default 19.2 is where `eps / S` passes 1e-8.

## Value

A list with `n_censored_floored`, `max_nlogS`, `threshold`,
`n_nonmonotone`, `scale` and `n_obs`, returned invisibly when nothing
was floored. The offending row indices are the `"rows"` attribute, a
list with elements `censored` and `nonmonotone`.

## The censored-row floor

frmtmb forms a right-censored contribution as `log(1 - F(y))` on the
probability scale (`R/objective.R:100`), so the scored `log S` carries
absolute error about `.Machine$double.eps / S`. Past `-log S` of about
19.2 that error passes 1e-8; past 30 the term is FLAT, its gradient
exactly zero, and the optimizer prices the row at a constant.

The size of the error is a property of the data rather than of the
family: a floored row contributes -35.127363 instead of its own
`-log S`, so the reported log likelihood is short by about `-log S - 35`
per floored row. Two runs of one 600-subject design differing only in
seed give 2.4e+03 and 2.166e+04, both converged without a warning and
both with the treatment coefficient out by tens of percent.

The quantity checked is `-log S` at the fitted parameters, on every
censored row. It is one quantity for all three scales: `exp(eta)` on
`"hazard"`, `log1p(exp(eta))` on `"odds"` and `-log(Phi(-eta))` on
`"normal"`. On the hazard scale it is the cumulative hazard `H`.

The real fix is a complementary log-CDF slot in core, so that a family
can hand back `log S` instead of `F`. See `dev/spline-seam-proposal.md`
in the package sources.

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

The refusal is POST-FIT.
[`logLik()`](https://rdrr.io/r/stats/logLik.html) reads
`object$opt$objective` directly (`R/methods-fit.R:233-240`) and the
family protocol has no hook that runs when a fit finishes, so nothing in
this package can make [`logLik()`](https://rdrr.io/r/stats/logLik.html)
or [`AIC()`](https://rdrr.io/r/stats/AIC.html) refuse on their own. The
optimizer may therefore have walked through, or stopped inside, the flat
region before this function is ever called. Call it on every
[`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
fit whose data carry censoring;
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
and its two companions call it for you.

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
