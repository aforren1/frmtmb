# The Royston and Parmar flexible parametric survival family

Models the log cumulative hazard as a natural cubic spline in log time.
It is the flexible parametric survival model of Royston and Parmar
(2002), and it is parameterized exactly as `flexsurv::flexsurvspline()`
parameterizes it, so a `gamma` vector means the same curve in both
packages and the two log likelihoods are the same number.

## Usage

``` r
royston_parmar(
  df = 3,
  knots = NULL,
  bknots = NULL,
  scale = c("hazard", "odds", "normal")
)
```

## Arguments

- df:

  Degrees of freedom: the number of interior knots plus one, which is
  Royston and Parmar's convention. `df = 1` is no interior knot at all
  and gives the Weibull, log-logistic or lognormal model that `scale`
  names. flexsurv counts the same spline with `k = df - 1`.

- knots:

  Interior knots on the log time scale, given explicitly. `df` is then
  read off their number and any `df` argument is checked against them.

- bknots:

  The two boundary knots, on the log time scale. Defaults to the range
  of the log uncensored times.

- scale:

  `"hazard"`, `"odds"` or `"normal"`.

## Value

A `frmtmb_family`.

## Details

Writing `x = log(t)`, the family fits

\$\$g(S(t)) = \gamma_0 + \gamma_1 x + \sum_j \gamma\_{j+1} v_j(x)\$\$

where `v_j` are the natural cubic spline basis functions of Royston and
Parmar, linear beyond the boundary knots, and `g` is chosen by `scale`:

- `"hazard"`:

  `g(S) = log(-log(S))`, the log cumulative hazard. Covariates on the
  first coefficient are PROPORTIONAL HAZARDS. With no interior knots
  this is a Weibull model.

- `"odds"`:

  `g(S) = log(1/S - 1)`, the log cumulative odds of failure. Covariates
  on the first coefficient are proportional odds. With no interior knots
  this is a log-logistic model.

- `"normal"`:

  `g(S) = -Phi^-1(S)`, the probit scale. With no interior knots this is
  a lognormal model.

## The distributional parameters

`mu` is `gamma_0` and every other coefficient is `gamma1`, `gamma2`, and
so on, all with identity links. There are `df + 1` of them.

A formula on `mu` is the proportional-hazards (or -odds, or -probit)
model, because `gamma_0` shifts the whole curve. A formula on any other
coefficient is a TIME-VARYING effect, since that coefficient multiplies
a function of log time; this is what `flexsurv` spells `anc =`.

    frm(bf(t | cens(censored) ~ trt, gamma1 ~ trt),
        family = royston_parmar(df = 3), data = dat)

## Knots

The spline needs `df + 1` knots: two boundary knots and `df - 1`
interior ones. By default they are placed at equally spaced quantiles of
the log UNCENSORED times, which is Royston and Parmar's own rule and
flexsurv's default; the boundary knots are then the smallest and largest
log uncensored time.

The quantiles cannot be taken until the response is in hand, so they are
taken at frame assembly through `family_finalize()`, and the family
object the fit carries has its knots baked in. `knots =` and `bknots =`
pin them instead, in which case they are taken as given and no data is
consulted. Both are on the LOG time scale, as flexsurv's are.

## Monotonicity is not enforced, and the floor is not free

The cumulative hazard has to increase, which means the spline's
derivative in `x` has to stay positive, and nothing in this
parameterization holds it there. flexsurv does not hold it either: the
model is fitted unconstrained and a fit whose spline turns over inside
the data range is over-parameterized, not something the software should
have prevented.

What this family does differently is refuse to return `NaN` for it.
Where the derivative goes non-positive there is no hazard and the true
log density is `-Inf`, so the log density becomes a large finite
negative number instead: `NaN` stops the optimizer, and inside a
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
one `NaN` component poisons the log-sum-exp of all of them.

That floor is NOT inert when it is used. It keeps the fit alive and it
makes [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`AIC()`](https://rdrr.io/r/stats/AIC.html) a pseudo-likelihood rather
than a density: a 60 percent cure-fraction dataset has been measured
converging, without a warning, with 6 floored rows and a reported log
likelihood 3952 units away from the density's. Nothing in the fitted
object says so, because the family protocol has no hook that runs when a
fit finishes.

[`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md)
is where that goes to be read, and it REFUSES by default. Call it on
every fit.
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
and its two companions call it for you, because the fitted curve of a
floored fit is a floor artifact too.

## Censoring, truncation, and the accuracy limit on log S

The family declares both a density and a distribution function, so
`cens()` and [`trunc()`](https://rdrr.io/r/base/Round.html) both work,
and right, left and interval censoring all reach the likelihood through
frmtmb's own machinery. `cens()` takes frmtmb's coding: `0` (or
`"none"`) is an observed event and `1` (or `"right"`) is right censored,
which is the OPPOSITE of the `status` column of a `Surv()` object. Pass
`1 - status`.

Both are CONDITIONAL rather than unqualified. frmtmb forms a
right-censored contribution as `log(1 - F(y))` on the PROBABILITY scale
(`R/objective.R:100`), and core offers a family no complementary log-CDF
slot to hand back `log S` directly, so the scored `log S` carries
absolute error about `.Machine$double.eps / S` whatever this family does
internally. Measured, by forming `F` for a given `-log S` and reading
`log(1 - F)` back (the three scales agree to every printed digit):

|          |           |                |
|----------|-----------|----------------|
| `-log S` | computed  | absolute error |
| 10       | -10       | 1.3e-13        |
| 19.2     | -19.2     | 2.4e-10        |
| 30       | -29.99983 | 1.7e-04        |
| 36       | -34.94504 | 1.05           |
| 40       | -35.12736 | 4.87           |

`log S` is floored at -35.127363 whatever the model says, and past
`-log S` of 30 the term is FLAT: its gradient is exactly zero, so the
optimizer prices that row at a constant and fits the rest as if it were
free.

How wrong the answer gets is a property of the DATA, not of the family.
A floored censored row contributes -35.127363 instead of its own
`-log S`, so the reported log likelihood is short by about `-log S - 35`
for each such row. That is thousands to tens of thousands as soon as one
censored time sits well past the event times: two runs of the same
600-subject design, differing only in seed, give 2.4e+03 and 2.166e+04.
Both converged without a warning and both put the treatment coefficient
out by tens of percent, on data flexsurv declines to fit at all.

`-log S` is the cumulative hazard `H` on the `"hazard"` scale,
`log(1 + exp(eta))` on `"odds"` and `-log(Phi(-eta))` on `"normal"`.
[`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md)
checks it on every censored row and refuses past 19.2, which is where
`eps / S` passes 1e-8. The real fix is a core one, an `lccdf` slot, and
it is in `dev/spline-seam-proposal.md`.

## References

Royston, P. and Parmar, M. K. B. (2002) Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data. *Statistics in Medicine* 21, 2175-2197.

## See also

[`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md),
which must be called on any fit whose data carry censoring;
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
for reading the fitted log cumulative hazard off with a band.

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
frmtmb::fixef(fit)$mu
#> (Intercept)         trt 
#>  -1.7692803   0.6424672 
```
