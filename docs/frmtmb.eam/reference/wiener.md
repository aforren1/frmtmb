# The Wiener first-passage time family

A drift-diffusion model for a two-choice decision: a noisy evidence
accumulator starts between two boundaries and the response time is the
first time it touches one of them. The family models the response time;
which boundary was touched is data, supplied through the `dec()`
addition term.

## Usage

``` r
wiener(
  max_ndt = NULL,
  variability = character(0),
  nodes = c(sz = 7L, st = 21L),
  allow_unreachable = FALSE,
  link = "identity"
)
```

## Arguments

- max_ndt:

  Upper bound for the non-decision time, in the units of the response,
  applied to every row. `NULL`, the default, takes the fastest response
  of each row's `ndt_group()`, or of the whole data set when the model
  has no `ndt_group()`. It cannot be combined with `ndt_group()`. Give
  it when a component of a
  [`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
  needs the bound up front.

- variability:

  Which across-trial variability parameters to estimate: any of `"sv"`
  (drift rate), `"sz"` (start point) and `"st"` (non-decision time). The
  default estimates none, which is the plain Wiener model.

- nodes:

  Gauss-Legendre node counts for the two quadratures, as a named vector.
  Only the entries for the `variability` parameters in use are read. The
  defaults are measured rather than chosen: `sz` reaches 1e-10 in the
  log density by 7 nodes everywhere it was probed, and `st` needs more
  because a response time below `ndt + st / 2` cuts the range and the
  integrand turns on sharply at the cut. See
  [`vignette("ddm")`](https://aforren1.github.io/frmtmb/frmtmb.eam/articles/ddm.md)
  for the measurement.

- allow_unreachable:

  Permit a `max_ndt` above the smallest response time. Only correct
  inside a mixture, where another component carries the rows the Wiener
  density cannot reach.

- link:

  Link for the drift rate. Identity by default, and there is rarely a
  reason to change it: the drift rate is signed.

## Value

A `frmtmb_family`.

## Details

The parameterization is brms's `wiener` family, name for name:

- `mu`:

  Drift rate, the mean rate of evidence accumulation (brms and the
  literature also call it `v`). Identity link, so it is signed: positive
  drift favors the upper boundary.

- `bs`:

  Boundary separation, the distance between the two boundaries (`a`).
  Log link.

- `ndt`:

  Non-decision time, the part of the response time spent encoding and
  moving rather than deciding (`t0` or `tau`). Bounded link; see
  Non-decision time below.

- `bias`:

  Relative start point in (0, 1), the fraction of the boundary
  separation the accumulator starts at (`w`). Logit link. 0.5 is
  unbiased.

`variability` adds Ratcliff's three across-trial variability parameters
to that set; see Across-trial variability below.

## The decision indicator

The boundary a trial ended at is data, and it reaches the density as an
addition term:

    frm(bf(rt | dec(response) ~ condition), family = wiener(), data = dat)

`dec()` is spelled as brms spells it and takes what brms takes: a factor
or character vector whose SECOND level is the upper boundary (so
`"lower"`/`"upper"` and `c(FALSE, TRUE)` both work as they read), or a
numeric 0/1 column. The package contributes the term to frmtmb's
addition-term registry when it loads.

`vint()` also still works, and carries the indicator as a plain 0/1
integer column:

    frm(bf(rt | vint(upper) ~ condition), family = wiener(), data = dat)

Supplying neither is refused with a message that says so, because the
failure is otherwise silent.

## The unit of the response

SECONDS. Every default in this package assumes it: the starting values,
the bound the `ndt` link is scaled onto, and `gddm_control(dt = )` and
the window it takes from the data. Nothing in any of the likelihoods
refuses milliseconds, and a fit to millisecond data converges and
reports a boundary separation three orders of magnitude out.

So every family here warns when the fastest response in the data is
above 20, names milliseconds as the likely cause and says to divide by
1000. It is a ceiling on the fastest response in the WHOLE data set
rather than on any one trial, so one slow trial does not reach it.

It is a warning rather than a refusal because it CAN fire on a correct
model: a slow task with a short session. Measured, the rate is zero over
192 designs at the parameters a two-choice task produces, and rises to
0.07 at 20 trials of a task whose median response is 39 seconds, 0.56 at
60 seconds and 1.00 at 100 seconds. A deliberation or matrix-reasoning
design is where that lands. The warning carries the class
`frmtmb_eam_units_warning` so that such a design can silence this one
condition and keep the rest. `NEWS.md` carries the full tables.

## Non-decision time, and what its coefficients mean

The density is zero for a response time at or below `ndt`, so the
likelihood has a hard edge at the fastest response and an ordinary log
link would let the optimizer walk straight over it. The `ndt` link makes
the constraint structural instead, and it does so in one of two ways.

**Without `ndt_group()`, `ndt` is a TIME**, on a logit scaled onto
`(0, ub)` with `ub` the fastest response in the whole data set, or
`max_ndt` when you give one. This is the parameterization the family has
always had. `predict(dpar = "ndt", type = "response")` reports seconds,
a `prior(class = "ndt")` is a density on those seconds, and a
`bf(ndt = 0.2)` constant is 0.2 seconds.

**With `ndt_group()`, `ndt` is a FRACTION of the row's own bound**, on a
plain logit, and the density multiplies it by that bound. The bound is
the fastest response of the row's group. Write

    frm(bf(rt | dec(response) + ndt_group(subject) ~ coherence,
           ndt ~ 1 + (1 | subject), bias = 0.5),
        family = wiener(), data = dat)

and each subject's non-decision time is bounded by its own fastest
response. That is what a random effect on `ndt` needs; the next section
is what one global bound does to it. The price is that `ndt` is on a
different scale: `predict(dpar = "ndt", type = "response")` reports the
fraction, and
[`ndt_time()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_time.md)
reports the time for either parameterization. A `prior(class = "ndt")`
and a `bf(ndt = )` constant are fractions under a grouping too.

Give the grouping as a factor, a character vector, a logical or integer
codes. It is keyed on the group's LABEL, so subsetting,
[`droplevels()`](https://rdrr.io/r/base/droplevels.html),
[`relevel()`](https://rdrr.io/r/stats/relevel.html) and a prediction
grid you build yourself all pair a row with the same bound the fit used.
`max_ndt` and `ndt_group()` cannot be combined, because they set the
same bound to different things, and an `ndt_group()` no family reads is
refused rather than carried into the fit unused.

A `max_ndt` above the smallest response time is refused, because for
this family alone it admits parameter values at which some observed row
has no likelihood. Inside a
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
that is exactly what the other component is for, so
`allow_unreachable = TRUE` lifts the refusal; see Mixtures. A mixture
never finalizes its components, so a component needs `max_ndt` given up
front and cannot use `ndt_group()`.

Whichever parameterization a model is in, the bound is fixed when the
model frame is assembled and is a property of the FITTED data, so a
prediction on new rows is scaled by the bound the fit used rather than
by one re-derived from the new rows. A group the fit never saw has no
bound and is refused.

Past a linear predictor of about 37 the logit saturates in double
precision. Nothing guards that, and nothing needs to: the density falls
off a cliff as the decision time goes to zero, so the log likelihood is
already unreachable long before the link runs out of digits.

## Why one bound is the wrong constraint under a random effect

A single bound is the GLOBAL fastest response, so a subject whose true
`ndt` is above it cannot be represented at any value of the random
effect. At 30 subjects by 400 trials with a between-subject spread of 26
ms on a mean of 250 ms, 20 of the 30 subjects are in that position and
NONE is inconsistent with its own data. Without `variability` such a fit
does not converge; with `variability = "sv"` it converges, reports
nothing from `diagnose()`, and returns a population `ndt` pinned at the
bound and wrong by ten percent with a standard error of 7.2e-06 on it.

What settles it is what happens as data accumulates. On that design, the
per-subject root mean squared error of the fitted non-decision times:

|        |                   |                  |
|--------|-------------------|------------------|
| trials | with ndt_group(s) | one global bound |
| 100    | 20.91 ms          | 27.23 ms         |
| 200    | 14.08 ms          | 20.51 ms         |
| 400    | 7.67 ms           | 30.37 ms         |

The per-group bound converges on the truth and the global bound does
not, because more data lowers the global minimum and tightens the
ceiling on every subject at once.

## What sd(ndt) does and does not tell you

A grouped model estimates `ndt` as a fraction of each group's own floor,
so the fitted per-subject times vary with those floors even when the
variance component is exactly zero. On the design above an estimator
with NO random effect on `ndt` returns a between-subject standard
deviation of 0.02748 against a truth of 0.02629, where the full model
returns 0.02543, and at 100 trials per subject the two are identical in
every digit.

The variance component is not empty: at 400 trials it buys 8.44
log-likelihood units and cuts the per-subject error from 11.90 ms to
7.67 ms. But a between-subject standard deviation is the wrong statistic
to read that off. Compare the per-subject non-decision times from
[`ndt_time()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_time.md)
against what you believe, and compare log-likelihoods, rather than
reading `VarCorr()`'s `ndt` row as evidence that the component was
estimated.

The bound's own quality is a function of the group's trial count. The
group's fastest response overshoots its true non-decision time by about
47 ms at 400 trials and 78 ms at 50, and a group of one trial gets that
trial's own response time as its bound. Nothing refuses a small group,
because any threshold would fire on correct models; the trial counts are
recorded on the fitted family, at `family(fit)$ndt_bound$sizes`.

## Across-trial variability

Ratcliff's full diffusion model draws three of the four parameters
afresh on every trial. `variability` names which of those to estimate,
and each one it names becomes an ordinary distributional parameter that
takes its own formula:

    frm(bf(rt | dec(response) ~ coherence, bias = 0.5),
        family = wiener(variability = c("sv", "sz", "st")), data = dat)

- `sv`:

  Standard deviation of a normal drift rate. Log link.

- `sz`:

  Width of a uniform relative start point, centered on `bias`, on the
  same (0, 1) scale as `bias`. Logit link, so the width is below 1 and
  the start point stays inside the boundaries whenever `bias` is 0.5.

- `st`:

  Width of a uniform non-decision time, centered on `ndt`, and on
  whichever scale `ndt` is on: a duration in the units of the response,
  on a logit scaled onto `(0, 2 * bound)`, or a FRACTION of twice the
  row's own bound under `ndt_group()`.

The likelihood is the analytic Wiener density averaged over those
distributions, and the three are done three different ways because they
are three different integrals. The drift integral is Gaussian against an
exponential-quadratic and is evaluated in CLOSED FORM: it is exact, it
takes no nodes, and there is nothing to tune. The other two are uniform
and are evaluated by fixed-node Gauss-Legendre quadrature, whose node
counts are the `nodes` argument.

The node positions and counts are decided when the family object is
built and are constants from then on, because a node count that moved
with a parameter would be a branch on a parameter and an
automatic-differentiation tape cannot record one. A parameter only
rescales the interval the fixed nodes are mapped onto.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) follow the
variability rather than ignoring it. Both condition on the boundary a
row ended at, and conditioning reweights which per-trial parameters that
row could have had: the trials that reached a boundary are not a fair
sample of the drift rates that could have produced them. So the fitted
mean is a ratio of two quadratures and the simulator accepts a drawn
drift rate and start point with the boundary probability they imply. The
plain closed-form mean is not a usable approximation here: at an
unbiased start point it returns the same number for both boundaries, and
the full model does not.

Two limits of the parameterization are worth knowing. The uniform start
point stays inside the boundaries by construction only when `bias` is
0.5, which is the usual case and the one the logit link on `sz` is
scaled for; at a strongly biased start a wide `sz` can push the range
past a boundary, where the density is near zero and the likelihood is a
barrier rather than a cliff. And nothing holds `ndt - st / 2` above
zero, so a fit is free to report a range that includes a negative
non-decision time. Neither can be made structural from outside frmtmb:
both are joint constraints on two distributional parameters, and a link
is a property of one.

## Mixtures

A contaminant component covers the trials the diffusion process cannot
produce, which is the standard treatment for fast guesses. The Wiener
component then wants a non-decision time that some observed rows fall
below, so pass the bound and lift the refusal:

    frm(bf(rt | dec(response) ~ 1, bias1 = 0.5),
        family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
                         lognormal()),
        data = dat)

A row below the non-decision time gets a log density of about
`-1 / delta` where `delta` is a billionth of the smallest response time:
it exponentiates to exactly zero, which is the right likelihood for the
component, and it differentiates to exactly zero, which a true `-Inf`
would not.

## Accuracy

The density is the Navarro and Fuss (2009) pair of series, both
evaluated at a fixed truncation and combined with a smooth weight,
because an automatic-differentiation tape cannot choose between them on
a parameter. It agrees with
[`RWiener::dwiener()`](https://rdrr.io/pkg/RWiener/man/wienerdist.html)
to better than 1e-12 relative on the log scale over normalized times
from 1e-3 to 50. See
[`vignette("ddm")`](https://aforren1.github.io/frmtmb/frmtmb.eam/articles/ddm.md).

## References

Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
for first-passage times in Wiener diffusion models. *Journal of
Mathematical Psychology*, 53(4), 222-230.

Ratcliff, R. and Tuerlinckx, F. (2002). Estimating parameters of the
diffusion model: approaches to dealing with contaminant reaction times
and parameter variability. *Psychonomic Bulletin & Review*, 9(3),
438-481.

## Examples

``` r
set.seed(1)
dat <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
           family = wiener(), data = dat)
fixef(fit)
#> $mu
#> (Intercept) 
#>   0.6956708 
#> 
#> $bs
#> (Intercept) 
#>   0.3349967 
#> 
#> $ndt
#> (Intercept) 
#>    1.935868 
#> 
#> $bias
#> (Intercept) 
#>           0 
#> 
```
