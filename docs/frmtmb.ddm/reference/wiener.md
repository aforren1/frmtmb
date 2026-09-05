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

  Upper bound for the non-decision time, in the units of the response.
  `NULL`, the default, takes it from the data.

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
  [`vignette("ddm")`](https://aforren1.github.io/frmtmb/frmtmb.ddm/articles/ddm.md)
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

## Non-decision time

The density is zero for a response time at or below `ndt`, so the
likelihood has a hard edge at `ndt = min(rt)` and an ordinary log link
would let the optimizer walk straight over it. The `ndt` link is a logit
scaled onto `(0, max_ndt)` instead, which makes the constraint
structural rather than a thing the optimizer has to discover.

`max_ndt` defaults to the smallest response time in the data, found when
the model frame is assembled. Give it explicitly to pin the bound, which
is worth doing when you will
[`predict()`](https://rdrr.io/r/stats/predict.html) on new data whose
minimum differs from the training minimum.

A `max_ndt` above the smallest response time is refused, because for
this family alone it admits parameter values at which some observed row
has no likelihood. Inside a
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
that is exactly what the other component is for, so
`allow_unreachable = TRUE` lifts the refusal; see Mixtures.

Past a linear predictor of about 37 the logit saturates in double
precision and `ndt` rounds to `max_ndt` exactly. Nothing guards that,
and nothing needs to: the density falls off a cliff as the decision time
goes to zero, so the log likelihood is already unreachable long before
the link runs out of digits.

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

  Width of a uniform non-decision time, centered on `ndt`, in the units
  of the response. Logit link scaled onto `(0, 2 * max_ndt)`.

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
[`vignette("ddm")`](https://aforren1.github.io/frmtmb/frmtmb.ddm/articles/ddm.md).

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
