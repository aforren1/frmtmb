# The Wiener first-passage time family

A drift-diffusion model for a two-choice decision: a noisy evidence
accumulator starts between two boundaries and the response time is the
first time it touches one of them. The family models the response time;
which boundary was touched is data, supplied through the `vint()`
addition term.

## Usage

``` r
wiener(max_ndt = NULL, link = "identity")
```

## Arguments

- max_ndt:

  Upper bound for the non-decision time, in the units of the response.
  `NULL`, the default, takes it from the data.

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

## The decision indicator

brms spells the boundary a trial ended at as `y | dec(decision)`.
frmtmb's formula grammar has no `dec()`, and its list of addition terms
is closed, so this family reads the indicator from `vint()` instead:

    frm(bf(rt | vint(upper) ~ condition), family = wiener(), data = dat)

`upper` must be 1 for a response at the upper boundary and 0 for the
lower one. Unlike brms's `dec()`, it will not accept a factor or the
strings `"upper"` and `"lower"`; code it yourself. Omitting `vint()`
altogether is refused with a message that says this, because the failure
is otherwise silent.

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

Past a linear predictor of about 37 the logit saturates in double
precision and `ndt` rounds to `max_ndt` exactly. Nothing guards that,
and nothing needs to: the density falls off a cliff as the decision time
goes to zero, so the log likelihood is already unreachable long before
the link runs out of digits.

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

## Examples

``` r
set.seed(1)
dat <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
fit <- frm(bf(rt | vint(upper) ~ 1, bias = 0.5),
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
