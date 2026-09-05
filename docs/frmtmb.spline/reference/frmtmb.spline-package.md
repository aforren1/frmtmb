# frmtmb.spline: Spline Curves and Curve Inference for 'frmtmb' Models

Two things a fitted curve needs after it is fitted. First, curve
inference: a smooth evaluated on a grid with pointwise and SIMULTANEOUS
confidence bands, its first and second derivatives with delta-method
standard errors, and the features movement papers report, the time of a
peak and the time a curve crosses a level, each with a standard error
from the implicit-function delta method. The simultaneous band is the
max-deviation simulation of Ruppert, Wand and Carroll (2003), the same
construction 'gratia' uses on 'mgcv' fits, and it is checked against
'gratia' inside its Monte Carlo error. Second, a Royston and Parmar
(2002) flexible parametric survival family, which writes the log
cumulative hazard as a natural cubic spline in log time. It is
parameterized exactly as 'flexsurv::flexsurvspline' parameterizes it, so
the two log likelihoods are the same number rather than two numbers that
ought to agree, and covariates reach the spline coefficients as
proportional hazards on the first and as time-varying effects on the
rest.

## What this package reads that frmtmb does not promise

Two of the three things this package reads off a fitted object are
documented seams and one is not. It is named here so that a reader does
not have to discover it from the sources.

- `fit$estimates` and `fit$obj` both have precedent: the
  `frmtmb::frmtmb-extension-api` example reads `fit$estimates`, and
  [`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
  reads `fit$obj`.

- **`fit$cache$Vjoint` has no precedent and is an internal.** It is the
  memo written by frmtmb's `get_joint_cov()`, which is not exported and
  not documented, and neither the slot nor its `list(V =, names =)`
  shape appears anywhere in frmtmb's own documentation.
  [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
  reads it because the alternative, recomputing the joint covariance
  here, was both far slower and WRONG on an autoscaled fit.

The consequences of that reach are bounded and are set out in full under
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md),
section "The one internal this reaches into": a change to the name or
the shape costs speed, a change to the meaning is caught by the
covariance check that every call makes, and absence is the ordinary case
on a fresh fit and is handled by warming the cache first. None of the
three can produce a wrong number.

`dev/spline-seam-proposal.md` Part 1a asks core for an exported
accessor, which would retire the reach entirely.

## Two limits that are core's to fix

- [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) on a
  [`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
  fit report whatever the optimizer reached, and neither can be gated
  from an extension: [`logLik()`](https://rdrr.io/r/stats/logLik.html)
  reads `object$opt$objective` and the family protocol has no hook that
  runs when a fit finishes. Where the likelihood was floored, those two
  numbers are wrong and say nothing.
  [`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md)
  refuses, and
  [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
  calls it, but a user who reads
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and nothing else gets no
  signal. The fix is core's: an `lccdf` slot, or a post-fit family hook.

- A mapped random-effect block is untested here, because
  [`frmtmb::frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.html)
  takes no `map` argument and there is therefore no supported route to
  one. The nearest reachable analogue is a distributional parameter held
  fixed, which sets `betad_fixed_idx` and takes the same index-remapping
  path: `bf(y ~ s(x, k = 8), sigma = 0.5)` works at a measured
  `cov_rel_error` of 1.55e-15.

## See also

Useful links:

- <https://aforren1.github.io/frmtmb/frmtmb.spline>

- <https://github.com/aforren1/frmtmb>

- Report bugs at <https://github.com/aforren1/frmtmb/issues>

## Author

**Maintainer**: Alex Forrence <alex.forrence@gmail.com>
([ORCID](https://orcid.org/0000-0002-9728-6337))

Authors:

- Alex Forrence <alex.forrence@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9728-6337))
