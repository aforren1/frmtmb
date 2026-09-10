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

Nothing, since frmtmb 0.52.0. Everything read off a fitted object is now
a documented seam.

- `fit$estimates` and `fit$obj` both have precedent: the
  `frmtmb::frmtmb-extension-api` example reads `fit$estimates`, and
  [`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
  reads `fit$obj`.

- The joint covariance of a grid prediction comes from
  [`frmtmb::frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.html),
  which is exported and documented. Up to frmtmb 0.51.0 there was no
  such seam and
  [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
  read `fit$cache$Vjoint`, the memo written by an unexported
  `get_joint_cov()`. That reach is gone, and with it the whole section
  this one replaces.

The covariance check survives the change and still runs on every call:
what it now verifies is that this package reads the seam correctly,
rather than that a reconstruction reproduced core's number.
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md),
section "The route to the covariance", sets it out.

## What core has since supplied, and what is still missing

Three things this package used to work around are seams in frmtmb
0.52.0.

- `lccdf`. A right-censored row is scored from `log S` directly, in
  closed form on all three
  [`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
  scales, so the floor at -35.127363 and the flat region past
  `-log S = 30` are both gone.

- `post$fit_check`, the fit-end family hook. A non-monotone fit warns as
  it is returned rather than only when someone calls
  [`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md).

- [`frmtmb::frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.html),
  above.

What is still missing is one thing, and it is the protocol's rather than
this family's. `frmtmb_structure(loglik =)` returns one AD scalar, so
core never sees the individual factors of a likelihood that factorizes.
A **per-row log-likelihood slot** would give `loo()` and `waic()` the
pointwise matrix they need, and a **per-group log-likelihood slot**
would give `frm(importance =)` the one value per group it corrects with.
Neither exists, so [`logLik()`](https://rdrr.io/r/stats/logLik.html) and
[`AIC()`](https://rdrr.io/r/stats/AIC.html) on a
[`royston_parmar()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/royston_parmar.md)
fit still report whatever the optimizer reached with no way to say which
rows carried it;
[`rp_floored()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/rp_floored.md)
is where that is read instead.

## One limit that is core's to fix

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
