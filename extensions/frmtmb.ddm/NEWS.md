# frmtmb.ddm 0.2.0

Ratcliff's full diffusion model, and the queued defects.

* `wiener(variability = )` adds across-trial variability to the family
  rather than forking it. Naming any of `"sv"` (drift rate), `"sz"`
  (start point) and `"st"` (non-decision time) turns that one into an
  ordinary distributional parameter, with its own link and its own
  formula. `wiener()` with no arguments is the model it always was: on
  a grid of 11520 parameter combinations its log density is
  byte-for-byte identical to the previous release everywhere the
  density is defined, and its tests are unchanged. The one difference
  is below the non-decision time, where it now returns `-Inf` instead
  of `NaN`, which is the deliberate fix below.
* The likelihood is the analytic Wiener density averaged over those
  distributions, and the three are done three different ways. The drift
  integral is Gaussian against an exponential-quadratic and is evaluated
  in CLOSED FORM: it agrees with adaptive quadrature of the same thing
  to better than 1e-13 relative and takes no nodes, so estimating `sv`
  is FREE relative to the plain density. Measured at the fit level over
  three seeds, `sv` costs 0.97 to 1.07 times plain, and at the density
  level both are 0.0002 s per call. The start-point and
  non-decision-time integrals are uniform and use fixed-node
  Gauss-Legendre quadrature, and those do cost: 0.0056 s per call for
  `sz` at 7 nodes and 0.0170 s for `st` at 21.
* Node counts are the `nodes` argument and the defaults are measured:
  `sz` reaches machine precision at 7 nodes, `st` reaches 1e-9 at 21.
  They differ because the non-decision-time range is cut by the response
  time on a fast trial, and the integrand turns on sharply at the cut.
  Node positions and counts are fixed when the family object is built,
  because an automatic-differentiation tape cannot record a branch on a
  parameter; a parameter only rescales the interval they map onto.
* Variability parameters at zero reproduce the plain Wiener density to
  better than 1e-13 in the log density, which is floating-point rounding
  on a differently associated sum rather than a quadrature error.
* Recovery on simulated data, 8 replicates of 1500 trials with a normal
  drift rate and a uniform non-decision time: every parameter's Monte
  Carlo mean is within 2.2 of its own standard errors of the value it
  was generated from.
* Across-trial variability is NOT frmtmb's `quadrature = TRUE`. That
  marginalizes random effects by Gauss-Kronrod, is wired to the
  random-effect coefficient vector by name, and refuses a model with no
  random-effect block. This integral shares nothing between trials, has
  no level to estimate, and exists in models with no grouping factor at
  all, so it lives inside the density.
* `ddm_simulate()` takes `sv`, `sz` and `st`, drawing each trial's own
  parameters and then running the ordinary process, so the simulator
  states the model independently of the density.
* `fitted()` and `simulate()` follow the variability rather than
  ignoring it. Conditioning on the boundary a row ended at reweights
  which per-trial parameters that row could have had, so the fitted mean
  is a ratio of two quadratures and the simulator accepts a drawn drift
  rate and start point with the boundary probability it implies. The
  plain closed form would not do: at an unbiased start point it returns
  the same mean for both boundaries, where 40000 simulated trials put
  them 0.06 s apart. Both are checked against simulated data, which
  knows nothing about how either is computed.

## Defects fixed

* **The density returns `-Inf` below the non-decision time, not `NaN`.**
  A likelihood of zero is the right answer there; `NaN` is not an answer
  and propagates through every component of a mixture. The normalized
  time is held at the smallest positive double, which is bit-for-bit
  inert everywhere the density was already right.
* **`wiener(max_ndt = )` above the smallest response time is usable in a
  mixture.** The refusal is correct for the family on its own and wrong
  inside a mixture, where the contaminant component is exactly what
  covers the trials the diffusion cannot produce, so
  `allow_unreachable = TRUE` lifts it and the refusal now names it. A
  trial below the non-decision time then gets a log density that
  exponentiates to zero AND differentiates to zero, which a true `-Inf`
  does not: `-Inf` yields a `NaN` gradient and stops the fit.

## What frmtmb 0.49.0 let this package delete

* **`dec()` is the spelling now.** frmtmb gained
  `frmtmb_register_aterm()`, this package registers `dec` when it loads,
  and `rt | dec(response) ~ x` works and takes a factor, a character
  vector or a logical the way brms does. `vint()` carries the same thing
  as a 0/1 integer and is unchanged.
* **The environment the link closures read is gone.** The bound on the
  non-decision time is a property of the response, and the family object
  is built before `frm()` has any; this package used to have `valid_y()`
  write the bound into an environment, which worked only for as long as
  an undocumented slot order held. `family_finalize()` is the documented
  slot for it and the family now derives itself from the data there.
* One hand-rolled check remains, and is not `required_aterms`'s fault:
  that argument names the terms a density needs ALL of, and this family
  needs EITHER `dec()` or `vint()`. See `dev-findings.md`.

# frmtmb.ddm 0.1.0

First release. A Wiener first-passage time family for two-choice
response times, written entirely against frmtmb's exported extension
API.

* `wiener()` is the drift-diffusion family, with brms's dpar names and
  links: `mu` (drift rate, identity), `bs` (boundary separation, log),
  `ndt` (non-decision time, bounded) and `bias` (relative start point,
  logit). A model written for `brms::wiener()` reads the same here, and
  the parameterization is pinned against brms in the test suite.
* The density is the Navarro and Fuss (2009) pair of series. Both are
  evaluated at a fixed truncation and their logs blended with a
  logistic weight in the normalized time, because an AD tape cannot
  choose between them on a parameter. It agrees with
  `RWiener::dwiener()` to better than 1e-11 relative over normalized
  times from below 1e-3 to 50, where a fixed truncation of the
  small-time series alone is wrong by tens of percent past about 8.
* The decision indicator reaches the density through `vint()`, because
  frmtmb's addition terms are a closed set and brms's `dec()` is not one
  of them. Omitting it is refused with a message naming both spellings,
  rather than silently producing a log likelihood over no rows.
* `ndt` uses a logit scaled onto `(0, max_ndt)`, so the support
  constraint is structural rather than something the optimizer has to
  discover. `max_ndt` defaults to the smallest response time in the
  data, found at frame assembly.
* `post$mean_fn` gives the mean response time in closed form,
  conditional on the boundary the row ended at, with the zero-drift
  limits handled explicitly. `sim` draws conditionally by inverse
  transform through `RWiener::qwiener()`, with a discretized forward
  simulation as the fallback when RWiener is absent.
* `ddm_simulate()` draws response times and boundary choices jointly,
  for building example and test data sets.
* The package registers its compatibility rows with the core at load
  through `frmtmb::frmtmb_register_compat()`, so `frm_compat("wiener")`
  states what was exercised and what was not.
* Deliberate omissions: no `lcdf`, so `cens()` and `trunc()` are
  refused; no variance function and no unit deviance, so
  `residuals()` answers `type = "response"` only.
* `dev-findings.md` in the package source records what building this
  from outside frmtmb cost, as an acceptance test of the extension API.
