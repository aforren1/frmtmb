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
