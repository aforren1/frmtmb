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

## The linear ballistic accumulator

* `lba(n)` is a race between `n` accumulators, each rising in a straight
  line from a uniform start point on `(0, A)` to a common threshold at a
  normally distributed rate. It is the family for choices with more than
  two alternatives: a diffusion between two absorbing boundaries admits
  exactly two responses, and no reparameterization of `wiener()` reaches
  a third. The likelihood is closed form for any `n`.
* Each accumulator's drift mean is its own distributional parameter
  (`v1`, ..., `vn`, identity link) and all of them are primary, so the
  main formula reaches every drift and each gets its own coefficients.
  Give one its own formula and a covariate moves that alternative alone,
  which is the capability neither two-choice family can offer.
* The remaining parameters are `A` (start-point range, log), `k`
  (threshold above the start-point range, log) and `ndt` (non-decision
  time, a logit scaled onto `(0, max_ndt)` as in `wiener()`). The
  threshold is `A + k` rather than a parameter of its own, so `b > A` is
  structural: a threshold inside the start-point range, where trials
  would begin already finished, is not a state the optimizer can reach.
* The drift standard deviation is fixed, because the model is identified
  only up to a common rescaling of `A`, the threshold, the drifts and
  the drift standard deviation. It is the family argument `sd_v` rather
  than a hidden default, and is carried on the family object.
* Drift rates are truncated at zero by default, following `rtdists` and
  its `posdrift = TRUE`. Every accumulator then arrives eventually and
  the choice probabilities sum to one. `posdrift = FALSE` gives the
  untruncated convention, under which the response distribution is
  defective; the two are different models, not a rescaling of each
  other, which matters when comparing against another package.
* The single-accumulator density agrees with `rtdists::dlba_norm()` to
  better than 1e-11 relative wherever `rtdists` is itself accurate, and
  the race with `rtdists::n1PDF()` for two, three and four
  accumulators. In the fast tail the two diverge by design and this one
  is the better: `rtdists` writes the normal difference
  `Phi(g) - Phi(h)` as a subtraction of two lower tails, which returns
  exactly zero once `pnorm(h)` saturates. Written in log space instead,
  as `exp(la) * -expm1(lb - la)` on the two upper-tail logs, nothing
  saturates and no comparison is needed, which matters because RTMB
  refuses comparison on AD types. Against a 200-bit Rmpfr reference on
  1144 points, the subtractive form returns exactly zero on 68 of the
  80 tail rows and is already 2.7e-3 wrong on 8 rows outside the tail;
  the log-space form holds 5.0e-14 in the bulk and 3.6e-4 in the tail.
  Four rows still exceed 1e-6, so it is an improvement, not a proof.
* That mattered beyond accuracy. The subtractive form errs in the value
  only, and a tape differentiates the function that was written, so
  value and gradient described different surfaces near the
  non-decision-time bound. It did not move the point estimate, which
  the likelihood keeps out of that region, but Hessian-based standard
  errors were up to 16.5 percent off before the change.
* The survival function is written out directly rather than as
  `1 - plba_norm()`, which returns exactly zero, and so a
  log-contribution of `-Inf`, while the true survival is still around
  1e-19; the direct form agrees with a quadrature of the density down
  to survivals of 1e-23. Under `posdrift = FALSE` the never-arriving
  mass it adds back is spelled `Phi(-v/s)` rather than `1 - Phi(v/s)`,
  which is the same saturation trap one line further on.
* The choice reaches the density through `vint()`, as the Wiener
  family's boundary indicator does, and is declared with
  `required_aterms`, so omitting it is refused by name. A choice outside
  `1..n`, a non-positive response time and a non-decision-time bound
  above the fastest response are each refused with their own message.
* The non-decision-time bound is derived from the response through
  `family_finalize()`, so the family object carries the link it will
  actually use rather than an environment filled in later.
* `lba_simulate()` draws choice and time jointly from the generative
  process. `simulate()` on a fitted object redraws times conditional on
  each row's observed choice, since the choice is data.
* Deliberate omissions: no `lcdf`, so `cens()` and `trunc()` are
  refused; no `post$mean_fn`, because the mean of the race has no closed
  form and a quadrature per row was not worth writing for `fitted()`.
