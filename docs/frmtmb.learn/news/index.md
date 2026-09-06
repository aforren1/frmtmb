# Changelog

## frmtmb.learn 0.1.0

First release. The value-learning models of the reinforcement-learning
and computational-psychiatry literature, written as frmtmb families.

### What the grammar buys

- Every parameter of every family is an ordinary distributional
  parameter with its own linear predictor. A learning rate takes a
  condition effect, a smooth term or a correlated per-subject random
  effect the way a mean does, and the effect comes back with a standard
  error. A reversal task fitted as `alpha ~ after_reversal + (1 | id)`
  gives the CHANGE in learning rate with an interval, which two
  separately fitted models do not. There is no separate reversal family
  here for that reason.

### One recursion, six families

- [`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md),
  the Rescorla-Wagner delta rule on two arms.

- [`bandit2arm_dual()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_dual.md),
  separate rates for gains and losses. `split` chooses whether the split
  is on the sign of the PREDICTION ERROR (the default) or of the
  OUTCOME, which is hBayesDM’s `prl_rp`. Measured: on binary payoffs the
  two are the SAME model, because the value store stays inside the
  payoff range and the two signs agree on every trial; they differ only
  on graded payoffs.

- [`prl_fictitious()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/prl_fictitious.md),
  counterfactual updating, where the option that was not chosen moves in
  the opposite direction.

- [`bandit4arm2_kalman_filter()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit4arm2_kalman_filter.md),
  the restless four-armed bandit of Daw and others (2006). Its learning
  rate is not a parameter: the Kalman gain is derived from a posterior
  variance the filter carries, so it falls as a subject learns.

- [`ts_par7()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/ts_par7.md),
  the two-stage model-based and model-free hybrid of Daw and others
  (2011), with `w`, `lambda`, two learning rates, two sensitivities and
  perseveration.

- [`igt_pvl_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_pvl_delta.md),
  prospect-valence learning for the Iowa gambling task. ORL is not here;
  [`?igt_pvl_delta`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_pvl_delta.md)
  says why.

- The engine is internal and is the point of the package. It walks
  trials once and updates every subject at each step, taking the
  learning rule and the choice rule as tape-safe functions, with padding
  and a mask for unequal trial counts. A new family is the model and
  nothing else. It runs on the tape and off it, so the likelihood, the
  fitted trajectory and the simulator are ONE recursion rather than
  three copies that have to be kept in step.

### Reading a fit

- [`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
  returns the per-trial value estimates each choice was made on, the
  prediction error the outcome produced, and the fitted probability of
  the choice that was made. These are the quantities papers plot.
- [`frm_task_design()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_design.md)
  builds a trial-level design for each task with the payoff schedule of
  every option fixed in advance, which is what makes a draw coherent: a
  simulated subject that takes the arm the real one did not still has to
  be paid.
- [`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md)
  draws whole datasets from the generative process, taking parameters
  directly rather than through a formula.
  [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.html)
  is the other route and goes through the fitted grammar; the two share
  the recursion and are checked against each other.
- [`frm_learn_families()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_learn_families.md)
  is the reference table, including the map from each parameter to
  hBayesDM’s spelling of it.

### What these families refuse, and why

- **[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  `predict(type = "response")` and every kind of residual refuse.** The
  response is the option a subject took, coded 1 to K, and it is
  nominal: arm 2 is not twice arm 1 and the Iowa gambling task’s four
  decks have no order at all. Core forms a fitted value as the
  conditional mean of the response and a residual as `y - mean`, so a
  mean declared here would make three methods return arithmetic on a
  category code.
  [`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
  returns everything a mean would have carried and more. Note that
  frmtmb’s own `rw_delta` example DOES have
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): it codes a
  two-armed choice 0 and 1 and returns P(arm 1), which does not survive
  a fourth arm.
- **`frm(importance =)` is refused by name for every family**, and the
  refusal is wider than the mathematics. The correction needs one
  log-likelihood value per group and these likelihoods have one, but
  `frmtmb_structure(loglik =)` returns a single AD scalar and there is
  nowhere to put the factors.
  `frm_compat("bandit2arm_delta", "importance")` names the seam at the
  console. The Laplace error is measured by simulation instead.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html) is refused for
  [`ts_par7()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/ts_par7.md)
  alone: one two-step trial’s draw is three numbers and a response
  vector holds one. Use
  [`frm_task_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_simulate.md).
- [`weights()`](https://rdrr.io/r/stats/weights.html), `cens()`,
  [`trunc()`](https://rdrr.io/r/base/Round.html) and `se()` refuse in
  the family’s own words: each reshapes a per-row likelihood
  contribution, and a trial’s contribution here is conditional on every
  earlier trial of the same subject.

### Validation

- Every family is checked by an IDENTITY against an independent Stan
  program of the same model, written from the published equations,
  evaluated at frmtmb’s own estimates. The programs declare exactly the
  parameters frmtmb estimates on the scales frmtmb estimates them on, so
  the map is the identity and the comparison carries no Jacobian. The
  tier is gated on `FRMTMB_BRMS_FIT_TESTS` and `NOT_CRAN`, like the brms
  comparison tier, because it compiles Stan programs.
- Parameter recovery at a realistic scale, with bias and coverage.
- The Laplace caveat is measured by simulation at few trials per
  subject, and the numbers are in
  [`vignette("learning")`](https://aforren1.github.io/frmtmb/frmtmb.learn/articles/learning.md)
  and in
  [`?frmtmb.learn`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frmtmb.learn-package.md).
- No hBayesDM code or data is used anywhere. Its models are published
  equations and the families here are written from them. A cross-check
  against it lives in `dev/`, is never run by the suite, and has not
  been run at all: hBayesDM 2.0.0 fits through cmdstanr, and CmdStan’s
  vendored TBB does not build under this machine’s compiler. The
  script’s header records the whole chain. So this release has no check
  against an outside implementation; what it has instead is an exact
  identity against its own independent Stan programs and a longhand
  reference per family.

### Deliberate omissions

- RLDDM, a learning rule feeding a drift-diffusion choice rule, is the
  next family and is not built: the choice rule belongs to `frmtmb.eam`,
  and building it here would make this package import a sibling
  extension. The engine already has the shape it needs.
- No exploration bonus on the Kalman filter, and no ORL for the Iowa
  gambling task. Both are short additions on this engine and both are
  left out for stated reasons in their families’ help.
