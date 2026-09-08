# frmtmb.learn 0.2.0

The three things 0.1.0 named as left out are in, and the seam it was
waiting on has closed. Requires frmtmb 0.53.0 for the two factorization
slots, and frmtmb.eam 0.4.1 for the Wiener density.

## frm(importance = ) works

* Every family declares `frmtmb_structure(loglik_row =)` and
  `(loglik_group =)` beside `loglik`, so `frm(importance = )` is
  admitted instead of refused by name. A subject is the group and a
  trial is the row: two subjects share nothing but the parameters, and
  a subject's likelihood is a product over trials given the value
  trajectory. Both facts are properties of the recursion rather than of
  any one learning rule, so both hold for all eight families.

  All three quantities come off ONE call to the recursion, which now
  yields per-trial factors rather than a running total. The change
  costs nothing: measured at 40 subjects by 100 trials, tape build
  0.250 s against 0.260 s and one gradient 98.5 ms against 103.0 ms,
  with the objective identical to full printed precision.

  `unit` is unchanged, "one subject's trial sequence". These slots say
  how finely the likelihood FACTORIZES; `unit` says what may honestly
  be left OUT, and dropping a trial changes every later trial's value
  store.

* **What the correction is worth, measured, and it depends on the
  design.** `bandit2arm_delta()` on a reversal task, 40 subjects, one
  scalar random intercept, truth `sd(alpha) = 0.5`, six datasets per
  cell. At 100 trials, on the three of six datasets whose variance
  component the Laplace fit had not collapsed, the correction moves
  `sd(alpha)` from 0.38-0.45 to 0.45-0.51, that is up and toward the
  truth, by 0.11 to 0.18 log units, with a Monte Carlo standard error
  of 0.05-0.06 and a smallest effective sample size per draw of
  0.87-0.93. Of the other three, two refused on the correction's own
  convergence guard and ONE HAD ALREADY COLLAPSED to 0.0004 and
  reported the same spurious +1.822 the short-session cell does, so a
  collapsed component is not a short-session problem and
  `sqrt(VarCorr(fit))` is worth checking at every trial count. The
  fixed effects move by less than 0.02.

  At 20 trials it is not usable, and the reason is worth knowing: five
  of six Laplace fits have already collapsed `sd(alpha)` to 1e-4, so
  there is nothing to reweight, and the correction then walks at its
  own step cap for every round while WARNING that it did not converge.
  A shift reported with `fit$importance$capped` true and every entry of
  `fit$importance$moves` equal is that step cap, not an estimate, and
  the warning's advice to raise `importance_rounds` enlarges it in
  exact proportion: five rounds of a 0.3645 cap give 1.822 and ten give
  3.645. Check `sqrt(VarCorr(fit))` before believing a correction.

  The model must be grouped on the family's own subject. `(1 | other)`
  gives a per-subject likelihood against a per-other proposal, and the
  core refuses it by name.

* `residuals(type = "deviance")` still refuses, and for a different
  reason than before. The magnitude is now available, because
  `loglik_row()` carries each trial's own log-likelihood and, for the
  seven option-code families, its saturated value of zero. What is
  missing is the SIGN: the response is nominal and has no mean to
  depart from. `rlddm()`'s row has no saturated value at all, because a
  density's supremum at a fixed response time is unbounded, and its
  refusal says so.

## Two new families

* `rlddm()`, a delta learning rule whose value difference drives the
  drift rate of a Wiener first-passage density, so choices and response
  times are one likelihood (Pedersen, Frank and Biele 2017). The
  response is the response TIME and the boundary reached travels in
  `dec()`. **hBayesDM does not carry this model**: it ships the two
  halves separately as `bandit2arm_delta` and `choiceRT_ddm`, and
  `?rlddm` maps the parameters onto the second of those rather than
  naming a counterpart that does not exist.

  The density comes from `frmtmb.eam::wiener_lpdf()`, newly exported
  for the purpose, rather than through a colon. This is the only
  dependency one extension of frmtmb has on another and it is one
  function.

  Recovered at 30 subjects by 100 trials, 60 replicates: every bias
  inside its Monte Carlo error, coverage 0.95 to 0.98 on `drift`, `bs`,
  `ndt` and `bias` and 0.87 on `alpha`, 60 of 60 usable. The pair that
  trades off is not the obvious one: `drift` with `bs` correlate at
  0.17 across replicates, while `bs` with `ndt` reach -0.57 and `drift`
  with `bias` +0.56.

  `simulate()` refuses, as `ts_par7()` does: one trial's draw is two
  numbers and a response vector holds one. `frm_task_simulate()`
  returns whole data frames with both columns.

* `igt_orl()`, outcome-representation learning for the Iowa gambling
  task (Haines, Vassileva and Ahn 2018). Three value stores per deck,
  separate learning rates for gains and losses that SWAP between the
  played deck and the three that were not, and a perseverance term.

  0.1.0 left it out on the grounds that it is weakly identified.
  Measured, that is half right. Every parameter recovers at 30 subjects
  by 100 trials: biases at or inside their Monte Carlo error, coverage
  0.92 to 0.98, 60 of 60 usable. What is weak is the SEPARATION of the
  three choice-rule parameters, visible in the correlation of the
  estimates rather than in their bias: `betaF` with `betaP` at -0.77,
  `k` with `betaP` at -0.53, `k` with `betaF` at +0.45. Report `betaF`
  and `betaP` together rather than one alone.

## The Kalman filter's exploration bonus

* `bandit4arm2_kalman_filter(bonus = TRUE)` adds Daw and others'
  exploration bonus, `tau * (mu + phi * sqrt(s))`, so an arm the
  subject is uncertain about is worth more than its posterior mean.
  `phi` is an ordinary distributional parameter and takes a formula
  like any other. `bonus = FALSE` is the default and fits exactly the
  model 0.1.0 fitted; the suite pins that by holding `phi` at zero with
  `bf(phi = 0)` and requiring the two log-likelihoods to agree.

## Checks

* Three new rows in the Stan identity tier, all written from the
  published equations: the Kalman filter with the bonus at -2.8e-13,
  `igt_orl` at -6.8e-13, and `rlddm` at 3.2e-09. **The last one is not
  exact and that makes it the strongest row.** Every other program in
  that file evaluates the same arithmetic as the engine in a different
  order; Stan implements the Wiener density itself, so this row
  compares two independent implementations of the density as well as
  two of the recursion. Which side the difference is on is settled by a
  third implementation rather than by the size of the residual:
  measured over 4000 rows spanning the region the fixture visits,
  `frmtmb.eam`'s density differs from `RWiener`'s by a median of 0 and
  a 99th percentile of 3.6e-15 per row, so a per-row difference of
  about 4e-12 against Stan sits on Stan's side.

* Longhand references for `igt_orl()` and `rlddm()` join the six that
  0.1.0 shipped, in `test-reference.R`: one subject and one trial at a
  time, `if` and `[` for selection, the opposite spelling to the
  engine's.

* `test-factorization.R` is new: the three likelihood slots agree with
  each other and with `logLik()`, a subject's value is the sum of its
  own rows in the order the core fixes, the slots honor the stacked
  design the correction calls them with, and a replicate whose
  parameters differ gives different values, without which the stacking
  check would pass on a family that ignored the stacking entirely.

# frmtmb.learn 0.1.0

First release. The value-learning models of the reinforcement-learning
and computational-psychiatry literature, written as frmtmb families.

## What the grammar buys

* Every parameter of every family is an ordinary distributional
  parameter with its own linear predictor. A learning rate takes a
  condition effect, a smooth term or a correlated per-subject random
  effect the way a mean does, and the effect comes back with a standard
  error. A reversal task fitted as
  `alpha ~ after_reversal + (1 | id)` gives the CHANGE in learning rate
  with an interval, which two separately fitted models do not. There is
  no separate reversal family here for that reason.

## One recursion, six families

* `bandit2arm_delta()`, the Rescorla-Wagner delta rule on two arms.
* `bandit2arm_dual()`, separate rates for gains and losses. `split`
  chooses whether the split is on the sign of the PREDICTION ERROR (the
  default) or of the OUTCOME, which is hBayesDM's `prl_rp`. Measured:
  on binary payoffs the two are the SAME model, because the value store
  stays inside the payoff range and the two signs agree on every trial;
  they differ only on graded payoffs.
* `prl_fictitious()`, counterfactual updating, where the option that was
  not chosen moves in the opposite direction.
* `bandit4arm2_kalman_filter()`, the restless four-armed bandit of Daw
  and others (2006). Its learning rate is not a parameter: the Kalman
  gain is derived from a posterior variance the filter carries, so it
  falls as a subject learns.
* `ts_par7()`, the two-stage model-based and model-free hybrid of Daw
  and others (2011), with `w`, `lambda`, two learning rates, two
  sensitivities and perseveration.
* `igt_pvl_delta()`, prospect-valence learning for the Iowa gambling
  task. ORL is not here; `?igt_pvl_delta` says why.

* The engine is internal and is the point of the package. It walks
  trials once and updates every subject at each step, taking the
  learning rule and the choice rule as tape-safe functions, with padding
  and a mask for unequal trial counts. A new family is the model and
  nothing else. It runs on the tape and off it, so the likelihood, the
  fitted trajectory and the simulator are ONE recursion rather than
  three copies that have to be kept in step.

## Reading a fit

* `frm_value_trace()` returns the per-trial value estimates each choice
  was made on, the prediction error the outcome produced, and the fitted
  probability of the choice that was made. These are the quantities
  papers plot.
* `frm_task_design()` builds a trial-level design for each task with the
  payoff schedule of every option fixed in advance, which is what makes
  a draw coherent: a simulated subject that takes the arm the real one
  did not still has to be paid.
* `frm_task_simulate()` draws whole datasets from the generative
  process, taking parameters directly rather than through a formula.
  `frm_simulate()` is the other route and goes through the fitted
  grammar; the two share the recursion and are checked against each
  other.
* `frm_learn_families()` is the reference table, including the map from
  each parameter to hBayesDM's spelling of it.

## What these families refuse, and why

* **`fitted()`, `predict(type = "response")` and every kind of residual
  refuse.** The response is the option a subject took, coded 1 to K, and
  it is nominal: arm 2 is not twice arm 1 and the Iowa gambling task's
  four decks have no order at all. Core forms a fitted value as the
  conditional mean of the response and a residual as `y - mean`, so a
  mean declared here would make three methods return arithmetic on a
  category code. `frm_value_trace()` returns everything a mean would
  have carried and more. Note that frmtmb's own `rw_delta` example DOES
  have `fitted()`: it codes a two-armed choice 0 and 1 and returns
  P(arm 1), which does not survive a fourth arm.
* **`frm(importance =)` is refused by name for every family**, and the
  refusal is wider than the mathematics. The correction needs one
  log-likelihood value per group and these likelihoods have one, but
  `frmtmb_structure(loglik =)` returns a single AD scalar and there is
  nowhere to put the factors. `frm_compat("bandit2arm_delta",
  "importance")` names the seam at the console. The Laplace error is
  measured by simulation instead.
* `simulate()` is refused for `ts_par7()` alone: one two-step trial's
  draw is three numbers and a response vector holds one. Use
  `frm_task_simulate()`.
* `weights()`, `cens()`, `trunc()` and `se()` refuse in the family's own
  words: each reshapes a per-row likelihood contribution, and a trial's
  contribution here is conditional on every earlier trial of the same
  subject.

## Validation

* Every family is checked by an IDENTITY against an independent Stan
  program of the same model, written from the published equations,
  evaluated at frmtmb's own estimates. The programs declare exactly the
  parameters frmtmb estimates on the scales frmtmb estimates them on, so
  the map is the identity and the comparison carries no Jacobian. The
  tier is gated on `FRMTMB_BRMS_FIT_TESTS` and `NOT_CRAN`, like the brms
  comparison tier, because it compiles Stan programs.
* Parameter recovery at a realistic scale, with bias and coverage.
* The Laplace caveat is measured by simulation at few trials per
  subject, and the numbers are in `vignette("learning")` and in
  `?frmtmb.learn`.
* No hBayesDM code or data is used anywhere. Its models are published
  equations and the families here are written from them. A cross-check
  against it lives in `dev/`, is never run by the suite, and has not
  been run at all: hBayesDM 2.0.0 fits through cmdstanr, and CmdStan's
  vendored TBB does not build under this machine's compiler. The
  script's header records the whole chain. So this release has no check
  against an outside implementation; what it has instead is an exact
  identity against its own independent Stan programs and a longhand
  reference per family.

## Deliberate omissions

* RLDDM, a learning rule feeding a drift-diffusion choice rule, is the
  next family and is not built: the choice rule belongs to
  `frmtmb.eam`, and building it here would make this package import a
  sibling extension. The engine already has the shape it needs.
* No exploration bonus on the Kalman filter, and no ORL for the Iowa
  gambling task. Both are short additions on this engine and both are
  left out for stated reasons in their families' help.
