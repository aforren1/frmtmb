# frmtmb.learn (development version)

* **`rlddm()` can bound the non-decision time PER GROUP.** Write
  `rt | dec(choice) + reward(pay1, pay2) + ndt_group(id) ~ ...` and
  each row's `ndt` is bounded by its own learner's fastest response
  instead of by the fastest response of every learner in the data set.
  That is what a random effect on `ndt` needs. On the scale tier's own
  design, 100 learners by 200 trials, the same data and the same model:

  | | 0.3.0, one global bound | with `ndt_group(id)` |
  |---|---|---|
  | convergence code | 1 | 0 |
  | maximum absolute gradient | 6.36e+09 | 3.12e-03 |
  | positive definite Hessian | no | yes |
  | standard errors that are `NaN` | 4 | 0 |
  | log-likelihood | -8333.1 | -7781.2 |
  | `sd(ndt)` on the link | 2.236 | 0.3261 |
  | learners below their own fastest response | (n/a) | 100 of 100 |
  | per-learner `ndt` error against the drawn truths | (n/a) | 14.0 ms |

  551.9 log-likelihood units at the same fourteen parameters. The `ndt`
  variance component is the clearest single symptom: under one bound it
  is a random effect running away on a scaled logit against a wall.

  Under the grouping `ndt` is a FRACTION of the row's own bound and
  `predict(dpar = "ndt", type = "response")` reports that fraction;
  `frmtmb.eam::ndt_time()` reports the non-decision time in the
  response's own units under either parameterization. `max_ndt` and
  `ndt_group()` together are refused, because they set one bound to two
  different things. A model that does not write `ndt_group()` is
  unchanged: the log-likelihood, every coefficient, every fitted value
  and the value trace are identical in every digit to 0.3.0's.

* The bound, its refusals and the `ndt_group()` coercion now come from
  `frmtmb.eam::ndt_bound()` rather than from a second copy of the
  scaled logit written here. The copy is what made this family inherit
  the one-global-bound defect in the first place.

* **BREAKING: `bf(ndt = )` on a bare `rlddm()` is refused and needs
  `rlddm(max_ndt = )`.** `bf()` transforms a pinned constant at PARSE
  time, before the response has been seen, so a family with no bound
  yet has no scale to transform it on. Pass the bound instead:
  `rlddm(subject = id, trial = trial, max_ndt = 0.26)`. This is the
  contract `frmtmb.eam::wiener()` and its siblings have had since
  frmtmb.eam 0.7.0, and `rlddm()` now shares it rather than being the
  exception.

  **The capability is respelled, not lost.** The same model under the
  new spelling is the same fit: `bf(ndt = 0.2)` with
  `rlddm(max_ndt = 0.26)` gives a log-likelihood of -86.3374111 on this
  version and -86.3374111 on 0.3.0, which is also what
  `bf(ndt = 0.2)` on a bare 0.3.0 family gave. Add the argument and the
  numbers do not move.

  **Nothing that worked was wrong, and no fitted model needs redoing.**
  Through 0.3.0 an in-range `bf(ndt = 0.2)` was fitted at exactly 0.2:
  the parse-time call only range-checks the constant, and the transform
  that reaches the parameter runs later, against the settled link. What
  the old behavior did NOT do is check the constant against the bound
  the fit would actually use. `bf(ndt = 0.3)` with `max_ndt = 0.26` was
  accepted, reached the objective as `NaN`, and surfaced as
  `NA/NaN gradient evaluation` from the optimizer, naming neither the
  parameter nor the constant. That case is now refused before the
  formula is parsed, by name: `Constant ndt = 0.3 is not in the range
  of the scaled_logit link`.

* `summary()` on an ungrouped `rlddm()` fit prints the `ndt` link as
  `scaled_logit` where 0.3.0 printed `scaled_logit(0, 0.2615)`. The
  link's arithmetic is unchanged and the bound is now on the fitted
  family at `family(fit)$ndt_bound$ub`, which 0.3.0 did not carry.

* THIS PACKAGE'S FLOOR ON frmtmb.eam HAS TO RISE. `rlddm()` now calls
  `ndt_bound()`, `ndt_bound_attach()`, `ndt_bound_of()`,
  `ndt_bound_pending()` and `ndt_apply()`, which frmtmb.eam gains in
  the same round. `DESCRIPTION` still says `frmtmb.eam (>= 0.6.0)` and
  is not edited here, because the version those exports ship in is not
  this lane's to choose. `library(frmtmb.learn)` against frmtmb.eam
  0.6.0 fails at namespace load, so the floor must rise in the same
  commit that lands the exports.

# frmtmb.learn 0.3.0

A duplicated reward schedule fitted correctly and then simulated a
different experiment. `reward(pay1, pay2)` and `reward(rec, rec)` give
a bitwise identical log-likelihood and identical `fixef()`, because the
density reads only the chosen arm, while the simulator drew a P(better
arm) of 0.5470 against 0.8616 with 0.8633 observed. Exact chance, from
a task nobody ran, with no error anywhere.

* The guard is DERIVED rather than declared: the schedule columns are
  read off the addition terms a family already names, so all eight
  families are covered and one that gains a `reward()` term later is
  covered without anyone remembering. What that costs is written down:
  a schedule under a third name is invisible, so a test asserts every
  multi-column term is classified or excluded with a reason.

* A compatibility row said `prl_fictitious()` reads both columns. It
  does not: its update takes the chosen payoff and flips the sign.
  That row has been false since 0.2.0.

* `rlddm()` inherits frmtmb.eam's non-decision-time bound and the
  defect that comes with it. See that package's 0.6.0 notes; do not
  put a random effect on `ndt` yet.

* **A draw needs a column a fit does not, and `simulate()` now refuses
  when the data does not carry it.** ALL EIGHT families read only the
  CHOSEN option's payoff. That is what lets a record of the received
  outcome alone be fitted, by passing the one column once per option,
  and it is exact rather than approximate: each option's prediction
  error is multiplied by a 0/1 chosen indicator, and zero times a finite
  number is exactly zero, so the objective is bitwise the same function
  of the parameters whatever the unchosen entry holds. Measured by
  replacing the unchosen entries with `N(100, 50)` noise and comparing
  the objective at one parameter vector: bitwise identical.

  A DRAW is a different question, and the difference had been silent. A
  simulated subject chooses for itself, so paying it needs the schedule
  of every option, and with the columns duplicated there is no
  schedule. Measured on the shipped two-armed design, 30 subjects by
  100 trials, where arm 1 pays with probability 0.7 and arm 2 with 0.3:

  | family, route | real | duplicated |
  | --- | --- | --- |
  | bandit2arm_delta, simulate | 0.7444 (0.0024) | 0.5004 (0.0035) |
  | prl_fictitious, simulate | 0.8616 (0.0027) | 0.5470 (0.0090) |
  | prl_fictitious, task_simulate | 0.8367 (0.0037) | 0.5072 (0.0083) |
  | rlddm, task_simulate | 0.1443 (0.0036) | 0.5040 (0.0071) |

  The statistic is the proportion of trials in the second half that took
  arm 1, over 40 draws (20 for `rlddm()`); the observed proportions in
  the data being fitted are 0.7593 and 0.8633. Every duplicated column
  gives exact chance, because with both arms paying the same on every
  trial there is nothing to learn, and each came back as a
  plausible-looking vector of option codes with no error and no
  warning. `ts_par7()` drew 800 rows the same way.

  The signature is EVERY column of `reward()` or `payoff()` identical on
  every row. Some columns identical is not it, because the data then
  still says what at least one option not taken would have paid, and the
  test file pins that case as accepted. `stage2()` is not read as a
  schedule at all: its two columns are the observed state and choice.

  The guard is DERIVED from the terms a family names rather than
  declared per family, so a ninth family cannot be added without the
  question being asked. `ln_family(counterfactual =)` takes `NULL` to
  derive, `identical(FALSE)` to opt out, or the columns themselves, and
  refuses anything else while the family is being built rather than
  when a draw is asked for. It does not take `TRUE`: deriving is
  already the default, so there is nothing for `TRUE` to mean, and
  reading it as an opt-in would make the one spelling that says yes the
  one that turned the guard off.

  What the derivation is weaker at than a declaration is a NAME it does
  not know. A schedule arriving through a term that is not `reward()`
  or `payoff()` is not guarded, and nothing says so; a forgotten
  declaration would at least sit in the family's source. The package's
  own registration table is the only place that is visible, so the test
  suite now asserts that every term it registers at arity 2 or more is
  classified as a payoff schedule or explicitly excluded with a reason.
  Registering a fourth term without answering the question turns the
  suite red. All four routes to a draw are covered:
  `simulate()`, `frm_simulate()` and `frmtmb.sample::posterior_predict()`
  through the family's `sim_ctx` slot, and `frm_task_simulate()` on the
  design it is handed, which is the only route `rlddm()` and
  `ts_par7()` have because `simulate()` already refuses them and their
  own refusal messages send users there.

  `newdata` is not a way round the refusal and the message says so: the
  formula names one column twice and `newdata` is read through that
  same formula, so a user whose record holds the received outcome alone
  has to refit with a column per option, or build a schedule with
  `frm_task_design()`.

* **A released compatibility row was wrong, and it is corrected.** Since
  0.1.0 the `prl_fictitious()` / `reward()` row has said "Both columns
  are READ here rather than only carried: counterfactual updating moves
  the unchosen option's value too, so the second column enters the
  likelihood and not just the simulator". It does not. The update forms
  the outcome as `c1 * reward1 + c2 * reward2` with the CHOSEN
  indicators and then FLIPS ITS SIGN, so what moves the unchosen value
  is the negative of the realized outcome and the second column never
  enters. Measured: replacing the unchosen entries with `N(100, 50)`
  noise leaves the log-likelihood bitwise unchanged at
  -937.55332137794574.

  `?prl_fictitious` has always said this correctly, that the
  counterfactual outcome is the negative of the realized one "rather
  than a reading of the second `reward()` column". Two documents in one
  package said opposite things about the same family, and a user
  choosing a family on that table would have chosen this one for a
  property it does not have. An earlier draft of the release above
  exempted `prl_fictitious()` from the new guard on the strength of the
  wrong one.

* **`reward(pay)` with one argument is still refused, and the refusal
  belongs to frmtmb rather than to this package.** An addition term's
  arity is fixed when it is registered: `frmtmb_register_aterm(arity =)`
  takes one whole number and the parser refuses any other argument
  count, so `reward(pay1)` stops at "`reward()` takes 2 arguments, not
  1" before this package sees anything. Registering `reward` at arity 1
  instead would refuse `reward(pay1, pay2)`, which every example,
  vignette and test writes.

  What the core seam buys is the SPELLING and not the capability. A
  one-column route needs no core change at all: a second term name
  registered at arity 1 keys its value at `aterms[["reward1"]]`, which
  is the key the two-column spelling's first argument already uses, so
  a family that reads `reward2` with a fallback takes both spellings.
  Verified as far as the parser: `reward1(rec)` is routed to that key
  and the shipped family then refuses on its own `required_aterms`
  declaration, which is this package's to change and not core's.
  `reward(rec, rec)` remains the spelling that works today and is
  exactly right for the likelihood.

# frmtmb.learn 0.2.1

Requires frmtmb 0.55.0. The hazard-container lint runs in this
package's own check, and `frm(importance = )` now says when the
correction has stalled rather than advising more rounds, which
this package's own documentation had to say in prose instead.

# frmtmb.learn 0.2.0

The three things 0.1.0 named as left out are in, and the seam it was
waiting on has closed. Requires frmtmb 0.53.0 for the two factorization
slots, and frmtmb.eam 0.5.0 for the Wiener density.

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
