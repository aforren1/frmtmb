#' @keywords internal
#'
#' @section What this package is for:
#' A learning model is usually fitted one subject at a time, or
#' hierarchically with a Bayesian sampler and a fixed prior. Written as
#' frmtmb families, the same models get the ordinary formula grammar:
#' every parameter is a distributional parameter with its own linear
#' predictor, so a learning rate takes a condition effect, a smooth term
#' or a correlated per-subject random effect the way a mean does, and
#' the effect comes back with a standard error. A reversal task fitted
#' as `alpha ~ after_reversal + (1 | id)` gives the CHANGE in learning
#' rate with an interval, which two separately fitted models do not.
#'
#' @section One recursion, eight families:
#' Every family here is the same walk: carry a value store per subject,
#' visit the trials in order, turn the store into a choice probability
#' and the outcome into a new store. The walk is written once. A family
#' supplies a choice rule and a learning rule, and is a few dozen lines.
#'
#' Seven of the eight turn the store into a SOFTMAX over options.
#' [rlddm()] turns it into the drift rate of a two-boundary diffusion,
#' so its trial contributes the joint density of which boundary was
#' reached and when. That is the one place the engine had to grow a
#' seam: a family may supply a log density in place of a choice rule,
#' and everything else about the walk, the mask and the three likelihood
#' slots is unchanged.
#'
#' The walk loops over TRIALS and vectorizes over SUBJECTS, which is not
#' a style preference. Written one row at a time with the store updated
#' by sub-assignment, RTMB pays for it at tape construction and the
#' penalty grows with the row count. The same code runs on the tape and
#' off it, so the likelihood, [frm_value_trace()] and the simulator are
#' one recursion rather than three copies.
#'
#' @section The Laplace caveat, and where the answers are safe:
#' frmtmb integrates the subject effects out with a Laplace
#' approximation, exact only when the conditional log-density is
#' quadratic. For binary choices it is not, and the fewer trials a
#' subject has the less quadratic it is. Measured by simulation on
#' [bandit2arm_delta()] at 40 subjects (`dev/learn-findings.md`, and the
#' Laplace section of `vignette("learning")`):
#'
#' | parameter | bias at 100 trials | coverage | bias at 20 trials | coverage |
#' |---|---|---|---|---|
#' | `alpha_(Intercept)` | +0.009 | 0.97 | -0.026 | 0.92 |
#' | `alpha_after` | -0.015 | 0.94 | +0.009 | 0.91 |
#' | `tau_(Intercept)` | +0.003 | 0.95 | +0.016 | 0.96 |
#' | `log sd(alpha)` | -0.266 | 0.99 | -1.596 | 0.95 |
#'
#' * the FIXED effects survive short sessions. At a fifth of the trials
#'   their bias is still inside Monte Carlo error and their Wald
#'   intervals still cover near the nominal rate.
#' * the variance component's POINT ESTIMATE does not. At 20 trials
#'   `log sd(alpha)` comes out 1.6 too low, a standard deviation about a
#'   fifth of the true one.
#' * its INTERVAL does survive, which is the part worth knowing. The
#'   spread of the estimates rises from 0.88 to 3.01 as the bias grows,
#'   so the interval widens at least as fast as the point estimate
#'   degrades and still covers 0.95. A reader of the point estimate is
#'   misled; a reader of the interval is not. Treat a subject-level
#'   standard deviation from twenty binary trials as a lower bound.
#'
#' Do not read the bias as Laplace error alone: a variance component
#' estimated by maximum likelihood from binary data with few levels is
#' biased downward whether or not the integral is approximated, and the
#' simulation does not separate the two causes. Separating them is what
#' `frm(importance =)` exists for, and it now runs; see the next
#' section for what it is worth.
#'
#' @section The seam that closed, and what the correction is worth:
#' `frm(importance =)` was refused for every family in the first release
#' and works now. The correction reweights draws from the Laplace
#' Gaussian and needs one log-likelihood value per GROUP; these
#' likelihoods factorize over subjects and again over trials, so the
#' values always existed, and what was missing was a slot to put them
#' in. `frmtmb_structure()` grew two, `loglik_row` and `loglik_group`,
#' and every family here declares both off the same recursion the
#' objective tapes. A subject is the group, a trial is the row, and
#' `unit` is unchanged at "one subject's trial sequence", because
#' dropping a trial changes every later trial's value store.
#'
#' What the correction is worth depends on the design rather than on the
#' family. Measured on [bandit2arm_delta()], 40 subjects, one scalar
#' random intercept, six datasets per cell, truth `sd(alpha) = 0.5`
#' (`dev/learn2-findings.md` has the per-dataset table):
#'
#' | trials | datasets | Laplace `sd(alpha)` | corrected | mcse | min ESS per draw |
#' |---|---|---|---|---|---|
#' | 100 | 3 of 6 informative | 0.38 to 0.45 | 0.45 to 0.51 | 0.05 to 0.06 | 0.87 to 0.93 |
#' | 100 | 1 of 6 collapsed | 0.0004 | 0.0025 | 0.000 | 1.000 |
#' | 100 | 2 of 6 refused | | | | |
#' | 20 | 5 of 6 collapsed | 0.0001 | unchanged | 0.000 | 1.000 |
#' | 20 | 1 of 6 informative | 0.371 | 0.729 | 0.103 | 0.743 |
#'
#' A COLLAPSED DATASET IS NOT A SHORT-SESSION PROBLEM. One of the six
#' at 100 trials has it too, and reports the same spurious +1.822 as
#' the five at 20 trials do. Check `VarCorr(fit)` on the Laplace
#' fit at every trial count, not only short ones.
#'
#' * at 100 trials, on the three of six datasets whose variance
#'   component is identified, the correction moves it UP by 0.11 to
#'   0.18 log units and toward the truth, its Monte Carlo error is a
#'   third of the shift, and the effective sample sizes are near 1. It
#'   costs 30 to 80 seconds against 2 for the Laplace fit. Of the other
#'   three, two refused on the correction's own convergence guard and
#'   one had already collapsed.
#' * at 20 trials it is not usable, and the reason is worth knowing.
#'   Five of six Laplace fits have already collapsed `sd(alpha)` to
#'   1e-4, so there is nothing left to reweight; the correction then
#'   walks at its own step cap for every round and WARNS that it did
#'   not converge. A shift reported with `fit$importance$capped` true
#'   and every entry of `fit$importance$moves` equal is that step cap,
#'   not an estimate. The one 20-trial dataset with a real variance
#'   component overshoots, 0.37 to 0.73 against a truth of 0.5.
#'
#' DO NOT FOLLOW THE WARNING'S ADVICE ON A CAPPED RUN. It says to raise
#' `frmtmb_control(importance_rounds =)`, and on a collapsed variance
#' component that makes the artifact bigger in exact proportion:
#' measured on one 20-trial dataset, five rounds of a 0.3645 cap give a
#' shift of 1.822 and ten rounds give 3.645, both exactly rounds times
#' the cap. The number is not an estimate that more iterations would
#' sharpen. Read `fit$importance$capped` and `$moves` first, and if the
#' moves are all equal, treat the correction as declining rather than
#' answering.
#'
#' So check `VarCorr(fit)` on the Laplace fit before believing a
#' correction, and read the warning if there is one.
#'
#' The correction requires the model to be grouped on the family's own
#' subject. `(1 | something_else)` gives a per-subject likelihood
#' against a per-something-else proposal, and the core refuses it by
#' name rather than adding up numbers that do not add up.
#'
#' THE SAME SLOTS ALSO TURN `loo()` ON, at subject granularity, and an
#' earlier draft of this page said the opposite. `frm()` itself is
#' maximum likelihood and has no draws to average over, so `loo()` on a
#' `frmtmb_fit` still refuses for that core-wide reason; the route that
#' changed is `frmtmb.sample`, which admits any structure declaring
#' `loglik` together with either factorization slot and then reads the
#' COARSEST one. Measured on a 6-subject, 180-row fit: the pointwise
#' matrix is admitted, it has 6 columns rather than 180, it carries
#' `attr(x, "unit")` of `"one subject's trial sequence"`, and it sums to
#' `logLik(fit)` to twelve digits.
#'
#' Read that literally. A column is a SUBJECT, so leaving one out drops
#' that subject's whole sequence, which is the only honest leave-one-out
#' for a recursion: a trial cannot be dropped because every later
#' trial's value store depends on it. `loo()` prints "Computed from N by
#' K" and says nothing about what a column is, so `frmtmb.sample`
#' messages the unit when the matrix carries one.
#'
#' Deviance residuals are
#' closer but still refused: the magnitude is now available from
#' `loglik_row()`, and what is missing is the SIGN, because seven of the
#' eight families have a nominal response with no mean to depart from
#' and the eighth has no exported route to the mean of a first-passage
#' time.
#'
#' @section Recovery tables for the three families measured this round:
#' `dev/learn-recovery-round2.R`, 60 replicates each, 30 subjects by 100
#' trials, one random intercept on the primary parameter. Everything on
#' the NATURAL scale, with each Wald interval's endpoints pushed through
#' the same monotone link, because [rlddm()]'s non-decision time is
#' bounded by its own dataset's fastest response and its link scale is
#' therefore not comparable across replicates.
#'
#' **[rlddm()].** Every bias inside its Monte Carlo error, 60 of 60
#' usable.
#'
#' | parameter | truth | bias | mc se | sd of est | coverage |
#' |---|---|---|---|---|---|
#' | `alpha` | 0.35 | -0.002 | 0.003 | 0.026 | 0.87 |
#' | `drift` | 3.00 | +0.001 | 0.009 | 0.069 | 0.98 |
#' | `bs` | 1.60 | 0.000 | 0.002 | 0.019 | 0.95 |
#' | `ndt` | 0.20 | +0.001 | 0.000 | 0.003 | 0.98 |
#' | `bias` | 0.50 | -0.002 | 0.001 | 0.007 | 0.95 |
#'
#' The pair that co-varies is `bs` with `ndt` at -0.57 and `drift` with
#' `bias` at +0.56. `drift` with `bs`, which is the pair a first-passage
#' density suggests, is 0.17.
#'
#' **[igt_orl()].** Every parameter recovers; coverage 0.92 to 0.98, 60
#' of 60 usable.
#'
#' | parameter | truth | bias | mc se | sd of est | coverage |
#' |---|---|---|---|---|---|
#' | `Arew` | 0.3 | -0.008 | 0.004 | 0.032 | 0.97 |
#' | `Apun` | 0.1 | 0.000 | 0.001 | 0.007 | 0.98 |
#' | `k` | 0.5 | -0.004 | 0.014 | 0.107 | 0.92 |
#' | `betaF` | 1.0 | -0.015 | 0.013 | 0.097 | 0.95 |
#' | `betaP` | 1.0 | +0.012 | 0.010 | 0.075 | 0.97 |
#'
#' What is weak is the SEPARATION of the choice-rule parameters rather
#' than any one estimate: `betaF` with `betaP` correlate at -0.77, `k`
#' with `betaP` at -0.53 and `k` with `betaF` at +0.45. The two learning
#' rates correlate at 0.38 and with none of the three above 0.21.
#'
#' **[bandit4arm2_kalman_filter()] with `bonus = TRUE`**, and with
#' `center`, `mu0` and `sigma0` held at the task's own values through
#' `bf(name = value)`, because those three do not separate on a session
#' this size whether or not the bonus is on.
#'
#' | parameter | truth | bias | mc se | sd of est | coverage |
#' |---|---|---|---|---|---|
#' | `tau` | 0.15 | +0.005 | 0.001 | 0.011 | 0.85 |
#' | `lambda` | 0.98 | 0.000 | 0.000 | 0.003 | 1.00 |
#' | `sigmaD` | 3.00 | -0.045 | 0.026 | 0.203 | 0.95 |
#' | `phi` | 1.50 | +0.059 | 0.018 | 0.141 | 0.93 |
#'
#' `phi` trades off with `sigmaD`, at -0.81, and not with `tau`, at
#' 0.12. Both scale the same term: `phi` multiplies the posterior
#' standard deviation and `sigmaD` sets how fast it grows.
#'
#' @section A correlated block over every parameter, at 100 by 200:
#' The tables above put ONE random intercept on the primary parameter.
#' This one puts `(1 | p | id)` on every parameter at once, at the scale
#' the field's designs sit at, and asks whether the CORRELATIONS come
#' back as well as the variances. `dev/learnhier-findings.md` has the
#' construction, the seeds and the scripts; the summary is here.
#'
#' The design draws a block with real off-diagonal structure, because a
#' design that draws independent deviations cannot say whether an
#' estimator recovers a correlation, only whether it invents one. Both
#' questions are answered by running the block twice.
#'
#' **[bandit2arm_delta()]**, `(1 | p | id)` on the learning rate and the
#' inverse temperature, 100 learners by 200 trials, 60 replicates per
#' arm. Everything on the link the family estimates it on.
#'
#' | component | truth | mean estimate | bias | mc se | coverage |
#' |---|---|---|---|---|---|
#' | `alpha_(Intercept)` | -0.619 | -0.6191 | -0.0001 | 0.0091 | 0.92 |
#' | `tau_(Intercept)` | 1.099 | 1.1019 | +0.0033 | 0.0041 | 1.00 |
#' | `sd(alpha)` | 0.5 | 0.4925 | -0.0075 | 0.0088 | 0.92 |
#' | `sd(tau)` | 0.3 | 0.2980 | -0.0020 | 0.0039 | 0.93 |
#' | `cor(alpha, tau)` | 0.5 | 0.5322 | +0.0322 | 0.0181 | 0.97 |
#'
#' The same design with a correlation of ZERO returns 0.0214 with a
#' Monte Carlo error of 0.0192 and coverage 0.93, so the estimator does
#' not manufacture a correlation that is not there. Both biases on the
#' correlation are inside twice their own Monte Carlo error.
#'
#' 60 of 60 fits in each arm reached convergence code 0 with a positive
#' definite Hessian and no `NaN` standard errors, and no component
#' collapsed: over the 120 fits the smallest `sd(alpha)` was 0.328 and
#' the smallest `sd(tau)` 0.194. Collapse is documented for this family
#' at 20 to 60 trials and is absent at 200.
#'
#' **READ THE INTERVAL, NOT THE POINT ESTIMATE, FOR A CORRELATION.** The
#' spread of the estimates across replicates is 0.140 in both arms and
#' the observed range is 0.172 to 0.804 against a truth of 0.5. One
#' dataset of this size locates a correlation to about half a unit, and
#' the Wald interval says so where the point estimate does not.
#'
#' At 60 replicates a coverage carries a standard error of 2.8 points at
#' the nominal 95, so nothing between about 89 and 100 percent in these
#' tables is distinguishable from nominal.
#'
#' **[rlddm()]** under `ndt_group(id)`, `(1 | p | id)` on `alpha`,
#' `drift`, `bs` and `ndt`, `bias` held at 0.5, 60 replicates. 60 of 60
#' fitted, none with a log-likelihood below the objective at its own
#' realized truths, no `NaN` standard errors.
#'
#' | component | truth | mean estimate | bias | mc se | coverage |
#' |---|---|---|---|---|---|
#' | `alpha_(Intercept)` | -0.619 | -0.6224 | -0.0033 | 0.0087 | 0.97 |
#' | `drift_(Intercept)` | 2.500 | 2.5075 | +0.0075 | 0.0142 | 0.95 |
#' | `bs_(Intercept)` | 0.4055 | 0.3993 | -0.0062 | 0.0023 | 0.98 |
#' | `sd(alpha)` | 0.5 | 0.4746 | -0.0254 | 0.0071 | 0.90 |
#' | `sd(drift)` | 1.0 | 0.9981 | -0.0019 | 0.0098 | 0.95 |
#' | `sd(bs)` | 0.2 | 0.1986 | -0.0014 | 0.0017 | 0.97 |
#' | `cor(alpha, drift)` | 0.4 | 0.3803 | -0.0197 | 0.0164 | 0.95 |
#' | `cor(alpha, bs)` | 0.0 | 0.0127 | +0.0127 | 0.0150 | 0.95 |
#' | `cor(drift, bs)` | -0.3 | -0.2836 | +0.0164 | 0.0130 | 0.93 |
#'
#' Every parameter with a population truth on the scale it is fitted on
#' recovers. The four that have none are the `ndt` row and column, and
#' two of those do NOT recover against their realized targets:
#' `sd(ndt)` at a rate of 0.28 and `cor(bs, ndt)` at 0.08, where every
#' other component on either target sits between 0.93 and 1.00.
#'
#' **That is a property of the parameterization and not of the
#' estimator**, and `?rlddm` has the one line of algebra: under
#' `ndt_group()` the fitted deviation is `log(ndt) - log(m)` where `m`
#' is an observed minimum, so about 94 percent of the link-scale
#' component's variance is a property of the data. Read the per-learner
#' non-decision times through [frmtmb.eam::ndt_time()] instead; they
#' never pass through that term.
#'
#' The same fits are checked against an independent Stan program of the
#' same model with the correlated block added, at frmtmb's own
#' estimates: agreement to 2.6e-15 relative for the delta learner's
#' two-by-two block, 4.9e-15 at the full 100 by 200 design, and 1.7e-11
#' for `rlddm()`'s four-by-four block, which is looser because Stan
#' implements the Wiener density itself and that row therefore checks
#' two independent implementations of the density as well as of the
#' recursion.
#'
#' @section What this package reads that frmtmb does not promise:
#' Nothing. Every accessor it uses is exported and documented:
#' `frmtmb_family()`, `frmtmb_structure()`, `frmtmb_register_aterm()`,
#' `frmtmb_register_compat()`, `compat_rule_builder()`,
#' `single_response()`, `eval_dpars()` and `frame_block_of()`. There is
#' one thing it wanted and could not have, and it is recorded rather
#' than worked around: `frm_compat_features()` shows that an addition
#' term is registered but not at what ARITY, so this package cannot
#' verify that a `reward` term another package registered is the
#' two-column one its families need. It declines to register over the
#' top of an existing term, and the mismatch surfaces one step later as
#' `frm()` refusing a missing `reward2`.
#'
#' @section Not built, and named:
#' The three things the first release left out are in: [rlddm()],
#' [igt_orl()] and the Kalman filter's exploration bonus. What is still
#' out:
#'
#' * **`fitted()` for [rlddm()].** Unlike the seven softmax families,
#'   this one's response is ordered and its conditional mean exists. The
#'   mean of a Wiener first-passage time conditional on the boundary
#'   reached belongs to `frmtmb.eam`, which exports the density and not
#'   the mean, so declaring `fitted_mean` here would mean deriving it
#'   again. That is a second seam and it is recorded rather than guessed
#'   at. `frm_value_trace()` gives the per-trial density and drift rate
#'   meanwhile.
#' * **A mixture over learning strategies.** A real model and a wanted
#'   one. Core's `mixture()` combines per-ROW densities and would need
#'   per-sequence ones; the families now declare their per-sequence
#'   values, so what is left is a change under core's `R/` rather than a
#'   declaration this package can make.
#' * **`accepts_aterms` on these families.** frmtmb gained an
#'   addition-term allow-list, and these families do not declare one, so
#'   a term none of them reads is still carried without effect. Every
#'   term that would reshape a per-row contribution is already refused
#'   by name with a better message than an allow-list would give, so
#'   what the declaration would add is narrow.
#' * **A recovery table for every family.** Four of the eight have one
#'   now ([bandit2arm_delta()], [rlddm()], [igt_orl()] and the Kalman
#'   filter with `bonus = TRUE`). The other four have the Stan identity,
#'   the longhand reference and a smoke-level recovery assertion, which
#'   establishes that they compute the model they claim; what they do
#'   not have is a published bias-and-coverage table. The reason is time
#'   rather than principle: each table is 60 to 200 fits.
#' * **[frm_task_design()] for a task the user brings.** The design
#'   helpers cover the tasks these families are written for. A user
#'   fitting their own data needs none of it: the family reads whatever
#'   columns the formula names.
#'
#' @examples
#' # a reversal task, and the change in learning rate with an interval
#' d <- frm_task_design("reversal", n_subject = 10, n_trial = 40, seed = 5)
#' d$choice <- frm_task_simulate(
#'   bandit2arm_delta(subject = id, trial = trial), d,
#'   pars = list(alpha = 0.35, tau = 3), seed = 5)[[1]]$choice
#' fit <- frmtmb::frm(
#'   frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal, tau ~ 1),
#'   family = bandit2arm_delta(subject = id, trial = trial), data = d)
#' frmtmb::fixef(fit)
"_PACKAGE"

# frmtmb is a Depends, so that one library(frmtmb.learn) call gives a
# user the formula grammar, the frm() a model is fitted with, and the
# family constructors. It matches every other family extension in the
# repository. The seams this package builds on are imported by name as
# well, because a namespace that is loaded and not attached reaches
# nothing through the search path.
#' @importFrom frmtmb frmtmb_family frmtmb_structure frmtmb_register_compat
#'   compat_rule_builder single_response eval_dpars frame_block_of
#' @importFrom stats setNames
NULL

# WHY frmtmb.eam IS IMPORTED BY NAME AND NOT ONLY THROUGH `::`.
#
# `rlddm()` requires the `dec()` addition term, and `dec()` reaches
# frmtmb's registry from frmtmb.eam's own `.onLoad()`. A package named
# in DESCRIPTION `Imports:` alone is only required to be INSTALLED; its
# namespace is not loaded when this package is attached, and
# `frmtmb.eam::wiener_lpdf()` loads it when the function is CALLED,
# which is long after `frm()` has parsed the formula and refused the
# term. So a clean session gave "Addition term `dec()` is not
# supported" on a user's own data, and every test and example in this
# package masked it by calling frm_task_simulate() first, which reaches
# frmtmb.eam::ddm_simulate() and loads the namespace as a side effect.
#
# An `importFrom` is what actually loads it, and R loads a package's
# imports before running its `.onLoad()`, so `dec()` is in the registry
# before this package's own compat rules name it.
# `test-clean-session.R` runs a fresh R process to keep this true.
# The same reasoning covers `ndt_group()`, which frmtmb.eam registers in
# the same `.onLoad()` and which `rlddm()` now reads through frmtmb.eam's
# own bound seam rather than through a second copy of the link.
#' @importFrom frmtmb.eam wiener_lpdf ddm_simulate ndt_bound
#'   ndt_bound_attach ndt_bound_of ndt_bound_pending ndt_apply
NULL
