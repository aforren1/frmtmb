# frmtmb.learn: Reinforcement-Learning Families for 'frmtmb' Models

The value-learning models of the reinforcement-learning and
computational-psychiatry literature, written as 'frmtmb' families so
that each of their parameters is an ordinary distributional parameter
with its own linear predictor. A learning rate then takes a condition
effect, a smooth term or a correlated per-subject random effect the way
a mean does, and it comes back with a standard error, which is what a
separate model per group does not give. Eight families share one
recursion that walks trials once and updates every subject at each step:
the two-armed delta learner and its dual-rate and counterfactual
variants, a Kalman filter over the restless four-armed bandit of Daw and
others (2006) with an optional exploration bonus, the two-stage
model-based and model-free hybrid of Daw and others (2011),
prospect-theory and outcome-representation learners for the Iowa
gambling task, and a joint model of choices and response times in which
the learned value difference drives the drift rate of a Wiener
diffusion. Families are named as 'hBayesDM' names them where a name
exists. Each one is checked by an identity against an independent 'Stan'
program of the same model at the same estimates, and by parameter
recovery at a realistic scale.

## What this package is for

A learning model is usually fitted one subject at a time, or
hierarchically with a Bayesian sampler and a fixed prior. Written as
frmtmb families, the same models get the ordinary formula grammar: every
parameter is a distributional parameter with its own linear predictor,
so a learning rate takes a condition effect, a smooth term or a
correlated per-subject random effect the way a mean does, and the effect
comes back with a standard error. A reversal task fitted as
`alpha ~ after_reversal + (1 | id)` gives the CHANGE in learning rate
with an interval, which two separately fitted models do not.

## One recursion, eight families

Every family here is the same walk: carry a value store per subject,
visit the trials in order, turn the store into a choice probability and
the outcome into a new store. The walk is written once. A family
supplies a choice rule and a learning rule, and is a few dozen lines.

Seven of the eight turn the store into a SOFTMAX over options.
[`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md)
turns it into the drift rate of a two-boundary diffusion, so its trial
contributes the joint density of which boundary was reached and when.
That is the one place the engine had to grow a seam: a family may supply
a log density in place of a choice rule, and everything else about the
walk, the mask and the three likelihood slots is unchanged.

The walk loops over TRIALS and vectorizes over SUBJECTS, which is not a
style preference. Written one row at a time with the store updated by
sub-assignment, RTMB pays for it at tape construction and the penalty
grows with the row count. The same code runs on the tape and off it, so
the likelihood,
[`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
and the simulator are one recursion rather than three copies.

## The Laplace caveat, and where the answers are safe

frmtmb integrates the subject effects out with a Laplace approximation,
exact only when the conditional log-density is quadratic. For binary
choices it is not, and the fewer trials a subject has the less quadratic
it is. Measured by simulation on
[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md)
at 40 subjects (`dev/learn-findings.md`, and the Laplace section of
[`vignette("learning")`](https://aforren1.github.io/frmtmb/frmtmb.learn/articles/learning.md)):

|                     |                    |          |                   |          |
|---------------------|--------------------|----------|-------------------|----------|
| parameter           | bias at 100 trials | coverage | bias at 20 trials | coverage |
| `alpha_(Intercept)` | +0.009             | 0.97     | -0.026            | 0.92     |
| `alpha_after`       | -0.015             | 0.94     | +0.009            | 0.91     |
| `tau_(Intercept)`   | +0.003             | 0.95     | +0.016            | 0.96     |
| `log sd(alpha)`     | -0.266             | 0.99     | -1.596            | 0.95     |

- the FIXED effects survive short sessions. At a fifth of the trials
  their bias is still inside Monte Carlo error and their Wald intervals
  still cover near the nominal rate.

- the variance component's POINT ESTIMATE does not. At 20 trials
  `log sd(alpha)` comes out 1.6 too low, a standard deviation about a
  fifth of the true one.

- its INTERVAL does survive, which is the part worth knowing. The spread
  of the estimates rises from 0.88 to 3.01 as the bias grows, so the
  interval widens at least as fast as the point estimate degrades and
  still covers 0.95. A reader of the point estimate is misled; a reader
  of the interval is not. Treat a subject-level standard deviation from
  twenty binary trials as a lower bound.

Do not read the bias as Laplace error alone: a variance component
estimated by maximum likelihood from binary data with few levels is
biased downward whether or not the integral is approximated, and the
simulation does not separate the two causes. Separating them is what
`frm(importance =)` exists for, and it now runs; see the next section
for what it is worth.

## The seam that closed, and what the correction is worth

`frm(importance =)` was refused for every family in the first release
and works now. The correction reweights draws from the Laplace Gaussian
and needs one log-likelihood value per GROUP; these likelihoods
factorize over subjects and again over trials, so the values always
existed, and what was missing was a slot to put them in.
`frmtmb_structure()` grew two, `loglik_row` and `loglik_group`, and
every family here declares both off the same recursion the objective
tapes. A subject is the group, a trial is the row, and `unit` is
unchanged at "one subject's trial sequence", because dropping a trial
changes every later trial's value store.

What the correction is worth depends on the design rather than on the
family. Measured on
[`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md),
40 subjects, one scalar random intercept, six datasets per cell, truth
`sd(alpha) = 0.5` (`dev/learn2-findings.md` has the per-dataset table):

|  |  |  |  |  |  |
|----|----|----|----|----|----|
| trials | datasets | Laplace `sd(alpha)` | corrected | mcse | min ESS per draw |
| 100 | 3 of 6 informative | 0.38 to 0.45 | 0.45 to 0.51 | 0.05 to 0.06 | 0.87 to 0.93 |
| 100 | 1 of 6 collapsed | 0.0004 | 0.0025 | 0.000 | 1.000 |
| 100 | 2 of 6 refused |  |  |  |  |
| 20 | 5 of 6 collapsed | 0.0001 | unchanged | 0.000 | 1.000 |
| 20 | 1 of 6 informative | 0.371 | 0.729 | 0.103 | 0.743 |

A COLLAPSED DATASET IS NOT A SHORT-SESSION PROBLEM. One of the six at
100 trials has it too, and reports the same spurious +1.822 as the five
at 20 trials do. Check `sqrt(VarCorr(fit))` on the Laplace fit at every
trial count, not only short ones.

- at 100 trials, on the three of six datasets whose variance component
  is identified, the correction moves it UP by 0.11 to 0.18 log units
  and toward the truth, its Monte Carlo error is a third of the shift,
  and the effective sample sizes are near 1. It costs 30 to 80 seconds
  against 2 for the Laplace fit. Of the other three, two refused on the
  correction's own convergence guard and one had already collapsed.

- at 20 trials it is not usable, and the reason is worth knowing. Five
  of six Laplace fits have already collapsed `sd(alpha)` to 1e-4, so
  there is nothing left to reweight; the correction then walks at its
  own step cap for every round and WARNS that it did not converge. A
  shift reported with `fit$importance$capped` true and every entry of
  `fit$importance$moves` equal is that step cap, not an estimate. The
  one 20-trial dataset with a real variance component overshoots, 0.37
  to 0.73 against a truth of 0.5.

DO NOT FOLLOW THE WARNING'S ADVICE ON A CAPPED RUN. It says to raise
`frmtmb_control(importance_rounds =)`, and on a collapsed variance
component that makes the artifact bigger in exact proportion: measured
on one 20-trial dataset, five rounds of a 0.3645 cap give a shift of
1.822 and ten rounds give 3.645, both exactly rounds times the cap. The
number is not an estimate that more iterations would sharpen. Read
`fit$importance$capped` and `$moves` first, and if the moves are all
equal, treat the correction as declining rather than answering.

So check `sqrt(VarCorr(fit))` on the Laplace fit before believing a
correction, and read the warning if there is one.

The correction requires the model to be grouped on the family's own
subject. `(1 | something_else)` gives a per-subject likelihood against a
per-something-else proposal, and the core refuses it by name rather than
adding up numbers that do not add up.

THE SAME SLOTS ALSO TURN `loo()` ON, at subject granularity, and an
earlier draft of this page said the opposite. `frm()` itself is maximum
likelihood and has no draws to average over, so `loo()` on a
`frmtmb_fit` still refuses for that core-wide reason; the route that
changed is `frmtmb.sample`, which admits any structure declaring
`loglik` together with either factorization slot and then reads the
COARSEST one. Measured on a 6-subject, 180-row fit: the pointwise matrix
is admitted, it has 6 columns rather than 180, it carries
`attr(x, "unit")` of `"one subject's trial sequence"`, and it sums to
`logLik(fit)` to twelve digits.

Read that literally. A column is a SUBJECT, so leaving one out drops
that subject's whole sequence, which is the only honest leave-one-out
for a recursion: a trial cannot be dropped because every later trial's
value store depends on it. `loo()` prints "Computed from N by K" and
says nothing about what a column is, so `frmtmb.sample` messages the
unit when the matrix carries one.

Deviance residuals are closer but still refused: the magnitude is now
available from `loglik_row()`, and what is missing is the SIGN, because
seven of the eight families have a nominal response with no mean to
depart from and the eighth has no exported route to the mean of a
first-passage time.

## Recovery tables for the three families measured this round

`dev/learn-recovery-round2.R`, 60 replicates each, 30 subjects by 100
trials, one random intercept on the primary parameter. Everything on the
NATURAL scale, with each Wald interval's endpoints pushed through the
same monotone link, because
[`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md)'s
non-decision time is bounded by its own dataset's fastest response and
its link scale is therefore not comparable across replicates.

**[`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md).**
Every bias inside its Monte Carlo error, 60 of 60 usable.

|           |       |        |       |           |          |
|-----------|-------|--------|-------|-----------|----------|
| parameter | truth | bias   | mc se | sd of est | coverage |
| `alpha`   | 0.35  | -0.002 | 0.003 | 0.026     | 0.87     |
| `drift`   | 3.00  | +0.001 | 0.009 | 0.069     | 0.98     |
| `bs`      | 1.60  | 0.000  | 0.002 | 0.019     | 0.95     |
| `ndt`     | 0.20  | +0.001 | 0.000 | 0.003     | 0.98     |
| `bias`    | 0.50  | -0.002 | 0.001 | 0.007     | 0.95     |

The pair that co-varies is `bs` with `ndt` at -0.57 and `drift` with
`bias` at +0.56. `drift` with `bs`, which is the pair a first-passage
density suggests, is 0.17.

**[`igt_orl()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_orl.md).**
Every parameter recovers; coverage 0.92 to 0.98, 60 of 60 usable.

|           |       |        |       |           |          |
|-----------|-------|--------|-------|-----------|----------|
| parameter | truth | bias   | mc se | sd of est | coverage |
| `Arew`    | 0.3   | -0.008 | 0.004 | 0.032     | 0.97     |
| `Apun`    | 0.1   | 0.000  | 0.001 | 0.007     | 0.98     |
| `k`       | 0.5   | -0.004 | 0.014 | 0.107     | 0.92     |
| `betaF`   | 1.0   | -0.015 | 0.013 | 0.097     | 0.95     |
| `betaP`   | 1.0   | +0.012 | 0.010 | 0.075     | 0.97     |

What is weak is the SEPARATION of the choice-rule parameters rather than
any one estimate: `betaF` with `betaP` correlate at -0.77, `k` with
`betaP` at -0.53 and `k` with `betaF` at +0.45. The two learning rates
correlate at 0.38 and with none of the three above 0.21.

**[`bandit4arm2_kalman_filter()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit4arm2_kalman_filter.md)
with `bonus = TRUE`**, and with `center`, `mu0` and `sigma0` held at the
task's own values through `bf(name = value)`, because those three do not
separate on a session this size whether or not the bonus is on.

|           |       |        |       |           |          |
|-----------|-------|--------|-------|-----------|----------|
| parameter | truth | bias   | mc se | sd of est | coverage |
| `tau`     | 0.15  | +0.005 | 0.001 | 0.011     | 0.85     |
| `lambda`  | 0.98  | 0.000  | 0.000 | 0.003     | 1.00     |
| `sigmaD`  | 3.00  | -0.045 | 0.026 | 0.203     | 0.95     |
| `phi`     | 1.50  | +0.059 | 0.018 | 0.141     | 0.93     |

`phi` trades off with `sigmaD`, at -0.81, and not with `tau`, at 0.12.
Both scale the same term: `phi` multiplies the posterior standard
deviation and `sigmaD` sets how fast it grows.

## What this package reads that frmtmb does not promise

Nothing. Every accessor it uses is exported and documented:
`frmtmb_family()`, `frmtmb_structure()`, `frmtmb_register_aterm()`,
`frmtmb_register_compat()`, `compat_rule_builder()`,
`single_response()`, `eval_dpars()` and `frame_block_of()`. There is one
thing it wanted and could not have, and it is recorded rather than
worked around: `frm_compat_features()` shows that an addition term is
registered but not at what ARITY, so this package cannot verify that a
`reward` term another package registered is the two-column one its
families need. It declines to register over the top of an existing term,
and the mismatch surfaces one step later as `frm()` refusing a missing
`reward2`.

## Not built, and named

The three things the first release left out are in:
[`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md),
[`igt_orl()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_orl.md)
and the Kalman filter's exploration bonus. What is still out:

- **[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) for
  [`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md).**
  Unlike the seven softmax families, this one's response is ordered and
  its conditional mean exists. The mean of a Wiener first-passage time
  conditional on the boundary reached belongs to `frmtmb.eam`, which
  exports the density and not the mean, so declaring `fitted_mean` here
  would mean deriving it again. That is a second seam and it is recorded
  rather than guessed at.
  [`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md)
  gives the per-trial density and drift rate meanwhile.

- **A mixture over learning strategies.** A real model and a wanted one.
  Core's `mixture()` combines per-ROW densities and would need
  per-sequence ones; the families now declare their per-sequence values,
  so what is left is a change under core's `R/` rather than a
  declaration this package can make.

- **`accepts_aterms` on these families.** frmtmb gained an addition-term
  allow-list, and these families do not declare one, so a term none of
  them reads is still carried without effect. Every term that would
  reshape a per-row contribution is already refused by name with a
  better message than an allow-list would give, so what the declaration
  would add is narrow.

- **A recovery table for every family.** Four of the eight have one now
  ([`bandit2arm_delta()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/bandit2arm_delta.md),
  [`rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.md),
  [`igt_orl()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/igt_orl.md)
  and the Kalman filter with `bonus = TRUE`). The other four have the
  Stan identity, the longhand reference and a smoke-level recovery
  assertion, which establishes that they compute the model they claim;
  what they do not have is a published bias-and-coverage table. The
  reason is time rather than principle: each table is 60 to 200 fits.

- **[`frm_task_design()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_task_design.md)
  for a task the user brings.** The design helpers cover the tasks these
  families are written for. A user fitting their own data needs none of
  it: the family reads whatever columns the formula names.

## See also

Useful links:

- <https://aforren1.github.io/frmtmb/frmtmb.learn>

- <https://github.com/aforren1/frmtmb>

- Report bugs at <https://github.com/aforren1/frmtmb/issues>

## Author

**Maintainer**: Alex Forrence <alex.forrence@gmail.com>
([ORCID](https://orcid.org/0000-0002-9728-6337))

Authors:

- Alex Forrence <alex.forrence@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9728-6337))

## Examples

``` r
# a reversal task, and the change in learning rate with an interval
d <- frm_task_design("reversal", n_subject = 10, n_trial = 40, seed = 5)
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = id, trial = trial), d,
  pars = list(alpha = 0.35, tau = 3), seed = 5)[[1]]$choice
fit <- frmtmb::frm(
  frmtmb::bf(choice | reward(pay1, pay2) ~ after_reversal, tau ~ 1),
  family = bandit2arm_delta(subject = id, trial = trial), data = d)
frmtmb::fixef(fit)
#> $alpha
#>         (Intercept) after_reversalafter 
#>          -0.9385871           0.1404640 
#> 
#> $tau
#> (Intercept) 
#>    1.118225 
#> 
```
