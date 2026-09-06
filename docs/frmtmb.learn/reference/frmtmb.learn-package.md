# frmtmb.learn: Reinforcement-Learning Families for 'frmtmb' Models

The value-learning models of the reinforcement-learning and
computational-psychiatry literature, written as 'frmtmb' families so
that each of their parameters is an ordinary distributional parameter
with its own linear predictor. A learning rate then takes a condition
effect, a smooth term or a correlated per-subject random effect the way
a mean does, and it comes back with a standard error, which is what a
separate model per group does not give. Six families share one recursion
that walks trials once and updates every subject at each step: the
two-armed delta learner and its dual-rate and counterfactual variants, a
Kalman filter over the restless four-armed bandit of Daw and others
(2006), the two-stage model-based and model-free hybrid of Daw and
others (2011), and a prospect-theory learner for the Iowa gambling task.
Families are named as 'hBayesDM' names them where a name exists. Each
one is checked by an identity against an independent 'Stan' program of
the same model at the same estimates, and by parameter recovery at a
realistic scale.

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

## One recursion, six families

Every family here is the same walk: carry a value store per subject,
visit the trials in order, turn the store into a choice probability and
the outcome into a new store. The walk is written once. A family
supplies a choice rule and a learning rule, and is a few dozen lines.

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
`frm(importance =)` exists for, and every family here refuses it.

## The seam this package is waiting on

`frm(importance =)` is refused for every family, and the refusal is
wider than the mathematics. The correction needs one log-likelihood
value per GROUP. These likelihoods factorize over subjects and again
over trials, so the values exist; `frmtmb_structure(loglik =)` returns
one AD scalar for the whole response and there is no slot to put them
in. The same gap costs deviance residuals, which need a per-row
saturated comparison. It does NOT cost `loo()` or `waic()`: those refuse
for every frmtmb fit, an ordinary gaussian one included, because `frm()`
is maximum likelihood and an elpd averages the likelihood over draws.
That refusal is core-wide, unrelated to this seam, and closing the seam
would not deliver it. The design is written up in frmtmb's
`dev/rl-findings.md` under "Protocol seams": one slot carrying the
finest factorization a family has, per row where one exists and per
group otherwise, with `unit` left as the separate declaration of the
leave-one-out granularity. This package computes a per-row conditional
likelihood already, for
[`frm_value_trace()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/frm_value_trace.md),
so it would fill a per-row slot the day one exists.
`frm_compat("bandit2arm_delta", "importance")` says so at the console.

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

RLDDM, a learning rule feeding a drift-diffusion choice rule, is the
next family and is deliberately absent. The choice rule belongs to
`frmtmb.eam`, and building it here would make this package import a
sibling extension. The engine already has the shape it needs: a choice
rule is a function from a value store to a per-option quantity, and the
Wiener first-passage density is one such function.

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
