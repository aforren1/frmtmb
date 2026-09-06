# frmtmb.learn: a learning-model extension

Lane: `wt-learn`, worktree `C:/Users/adf44/source/r/frmtmb-wt-learn`,
branch point `de9d639` (frmtmb 0.52.0 plus two CI fixes). Findings are
written as they are measured; a number in this file was produced by a
script in this lane, not remembered.

## What this package is

The repository already carries a reinforcement-learning worked example:
`inst/rl/rw-delta.R`, promoted through `vignette("reinforcement-learning")`.
It is one family written out longhand to show the structured-family
protocol. This package promotes it to a package of families, and the
thing it adds is the part the example says a real package would need:

> **More of the literature.** A real package would carry the family of
> models this one is the simplest member of: separate learning rates
> for gains and losses, a decay term, more than two arms. All of them
> share this block and this loop shape, and differ only in the body.
> (`vignettes/reinforcement-learning.Rmd`)

So the deliverable is one recursion engine and a set of family bodies
on top of it.

## Toolchain, established before any code

| fact | value |
|---|---|
| R | 4.6.1 (2026-06-24 ucrt), Windows 11 |
| RTMB | 1.9 |
| core installed into ln-lib | frmtmb 0.52.0 from this worktree |
| rstan / StanHeaders | 2.32.7 / 2.39.1, user library, read fallback |
| hBayesDM | NOT installed at start; see the cross-check section |

## Probes that settled the engine's shape

Measured, not assumed, before the engine was written.

* **`sign()` tapes, and its branch derivatives are the right ones.**
  `MakeTape(function(x) {s <- sign(x - 0.5); sum(0.5*(1+s)*x + 0.5*(1-s)*2*x)})`
  returns the correct value and the jacobian `(2, 1, 1)` at
  `(0.1, 0.7, 0.9)`. This is what a dual-learning-rate family needs: the
  sign of a PREDICTION ERROR is a comparison on an AD quantity, and RTMB
  has neither a comparison operator nor `CondExp`.
* **`RTMB::logspace_add()` works on plain doubles as well as on
  advectors**, returning `2.126928` for `logspace_add(0, 2)`. So does
  `sign()`. That is what lets ONE recursion serve the taped likelihood,
  `fitted()` and the simulator. The worked example needed two copies of
  its loop (`rw_loglik` and `rw_replay`); this package needs one.
* **An aterm is visible through `frm_compat_features()`** as
  `name = "reward()"`, `key = "reward"`, `kind = "aterm"`, so a second
  registration of a name can be DETECTED. Its arity cannot: the table
  carries name, key and kind only. Recorded as a seam below.

## The engine

`extensions/frmtmb.learn/R/engine.R`. One `ln_recurse()` walks trials
once and updates every subject at each step. A family supplies three
tape-safe functions through `ln_spec()`: `init(ns, d1)`,
`choice(state, d, j)` returning per-option utilities, and
`update(state, d, ch)` returning the new store plus anything worth
recording. A family file is then the model and nothing else, which is
what makes the six of them 40 to 120 lines each.

Four decisions in it that are not obvious.

**One recursion serves three jobs.** `mode` is "loglik" (an AD scalar),
"trace" (per-row vectors) or "simulate" (a draw). The worked example
this generalizes needs two copies of its loop, `rw_loglik` and
`rw_replay`, kept in step by hand. One copy is possible only because
`RTMB::logspace_add()` and `sign()` both accept doubles as well as
advectors, which was measured before the engine was written rather than
assumed.

**The mask is applied by the engine, not by the rule.** Padded cells are
frozen with `ln_blend()`, `state + m * (new - state)`, so a learning
rule never has to thread the mask correctly. The worked example folds
the mask into its learning rate instead, which is one operation cheaper
per state slot and has to be got right again in every new rule.

**A length-one parameter is not broadcast.** `ln_at()` indexes only
entries longer than 1, so a distributional parameter with no formula
stays one number instead of becoming `n` tape nodes that are then read
back `n_subj` times per trial. The worked example broadcasts
unconditionally.

**Utilities, not probabilities.** `choice()` returns one unnormalized
log-probability per option and the engine does the softmax through
`ln_lse()`. For two options this is the same arithmetic as the worked
example's `logspace_add(0 * eta, eta)`; for four it is the only spelling
that stays branch-free. A trial may hold more than one decision
(`n_option = c(2, 2)`), which is what lets the two-step task be the same
engine.

## Two AD facts the families needed

* **`sign()` selects a branch that a comparison cannot.** RTMB refuses a
  comparison on an advector and has no `CondExp`. `bandit2arm_dual(split
  = "pe")` needs the sign of a PREDICTION ERROR, which is an AD
  quantity; `0.5 * (1 + sign(pe))` selects the rate exactly, and the
  kink at `pe == 0` is the model's own, because the update is zero there
  under both rates so the branches agree at the crossing. `split =
  "outcome"` needs the sign of the OUTCOME, which is data, so it selects
  with a column of zeros and ones and adds nothing to the tape.
* **`ln_max2(a, b) = 0.5 * (a + b + abs(a - b))`** gives `ts_par7()` the
  best stage-two option without a branch. `abs()` is defined on an
  advector, and the kink is the maximum's own.

## Measured: the two dual-rate splits are ONE model on binary payoffs

Not a coincidence, and pinned by a test. With payoffs in 0/1 the value
store stays in [0, 1), so `pe = pay - Q` is positive exactly when `pay`
is 1, and the two splits select the same rate on every trial. Measured:
identical log-likelihoods to 1e-8 on the same data. They are different
models only when the payoff is GRADED, which the same test shows by
drawing payoffs from `runif(-1, 1)` and getting two different optima. A
reader who codes a reversal task 0/1 and expects `split` to matter is
wrong, and the help says so.

## The decision this package is most likely to be questioned on

**No `fitted()`, no `predict(type = "response")`, no residuals of any
type.** Every family declines to declare `fitted_mean`.

Core forms a fitted value as the conditional mean of the response and a
residual as `y - mean` (`R/predict.R:2399`), and the response here is
the option a subject took, coded 1 to K. It is NOMINAL. Arm 2 is not
twice arm 1 and the Iowa gambling task's four decks have no order at
all, so there is no number `mu` for which `y - mu` is a residual.
Supplying one would make three methods return arithmetic on a category
code, silently.

This was found by measurement, not by reading. The first draft DID
supply `fitted_mean` as the probability of the observed choice, and the
suite caught pearson residuals coming back as `(y - p) / sqrt(p(1-p))`
with `y` in 1/2: 3.5 where the binary residual is 1.4.

What replaces it is `frm_value_trace()`, which returns more than a mean
could: the per-option value estimates each choice was made on, the
prediction error, and the fitted probability of the choice that was
made. `core::cox()` and `frmtmb.spline::royston_parmar()` decline to
invent a mean for the same kind of reason, and both point at a richer
accessor instead.

The cost is real and is stated in the help rather than left to be
discovered: frmtmb's own `rw_delta` example DOES have `fitted()`,
because it codes a two-armed choice 0/1 and returns P(arm 1), so
`y - p` is the ordinary binary residual. That does not survive a fourth
arm. One engine over two and four options was judged worth more than
`fitted()` on the two-option half of it.

## Measured: what sum(log(p)) is, and what it is not

`frm_value_trace()$p` is the per-trial factor of the likelihood, so
`sum(log(p))` is the CONDITIONAL data log-likelihood given the fitted
subject effects.

| fit | sum(log(p)) | logLik(fit) | difference |
|---|---|---|---|
| no random effects | -159.7215 | -159.7215 | 2.8e-14 |
| one random intercept, sd collapsed to 9.9e-05 | -159.7215 | -159.7215 | 1.8e-07 |
| one random intercept, 30 subjects, sd = 1.07 | -1507.816 | -1542.488 | 34.672 |

The first row is the identity and the suite asserts it. The third is the
honest general statement: `logLik()` reports the Laplace-approximated
MARGINAL likelihood, which adds the random-effect density and the
curvature term a conditional quantity does not carry. The middle row is
the trap: on a fit whose variance component collapses, the two agree to
1.8e-07, and a test written against that fixture would assert a false
identity and pass. This was caught because the first draft asserted the
identity on a hierarchical fixture whose sd happened to collapse.

## Two seams recorded rather than worked around

* **An addition term's ARITY is not readable through any exported
  route.** `frm_compat_features()` shows a registered term's name, key
  and kind, so this package can detect that `reward` is already
  registered and decline to register over the top of it, which is what
  `ln_register_aterms()` does. That is a SILENT SKIP rather than a
  refusal, and it is the right behaviour: reloading a package re-runs
  `.onLoad()` and must stay harmless. An earlier summary of this lane
  called it a refusal; it is not, and `zzz.R` never claimed otherwise. It cannot check that the existing term
  has the two columns its families need.
  `frmtmb_register_aterm()` REPLACES a same-named entry silently, which
  is the sharp edge frmtmb's own `dev/rl-findings.md` records under "Two
  core sharp edges". The consequence here is bounded: a family declares
  `required_aterms` by indexed name, so a `reward` registered at arity 1
  makes `frm()` refuse the fit for a missing `reward2` rather than fit
  something wrong.
* **`frm(importance =)` is refused for every family, by name**, and the
  refusal is wider than the mathematics. The correction needs one
  log-likelihood value per GROUP; these likelihoods factorize over
  subjects and again over trials, so the values exist, but
  `frmtmb_structure(loglik =)` returns one AD scalar and there is no
  slot to put them in. This package already computes a per-ROW
  conditional likelihood, for `frm_value_trace()`, so it would fill a
  per-row slot the day one exists. The design is frmtmb's own
  `dev/rl-findings.md` "Protocol seams" section: the slot carries the
  finest factorization the family HAS, with `unit` left as the separate
  declaration of leave-one-out granularity. Not implemented here: it is
  a change under core's `R/`, which this lane has no mandate for.

## An unreachable guard, kept and labelled

`ln_block()` refuses a trial column with missing values. It is not
reachable through `frm()`: the trial variable is in the model frame,
because the family puts it there through `frame_vars`, so `na.action`
drops the row and reports "1 row removed because of missing values"
before the block runs. The guard is kept because a trial expression that
COMPUTES an NA from non-missing columns would still meet it, and the
test asserts the na.action behaviour that actually happens.

## Measured: the Iowa gambling deck-B effect

`frm_task_design("igt")` uses the classic contingencies, two bad decks
(1 and 2) with large gains and larger losses and two good decks (3 and
4). Under the PVL utility at `shape = 0.4` and `lambda = 2` the model
does NOT sort them that way, and the first version of the test asserted
that it did and failed.

| deck | mean payoff | mean subjective utility |
|---|---|---|
| 1 | -0.306 | -0.752 |
| 2 | -0.359 | +0.314 |
| 3 | +0.251 | +0.381 |
| 4 | +0.250 | +0.418 |

Deck 2's loss is rare and large, and the compressive exponent shrinks it
to a fraction of its size, so decks 2, 3 and 4 come out within 0.11 of
each other while deck 1, whose loss is frequent, sits far below. That is
the well-known deck-B effect: a property of the model rather than a
fault in the design. Measured choice share of deck 1: 0.270 over the
first ten trials, which is chance, and 0.100 over the last twenty. The
test follows the model.

## Suite status

Ungated suite, one process: **175 passing, 0 failing, 6 skipped**. The
six skips are the Stan identity tier, which is gated on
`FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true` like the brms tier.

## The bug the Stan tier caught, and what it says about the engine

The Kalman family failed its Stan identity by 712.6 log-likelihood
units with an inner gradient of 149.9, on a run where the other five
families agreed to between 8.5e-14 and 2.3e-13. It was a real bug in
this package.

`ln_blend()`, which freezes a padded cell by mixing the old value store
with the new one, matched slots BY POSITION. That is correct only while
a learning rule happens to rebuild its store in the order `init()`
declared it. Five of the six rules do. The Kalman filter does not: it
fills one loop over the four arms and so returns `mu1, s1, mu2, s2, ...`
where `init()` returned `mu1..mu4, s1..s4`. Every posterior mean was
therefore blended with a posterior variance on every trial.

What makes it worth recording is how it behaved. The fit CONVERGED. It
recovered `tau` to 0.150 against a truth of 0.15, `lambda` to 0.982
against 0.98 and `sigmaD` to 3.04 against 3, which is what a smoke test
and a recovery table both call a pass. The likelihood it reported was
wrong.

**By how much depends on which two numbers you compare, so the fixture
and the comparison are named here rather than a bare figure.** Measured
by reverting `ln_blend()` to positional matching in a scratch copy of
the package, installing it to a second library, and running BOTH
versions on the Kalman fixture of `tests/testthat/test-reference.R`
(6 subjects, 40 trials, seed 73):

| version | engine `sum(log(p))` | longhand reference at the same estimates | gap |
|---|---|---|---|
| positional (buggy) | -163.0275 | -368.9654 | **205.94** |
| by name (shipped) | -211.8118 | -211.8118 | 0.0000 |

That 205.94 is the figure this document used to quote as "206", and it
is the gap between what the buggy engine REPORTED and what its own
parameters were actually worth under the model. It is the comparison
`test-reference.R` makes, which is why it is the one recorded.

Two other comparisons on a different fixture give different numbers, and
all three are correct for their own question. On the identity tier's
Kalman fixture (10 subjects, 500 rows, seed 104) the reviewer measured
70.6 units comparing the correct engine at the buggy fit's parameter
vector against the buggy engine's report there, and 164.3 units
comparing converged `logLik` against converged `logLik`. A bug that
changes the objective changes where the optimizer stops as well as what
it reports, so "how wrong was it" has no single answer; what all three
comparisons agree on is that the answer is of order a hundred
log-likelihood units on a few hundred rows, while the parameter
estimates look fine.

Two checks were added, and neither is the Stan tier:

* `ln_blend()` now matches by NAME and refuses a rule that returns a
  different set of slot names than `init()` declared. The refusal is
  cheap: it runs once per trial at tape build, not per gradient.
* `tests/testthat/test-reference.R` writes every family out longhand,
  one subject and one trial at a time with `if` and `[` for selection,
  and compares against `sum(log(frm_value_trace(fit)$p))` on fits with
  no random effects. It is UNGATED, so it runs in CI where the Stan tier
  does not, and it is what caught the bug once the Stan tier had pointed
  at the family. Six families, six references.

The lesson for the engine's design is the one the protocol already
implies and this made concrete: a value store is a NAMED collection and
every access to it, including the engine's own, must go through the
name. The repository's RTMB notes record `$` partial matching as the
same class of hazard.

## Identity against Stan, per family

Method: each program declares exactly the parameters frmtmb estimates on
the scales frmtmb estimates them on, so the map is the identity and the
comparison carries no Jacobian; the random-effect standard deviation
rides along as data. `frm()` carries no default prior, the programs use
`target +=` with `_lpmf` and `_lpdf` so they drop no normalizing
constant, and every declared parameter is unbounded so
`adjust_transform = FALSE` has no Jacobian to discard. The difference is
expected to be exactly zero rather than merely small.

The gradient column is the check that catches a wrong map: frmtmb puts
the subject effects at their conditional modes, so Stan's gradient with
respect to THEM must vanish. The gradient with respect to the fixed
effects is not zero there and is not asserted on, because maximum
likelihood makes the marginal likelihood stationary, not the joint one.

## Identity against Stan, measured

| family | subjects, rows | Stan log_prob | frmtmb joint | difference | max grad on subject effects | max grad overall | population grad vs frmtmb's |
|---|---|---|---|---|---|---|---|
| `bandit2arm_delta` | 12, 720 | -331.1511111995 | -331.1511111995 | 1.705e-13 | 2.331e-15 | 2.760e+00 | 1.332e-15 |
| `bandit2arm_delta (displaced)` | 12, 720 | -333.4730381567 | -333.4730381567 | -2.842e-13 | NA | NA | NA |
| `bandit2arm_dual (pe)` | 12, 720 | -404.0743384819 | -404.0743384819 | 5.684e-14 | 3.974e-02 | 2.013e+00 | 8.882e-15 |
| `bandit2arm_dual (outcome)` | 12, 720 | -312.2230591683 | -312.2230591683 | -2.842e-13 | 1.443e-15 | 3.462e-04 | 9.864e-15 |
| `prl_fictitious` | 12, 720 | -156.7472771410 | -156.7472771410 | 8.527e-14 | 8.983e-11 | 5.197e+00 | 1.155e-14 |
| `bandit4arm2_kalman_filter` | 10, 500 | -141.7830467812 | -141.7830467812 | 0.000e+00 | 9.544e-14 | 1.588e-04 | 1.323e-13 |
| `igt_pvl_delta` | 12, 720 | -786.0320609229 | -786.0320609229 | 5.684e-13 | 1.176e-14 | 1.418e-06 | 1.029e-14 |
| `ts_par7` | 12, 720 | -741.8237574444 | -741.8237574444 | 2.274e-13 | 2.102e-15 | 4.518e-05 | 6.795e-14 |

Every difference is at the level of double-precision accumulation over
several hundred trials, which is what "exactly zero" looks like when the
two sides sum in different orders. The displaced row is a deliberately
non-stationary point, its subject effects shifted by 0.3, so the
agreement cannot be an artifact of both sides sitting at an optimum: the
objective moves from -331.15 to -333.47 and the two implementations move
together.

The gradient columns are the check that catches a wrong parameter map,
and the two are meant to differ. The one on the subject effects vanishes
because frmtmb puts them at their conditional modes. The overall one
does not, and is not asserted on, because maximum likelihood makes the
MARGINAL likelihood stationary rather than the joint one.

One more thing this tier caught, and it errored rather than
disagreeing. `bandit2arm_delta()` is the only family here with a single
non-primary parameter, so its program declares `vector[1] bd`, and a
length-one R numeric reaches rstan with no `dim`: "dims declared=(1);
dims found=()". Its row was simply MISSING from the table for two runs
rather than showing a bad residual, which is worth knowing about a
harness that writes its results to a log. `as.array()` on every
parameter block fixes it.

### The population-parameter gradient, added in the punch round

The last column is a second identity and it is not against zero. Stan's
gradient of the log density on the `b` and `bd` positions is compared
with frmtmb's own joint gradient there, `-obj$env$f(par, order = 1)` at
the matching names. It agrees to between 1.3e-15 and 1.3e-13.

It has to be compared against frmtmb rather than against zero, and the
"max grad overall" column is why: those gradients are 2.76, 2.01 and
5.20 on three of the rows. Maximum likelihood makes the MARGINAL
likelihood stationary, not the joint one, so the joint gradient with
respect to the population parameters is far from zero at the estimates.
The earlier version of this harness computed that number, logged it and
asserted nothing on it. It now carries an identity that would catch a
wrong fixed-effect map, which a single-point value check can miss.

### One row that used to assert nothing new

`bandit2arm_dual`'s two splits were previously both fitted on BINARY
payoffs, where this package proves they are the same model, so the
`(outcome)` row reproduced the `(pe)` row digit for digit and the
`outcome` branch of the Stan program was never exercised where it
differs. The fixture now draws payoffs from `runif(-1, 1)`: the two rows
are -404.0743384819 and -312.2230591683, and the test asserts they
differ.

That change surfaced a real property of the model, now recorded in
`?bandit2arm_dual`. On graded payoffs the `pe` split selects its rate
with `sign(pe)`, so the joint log density has a kink in the random
effects wherever a prediction error crosses zero and TMB's inner Newton
solve stops short of the conditional mode: the subject-effect gradient
is 3.97e-02 under `"pe"` against 1.44e-15 under `"outcome"` on the same
data, whose selector is data and therefore smooth. The VALUE identity is
untouched at 5.7e-14, because both sides evaluate at the same point, so
the fixture keeps the full tolerance on the value and relaxes only the
mode-convergence assertion, with the reason in the test.

Suite: 39 assertions in the gated tier, 0 failures.

## Timings, measured

Method: R 4.6.1, Windows 11, one core, warm session, designs and shapes
INTERLEAVED (one measurement of each per round, not all of one and then
all of the next), `gc(FALSE)` before every measurement,
`frm(dry_run = "objective")` with the `dry_run = "frame"` time
subtracted, median of 7 builds, and gradients timed over batches of 50
`obj$gr()` calls with the best batch of 3 reported. Script:
`extensions/frmtmb.learn/dev/learn-timing.R`.

### Scaling, bandit2arm_delta with one random intercept

| subjects | trials | rows | tape build | one gradient |
|---|---|---|---|---|
| 10 | 100 | 1000 | 0.14 s | 18.8 ms |
| 40 | 25 | 1000 | 0.05 s | 13.2 ms |
| 40 | 100 | 4000 | 0.17 s | 80.0 ms |
| 40 | 400 | 16000 | 0.66 s | 520.6 ms |
| 160 | 100 | 16000 | 0.41 s | 382.4 ms |
| 40 | 800 | 32000 | 1.42 s | 906.4 ms |

Both costs track the ROW count, which is the tape's node count. Read the
two 16000-row rows against each other for the part that does not. The
BUILD is 0.66 s against 0.41 s: more trials at a fixed row count costs
more to tape, because the R loop runs four times as many iterations, and
that is the expected direction.

The GRADIENT column at fixed rows is 520.6 ms against 382.4 ms, and this
lane does not offer an explanation for it. The tape is the same size
either way, so tape size does not predict it; the two designs also
differ in the size of the inner Laplace problem (40 subject effects
against 160), and the direction is the opposite of what more inner
parameters would suggest. Do not read a rule from that pair.

### Neither timing column is portable between machines

The BUILD column is quoted to two significant figures and should be read
for its SHAPE, not its level. A reviewer reran two rows on a quiet
machine with the same interleaved method and medians of 5 and got 0.27 s
against this document's 0.14 s at 1000 rows, and 1.55 s against 0.66 s
at 16000: about 2x at both sizes. The scaling is the same either way
(5.7x their build for 16x the rows, against 4.7x here), which is the
part that transfers. The same caveat applies to the elementwise factors
below, where the reviewer measured 7.1x and 22.1x against 3.6x and
14.2x here: the growth reproduces, the multipliers do not.

### The gradient column is not stable to better than about 1.5x

Worth stating because it bounds every gradient number above. The SAME
design size, 40 subjects by 100 trials, measured three times on three
simulated datasets, gave 80.0, 108.6 and 68.2 ms. The measurement is
batched and takes the best of three batches, so this is not clock noise:
`obj$gr()` on a random-effects object re-solves the inner problem, and
how many Newton steps that takes depends on the DATA. Read the gradient
column for its order of magnitude and its scaling, not to three
significant figures. This is the mistake frmtmb's own `dev/rl-findings.md`
had to retract, in a different form.

### The elementwise penalty, paid at tape build

The engine loops over TRIALS and vectorizes over SUBJECTS. The
alternative is one iteration per ROW with the value store updated by
sub-assignment, which is how the model is usually stated. The family
measured here differs from `bandit2arm_delta()` in the loop shape and in
nothing else.

| rows | vectorized | elementwise | factor | objective gap |
|---|---|---|---|---|
| 1000 | 0.11 s | 0.40 s | 3.6 | 3.4e-13 |
| 4000 | 0.16 s | 1.52 s | 9.5 | 3.6e-12 |
| 10000 | 0.35 s | 4.53 s | 12.9 | 3.5e-11 |
| 20000 | 0.74 s | 10.51 s | 14.2 | 1.3e-11 |

The two spellings agree on the objective throughout, so this is a cost
comparison of one model rather than of two. The factor GROWS with the
row count, from 3.6x to 14.2x, because RTMB's replacement operator
copies the vector it writes into. That reproduces the shape frmtmb's
`dev/rl-findings.md` measured (5.3x at 1000 rows to 12.7x at 20000) on a
different family, and it is the shape rather than any single multiplier
that is the finding.

### And not paid per gradient

| loglik shape | tape build, 4000 rows | one gradient |
|---|---|---|
| vectorized over subjects | 0.16 s | 108.6 ms |
| row loop, sub-assigned store | 1.52 s | 90.4 ms |

The build differs by 9.5x and the gradient does not differ at all: 0.83x
is inside the 1.5x spread the section above establishes for this
measurement. That is what it has to be. Once the tape exists it is a
node list and the R code that built it is gone. An optimization runs
hundreds of gradients against one build, so the penalty is real but it
is paid once.

### Every family, one size

40 subjects by 100 trials, 4000 rows, one random intercept on the
primary parameter.

| family | value stores | decisions per trial | tape build | one gradient |
|---|---|---|---|---|
| `bandit2arm_delta` | 2 | 1 | 0.13 s | 68.2 ms |
| `bandit2arm_dual` | 2 | 1 | 0.21 s | 125.0 ms |
| `prl_fictitious` | 2 | 1 | 0.18 s | 48.8 ms |
| `bandit4arm2_kalman_filter` | 8 | 1 | 0.58 s | 109.2 ms |
| `ts_par7` | 8 | 2 | 0.46 s | 55.0 ms |
| `igt_pvl_delta` | 4 | 1 | 0.29 s | 94.6 ms |

The BUILD column is the one to read, and it tracks the number of value
stores the rule carries: two stores build in 0.13 to 0.21 s, four in
0.29 s, eight in 0.46 to 0.58 s. That is the engine's own arithmetic,
because `ln_blend()` touches every store on every trial. The gradient
column spans 48.8 to 125.0 ms with no ordering that survives the
stability bound above, so nothing should be read from it beyond "all six
are the same order of magnitude at this size".

## Parameter recovery, and the Laplace caveat, measured

Design: `bandit2arm_delta()` on a reversal task, 40 subjects, a learning
rate of 0.30 before the reversal and a `+1.0` logit effect after it, a
softmax sensitivity of 3, and a subject-level standard deviation of 0.5
on the learning rate's logit. 100 replicates per row, each drawn through
the family's own generative simulator and refitted. Script:
`extensions/frmtmb.learn/dev/learn-recovery.R`.

### 100 trials per subject

| parameter | truth | bias | Monte Carlo se | sd of estimates | coverage |
|---|---|---|---|---|---|
| `alpha_(Intercept)` | -0.847 | +0.009 | 0.015 | 0.150 | 0.97 |
| `alpha_after` | 1.000 | -0.015 | 0.017 | 0.173 | 0.94 |
| `tau_(Intercept)` | 1.099 | +0.003 | 0.003 | 0.034 | 0.95 |
| `log sd(alpha)` | -0.693 | -0.266 | 0.088 | 0.880 | 0.99 |

All three fixed effects are unbiased to within their Monte Carlo error
and cover at the nominal rate. 100 of 100 replicates produced a usable
interval.

### 20 trials per subject, the same everything else

| parameter | truth | bias | Monte Carlo se | sd of estimates | coverage |
|---|---|---|---|---|---|
| `alpha_(Intercept)` | -0.847 | -0.026 | 0.029 | 0.286 | 0.92 |
| `alpha_after` | 1.000 | +0.009 | 0.032 | 0.322 | 0.91 |
| `tau_(Intercept)` | 1.099 | +0.016 | 0.009 | 0.092 | 0.96 |
| `log sd(alpha)` | -0.693 | -1.596 | 0.301 | 3.012 | 0.95 |

98 or 99 of 100 replicates produced a usable interval, and 27 fits
warned about convergence.

### What this says, and one thing it corrects

The FIXED effects survive short sessions. At a fifth of the trials their
bias is still inside Monte Carlo error and their intervals still cover
near the nominal rate (0.91 to 0.96 against 0.95). A study that wants
the condition effect on a learning rate can work with short sessions.

The variance component's POINT ESTIMATE does not survive. `log sd(alpha)`
comes out 1.60 too low at 20 trials, which is a standard deviation about
a fifth of the true one, against 0.27 too low at 100 trials.

**Its INTERVAL does survive, and this lane wrote the opposite before
measuring it.** The first draft of `?frmtmb.learn` said the interval
"undercovers", by analogy with frmtmb's own RL vignette, which measured
0.73 at 20 trials on a different design. Measured here it covers 0.95,
and the reason is visible in the neighbouring column: the spread of the
estimates rises from 0.880 to 3.012, so the interval widens at least as
fast as the estimate degrades. The honest statement is therefore
narrower than the one it replaces, and more useful:

* a subject-level standard deviation estimated from twenty binary trials
  is a LOWER BOUND, and its point estimate should not be reported;
* its interval is not misleading, because it is wide enough to know it;
* a user who reads the point estimate is misled and a user who reads the
  interval is not.

Do not read the bias as Laplace error alone. A variance component
estimated by maximum likelihood from binary data with 40 levels is
biased downward whether or not the integral is approximated, and this
study does not separate the two causes. Separating them is what
`frm(importance =)` exists for, and every family here refuses it: see
the seam section above.

## Left out, and why

* **RLDDM.** A learning rule feeding a drift-diffusion choice rule, and
  the family this package would add next. Not built, deliberately: the
  Wiener first-passage density belongs to `frmtmb.eam`, and building it
  here would make this package import a sibling extension, which the
  repository's dependency shape does not allow. The engine already has
  the seam for it: a choice rule is a function from a value store to a
  per-option quantity, and a first-passage density is one such function.
  What it would need beyond that is a second response column, the
  reaction time, alongside the choice.
* **ORL for the Iowa gambling task.** PVL-delta ships instead.
  PVL-delta's four parameters recover at a realistic scale, which is the
  standard this package holds a family to; ORL carries two further value
  stores and two frequency-weighting parameters that are known to be
  weakly identified, so shipping it would mean shipping a recovery table
  that says so. It is about forty lines on this engine when wanted.
* **An exploration bonus on the Kalman filter.** Daw and others compare
  softmax choice against rules that add a bonus for uncertainty, and the
  filter here already carries the posterior variance such a rule needs
  (`frm_value_trace()` returns it as `s1` to `s4`). Left out because a
  bonus and `tau` are hard to separate at the data sizes these studies
  collect, so it would ship with a recovery table saying the two trade
  off.
* **The per-group log-likelihood slot.** The seam `frm(importance =)`,
  `loo()` and `waic()` all wait on. It is a change under core's `R/`
  (structure.R, objective.R, importance.R, predict.R), which this lane
  has no mandate for and which sibling lanes are editing this round. It
  is refused by name instead, with the seam named in the refusal.
* **A recovery table per family.** One family is measured, at two
  session lengths, 100 replicates each. The other five have the Stan
  identity, the longhand reference and a smoke-level recovery assertion
  in the suite, which establishes that they compute the model they
  claim; what they do not have is a published bias-and-coverage table.
  `dev/learn-recovery.R` is written against `bandit2arm_delta()` and
  would take a `family =` argument to generalize. The reason is time
  rather than principle: each table is 200 fits.
* **`frm_task_design()` for a task the user brings.** The design
  helpers cover the five tasks the six families are written for, with
  canonical column names. A user fitting their own data needs none of
  it: the family reads whatever columns the formula names.

## Two observations for the coordinator, outside this lane's ownership

* `dev/build-docs.R` prints "installed five packages" and its header
  says "the four extension subsites". Both counts were already wrong
  before this lane (six packages were installed, not five), and adding
  `frmtmb.learn` makes them wronger. Not touched: the brief allows one
  subsite entry in that file, which is what was added.
* `frmtmb_register_aterm()` still accepts a re-registration of an
  existing name at a DIFFERENT arity and clobbers the first entry
  silently, and no exported accessor reports a registered term's arity.
  frmtmb's own `dev/rl-findings.md` records this under "Two core sharp
  edges"; this package works around it by declining to register a name
  that is already present.

## Verification run

| check | result |
|---|---|
| package suite, one process, ungated | 184 passing, 0 failing, 6 skipped |
| the 6 skips | the Stan identity tier, gated on `FRMTMB_BRMS_FIT_TESTS` and `NOT_CRAN` |
| Stan identity tier, gated on | 39 assertions, 0 failing |
| `R CMD check --as-cran` | **Status: OK**, 0 errors, 0 warnings, 0 notes |
| roxygen idempotence | `man/` and `NAMESPACE` byte identical after a second run |
| vignette knit, pandoc 3.8.3 | renders, 48 KB |
| workflow against the spline one | structurally identical apart from two named skips |
| core `R/` changed | none |
| core `tests/` changed | none |

`R CMD check` was run against the built tarball with the core installed
in the lane's private library, `_R_CHECK_CRAN_INCOMING_=false`,
`_R_CHECK_FORCE_SUGGESTS_=false` and `NOT_CRAN=true`, with pandoc and
TinyTeX on PATH, so the run includes `checking PDF version of manual
... OK` rather than skipping it under `--no-manual`. The first run had one NOTE, "Non-standard
files/directories found at top level: '_pkgdown.yml' 'dev'", because
this package had no `.Rbuildignore`; every sibling extension has one and
this one now carries a byte-identical copy of `frmtmb.spline`'s. The
second run is clean.

Inside `R CMD check` the suite reports fewer assertions than it does
from sources, because `test-message-uniqueness.R` skips itself when
there are no sources to scan. That is the same behaviour the spline
package documents in its own workflow.

### The workflow, compared line by line

`diff` of the two files with comments and blank lines stripped and
`spline` rewritten to `learn` gives exactly one hunk: the dependency
step drops `rstan` and `hBayesDM` by name before installing. Everything
else, including the `NOT_CRAN` repetition at both job and step level, the
path filters, the core install from the checkout and the
`check-r-package` arguments, is the same. The two skips are argued in
the workflow header: `rstan` powers a tier that is gated on
`FRMTMB_BRMS_FIT_TESTS` anyway and would cost tens of minutes of Stan
compilation to assert nothing, and `hBayesDM` is used by nothing that
ships. This is the one place this package deliberately differs from the
spline job, which installs every one of its suggests because its two
headline comparisons depend on them.

### Scope

`git diff --name-only` against the branch point `de9d639`, plus
untracked files:

* changed: `README.md` (one table row, one install line), `_pkgdown.yml`
  (one Extensions dropdown entry), `dev/build-docs.R` (one subsite
  entry). Three files, six inserted lines, one deleted.
* added: `.github/workflows/check-frmtmb-learn.yaml`,
  `dev/learn-findings.md`, and everything under
  `extensions/frmtmb.learn/`.
* `R/` and `tests/` under the repository root: untouched, verified by
  `git diff --name-only de9d639 | grep -E "^(R|tests)/"` returning
  nothing.

## The hBayesDM cross-check: written, not run, and why

The script is `extensions/frmtmb.learn/dev/hbayesdm-crosscheck.R`. It
covers three families against their hBayesDM counterparts on data drawn
from this package's own simulators: `bandit2arm_delta`, `prl_rp` against
`bandit2arm_dual(split = "outcome")`, and `igt_pvl_delta` with the
`tau = 3^cons - 1` reparameterization applied rather than hidden.

It did not run here, and the blocker is in the toolchain rather than in
either package. Each step was checked, not assumed:

1. hBayesDM 2.0.0 (CRAN, dated 2026-09-01) installs cleanly as a binary
   into the lane's private library.
2. It no longer fits through rstan. A fit refuses with "Model fitting
   requires the 'cmdstanr' package"; its DESCRIPTION suggests
   `cmdstanr (>= 0.8.1)` and carries an `Additional_repositories` entry
   for the stan-dev r-universe.
3. cmdstanr installs from that r-universe, and `install_cmdstan()`
   fetches and largely builds CmdStan 2.39.0. `stanc.exe` and
   `make/local` end up in place and `cmdstan_version()` resolves.
4. CmdStan 2.39.0 vendors TBB 2020.3, which does not compile under this
   machine's GCC 14.3.0:

       tbb_2020.3/include/tbb/internal/../atomic.h:17:10: fatal error:
       internal/_deprecated_header_message_guard.h: No such file or directory

   Adding `TBB_INTERFACE_NEW=true` to `make/local` does not avoid it,
   because the vendored copy is what gets built either way. Tried once
   and reverted.

**The 1.x fallback is dead too, and for a more fundamental reason than
time.** hBayesDM 1.2.1 (2022-09-23) is the last rstan-backed release,
but 64 of its 65 Stan programs declare data in the pre-2.33 array
syntax (`int<lower=1, upper=T> Tsubj[N];`), which the parser in this
machine's rstan 2.32.7 with StanHeaders 2.39.1 rejects outright: "Ill-
formed declaration ... It looks like you are trying to use the old array
syntax." Since 1.x precompiles every model at install time, it cannot be
installed here at all, and no amount of build time changes that. So
neither the 2.x route nor the 1.x route is available on this machine,
and the question is closed rather than left open.

What that costs, stated plainly: this package has no check against an
outside implementation of these models. What it has instead is stronger
on every axis except independence of authorship, and that is worth
saying rather than glossing:

* an EXACT identity against Stan programs written from the published
  equations, through rstan, which does work on this machine, for all six
  families, at 0.0e+00 to 5.7e-13;
* a longhand reference per family in the ungated suite, written in the
  opposite style to the engine (one subject and one trial at a time,
  `if` and `[` for selection), which is what caught the one real bug
  this lane shipped and fixed;
* parameter recovery at a realistic scale.

The cross-check would add one thing those do not: confirmation that this
package's READING of each published model agrees with the reading the
field's reference implementation made. It needs a newer CmdStan whose
vendored TBB builds under GCC 14, or a machine with an older GCC. The
script's header records the whole chain so it is a one-step pickup.

### Suite counts, audited by name

| file | passing | failing | skipped |
|---|---|---|---|
| `test-engine.R` | 19 | 0 | 0 |
| `test-families.R` | 36 | 0 | 0 |
| `test-message-uniqueness.R` | 4 | 0 | 0 |
| `test-reference.R` | 6 | 0 | 0 |
| `test-stan-identity.R` | 0 | 0 | 6 |
| `test-surface.R` | 119 | 0 | 0 |
| **total, ungated** | **184** | **0** | **6** |

42 `test_that` blocks. The six skips are the whole of
`test-stan-identity.R`, which is gated on `FRMTMB_BRMS_FIT_TESTS=true`
and `NOT_CRAN=true`; run with both set, that file contributes 24 passing
assertions and 0 failures, and the file count above becomes 208 with 0
skips.

## The compatibility table this package registered

154 rows, through `frmtmb_register_compat()` from `.onLoad()`: 25 per
family for six families, plus 26 for `ts_par7` (it carries an extra
`stage2()` row), plus 3 for the two methods.

| status | rows | what they are |
|---|---|---|
| works | 44 | the grammar (`s()`, `smooth`, `us`, `|ID|`, `prior`), the addition terms, `simulate` for five families |
| refused | 85 | everything that follows from the likelihood not factorizing over rows, plus the missing mean |
| untested | 18 | `autoscale`, `mi()`, `nl`, three per family, each with a reason |
| conditional | 7 | `predict`, where the link scale works and the response scale does not |

The 25 rows are the same 25 for every family, because the FACT is the
same for every family: each is a consequence of the likelihood not
factorizing over rows or of the response being a nominal option code,
and neither depends on which learning rule the family carries. They are
written once in `ln_common_rules()` with the family name filled in at run
time, which is also what keeps the message-uniqueness property true.

Exactly ONE feature differs across the six families:

| feature | bandit2arm_delta | ts_par7 |
|---|---|---|
| `simulate` | works | refused |

`ts_par7` refuses because one two-step trial's draw is three numbers and
a response vector holds one. Everything else about the six is the same,
which is what one engine buys.

`frm_compat("bandit2arm_delta", "importance")` returns `refused` with
the seam named in the note, and `test-surface.R` asserts both the status
and that the note names the finest-factorization design.
