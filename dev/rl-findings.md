# A reinforcement-learning example on the structured-family protocol

Lane notes for the worked example that answers "it'd be neat to see a
reinforcement learning example". The deliverable is a vignette, a
shared family source and a test, not a new exported family.

## The model

Rescorla-Wagner delta learning on a two-armed bandit with a softmax
choice rule, the model hBayesDM calls `bandit2arm_delta` (Ahn, Haines
and Zhang 2017). Per subject, Q starts at 0 for both arms, after each
choice `Q[chosen] += alpha * (reward - Q[chosen])`, and
`P(arm 1) = plogis(beta * (Q1 - Q2))`. `alpha` is a dpar on (0, 1)
through a logit link and `beta` a dpar on (0, Inf) through a log link,
each with its own linear predictor.

No hBayesDM code is copied. hBayesDM is GPL-3; the equations are the
published model. The Stan program used for the identity check is
written here from those equations.

## Where the code lives

`inst/rl/rw-delta.R`, read by both the vignette and the test through
`system.file("rl", "rw-delta.R", package = "frmtmb")`.

A test helper cannot be the shared source, because a vignette has no
reliable path to `tests/testthat/`: `R CMD build` evaluates vignettes
with the working directory inside `vignettes/`, and pkgdown evaluates
them from the package root, so no one relative path reaches the test
directory in both. `inst/` is the one directory both an installed
package and a source build resolve the same way. The vignette does not
duplicate the code: it pulls the same file in with `knitr::read_chunk()`
so every chunk shown is the chunk the test runs.

## The spellings chosen, and why

### `reward(pay1, pay2)`, a registered addition term with arity 2

Registered with `frmtmb_register_aterm("reward", arity = 2L)` inside the
example. `vreal(pay1, pay2)` carries the same two numeric columns and
would have cost nothing to use. The registry is chosen because it is
the seam built for this exact case (the roxygen on
`frmtmb_register_aterm()` says so: a family outside frmtmb gives its
per-row data "the spelling its literature uses, instead of asking users
for `vint()`"), and because a worked example of a contributed family
should show the contributed spelling.

Arity 2 rather than 1 is a modeling decision, not a cosmetic one. The
term carries what EACH arm would have paid on this trial, not the
payoff received. The likelihood only ever reads the chosen arm's entry,
so the two spellings agree on observed data, and data that records the
received outcome alone passes it twice. The second column is what makes
the simulator coherent: a simulated choice needs the payoff of the arm
the subject did not take in the data, and a one-column `reward()` would
have to invent it.

`trial()` was considered as a third registered term and REJECTED. The
core already has `trials()`, and `dev/bracket-sweep.md` records aterm
containers as a CONFIRMED partial-matching hazard with three colliding
pairs already in the tree. Adding `trial` one character from `trials`
puts a fourth one there for the sake of a spelling.

### `subject` and `trial` on the family, not in the formula

`rw_delta(subject = id, trial = trial)`, reaching the frame through the
structure's `frame_vars` slot. This is the protocol's own answer: the
`frame_vars` documentation names "a grouping column, a time column, a
sequence id" as exactly what the slot is for, and `hmm()` spells its
`group =` and `time =` the same way. Putting them in the formula would
work (an aterm value reaches `frame_block` through `av`) but would say
that a sequence id is per-row data, which it is not.

### `primary_dpars = "alpha"`

The main right-hand side of `bf()` is the learning rate's predictor.
The explicit spelling `bf(choice | reward(..) ~ 1, alpha ~ x)` is
REFUSED by the core ("dpar(s) not available for family 'rw_delta':
alpha"), the same way brms refuses `mu ~ x` beside a main formula. So
the example writes `bf(choice | reward(pay1, pay2) ~ condition +
(1 | p | id), beta ~ 1 + (1 | p | id))` and says which predictor is
which.

## Padding, and what a mask means here

The recursion is sequential in trials and vectorized across subjects.
`frame_block` returns a `n_subj` by `n_trial` integer matrix `idx` of
row numbers in trial order, padded on the right to the longest subject
with a repeat of that subject's own first row, and a 0/1 `mask` of the
same shape. The mask multiplies the log-likelihood (a padded cell
contributes nothing) AND the learning rate (a padded cell updates
nothing), so unequal trial counts need no branch on the tape.

Measured waste, 40 subjects with trial counts drawn uniformly on
30 to 100: 2359 rows in a 40 by 100 grid, 41 percent of cells padded,
and the tape build cost 0.09s against 0.08s for the balanced design of
the same shape. The padding is paid for in tape nodes that the mask
multiplies by zero, and at this imbalance it is not measurable against
run-to-run noise. `hmm()` avoids padding entirely by sorting sequences
by decreasing length, which makes the sequences still running at step
`s` a prefix of the order; that is strictly better and strictly less
readable, and the example says so rather than doing it.

`weights()` is a different thing and is REFUSED in `check_spec`. A
trial's factor cannot be reweighted on its own, because its value
depends on every earlier trial of the same subject.

## NA rows: `keep_na = FALSE`, unlike `hmm()`

An unanswered trial produces no prediction error, so `Q` does not move
and dropping the row IS the correct recursion. This is the opposite of
`hmm()`, which sets `keep_na = TRUE` because the chain still
transitions through a time point that emits nothing. The contrast is
worth the vignette paragraph it gets: the protocol flag is a modeling
question, not a plumbing one.

## Timings

RE-MEASURED in the punch round; the first-pass table was single-shot
`system.time()` against Windows' ~15 ms clock granularity and should
not be quoted. Method now: R 4.6.1, Windows 11, one core, warm session,
designs INTERLEAVED (one build of each per round, not all of one
design), `gc(FALSE)` before every measurement, `frm(dry_run =
"objective")` with the `dry_run = "frame"` time subtracted, median of 7
builds, and gradients timed over batches of 50 `obj$gr()` calls with
the best batch of 3 reported.

| subjects | trials | rows | tape build | one gradient |
|---|---|---|---|---|
| 10 | 100 | 1000 | 0.04 s | 10.0 ms |
| 40 | 25 | 1000 | 0.03 s | 8.8 ms |
| 40 | 100 | 4000 | 0.08 s | 39.6 ms |
| 40 | 400 | 16000 | 0.28 s | 153.4 ms |
| 160 | 100 | 16000 | 0.22 s | 160.6 ms |
| 40 | 800 | 32000 | 0.61 s | 296.2 ms |

Both columns track ROWS, which is the tape's node count. The 16000-row
pair is where they part: 40x400 costs more to build than 160x100
because the R loop runs four times as many iterations, and costs the
SAME per gradient (ratio 0.96) because the tape is the same size.

Two things the first pass got wrong here.

The gradient ordering does NOT reverse at fixed rows. The first pass
reported 274 ms for 160x100 against 196 ms for 40x400 and explained it
by the random-effect count; interleaved and batched, the two are
153.4 ms and 160.6 ms, which is the same number. The explanation was
built on noise.

The build gap's SIZE is not reproducible to the precision it was
quoted at. This lane measured 27 percent twice (first pass and the
interleaved re-measurement, 1.27x); the reviewer measured 70 percent on
the same machine with the same interleaving. Two runs a factor of two
apart means the direction is the finding and the percentage is not, so
neither document quotes one any more.

What survives unchanged is the qualitative conclusion, and it is the
part worth keeping: "the tape grows with trials, not rows" is too
strong. The tape's SIZE is a function of rows either way. What
vectorization removes is the R-level operation count, and that is what
the next section prices.

## The elementwise alternative, measured

RE-MEASURED in the punch round. Same data (40 subjects, 100 trials,
4000 rows), same model, three spellings of `loglik`, interleaved,
median of 7 builds, gradients over batches of 100 calls with the best
of 3 batches reported. All three return the objective 1982.6691425636,
identical to 14 significant digits.

| loglik shape | tape build | one gradient |
|---|---|---|
| vectorized over subjects | 0.08 s | 7.3 ms |
| row loop, scalar Q per subject | 0.57 s | 7.2 ms |
| row loop, length-n log-density vector | 0.65 s | 7.2 ms |

CORRECTION 1, and the reason this table was on the punch list. The
first pass reported 50 / 72 / 95 ms in the gradient column and
concluded "the gradient follows at 1.4x and 1.9x". That is wrong. The
three shapes cost the SAME per evaluation, ratios 0.99x and 0.99x here
and 0.66x to 1.20x in the reviewer's independent batches. It has to be
that way: once the tape exists it is a node list and the R code that
built it is gone. The first pass timed one `gr()` call at a time
against a ~15 ms clock, so it was reading the clock, not the gradient.
The penalty is paid ENTIRELY at tape construction.

CORRECTION 2, on the shape of the build claim. The penalty is not a
constant multiplier at any scale. RTMB's `[<-` copies the vector it
writes into, so the sub-assigning form gets relatively worse as `n`
grows. Measured on this family, vectorized against the length-n
sub-assignment, values agreeing to 1e-11 throughout:

| rows | vectorized | elementwise | factor |
|---|---|---|---|
| 1000 | 0.03 s | 0.16 s | 5.3 |
| 4000 | 0.08 s | 0.60 s | 7.5 |
| 10000 | 0.14 s | 1.67 s | 11.9 |
| 20000 | 0.34 s | 4.33 s | 12.7 |

The reviewer measured the same growth further out on the canary's own
500-group Poisson GLMM: 1.3x at n = 1000, 3.2x at 5000, 8.6x at 20000,
26.5x at 50000, 25.1x at 100000.

## The ">1000x" note, corrected properly

The repo records elementwise `a[i] <- x` in a taping loop as ">1000x
slow" (RTMB memory item 6, sourced from a tmb-users thread). Three
things to say about it, and the first-pass version of this section got
one of them and missed two.

FIRST, THE REPO ALREADY KNEW. `dev/hmm-feasibility.md:139` records that
"the `[<-` gotcha is much milder here than the >1000x folklore",
measures roughly 1.4 to 1.5x on the tape build, and adds the
distinction this lane's note omitted: "the >1000x figure is about
assigning into n-length vectors; here the assigned vector is length K".
hmm measured the length-K case; this lane measured the length-n case,
which is the one the folklore is about. A consolidation should MERGE
the two notes, not land a second differently worded correction.

SECOND, NO SINGLE MULTIPLIER IS RIGHT. The factor is a function of `n`,
because sub-assignment copies the vector: about 5x at 1000 rows on this
family, about 13x at 20000, and about 25x by 50000 on the canary's
model, levelling there. So "over 1000x" and "15x" are both the wrong
SHAPE of claim, independently of which is numerically closer. The note
that replaces memory item 6 should say the penalty grows with `n`
because sub-assignment copies the vector, and quote a scale with each
number it gives. Nothing measured at any size the package tests comes
within two orders of magnitude of 1000x.

THIRD, `test-perf.R`'s canary needs no change, and the first-pass
suggestion that it did was misdirected. It quotes no number: it asserts
a tape build under 20 s at n = 100000 and says an elementwise loop
"would blow past the bound by orders of magnitude". At exactly that n
the elementwise build takes 78.5 s against the 20 s bound, so the
canary catches what it exists to catch. Only the memory note is wrong.

## Protocol seams: one gap, two granularities, refused by name

REVISED in the punch round. The first-pass version had the diagnosis
right and the proposed fix one slot too coarse.

THE GAP. `frmtmb_structure(loglik =)` returns one AD scalar for the
whole response. This family's likelihood DOES factorize: given the
parameters, `Q` is a deterministic function of the subject's earlier
data, so each trial contributes one Bernoulli factor and each subject
one product of them. What the likelihood is not is ROWWISE in the
protocol's sense, because a factor cannot be computed from its own
row's dpars and response alone. So the factors exist and the slot has
nowhere to put them.

CONSEQUENCE 1, per group. `frm(importance =)` is refused outright:

    `importance` cannot correct the 'rw_delta' family: it supplies its
    own log-likelihood, which does not factorize over rows, so a
    group's rows have no separable integrand to resample. This is the
    same restriction quadrature has. Use importance = 0

The refusal is correct as an implementation fact and wider than the
mathematics, and the CODE proves it rather than merely allowing it.
`R/importance.R:660-661` computes per-row densities and collapses them
on the very next line with a sparse `ng`-by-`n` indicator
(`smat`, built at `R/importance.R:348`), and `imp_group_map()` at
`R/importance.R:333-343` already REFUSES any design where a row reaches
more than one group. Nothing downstream of that product ever sees a row
again: `fn()` does a per-group log-mean-exp over draws and `imp_ess()`
reads the same `ng`-by-`nd` matrix. So group separability is a
precondition the correction ENFORCES today, not a property a new slot
would have to introduce, and a per-group slot drops straight into those
two lines. `quadrature = TRUE` is a genuine refusal by contrast: it
integrates one scalar random effect against a product of per-ROW
densities.

CONSEQUENCE 2, per row. No pointwise log-likelihood matrix, so `loo()`
and `waic()` cannot work, and deviance residuals have no per-row
saturated fit to compare against (`R/predict.R:2123` defines the unit
deviance as `2 * (saturated loglik - fitted loglik)`, which is
irreducibly per row). In the CORE suite the refusal a user meets for
`loo()` is one step earlier and unrelated ("loo() is a posterior
quantity and this is a maximum-likelihood fit"), because the elpd
machinery lives in `frmtmb.sample`.

THE SEAM, corrected. One slot is right; keying it on `unit` is not.

The first-pass proposal was "a `loglik` allowed to return a vector over
the structure's own `unit`", and claimed that one addition unlocks
importance, `loo()`, `waic()` AND deviance residuals. The last does not
follow, and this family is the counterexample to its own proposal: its
declared `unit` is "one subject's trial sequence", but its FINEST
factorization is the ROW. `rw_fitted()` already computes the per-trial
conditional choice probability, and `test-rl-example.R:124` asserts
that `sum(dbinom(y, 1, fitted(fit), log = TRUE))` equals the reference
data log-likelihood. A per-subject value could never produce a deviance
residual, so a slot keyed on `unit` would have the protocol promising a
per-row quantity out of per-group data.

`unit` and the factorization granularity are two different questions:
`unit` declares what may honestly be LEFT OUT (the leave-one-out unit),
and the slot should carry the finest factorization the family HAS.

The design that follows:

- the slot returns the finest factorization available, per row where
  one exists and per group otherwise, and says which;
- `unit` stays exactly as it is, the declaration of the leave-one-out
  unit;
- core sums rows into groups with the `smat` it already builds for the
  importance correction, or takes group values directly;
- `loo()` and `waic()` key their columns on `unit` and keep refusing
  where the unit is a group;
- deviance residuals require row granularity and refuse otherwise.

That is one slot plus one existing declaration, no more machinery than
the first-pass proposal, and it stops the protocol over-promising.

Two implementation facts belong in the design note.

STACKING. The correction evaluates at `ridx <- rep.int(seq_len(n), nd)`,
the whole design repeated once per draw, so a per-group slot must
return `ng` by `nd` and a sequential family must run its recursion over
subject-crossed-with-draw rather than subject. For `rw_delta` that is
benign and worth saying out loud: the loop is already vectorized across
subjects, so widening it to subject-by-draw is the same loop over a
longer vector. For a family not already vectorized, that is where the
cost lands.

GROUPING ALIGNMENT. The family's unit and the correction's group are
different objects that merely coincide here. `rw_delta(subject = id)`
with `(1 | p | id)` gives one family unit per importance group, but the
same family under `(1 | item)` would give a per-subject log-likelihood
against a per-item proposal, and summing the first into the second is
not defined. Core has what it needs to check this (`row_level` from
`imp_group_map()`), and the check must be required, because getting it
wrong is silent.

NOT IMPLEMENTED, and refused by name rather than worked around. It is a
change under `R/` (structure.R, objective.R, importance.R, loo.R,
predict.R) that this lane has no mandate for and that sibling lanes are
editing this round. What the example does instead: it refuses
`residuals(type = "deviance")` in its own words, sets `unit` to "one
subject's trial sequence", shows the `importance` refusal in the
vignette as output rather than as prose, and measures the Laplace
caveat by simulation instead.

### Two core sharp edges, for the consolidation, not this lane

`structure_unit()` (`R/structure.R:669`) is dead code in core: nothing
in `R/` calls it, and its only reader is `test-structure.R:47`. The
first-pass version of this document said its "only reader" was
`extensions/frmtmb.sample/R/loo.R:72`, which inlines the same
expression. That inlining is deliberate and commented at `loo.R:64-67`
(an out-of-tree extension must not reach for a `@noRd` internal), so
the actionable fact is narrower than the first pass implied: core's
accessor should be exported or deleted.

`frmtmb_register_aterm()` silently accepts a re-registration of the
same name at a DIFFERENT arity, clobbering the first entry. The
documented no-op for re-registration (R/parse.R:42-44) covers a package
reload; a changed arity is not that, and a stale arity-1 `reward` would
win silently.

## What the family answers

Implemented: `fitted()` and `predict(type = "response")` on the
training data (the per-trial choice probability, the natural fitted
value), pearson residuals through `fitted_var` (the binomial variance
of that probability), response residuals, and one structured simulator
serving `simulate()`, `posterior_predict()` and `frm_simulate()`.

Refused by name, each with its own sentence: `newdata` on the response
scale, `conditional_effects()`, `residuals(type = "osa")`,
`residuals(type = "deviance")`. Left at the protocol's conservative
default with a generic message: REML, quadrature, profiling,
multivariate, `cens()`/`trunc()`, `mi()`, `re.form`, cluster-robust
variance.

## Validation

### 1. Identity against a Stan program

`tests/testthat/helper-rl.R` carries a 30-line Stan program for the
same model, written from the published equations. The trick that keeps
it a clean identity: the Stan program declares exactly the parameters
frmtmb estimates, on the scales frmtmb estimates them on, including the
subject deviations themselves rather than a standardized version of
them. The covariance factor rides in as DATA. So every parameter is
unbounded, the map from frmtmb's estimates is the identity, and the
comparison carries no Jacobian at all, where the brms tier
(`brms_lp_check()`) has to accumulate one.

Measured, 30 subjects by 100 trials:

RE-DERIVED in the punch round. The first pass recorded -392.4486089
against the label "30 subjects by 100 trials". That value belongs to no
design shipped here: 3000 Bernoulli trials at this bandit run about
-0.487 per trial, and -392 corresponds to roughly 850 rows, so the
label and the number came from different runs. The residual, which is
the actual claim, was never in doubt and reproduces. Current values, 30
subjects by 100 trials, seed 5, `n = 3000`, `logLik = -1461.882`:

| quantity | at the estimates | at a displaced point |
|---|---|---|
| Stan `log_prob` | -1427.8453940 | -1460.3042080 |
| frmtmb joint log density | -1427.8453940 | -1460.3042080 |
| residual | -4.547e-13 | -2.501e-12 |
| max abs gradient on the 60 subject effects | 1.199e-14 | not asserted |
| max abs gradient overall | 5.211 | 116.7 |

Independently confirmed: the reviewer measured the same five numbers.

The residual is expected to be exactly zero, not merely small, and the
reason is checkable rather than hopeful. `fit$frame[["priors"]]` is
empty (core `frm()` carries no default prior; the sampling default-prior
machinery left with `frmtmb.sample`), the Stan program uses `target +=`
with `_lpmf` and `_lpdf` so it keeps every normalizing constant, and
every declared parameter is unbounded so `adjust_transform = FALSE` has
no Jacobian to drop.

The gradient rows are the point of the second check. frmtmb puts the
subject effects at their conditional modes, so Stan's gradient there
must vanish; the outer gradient is NOT zero at the marginal optimum and
is correctly not asserted on. The second column is a deliberately
non-stationary point (fixed effects and subject effects displaced by
0.25 cos(k)), so the agreement cannot be an artifact of both sides
sitting at a stationary point.

The compile is cached through `brms_stan_model()` from
`helper-brms.R`, which needs rstan but not brms. The gate is this
lane's own `skip_unless_rl_stan()`, on `FRMTMB_BRMS_FIT_TESTS` and
`NOT_CRAN`, matching the brms fit tier. hBayesDM was NOT installed or
consulted.

### 2. Parameter recovery

100 datasets per design, simulated through the family's own simulator
(`frm_simulate()` draws fresh subject effects from `rl_truth` and calls
`sim_ctx`), each refitted. Coverage is over the replicates that
produced a finite Wald interval; `n_ci` counts them.

40 subjects, 100 trials:

| parameter | truth | bias | MC se | sd of est | coverage | n_ci |
|---|---|---|---|---|---|---|
| `alpha_(Intercept)` | -0.619 | +0.0087 | 0.0178 | 0.178 | 0.940 | 100 |
| `alpha_conditiontrt` | 0.800 | -0.0181 | 0.0248 | 0.248 | 0.940 | 100 |
| `beta_(Intercept)` | 1.099 | +0.0021 | 0.0070 | 0.070 | 0.960 | 100 |
| `log sd(alpha)` | -0.693 | -0.2107 | 0.0573 | 0.573 | 0.929 | 98 |
| `log sd(beta)` | -0.916 | -0.0619 | 0.0167 | 0.167 | 0.950 | 100 |

Zero fits failed. The three fixed effects are unbiased to within Monte
Carlo error and cover at the nominal rate. The learning rate's
subject-level standard deviation is biased low.

The unstructured correlation is the least well behaved quantity in the
model: truth 0.30, mean estimate 0.482, sd across replicates 0.337.
Binary data with 40 levels does not identify it well.

### 3. The Laplace caveat, measured

`frm(importance =)` is refused (see the seam section), so the caveat is
measured by repeating the recovery study at 20 trials per subject and
reading what moves.

| parameter | bias @100 | cov @100 | bias @20 | cov @20 | n_ci @20 |
|---|---|---|---|---|---|
| `alpha_(Intercept)` | +0.009 | 0.94 | +0.064 | 0.949 | 97 |
| `alpha_conditiontrt` | -0.018 | 0.94 | -0.005 | 0.949 | 98 |
| `beta_(Intercept)` | +0.002 | 0.96 | -0.000 | 0.917 | 96 |
| `log sd(alpha)` | -0.211 | 0.93 | -0.628 | 0.726 | 84 |
| `log sd(beta)` | -0.062 | 0.95 | -0.247 | 0.885 | 96 |

The fixed effects hold up at 20 trials: bias within Monte Carlo error,
coverage near nominal. The variance components do not. `log sd(alpha)`
loses 0.63, which is a standard deviation about half the true one, its
coverage falls to 0.73, and 16 of 100 fits produce no usable interval
for it at all.

HONEST LIMIT of this measurement. It does not separate Laplace error
from the small-sample maximum-likelihood bias a variance component
carries anyway with binary data and 40 levels. Separating them is
precisely what `frm(importance =)` is for, and this family cannot ask
it. So the vignette states the operational conclusion instead: report
the fixed effects from short sessions, and treat a subject-level
standard deviation from 20 binary trials as a lower bound.

### 4. fitted() and residuals: implemented, not refused

`fitted_mean` returns the per-trial choice probability, replayed
numerically at the estimates, and `fitted_var` its binomial variance.
That makes `fitted()`, `predict(type = "response")` on the training
data, response residuals and pearson residuals all work. The test
asserts each of them against the independent scalar reference:
`sum(dbinom(y, 1, fitted(fit), log = TRUE))` equals the reference's
data log-likelihood to 1e-8, and the pearson residuals equal
`(y - p) / sqrt(p * (1 - p))`.

`residuals(type = "osa")` and `residuals(type = "deviance")` are
refused by name with their own sentences.

## Files

| file | what |
|---|---|
| `inst/rl/rw-delta.R` | the family, read by the vignette and the test |
| `vignettes/reinforcement-learning.Rmd` | the tutorial |
| `tests/testthat/helper-rl.R` | loads the family, holds the Stan program and gate |
| `tests/testthat/test-rl-example.R` | 8 tests, 2 of them the Stan tier |
| `_pkgdown.yml` | one line under the Guides articles |
| `NEWS.md` | one bullet |

Nothing under `R/` changed.

## Left out, and why

**hBayesDM as a third implementation.** The brief allowed installing it
into the private library as a Suggests-level cross-check. Not done, and
not because of the GPL: an identity check against it is not available.
hBayesDM's program carries hard-coded priors as `~` statements and a
non-centered hierarchical parameterization, so its `log_prob` is not
the quantity frmtmb maximizes and cannot be made into it without
editing its source. What would remain is an estimate-against-estimate
comparison between two different inferential procedures, which is a
much weaker statement than the exact identity the purpose-written Stan
program already gives (residual -4.5e-13). Two independent
implementations are checked against the tape as it is: a scalar R
reference in the test suite, and the Stan program.

**`check_laplace()` from frmtmb.sample.** It would sample the joint
with tmbstan and compare against the mode plus Wald errors, which is a
genuine Laplace check. It costs a dependency the vignette does not
otherwise need and minutes of NUTS per design, and it answers a
narrower question than the recovery study does (one dataset rather than
the sampling distribution). The recovery-at-20-trials measurement is
what the vignette carries instead.

**The prefix-ordered block.** `hmm()`'s no-padding layout is described
and priced in the vignette, not implemented, because the measurement
says the padding costs nothing at a realistic imbalance and the mask is
what a reader can follow.

**A `pointwise_loglik` slot.** The seam section above. Refused by name.

## Verification run

First pass, corrected in the punch round: the gate-off row said 82,
which was the gate-ON count. The counts below are after the punch
round, which added one test (three assertions) for the one-trial-per-
subject block.

| check | result |
|---|---|
| `test-rl-example.R`, gate off | PASS 79, FAIL 0, WARN 0, SKIP 2 (Stan tier) |
| `test-rl-example.R`, gate on | PASS 85, FAIL 0, WARN 0, SKIP 0 |
| `test-message-uniqueness.R`, own process | PASS 6, FAIL 0 |
| `test-bracket-access.R`, own process | PASS 8, FAIL 0 |
| vignette knit with pandoc on PATH | 8.9 s, no warnings |
| `R CMD check --as-cran` | Status: OK, zero errors, warnings and notes |

Inside that check the suite reports `FAIL 0 | WARN 0 | SKIP 176 | PASS
4029`. Three of this lane's tests are among the skips and all three are
correct: the recovery test carries `skip_on_cran()` and the run had
`NOT_CRAN=false`, and the two Stan tests are behind the opt-in gate. The
other six ran.

ONE TRAP WORTH RECORDING. A confirmation run of the test file reported
every test skipped with "inst/rl/rw-delta.R is not installed", and the
cause was the shell, not the code: the private library had been spelled
`/c/Users/...` (MSYS style) in that one command instead of
`C:/Users/...`, so R silently dropped it from `.libPaths()`, loaded the
user library's older frmtmb, and `system.file()` returned "". A
`skip_if_not()` guard turns a wrong library path into a green run with
nothing in it. Any run of these tests should be read together with its
skip count.

`R CMD check --as-cran --no-manual` with `_R_CHECK_CRAN_INCOMING_=false`
on the built tarball, private library, `NOT_CRAN=false`. Tests 200s,
vignette rebuild 189s, both OK. The brief expected the V8 NOTE; this
run produced NO note at all, so either the note is incoming-check only
or it depends on a package this library does not carry. Nothing in this
lane's files suppressed it.

The suite run inside that check is the one that matters for a side
effect this lane has: `helper-rl.R` registers the `reward()` addition
term when testthat loads helpers, which is BEFORE every test file, so
`test-compat.R` and `test-compat-register.R` see one extra entry in the
compatibility vocabulary. They pass, because they test membership and
before-and-after identity rather than an exact feature set. A future
test that pins the vocabulary exactly would have to account for it.

## Punch round, 2026-09-05

Against the review at `dev/review-rl.md` (GO-WITH-FIXES, 8 items). Every
correction is folded into the sections above; this is the record of what
changed and why. Three of the eight were published numbers that did not
reproduce, and all three came from one cause: single-shot
`system.time()` against Windows' ~15 ms clock granularity, plus one
value never refreshed after a design change.

### The measurement method, changed once for all three tables

Every timing in this document and in the vignette is now: warm session,
designs and shapes INTERLEAVED (one measurement of each per round, not
all of one then all of the next), `gc(FALSE)` before each measurement,
median of 7 tape builds, and gradients timed over a BATCH of 50 or 100
`obj$gr()` calls with the best batch of 3 reported. Script:
`rl-timing2.R` in the lane scratchpad.

| # | verdict | what changed |
|---|---|---|
| P1 | fixed | The three loglik shapes' gradient column was 50 / 72 / 95 ms and read as "the gradient follows at 1.4x and 1.9x". Batched, the three are 7.3 / 7.2 / 7.2 ms, ratios 0.99. The column is now reported as flat and the prose says the penalty is paid entirely at tape construction. Build column re-measured to 0.08 / 0.57 / 0.65 s. Both the vignette table and the findings table. |
| P2 | fixed | The scaling table's "the same pair reverses for the gradient" does not reproduce: 153.4 ms against 160.6 ms at 16000 rows is the same number, ratio 0.96, and the random-effect-count explanation was built on noise. Removed. The "27 percent" build gap is NOT reproducible to that precision either: this lane measures 1.27x on both its passes and the reviewer measures 1.70x on the same machine, so both documents now give the direction and say the magnitude is unstable, quoting both figures. |
| P3 | fixed | The identity value -392.4486089 labelled "30 subjects by 100 trials" was from a different design. Re-derived at the shipped fixture: -1427.8453940 at the estimates and -1460.3042080 displaced, residuals -4.5e-13 and -2.5e-12, inner gradient 1.199e-14, overall 5.211. Matches the reviewer's five numbers independently. Also recorded WHY the constant is exactly zero, now checked rather than assumed: `fit$frame[["priors"]]` is empty. The stale value was in the vignette too and is fixed there. |
| P4 | fixed | The ">1000x" correction now cites `dev/hmm-feasibility.md:139`, which made it first and recorded the length-K against length-n distinction this lane's note had omitted. It no longer quotes a single multiplier: the factor grows with `n` because sub-assignment copies the vector, measured here at 5.3x (1000 rows) to 12.7x (20000) on this family and by the reviewer at 1.3x to 26.5x on the canary's own model. The suggestion to re-measure `test-perf.R`'s canary is withdrawn: it quotes no number, and at its own n = 100000 the elementwise build takes 78.5 s against a 20 s bound, so it works. |
| P5 | fixed | The seam proposal is rewritten. Per-group is right for the importance correction and the code proves it (`R/importance.R:660-661` collapses per-row densities with `smat`; `imp_group_map()` already enforces group separability), but one slot keyed on `unit` is too coarse: deviance residuals need a per-row saturated comparison (`R/predict.R:2123`), and `rw_delta` already computes a per-row conditional log-likelihood for `fitted()`, so its finest factorization is the ROW while its `unit` describes the leave-one-out granularity. The slot now carries the finest factorization and `unit` stays the separate declaration. Added the `n * nd` stacking requirement and the grouping-alignment check. Corrected the `structure_unit()` claim: `frmtmb.sample`'s inlining is deliberate and commented, so the actionable fact is only that core's accessor is dead code. |
| P6 | fixed | Verification table said 82 assertions gate-off; that was the gate-on count. Re-run and re-recorded. |
| P7 | fixed | The vignette's `"../inst/rl/rw-delta.R"` fallback resolved only from `vignettes/`, contradicting the prose beside it. Removed; the chunk now uses `system.file(..., mustWork = TRUE)`, so a missing install fails loudly instead of silently reading a path that may not be there. |
| P8 | recorded, not fixed | Both are core, not this lane: `structure_unit()` (`R/structure.R:669`) is dead code, and `frmtmb_register_aterm()` silently accepts a re-registration at a different arity. Written up under "Two core sharp edges". |

### The one code defect

`inst/rl/rw-delta.R:67-74`, the reviewer's fix, kept as written. When
every subject has one trial, `vapply()` returns a plain vector rather
than a 1-by-`n_subj` matrix, `t()` of a vector is 1-by-`n_subj`, and
the block came back transposed: six subjects were read as one learner
with six trials and Q was carried across subject boundaries. `frm()`
fitted it without complaint, so the answer was silently wrong. The
`matrix(..., nrow = nt)` before each transpose is a no-op reshape when
`nt > 1` and is what fixes `nt == 1`.

New test at `tests/testthat/test-rl-example.R:83`. It asserts the block
is 6 by 1 rather than 1 by 6, and pins the likelihood at
`6 * log(0.5) = -4.158883083`, which is the value for ANY parameter
setting because every subject's first trial starts from Q1 = Q2 = 0.
The transposed block gave -4.356251586. A test that only checked the
dimensions would pass on a block that was right by accident; this one
fails on the arithmetic.
