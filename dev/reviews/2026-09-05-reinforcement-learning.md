# Review: the reinforcement-learning lane (wt-rl)

Reviewer's notes, written incrementally while reproducing. Lane branch
`wt-rl` at `C:/Users/adf44/source/r/frmtmb-wt-rl`, branched at 316a28b.
Every number below is measured on this machine unless it is quoted from
`dev/rl-findings.md`.

## Scope, confirmed

`git diff --name-only 316a28b` plus untracked is exactly:

| file | state |
|---|---|
| `NEWS.md` | modified, one bullet under a new development heading |
| `_pkgdown.yml` | modified, one line |
| `dev/rl-findings.md` | new |
| `inst/rl/rw-delta.R` | new, 335 lines |
| `tests/testthat/helper-rl.R` | new |
| `tests/testthat/test-rl-example.R` | new |
| `vignettes/reinforcement-learning.Rmd` | new |

`git diff --name-only 316a28b -- R/` is EMPTY. The claim "nothing under
`R/` changed" holds.

The main checkout was already dirty and at 9a1d3f3 when this review
started (a sibling lane is working there). Its `git status --porcelain`
was recorded before any work and is re-checked at the end.

## Environment

Private library `.../scratchpad/rrl-lib`, holding the worktree's core
(frmtmb 0.51.0) and `frmtmb.latent` built from the worktree. The user
library is a read-only fallback for rstan and StanHeaders. Every run
below loads frmtmb from `rrl-lib` and prints the resolved path, because
the lane's own notes record a run that silently loaded the wrong
frmtmb and skipped everything green.

## 1. The identity against Stan

### It reproduces, exactly, at both points

Run from the cached program with the worktree core, 30 subjects by 100
trials, seed 5:

| quantity | at the estimates | at the non-stationary point |
|---|---|---|
| Stan `log_prob` | -1427.8453940 | -1460.3042080 |
| frmtmb joint | -1427.8453940 | -1460.3042080 |
| residual | -4.547e-13 | -2.501e-12 |
| inner grad, max abs over the 60 subject effects | 1.199e-14 | not asserted |
| overall max abs grad | 5.211 | not asserted |

This is the lane's central validation claim and it holds completely.
The constant is zero to rounding, and nothing leaks into the identity:
`frame[["priors"]]` is empty (core `frm()` carries no default prior,
the sampling default-prior machinery having left with
`frmtmb.sample`), the Stan program uses `target +=` with `_lpmf` and
`_lpdf` so it keeps every normalizing constant, and every declared
parameter is unbounded so `adjust_transform = FALSE` has no Jacobian to
drop. The inner gradient vanishing at 1.2e-14 while the overall
gradient is 5.2 is exactly the right pair of facts: frmtmb puts the
subject effects at their conditional modes, and the outer gradient at a
marginal optimum is not the joint gradient. The test asserts only on
the inner block, which is correct.

### The Stan program is the same model

Read against `rw_loglik()` line by line, `helper-rl.R:42-88` matches on
every point: `q1 = q2 = rep_vector(0, S)` (Q starts at 0), the
chosen-arm-only update `q1[s] += alpha[i] * (pay1[i] - q1[s])` under
`if (choice[i] == 1)`, `e = beta[i] * (q1[s] - q2[s])` scored with
`bernoulli_logit_lpmf(choice[i] | e)` (softmax with beta multiplying
the Q difference), `alpha = inv_logit(eta_a)` and `beta = exp(eta_b)`
(logit and log links), and `multi_normal_cholesky_lpdf(u[s] | 0, L)`
with the covariance factor arriving as DATA. The padding is handled the
other way round from frmtmb (`if (mask[s, t] == 1)` skips where frmtmb
multiplies by the mask), which makes the agreement a stronger check
than a transcription would be.

The reshape the map depends on is asserted rather than assumed:
`test-rl-example.R:200` checks `rl_stan_pars(par)$u` equals
`ranef(fit)[[1]]`, so the level-major `b` ordering (RTMB gotcha 5)
cannot silently permute. Reproduced: TRUE.

### DEFECT 2: the recorded log_prob is stale

`dev/rl-findings.md` reports the identity as -392.4486089 at "30
subjects by 100 trials". At that design I measure -1427.8453940. The
number is not merely different from mine, it is implausible for the
design: -392.45 over 3000 Bernoulli trials is -0.13 per trial, or about
88 percent accuracy, and this bandit at the simulated parameters
achieves 76 percent (measured per-trial log-likelihood -0.472).
Sweeping designs, -392 corresponds to roughly 850 rows, not 3000:

| design | rows | joint |
|---|---|---|
| 30 x 100 | 3000 | -1427.85 |
| 20 x 60 | 1200 | -377.50 |
| 10 x 100 | 1000 | -472.70 |
| 30 x 30 | 900 | -284.70 |

So the findings table quotes a value from some earlier configuration
(20 x 60 is the fixture default and is the closest) against a design
label that no longer matches. The SCIENCE is unaffected - the residual,
which is the actual claim, reproduces at 1e-13 - but the table is wrong
and would mislead anyone trying to reproduce it. Punch list item.

## 2. The family

### The recursion is correct

`rw_loglik()` (inst/rl/rw-delta.R:88) walks one iteration per TRIAL,
each a handful of vector operations over all subjects at once. Read
against the stated model it is right on every point:

- `q1 <- rep(0, nrow(idx)); q2 <- q1` (:99): a subject's first trial
  starts from Q = 0 for both arms.
- `ll + sum(m * (c1 * eta - logspace_add(0 * eta, eta)))` (:109) is
  `log P(choice)` for `P(arm 1) = plogis(eta)`, `eta = beta * (Q1 - Q2)`.
- `a <- alpha[i] * m` (:112) then `q1 + (a * c1) * (pay1[i] - q1)` and
  `q2 + (a * (1 - c1)) * (pay2[i] - q2)`: only the chosen arm updates,
  and the mask zeroes both the contribution and the update, so a padded
  cell contributes nothing and learns nothing with no branch on the
  tape. This is the claim the lane makes and it holds.

`rw_replay()` (:125) is the same recursion off the tape and shares the
update lines verbatim, so the simulator draws from the process the
likelihood scores. One detail that had to be right and is: at a padded
cell `y[i[keep]] <- c1[keep]` (:143) declines to write, so a pad does
not overwrite the subject's real first-trial choice that `idx` repeats.

### DEFECT 1: the block transposes when every subject has one trial

`rw_block()` (inst/rl/rw-delta.R:67-70) builds `idx` and `mask` with
`t(vapply(..., integer(nt)))`. When `nt == 1`, `vapply` returns a plain
vector rather than an `nt`-by-`n_subj` matrix, and `t()` of a vector is
a 1-by-`n_subj` matrix. So `idx` comes back 1 by `n_subj` instead of
`n_subj` by 1, while the `n_subj` and `n_trial` slots beside it stay
correct.

Reproduced, six subjects with one trial each:

    block dim  : 1 x 6      (expected 6 x 1)
    rw_loglik  : -4.356251586
    correct    : -4.158883083   ( = 6 * log(0.5) )

The recursion then treats six subjects as ONE subject with six trials
and carries Q across subject boundaries. `frm()` fits it without
complaint, so the failure is silent and the answer is wrong, not an
error. The trigger is exactly `max(trial count) == 1`; a single subject
with many trials is fine, and so is any ragged design, because `vapply`
returns a matrix as soon as `nt > 1`.

Impact is low (a delta-rule dataset with one trial per subject is
degenerate) but the family validates other degenerate inputs loudly,
and this one is a one-line fix. Punch list item.

### The reward aterm's registration is sane

`frmtmb_register_aterm()` is documented (R/parse.R:42-44) as replacing
an earlier entry so that reloading a contributing package is not an
error, and the registry is a package-level environment, not persisted.
Reproduced on a fresh session: re-registering `reward` at arity 2 is a
quiet no-op, and re-sourcing the whole family file is fine. So the
vignette registers nothing that outlives the session, and the
`aterms[["reward1"]]` / `[["reward2"]]` spelling the family reads is
exactly the arity-above-one convention the roxygen documents.

One sharp edge, in core rather than in the lane: registering the same
name at a DIFFERENT arity is also silently accepted and clobbers the
first entry. Not this lane's to fix; noted because a stale arity-1
`reward` would silently win.

## 3. Timings

### The three loglik shapes reproduce on the tape build, NOT on the gradient

40 subjects x 100 trials, 4000 rows, all three returning the same
objective (1967.6888631596 for my seed, agreeing to 10 significant
digits across shapes, which is the lane's check and it passes).

| loglik shape | lane tape | mine tape | lane grad | mine grad |
|---|---|---|---|---|
| vectorized over subjects | 0.10 s | 0.09 s | 50.0 ms | ~12 ms |
| row loop, scalar Q | 0.73 s | 1.00 s | 72.0 ms | ~12 ms |
| row loop, length-n vector | 1.47 s | 1.29 s | 95.0 ms | ~12 ms |

The tape-build column reproduces: the vectorized form is 11x cheaper
than the scalar row loop and 14x cheaper than the sub-assigning one,
against the lane's 7x and 15x. The lane's headline "worst 15x on the
build" stands.

The GRADIENT column does not reproduce, and the lane's claim that "the
gradient follows at 1.4x and 1.9x" is WRONG. Measured with 100 `gr()`
calls per batch and three batches:

    round 1:  vec 11.70   scalarQ 11.00 (0.94x)   len-n 10.60 (0.91x)
    round 2:  vec 13.30   scalarQ 15.90 (1.20x)   len-n 13.00 (0.98x)
    round 3:  vec 16.70   scalarQ 11.00 (0.66x)   len-n 11.40 (0.68x)

and `fn()` alone over 200 calls gives ratios 0.91x to 1.07x. The three
shapes have the SAME run-time cost, which is what they must have: once
taped, the tape is a node list and the R code that built it is gone.
The lane's 50/72/95 ms are single-shot `system.time()` means at
Windows' ~15 ms clock granularity, which is the same artifact I hit on
my first pass (I measured 18/13/15 ms before batching). The correct
statement is that the elementwise penalty is paid ENTIRELY at tape
construction and not at all per evaluation.

### The "tape grows with trials, not rows" re-measurement

Interleaved, gc between, 7 build reps and 3 gradient batches:

| design | rows | tape (median) | grad (min) |
|---|---|---|---|
| 40 x 400 | 16000 | 0.56 s | 45.5 ms |
| 160 x 100 | 16000 | 0.33 s | 36.0 ms |
| 40 x 800 | 32000 | 1.17 s | 87.0 ms |

The lane's DIRECTION is right and its magnitude is understated: at
fixed rows, quadrupling the loop iterations costs 70 percent more tape
build here, against the lane's 27 percent. Its secondary claim that the
gradient reverses the ordering (274 ms for 160x100 against 196 ms for
40x400, blamed on 320 random effects against 80) does NOT reproduce: I
get 40x400 slower on the gradient too, 1.26x, so loop count costs on
both and the random-effect count does not dominate. On a first,
non-interleaved pass I got 32000 rows apparently CHEAPER than 16000,
which is how noisy this measurement is on this machine; the lane's
single-pass table should not be quoted as reproducible at the
two-significant-digit level.

The qualitative conclusion the lane draws is nonetheless correct and is
the part worth keeping: the tape's SIZE is a function of rows, so
"the tape grows with trials, not rows" is too strong, and what
vectorization removes is the R-level operation count.

### The ">1000x" note: the number is wrong, but so is having a number

Memory item 6 says elementwise `a[i] <- x` in taping loops is ">1000x
slow (tmb-users benchmark thread); observation-length loops are
forbidden in objective.R/covstruct.R, enforced by
tests/testthat/test-perf.R canary."

Two corrections, one of them the lane's and one the lane missed.

FIRST, the repo already knew. `dev/hmm-feasibility.md:139` records
"The `[<-` gotcha is much milder here than the >1000x folklore ...
roughly 1.4-1.5x on the tape build", and adds the nuance the RL lane's
note omits: "The >1000x figure is about assigning into n-length
vectors; here the assigned vector is length K". So hmm measured the
length-K case and RL measured the length-n case, and RL's is the one
that actually tests the claim. The RL lane presents its finding as a
fresh CORRECTION without citing the earlier one. Consolidation should
merge them rather than land a second, differently-worded correction.

SECOND, and this is the substantive finding: the penalty is not a
constant factor at all, so no single multiplier is right. Measured on
the CANARY's own model (the 500-group Poisson GLMM of
`tests/testthat/test-perf.R`), vectorized against a length-n
sub-assignment loop, values agreeing to 1e-12 throughout:

|  n | vectorized | elementwise | factor |
|---|---|---|---|
| 1000 | 0.06 s | 0.08 s | 1.3x |
| 5000 | 0.14 s | 0.45 s | 3.2x |
| 20000 | 0.49 s | 4.20 s | 8.6x |
| 50000 | 1.12 s | 29.70 s | 26.5x |
| 100000 | 3.13 s | 78.53 s | 25.1x |

The elementwise column grows superlinearly, because each `lp[i] <- ...`
copies the whole advector. So the FACTOR is a function of n, not a
constant: 1.3x at n = 1000, rising to about 25x by n = 50000 and
levelling there. That is why both "1000x" and "15x" are the wrong SHAPE
of claim, independently of which number is closer. The lane's 15x is a
true measurement at 4000 rows and its "at this scale" hedge is honest,
but the note that replaces memory item 6 should say "grows with n
because sub-assignment copies the vector; about 25x at n = 100000" and
not quote a single multiplier. Nothing at any size up to the largest
the package tests comes within two orders of magnitude of 1000x.

That also means the lane's suggestion to "re-measure the canary in
test-perf.R before quoting the number again" is misdirected: the canary
quotes no number. It asserts tape build under 20 s at n = 100000 and
says an elementwise loop "would blow past the bound by orders of
magnitude". At exactly that n the elementwise build takes 78.5 s
against a 20 s bound, so the canary does catch the defect it is there
to catch. It is sound and needs no change; only the memory note does.

## 4. Recovery, the Laplace caveat, and the pointwise seam

### Recovery reproduces

My own seeds (20260906..20260955), 50 replicates per design, the lane's
own simulator and truth. Zero fit failures at either design.

40 subjects x 20 trials:

| parameter | truth | bias | MC se | sd(est) | coverage | n_ci |
|---|---|---|---|---|---|---|
| `alpha_(Intercept)` | -0.619 | -0.001 | 0.044 | 0.312 | 0.958 | 48 |
| `alpha_conditiontrt` | 0.800 | +0.006 | 0.073 | 0.512 | 0.918 | 49 |
| `beta_(Intercept)` | 1.099 | -0.003 | 0.018 | 0.127 | 0.936 | 47 |
| `log sd(alpha)` | -0.693 | -0.764 | 0.145 | 1.023 | 0.705 | 44 |

40 subjects x 100 trials:

| parameter | truth | bias | MC se | sd(est) | coverage | n_ci |
|---|---|---|---|---|---|---|
| `alpha_(Intercept)` | -0.619 | +0.065 | 0.028 | 0.198 | 0.940 | 50 |
| `alpha_conditiontrt` | 0.800 | -0.054 | 0.042 | 0.294 | 0.900 | 50 |
| `beta_(Intercept)` | 1.099 | +0.005 | 0.011 | 0.077 | 0.940 | 50 |
| `log sd(alpha)` | -0.693 | -0.196 | 0.050 | 0.350 | 0.979 | 47 |

Every qualitative claim the lane makes holds on independent seeds:

- the three fixed effects are unbiased to within Monte Carlo error at
  both designs and cover near the nominal rate (0.90 to 0.96);
- `log sd(alpha)` is biased low at 100 trials, -0.196 against the
  lane's -0.211, which is close agreement;
- it collapses at 20 trials, -0.764 against the lane's -0.628, with
  coverage 0.705 against their 0.726, and 6 of 50 fits (12 percent)
  give no usable interval against their 16 of 100.

So the operational conclusion the vignette states - report fixed
effects from short sessions, treat a subject-level standard deviation
from 20 binary trials as a lower bound - is supported. The lane's
refusal to attribute all of that to Laplace error rather than
small-sample ML bias is the right call and is stated honestly.

### The importance refusal, quoted live

Reproduced verbatim, exactly as the lane records it:

    `importance` cannot correct the 'rw_delta' family: it supplies its
    own log-likelihood, which does not factorize over rows, so a
    group's rows have no separable integrand to resample. This is the
    same restriction quadrature has. Use importance = 0

Its trigger is `R/importance.R:126-133`, a loop over responses that
fires on any structure whose `loglik` slot is non-NULL. So the refusal
is on the PRESENCE of a whole-response likelihood, not on any property
of it. The lane's reading is exactly right.

### The pointwise_loglik proposal: right diagnosis, one slot too coarse

The lane's claim that the importance correction needs separability over
GROUPS and not over ROWS is not merely defensible, it is visible in the
code. `R/importance.R:659-662` is:

    ll  <- wts * row_lpdf(fam, yraw, yraw, dpv, av, extra)
    agg <- smat %*% RTMB::matrix(ll, n, nd)

`smat` is `gmap[["S"]]`, an `ng` by `n` indicator built at
`R/importance.R:348` that maps each row to its single grouping level,
and `R/importance.R:333-343` already REFUSES any design where a row
reaches several groups. So the per-row values are collapsed into
per-group sums on the line after they are produced, and nothing
downstream ever sees a row again: `fn()` does a per-group
log-mean-exp over draws, and `imp_ess()` reads the same `ng` by `nd`
matrix. A slot returning per-group log-likelihoods would replace those
two lines directly. The lane is right, and for a stronger reason than
it gives.

Three things the proposal as written does not account for.

FIRST, the stacking. The correction evaluates at
`ridx <- rep.int(seq_len(n), nd)`, the whole design repeated once per
draw, so a per-group slot must return `ng` by `nd` and a sequential
family must run its recursion over subject-crossed-with-draw, not
subject. For `rw_delta` that is benign and worth saying out loud: the
loop is already vectorized across subjects, so widening it to
subject-by-draw is the same loop over a longer vector, not a new
algorithm. For a family not already vectorized, this is where the cost
would land.

SECOND, grouping alignment. The family's unit and the correction's
group are different objects that merely coincide here.
`rw_delta(subject = id)` with `(1 | p | id)` gives one family unit per
importance group, but the same family under `(1 | item)` would give a
per-subject log-likelihood against a per-item proposal, and summing the
first into the second is not defined. Core has what it needs to check
this (`row_level` from `imp_group_map()`); the proposal should say the
check is required, because getting it wrong is silent.

THIRD, and this is the real over-claim: one slot keyed on `unit` cannot
serve deviance residuals. `R/predict.R:2123` defines the deviance
residual as `sign(y - E[Y]) * sqrt(w * d)` with unit deviance
`d = 2 * (saturated loglik - fitted loglik)`, which is irreducibly PER
ROW. A per-subject value cannot produce it. The findings doc says the
one addition "unlocks `frm(importance =)`, `loo()`, `waic()` and
deviance residuals"; for this family, whose declared `unit` is "one
subject's trial sequence", the last does not follow.

What makes this worth fixing rather than noting: `rw_delta` DOES have a
per-row conditional log-likelihood and already computes it.
`rw_fitted()` returns the per-trial choice probability, and
`test-rl-example.R:102` asserts that
`sum(dbinom(y, 1, fitted(fit), log = TRUE))` equals the reference data
log-likelihood. So the finest factorization available here is the ROW,
not the subject, and `unit = "one subject's trial sequence"` describes
the leave-one-out unit rather than the factorization unit. Those are
two different questions and one slot conflates them.

The design that follows: let the slot return the FINEST factorization
the family has (per row where one exists, per group otherwise), keep
`unit` as the declaration of what a column may honestly be left out at,
and let core dispatch. Importance sums rows into groups with the `smat`
it already builds, or takes group values directly. `loo()` and `waic()`
key their columns on `unit` and keep refusing where the unit is a
group. Deviance residuals require row granularity and refuse otherwise.
That is one slot plus one existing declaration, no more machinery than
the lane proposes, and it does not promise a per-row quantity from
per-group data.

One correction to the seam section's supporting detail. The findings
doc says `structure_unit()` "is defined in R/structure.R and called by
NOTHING in core; its only reader is
extensions/frmtmb.sample/R/loo.R:72, which inlines the same expression
rather than calling it". The first half is right - `R/structure.R:669`
is reached only by `test-structure.R:47` - but the inlining is
deliberate and commented (`loo.R:64-67`: "`structure` and its `unit`
field are documented parts of the frmtmb_structure() protocol, so they
are read directly rather than through core's one-line accessors"),
because an out-of-tree extension must not reach for a `@noRd` internal.
The actionable fact is narrower than the lane implies: core's
`structure_unit()` is dead code that consolidation should export or
delete.

## 5. Packaging

- `inst/rl/rw-delta.R` is read by both
  `vignettes/reinforcement-learning.Rmd:18` (through
  `knitr::read_chunk()`) and `tests/testthat/helper-rl.R:12-13`
  (through `system.file()`), so the page and the test genuinely cannot
  drift. The argument for `inst/` over a test helper is sound: a
  vignette has no one relative path reaching `tests/testthat/` under
  both `R CMD build` (cwd `vignettes/`) and pkgdown (cwd package root).
  There is precedent in the tree, `inst/extdata/habit-hierarchical.csv`
  already sharing material between a vignette and the suite, and the
  cost is 13 KB shipped to users who also gain a sourceable worked
  example. Right call.
- The vignette's fallback
  `if (!nzchar(rl_src)) rl_src <- "../inst/rl/rw-delta.R"`
  (`reinforcement-learning.Rmd:17`) resolves only when the cwd is
  `vignettes/`, which is the very asymmetry the findings doc argues
  makes a relative path unusable. It is harmless, because the
  `system.file()` path wins whenever the package is installed and it is
  under `R CMD build`, but it is inconsistent with the prose beside it.
- The vignette installs nothing and registers nothing that outlives the
  session: `frmtmb_register_aterm()` writes to a package-level
  environment and re-registration is a documented no-op.
- hBayesDM appears in no DESCRIPTION field and is not installed in
  either library. The claim that it was neither installed nor consulted
  is consistent with everything visible.
- `_pkgdown.yml` lists `reinforcement-learning` under Guides after
  `case-studies`. All eight vignettes in `vignettes/` are listed and
  none is orphaned, so the entry is complete and correctly placed.
- One real side effect, which the lane records itself: `helper-rl.R`
  registers `reward()` at helper-load time, so every test file in the
  suite runs with one extra entry in the compatibility vocabulary.
  `test-compat.R` and `test-compat-register.R` assert membership rather
  than an exact vocabulary, so they pass; a future test that pins the
  vocabulary exactly would break. Worth keeping in the findings doc.

## 6. Runs

One process each, private library, resolved frmtmb path printed at the
top of every run.

| run | result | lane claimed |
|---|---|---|
| `test-rl-example.R`, gate off (`NOT_CRAN=true`) | PASS 76, FAIL 0, SKIP 2 | 82 / 0 / 2 |
| `test-rl-example.R`, gate on | PASS 82, FAIL 0, SKIP 0 | 82 / 0 / 0 |
| `test-rl-example.R`, no `NOT_CRAN` at all | PASS 70, FAIL 0, SKIP 3 | not reported |
| `test-message-uniqueness.R` | PASS 6, FAIL 0, SKIP 0 | pass |
| `test-bracket-access.R` | PASS 8, FAIL 0, SKIP 0 | pass |
| vignette knit, pandoc 3.8.3 | 6.6 s, no warnings | 11.5 s, no warnings |
| `R CMD check --as-cran`, suite | FAIL 0, WARN 0, SKIP 176, PASS 4029 | identical |

Everything passes. Two small bookkeeping corrections.

DEFECT 3, minor. The gate-off assertion count is 76, not 82. 82 is the
gate-ON count, and the findings doc reports 82 for both rows
("82 pass with 2 skips gate off and 82 with 0 gate on"). The six
assertions that separate them are the two Stan tests. The lane appears
to have copied the gate-on total into the gate-off row.

The skip counts are right and, more importantly, correctly EXPLAINED:
with `NOT_CRAN` unset entirely the recovery test also skips (3 skips,
70 assertions), and the lane's own trap note about reading a run
together with its skip count is exactly the right instinct. In the
`--as-cran` run the three RL skips are `test-rl-example.R:178` (the
recovery test, `skip_on_cran()` with `NOT_CRAN=false`) and `:195` and
`:207` (the Stan tier), which is what the findings doc says.

The vignette knit renders the `importance` refusal into the page as
intended (the `error = TRUE` chunk is the only `Error:` in the output)
and produces no warnings.

`R CMD check --as-cran --no-manual` on the built tarball, private
library, `_R_CHECK_CRAN_INCOMING_=false`, `NOT_CRAN=false`:
**Status: OK**, with zero NOTE, WARNING or ERROR lines anywhere in the
log. Tests 124 s, vignette rebuild 123 s. The lane's report that the
expected V8 NOTE did not appear is confirmed; nothing in this lane
suppresses it.

### Main checkout

Main was recorded at 9a1d3f3 with a 77-line dirty status before any
work here. It is now at 528aad7 with a clean tree, because a SIBLING
lane committed during this review: 528aad7 is "Merge branch
'wt-spline-core'", authored by the repository owner, and it touches no
RL file. Every git command this review ran against the main checkout
was `status`, `rev-parse`, `log` or `show`, and no file under it was
written.

## Edits I made

Two, both in the worktree, both reported in full.

1. `inst/rl/rw-delta.R:67-74` - the DEFECT 1 fix. The two block
   constructors now read

       idx <- t(matrix(vapply(rows, function(r) c(r, rep(r[1L], nt - length(r))),
                              integer(nt)), nrow = nt))
       mask <- t(matrix(vapply(len, function(l) as.numeric(seq_len(nt) <= l),
                               numeric(nt)), nrow = nt))

   with a four-line comment saying why the `matrix()` is load-bearing.
   For `nt > 1` this is a no-op reshape of what `vapply()` already
   returns, so nothing changes; for `nt == 1` it stops `t()` from
   turning a dropped vector into a 1-by-`n_subj` block. Verified: the
   six-subject one-trial case now gives -4.158883083, exactly
   `6 * log(0.5)`, against -4.356251586 before. Re-ran the test file at
   both gates after installing the fix: 76/0/2 and 82/0/0, unchanged.

2. `dev/review-rl.md` - this file, which did not exist before.

No other file in the worktree was touched, and nothing was committed.

## Punch list

Blocking, in the sense that the numbers are published and wrong.

| # | file and line | what |
|---|---|---|
| P1 | `vignettes/reinforcement-learning.Rmd:396-399` | The "one gradient" column (50 / 72 / 95 ms) and the sentence "the gradient follows at 1.4x and 1.9x". All three shapes have the same run-time cost, ratios 0.66x to 1.20x over 100-call batches. Drop the column and say the penalty is entirely at tape construction. Same table in `dev/rl-findings.md` under "The elementwise alternative, measured". |
| P2 | `vignettes/reinforcement-learning.Rmd:409-421` | The scaling table's gradient column and the paragraph "The same pair reverses for the gradient". Interleaved, it does not reverse: 40x400 is slower than 160x100 on the gradient too (1.26x). The "27 percent" is 70 percent here. Either re-measure interleaved with gc, or reduce the claim to a direction with no percentage. |
| P3 | `dev/rl-findings.md`, Validation section 1 table | `log_prob` -392.4486089 labelled "30 subjects by 100 trials" is stale; that design gives -1427.8453940, and -392 is not plausible for 3000 Bernoulli trials at this bandit. Fix the value or the label. |

Non-blocking, but should land with the lane.

| # | file and line | what |
|---|---|---|
| P4 | `dev/rl-findings.md:165-175` | The ">1000x" correction should cite `dev/hmm-feasibility.md:139`, which already made it and already recorded the length-K against length-n distinction. It should also stop quoting a single multiplier: the factor grows with n (1.3x at n=1000, about 25x from n=50000). And the closing suggestion to re-measure the `test-perf.R` canary is misdirected - the canary quotes no number, and at its own n=100000 the elementwise build takes 78.5 s against its 20 s bound, so it works. |
| P5 | `dev/rl-findings.md`, seam section | "unlocks `frm(importance =)`, `loo()`, `waic()` and deviance residuals" over-claims the last one: deviance residuals need a per-ROW saturated comparison (`R/predict.R:2123`) that a per-group slot cannot give. The proposal should also name the `n * nd` stacking and the grouping-alignment check. See section 4 above. |
| P6 | `dev/rl-findings.md`, Verification run table | Gate-off row says 82 assertions; it is 76. 82 is the gate-on count. |
| P7 | `vignettes/reinforcement-learning.Rmd:17` | The fallback `"../inst/rl/rw-delta.R"` resolves only from `vignettes/`, which is the asymmetry the findings doc argues makes a relative path unusable. Harmless, but it contradicts the prose beside it. |
| P8 | `R/structure.R:669` (core, not this lane) | `structure_unit()` is dead code, reached only by `test-structure.R:47`. Consolidation should export it or delete it. Related: `frmtmb_register_aterm()` silently accepts a re-registration at a DIFFERENT arity, clobbering the first entry. |

## Verdict: GO-WITH-FIXES

The engineering is sound and the central claims are real. The recursion
is correct on every point I checked, the Stan identity is exact at
1e-13 at two different points with a provably zero constant, the
recovery study reproduces closely on independent seeds including the
awkward part (`log sd(alpha)` collapsing at 20 trials), the family
refuses what it cannot do in its own words, `R CMD check --as-cran` is
OK with zero notes, and the whole thing lands without touching `R/`.
The lane's habit of naming what it could not do rather than working
around it is the best thing about it, and the `pointwise_loglik`
diagnosis is not only correct but better supported by the code than the
lane realised.

What holds it back from a plain GO is that three published tables carry
numbers that do not reproduce (P1, P2, P3), one of them in a
user-facing vignette. All three come from the same cause: single-shot
`system.time()` on Windows, whose ~15 ms granularity swamps the
quantities being compared, plus one stale value never refreshed after a
design change. None of them changes a conclusion. All are edits to
prose and tables, not to code.

The one code defect (silent wrong likelihood when every subject has one
trial) is fixed above and verified.

### The pointwise_loglik seam, for the consolidation

A per-group log-likelihood is the right shape for the importance
correction and the code proves it: `R/importance.R:659-662` computes
per-row densities and immediately collapses them with a sparse
`ng`-by-`n` indicator, and `imp_group_map()` already refuses any design
where a row reaches two groups, so group separability is a
precondition the correction enforces today rather than a property a new
slot would have to introduce. A slot returning per-group values drops
straight into those two lines. But one slot keyed on `unit` is too
coarse to carry everything the lane hangs on it: deviance residuals
need a per-row saturated comparison that no per-group value can supply,
and `rw_delta` itself has a per-row conditional log-likelihood it
already computes for `fitted()`, so its declared `unit` of "one
subject's trial sequence" is describing the leave-one-out unit, not the
factorization unit. The seam that actually pays for itself is a slot
returning the FINEST factorization a family has, with `unit` kept as
the separate declaration of what may honestly be left out; core then
sums rows to groups with the `smat` it already builds for importance,
keys `loo()` and `waic()` columns on `unit`, and requires row
granularity for deviance. That is the same amount of machinery the lane
proposes and it stops the protocol promising a per-row quantity from
per-group data. Two implementation facts belong in the design note: the
correction evaluates a design stacked `n * nd` times, which a
sequential family must vectorize across subject-by-draw (cheap for this
family, since its loop is already vectorized across subjects), and the
family's grouping must be checked against `row_level` or a
`(1 | item)` model would silently sum the wrong things.

### Should this be an exported family in an extension next?

Yes, but after the pointwise seam lands and not before, and it belongs
beside the evidence-accumulation families in the cognitive-model
extension rather than in a package of its own, because RL-DDM hybrids
need both halves in one place.

# Punch re-check, 2026-09-05

Re-run against the lane's punch round (its record is at
`dev/rl-findings.md:592`). Same private library, same cache, same
worktree. Main is now at 21f947b and the lane merges onto it after this
verdict; nothing here was committed and the main checkout was not
touched.

## Scope is unchanged

`git diff --name-only 316a28b -- R/` is still EMPTY. The tracked diff is
still `NEWS.md` and `_pkgdown.yml`; the untracked set is the same five
lane files plus this review. Nothing crept in with the punch round.

## Reproduced

**P1, the gradient-cost equality.** Re-measured interleaved with
`gc(FALSE)`, median of 7 builds, gradients over batches of 100, all
three shapes on one dataset:

| shape | build (median) | gradient (best of 3) | objective |
|---|---|---|---|
| vectorized | 0.06 s | 12.30 ms | -1967.6888631596 |
| row loop, scalar Q | 0.60 s | 12.10 ms | -1967.6888631596 |
| row loop, length-n vector | 0.68 s | 12.20 ms | -1967.6888631596 |

Gradient ratios 0.98 and 0.99 against the lane's 0.99. The gradient
column is flat, which is the corrected claim, and it reproduces. The
build column agrees closely with the lane's re-measurement (0.08 / 0.57
/ 0.65 against my 0.06 / 0.60 / 0.68); my absolute gradient is 12.2 ms
against their 7.2 ms, which is machine state and not a disagreement,
because the claim is the ratio. Objective identical across shapes to
14 digits on my seed, as on theirs.

**The one-trial regression test.** `tests/testthat/test-rl-example.R:83`
runs and passes. It is a real test of the defect, not a vacuous one: it
asserts the block is 6-by-1 rather than 1-by-6 AND pins the likelihood
at `6 * log(0.5) = -4.158883083`, which I measured at -4.356251586
against the transposed block. The choice of quantity is good - every
subject's only trial starts from Q1 = Q2 = 0, so the value is
parameter-free and the test cannot drift with the fixture. My fix is
kept verbatim at `inst/rl/rw-delta.R:67-74`.

**Test and knit counts**, each in its own process, all matching the
lane's report exactly:

| run | measured | lane reported |
|---|---|---|
| `test-rl-example.R` gate off | PASS 79, FAIL 0, SKIP 2 | 79 / 0 / 2 |
| `test-rl-example.R` gate on | PASS 85, FAIL 0, SKIP 0 | 85 / 0 / 0 |
| `test-message-uniqueness.R` | PASS 6, FAIL 0, SKIP 0 | 6 |
| `test-bracket-access.R` | PASS 8, FAIL 0, SKIP 0 | 8 |
| vignette knit | 2.9 s, no warnings | 8.9 s |

The gate-on run executes the Stan tier with zero skips, so P3's
re-derived identity is asserted by the suite itself rather than only
recorded. The knit time differs only because my Stan and tape caches
were warm.

**P3, the identity.** The five numbers the lane re-derived
(-1427.8453940, -1460.3042080, -4.5e-13, -2.5e-12, 1.199e-14, overall
5.211) are the five I measured independently in the first pass. The
rendered page now carries -1427.8453940 in both the Stan and the frmtmb
row and contains no occurrence of the stale -392.44. `fit$frame[["priors"]]`
being empty is now recorded as a checked fact rather than an
assumption, which is the right strengthening.

**The vignette prose no longer claims a run-time penalty anywhere.**
Swept the whole file for the old language. There is no surviving
"1.4x", "1.9x", "gradient follows", "reverses", "50 / 72 / 95 ms" or
"392.44". What the page says now is correct: "The whole penalty is paid
at TAPE CONSTRUCTION, and none of it per evaluation", "The gradient
column is flat to within measurement noise", and, for the scaling
table, more trials at fixed rows "costs the same per gradient, because
the tape is the same size either way". The remaining hits on those
strings anywhere in the lane are all inside `dev/rl-findings.md`
passages that exist to record what was wrong, which is where they
belong.

**P2, the build gap.** Handled the way a disagreement should be. Both
documents now give the direction only and name the disagreement rather
than picking a winner: the vignette says "The size of that build gap is
not stable enough to quote. Two interleaved runs of this table on one
machine put it at 1.27x and 1.70x." That is the honest resolution, and
it is better than either of us quoting a single figure.

**P4, P5, P7.** P4 now cites `dev/hmm-feasibility.md:139`, quotes a
sweep instead of a multiplier, and withdraws the canary suggestion. P5
is rewritten to exactly the shape I argued for: the slot carries the
finest factorization available, `unit` stays as the separate
leave-one-out declaration, and both the `n * nd` stacking and the
`row_level` grouping-alignment check are named. P7 replaced the
relative fallback with `system.file(..., mustWork = TRUE)`, so a
missing install now fails loudly. P8 is correctly recorded as core work
rather than fixed here.

## One residual nit, non-blocking

`vignettes/reinforcement-learning.Rmd`, in "What a built-in family
would add", still carries the pre-P5 framing: "a pointwise
log-likelihood matrix, which `loo()` and `waic()` need ... A future
`pointwise_loglik` slot, or a `loglik` allowed to return a vector over
the structure's own `unit`, would serve both." For this family `unit`
is the subject sequence, so a vector over `unit` serves GROUP-level
`loo()`, not the row-level pointwise matrix the sentence just asked
for. The findings doc's rewritten seam section draws that distinction
correctly; the vignette paragraph did not get the same pass. It is a
loose sentence in a forward-looking section, not a wrong measurement
and not a user-facing number, so it should not hold the merge. Worth a
line when the seam actually lands.

## Verdict: GO, merge it

Every blocking item from the first review is fixed and independently
reproduced. The three non-reproducing tables are corrected at the
source and in the vignette, and the one that could not be settled is
now reported as unsettled with both figures named, which is the right
answer rather than a compromise. The code defect is fixed and pinned by
a regression test that pins arithmetic rather than shape. The suite,
the two boundary tests and the knit all match the lane's report exactly
on my machine, and the scope is still zero lines under `R/`.

The lane's handling of this round is worth recording: it did not defend
any of the three bad tables, it changed the measurement METHOD once for
all of them rather than patching numbers one at a time, and where our
two machines disagreed it published the disagreement. That is the
behavior that makes the rest of the document trustworthy.

No further edits by me in this round. My total edit footprint across
both rounds remains two files: the fix at `inst/rl/rw-delta.R:67-74`
and this review.
