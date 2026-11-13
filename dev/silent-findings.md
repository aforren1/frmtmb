# Lane silent: seven small silent defects (round of 2026-11-11)

Base: e8abed0a (frmtmb 0.81.0), measured against
`C:/Users/adf44/source/r/rellib-r20`. Lane library:
`C:/Users/adf44/source/r/wt-silent-lib`. Every script named here is in
`dev/` of this worktree and writes its log beside it
(`dev/silent-*.log` and `dev/silent-*-log/`, gitignored); the numbers
below are pasted from those logs. "base" is a run against rellib-r20,
"lane" against the lane library first.

Status: DONE. All seven items fixed, each with a test seen to fail on
rellib-r20; every suite green one file per process and in one session
three ways; changed tests green under OpenBLAS 0.3.26 LAPACK; checks
1 NOTE (V8) each for core, latent and survival.

## Summary

| Item | What was silent | Fix | Test (fails on base) |
|---|---|---|---|
| 650 | `dharma_residuals()` scored a censored row as an exact value | censored randomized quantile residual | core `test-dharma-cens.R` |
| 857 | `hmm()` dropped the component's `post$stationary` | copy per state | latent `test-hmm-stationary.R` |
| 1562 | mv `logLik()` `nobs` counted rows no response uses | count rows a response keeps; also `df.residual()`, `anova()` | core `test-subset-rate.R` |
| 4235 | `drop1()` left out `mo()`, `me()`, smooth terms | scope from the whole formula; `add1()` too | core `test-drop1-specials.R` |
| 3874 | cox knots read `range(y)` only | upper knot covers every upper bound read | survival `test-cox-bounds.R` |
| 4675 | two aterms could share a key | refused at registration | core `test-parse.R` |
| 720 | `ce_aterms()` dropped a payload | refuse, as `aterms_for_newdata()` | core `test-newdata-aterms.R` |

Failing output on base and passing output on the lane, one file per R
process (`dev/silent-run-test.R`, logs `dev/silent-base-log/`,
`dev/silent-lane-log/`, `dev/silent-d1-log/`):

```
base test-dharma-cens.R      pass=10 fail=5  error=0
base test-subset-rate.R      pass=89 fail=4  error=0
base test-drop1-specials.R   pass=8  fail=14 error=1   (final version)
base test-parse.R            pass=48 fail=9  error=0
base test-newdata-aterms.R   pass=56 fail=1  error=0 warn=1
base test-hmm-stationary.R   pass=2  fail=5  error=1 warn=1
base test-cox-bounds.R       pass=9  fail=7  error=0   (first version)
lane test-dharma-cens.R      pass=15 fail=0
lane test-subset-rate.R      pass=93 fail=0
lane test-drop1-specials.R   pass=24 fail=0
lane test-parse.R            pass=57 fail=0
lane test-newdata-aterms.R   pass=57 fail=0
lane test-hmm-stationary.R   pass=10 fail=0
lane test-cox-bounds.R       pass=16 fail=0
```

The base failures were read one by one; each is the defect, not the
test (excerpts under each item). Two base failures of the first
`test-cox-bounds.R` were the test's own: its truncation case had times
from 0.0004, so brms's lower knot was already 0 and no entry time could
lie below it. The case now shifts the times by 2.

## 650: `dharma_residuals()` on a `cens()` fit

Decision: compute the censoring-aware residual, not refuse.

What DHARMa does. DHARMa 0.5.0, the installed version, documents no
censoring support: `dev/silent-dharma-docs.R` searched every Rd page
for "censor" and found none, and neither vignette mentions it. A web
search found no DHARMa page on it either. Its default residual
(`method = "PIT"`, `getQuantile()`, printed by
`dev/silent-dharma-src.R`) is `runif(lower, lowerOrEqual)` with
`lower = mean(sim < obs)` and `lowerOrEqual = mean(sim <= obs)`: it
randomizes over ANY tie, which its Rd gives as the reason it replaced
the old integer method (families "partly continuous, partly discrete",
such as the Tweedie). So the backlog's "the point mass at each
censoring point is not a distribution DHARMa's rank transform can use"
has been wrong since DHARMa 0.3.1.

What is correct. The residual of a row that observes `a < Y <= b` is
`U ~ Uniform(F(a), F(b))`, Dunn and Smyth's randomized quantile
residual extended to censoring (the "randomized survival probability"
of the survival literature). It is uniform given each row's censoring
time c, so it needs no model of the censoring and allows row-varying
times: for right censoring, P(U <= u) = u for u <= F(c) from the
events, and F(c) + (1 - F(c)) (u - F(c)) / (1 - F(c)) = u above it.
This is an identity. An uncensored row keeps `F(y)`.

Construction (`dharma_censor()`, `R/diagnostics.R`): a censored row's
draws that fall inside its event are set to one point `s` of the event,
and its observation to `s`. DHARMa's PIT then draws over
`(mean(draws <= a), mean(draws <= b))`. Right: `s` a double just above
`y`, or `y + 1` on the integer lattice, where "censored at k" is
`Y > k`. Left: `s = y`. Interval: `s = y2`. Draws come from
`simulate()`, so a `trunc()` window already holds. On a censored row
the DHARMa object's `observedResponse` and draws hold `s`, which the
new Rd section says.

Base failure (`dev/silent-base-log/base-frmtmb-test-dharma-cens.log`):
on base a right-censored row's residual equals the bottom of its
interval, `mean(draws < y)`, and a left-censored row's the top, so
`all(r > lo)` and `all(r < hi)` fail; a right-censored count read
`(F(k - 1), F(k))` instead of `(F(k), 1)`.

Uniformity, 20 seeds each, n = 300, nsim = 250
(`dev/silent-dharma-unif.R`, logs `dev/silent-dharma-unif-{base,lane}.log`):
gaussian with row-varying right and left censoring and intervals
(censored share about a third), poisson right-censored at a row-varying
k in 1..4, and the gaussian model fitted to skewed errors (the wrong
model):

```
SUMMARY base: p<0.01 on gaussian 20/20, on poisson 20/20, on the wrong (skewed) model 20/20; median p gaussian 1.74e-06 poisson 1.33e-20 wrong 7.57e-07
SUMMARY lane: p<0.01 on gaussian 0/20, on poisson 0/20, on the wrong (skewed) model 14/20; median p gaussian 0.838 poisson 0.855 wrong 0.00567
```

So on base every correct censored model was rejected (a false alarm in
20 of 20), and on the lane none was, while the wrong model is still
rejected in 14 of 20.

Messages and docs changed with it: the two `residuals(type = "osa")`
refusals in `R/predict.R` that said `dharma_residuals()` "covers the
same ground through simulate(censored = TRUE)" and "is not an
alternative" now say it scores each censored row on its event; the
`?residuals.frmtmb_fit` and `?simulate.frmtmb_fit` paragraphs likewise.

## 857: `hmm()` drops the component's `post$stationary`

Confirmed on 0.81.0 (`dev/silent-hmm-probe.R base`, log
`dev/silent-hmm-probe-base.log`): `hmm(2, skew_normal())` builds, its
`post$stationary` is empty, and with two states skewed in opposite
directions (alpha = 5 and -5, 20 sequences of 40) it stalls:

```
base:
seed 1 free logLik -1019.4522 alpha -0.00122 -5.607 | from alpha=0 logLik -1056.9509 alpha 1.189e-14 -7.781e-15 escape NULL
seed 2 free logLik -985.6858 alpha 4.468 0.001169 | from alpha=0 logLik -1014.4305 alpha 2.946e-14 -4.916e-15 escape NULL
seed 3 free logLik -953.5887 alpha 4.86 -4.931 | from alpha=0 logLik -1007.4344 alpha -1.799e-14 -1.713e-14 escape NULL
lane:
seed 1 free logLik -998.2551 alpha 4.18 -5.451 | from alpha=0 logLik -1019.4522 alpha -8.365e-14 -5.607 escape 4.00000 37.49867
seed 2 free logLik -960.6004 alpha 4.332 -5.218 | from alpha=0 logLik -985.6858 alpha 4.468 -1.449e-14 escape 4.00000 28.74473
seed 3 free logLik -953.5887 alpha 4.86 -4.931 | from alpha=0 logLik -976.1491 alpha 5.126 -5.035e-14 escape 4.00000 31.28532
```

With the DEFAULT start ("free"), one state stops at alpha = 0 on seeds
1 and 2: every state's alpha starts from one pooled value (-2.01 on seed
1), so the state skewed against the pooled residual starts on the wrong
side. The only sign is the SE check's warning that the alpha SE is
unavailable. On the lane the default fit escapes: 21.20 (seed 1) and
25.09 (seed 2) log-likelihood units higher; seed 3 did not stall and is
unchanged to the printed digit.

Fix: `hmm()` builds the declaration alongside `init_fns`, keyed
`paste0(dp, k)`, as core's `mixture()` does, and passes it as
`frmtmb_family(post = )`. No core change.

Defect found, not fixed (filed): the escape is ONE pass. With both
states started at 0 ("from alpha=0" above) the escape ran 4 restarts,
each moving one state, and kept the best, which leaves the other state
at 0 (seeds 1 to 3: -1019.45, -985.69, -976.15 against -998.26,
-960.60, -953.59). It is core's `escape_stationary()` in `R/fit.R`, a
file five other lanes touch, and it changes `mixture()` too, so it is
filed rather than fixed.

The test pins the one-stalled-state case from an explicit start (state
1 at 0, state 2 at its default start), so it does not depend on the
default start's path. Rounding: see "Rounding" below.

## 1562: multivariate `logLik()` `nobs`

brms 2.23.0 on the same mvbf (`dev/silent-nobs-brms.R`, an empty
`brm()` fit, which holds the data `nobs.brmsfit()` reads):

```
A nobs() 80  nobs(resp = 'y1') 72
B nobs() 80  nobs(resp = 'y2') 60
```

`nobs.brmsfit()` is `nrow(model.frame(object))` without `resp`, so
frmtmb's `nobs()` (80) already matches brms and stays. brms has no
`logLik()` `nobs` or `BIC()`; those are frmtmb's and must count the
rows that carry a density. The count is now `loglik_nobs()`
(`R/frame.R`, replacing `mi_rows_gone()`): the union over responses of
the rows each one reads, so a row counts once if any response reads it.
It is the old count wherever every row is used.

Every reader of the count (`dev/silent-nobs-probe.R`, logs
`dev/silent-nobs-probe-{base,lane}.log`; A: both responses on rows 1 to
72 of 80; B: rows 1 to 72 and 1 to 60; C: anova against an 80-row fit;
D: no subset):

```
base A: logLik nobs attr 80  nobs() 80  BIC 441.137756435 (= -2ll + log(80) df)  df.residual 74
lane A: logLik nobs attr 72  nobs() 80  BIC 440.505593341 (= -2ll + log(72) df)  df.residual 66
base B: logLik nobs attr 80  df.residual 74
lane B: logLik nobs attr 72  df.residual 66
base C: anova(80-row fit, 72-row fit): ran
lane C: anova() needs fits with the same number of observations (got 72, 80); ...
both D: logLik nobs attr 80  nobs() 80  df.residual 74
```

- `AIC()` reads no n. `BIC()` reads the attribute.
- `df.residual()` read `n_obs`; it now reads the same count, so it and
  `BIC()` never disagree (it fed Rubin's complete-data df in
  `R/multiple.R`, where fits have no subset rows).
- `anova()` compared `n_obs`, so a 72-row fit and an 80-row fit were
  tested against each other (base C). It now compares the same count,
  and so does the pooled `anova()` of `R/multiple.R`.
- `loo()` on an ML fit refuses; on draws (frmtmb.sample) a subset model
  needs `resp =`, and then has one column per row of that response,
  which is `nobs(resp = )`. Nothing to change.

## 4235: `drop1()` leaves out `mo()` terms

Base (`dev/silent-drop1-probe.R base`, log
`dev/silent-drop1-probe-base.log`; n = 200, seed 11): `drop1()` of
`y ~ mo(k) + z`, `y ~ s(x) + z` and `y ~ me(xo, sdx) + z` each listed
`z` alone, because `terms(fit)` is the fixed design's terms, built after
the parser takes the specials out. The documentation said so ("Random-
effect, smooth, and mo()/mi() terms are not part of the deletion
scope"), the table did not. `add1()` worked for `mo(k)` and `s(x)` (LRT
68.737 against a hand refit 68.73707894) but:

- `add1(y ~ z, ~ . + me(xo, sdx))` reported LRT -543.37: the larger
  model's likelihood holds the measurements of `xo` (logLik -617.74
  against -346.06).
- `add1()` read the fixed design's terms too, so on
  `y ~ me(xo, sdx) + z` the scope `~ . + me(xo, sdx):z` offered
  `me(xo, sdx)` again with Df 0 and never the interaction
  (`dev/silent-add1-probe.R`), and on `y ~ x + (1 | g)` the scope
  `~ . + z + (1 | g)` offered `1 | g` as a fixed term and errored
  ("'|' not meaningful for factors", then "No complete observations").

Fix (`R/confint.R`): `drop1_terms()` takes the scope from the main
formula without its bar terms and without its autocorrelation terms;
`add1.frmtmb_fit()` builds a formula scope against the whole formula
and then calls `add1.default()`. A row whose refit changes the `me()`
calls (compared as `anova()` compares them) has `NA` statistics and a
message, in both. Lane (`dev/silent-drop1-probe-lane.log`): `drop1()`
lists `mo(k)` (Df 3, LRT 68.737, the hand refit's 68.73707894), `s(x)`
(Df 2, LRT 49.875, the same as `add1()` from `y ~ z`) and
`me(xo, sdx)` (NA, with the message), and `z` as before; `add1()` of
`me(xo, sdx)` is NA with the message. mi() predictors need a multivariate model, which
`drop1()` refuses and `add1()`'s update refuses, so they cannot reach
this.

## 3874: survival boundary knots read `range(y)` only

`bhaz_spec()` built the default knots from `range(y)` and checked given
knots against `range(y)`, while `cox_cb_entry()` clamps every other
bound the likelihood reads (interval upper bound, truncation bounds) to
the knots. So both the default and a given basis silently read zero
hazard past the last knot. Base failure (first version of
`test-cox-bounds.R`): default upper knot 8.9 against an interval bound
of 17.5, and every I-spline row of a bound past the knot all ones.

Fix (as revised in the fix round, below): `cox_read_bounds()` lists the
interval upper bounds (code 2 rows only, as the frame normalizes the
others to `y`) and the finite truncation upper bounds. The default upper
knot widens to them; given knots whose upper knot misses one are
refused with their own message; a response outside given knots keeps
the old message. Lower truncation bounds do not move the knots. A fit
without an upper bound past the data keeps brms's default exactly
(`identical()` in the test). The first version also widened the lower
knot to the truncation lower bounds, which was wrong (MAJOR 1 of the
review).

brms parity: brms's default reads the response alone, but brms cannot
fit an interval-censored Cox model (its `cox_log_lcdf()` reads `cbhaz`
at `Y`, so the interval has probability 0: `dev/upstream-bugs.md`,
brms-36), and the same function serves its truncation, so there is no
working brms behavior to match there.

## 4675: two registered addition terms can claim one key

Base (`dev/silent-aterm-probe.R base`):

```
register x (arity 3): registered
register x2: registered
   y | x2(z) + x(u, v, w) ~ 1 -> x2=v  x1=u  x3=w
   y | x(u, v, w) + x2(z) ~ 1 -> ERROR: Duplicated addition term `x2()`
register vint1: registered
register cens_y (arity 2): registered
```

`z` vanished without a word in the first spelling. Lane: every one of
these registrations is refused, naming the key and its owner. Fix
(`R/parse.R`): `aterm_key_clash()` runs before the compatibility
registration, the last step that can refuse, so a refusal leaves both
registries as they were (the test asserts it). Keys: `name` (arity 1 or
an option term, which `registered_aterm_of()` finds by name),
`name1..k` otherwise; core keys are the core term names, `vint<k>`,
`vreal<k>`, `trunc_lb`, `trunc_ub`, `cens_y2`, `se_sigma`, `mi_sd`,
`thres_gr`. Re-registering the same name is still a replacement. No
shipped extension is refused: eam's `dec` and `ndt_group`, learn's
`reward` (2), `payoff` (4), `stage2` (2) and survival's `bhaz` share no
key (the full suites load them; see below).

## 720: `ce_aterms()` drops a `vint()`/`vreal()` payload

Where it reaches. The grid holds every data variable at its mean, so a
payload fails on the grid only when its expression does not follow the
rows (`vreal(rep_len(c(1, 3), 120))` on a 120-row fit, 100-row grid).
In core's own `conditional_effects()` the drop was masked: the mean
path and `method = "predict"` both reach `aterms_for_newdata()`
(through `predict_mean_se()` and `predict_mean_response()`), which
refuses (`dev/silent-ce-trace.R`). The unmasked path is `ce_sim_one()`,
the seam frmtmb.sample's `method = "predict"` draws through, one call
per draw:

```
base: == ce_sim_one on a 100-row grid: class numeric length 100  NA 100
      warning: NAs produced
      ce_aterms names:  (empty means dropped)
```

A family that declares the term in `required_aterms` was already named
("did not reach the family"); one that does not drew `NA`.

Fix: the strict set of `ce_aterms()` adds every custom data term
(`is_custom_data_aterm()`: `vint<k>`, `vreal<k>`, registered terms) and
`rate()`, which `aterms_for_newdata()` also requires. Lane
(`dev/silent-ce-probe-lane.log`): the mean path now stops in
`ce_aterms()` ("conditional_effects() could not evaluate
vreal(rep_len(c(1, 3), 120)) on the effect grid"); the column spelling
`vreal(s)` is unchanged (rows 100, estimate range -0.233 to 3.70 on
both arms).

Defect found, not fixed (filed): at `resolution = 120` the same
expression has the grid's length and is read row by row on both arms
(estimate range -0.317 to 5.57), pairing the fitted rows' values with
grid rows.

## Rounding

The only new assertion that rests on an optimizer path is the hmm
escape test. `dev/silent-hmm-ulp.R` (log `dev/silent-hmm-ulp.log`)
perturbs the response to `y * (1 + 2e-15 k)`, k = 0..20, and repeats
the test's checks: `SUMMARY all four checks hold on 21 of 21 draws`
(2 restarts, gain 0.02123 of |logLik|, alpha1 4.18, logLik within
1.4e-15 of the reference, relative, on every draw). The other new
assertions are exact relations (an interval each residual must lie in,
counts, names, knots) or identities between a table row and the same
refit written by hand.

## Suites

Every suite, one test file per R process, 3 at a time, NOT_CRAN=true,
lane library first (`dev/release/run-par.sh silent-full
dev/silent-all.jobs 3` with `FRMTMB_RELLIB` set; log
`dev/release/silent-full.log`, per-file logs
`dev/release/silent-full-files/`). Every per-file log's `lib:` line
names `wt-silent-lib` for frmtmb (and for latent and survival).

```
SILENT-FULL ran 514 of 514
frmtmb files=319 pass=23546 fail=0 err=1 skip=239 warn=0
frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=33 pass=1848 fail=0 err=0 skip=5 warn=0
frmtmb.latent files=13 pass=697 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=87 pass=4050 fail=0 err=0 skip=18 warn=0
frmtmb.spatial files=8 pass=478 fail=0 err=0 skip=6 warn=0
frmtmb.survival files=16 pass=572 fail=0 err=0 skip=0 warn=0
```

The one error was `test-osa-inference.R`, which matched the old
refusal text "LATENT uncensored" that this lane rewrote (item 650). The
assertion now matches the new text; rerun alone:
`RESULT test-osa-inference.R pass=34 fail=0 error=0 skip=0 warn=0`
(`dev/silent-lane-log/lane-frmtmb-test-osa-inference.log`). The skips
are the ungated tier's gates (brms fits, Stan); no changed path reaches
a brms fit, so the gated tier was not run.

Sessions (each suite in one session, forward, reverse and preload;
`dev/release/run-sessions.sh silent-sess 3` with `FRMTMB_RELLIB` set;
log `dev/release/silent-sess.log`): `SILENT-SESS ran 27 of 27`, fail 0
and err 0 in all 27, and every LEAK line is tagged `[with load of
...]` (an extension or a package an extension loads), so none is a
test's own.

```
SESSION frmtmb forward files=306/319 pass=23549 fail=0 err=0 skip=239 warn=0
SESSION frmtmb reverse files=306/319 pass=23714 fail=0 err=0 skip=239 warn=0
SESSION frmtmb preload files=306/319 pass=23817 fail=0 err=0 skip=241 warn=0
SESSION frmtmb.coupling forward|reverse|preload files=11/11 pass=551 fail=0 err=0 skip=5
SESSION frmtmb.eam forward|reverse|preload files=33/33 pass=1848 fail=0 err=0 skip=5
SESSION frmtmb.latent forward|reverse|preload files=13/13 pass=697 fail=0 err=0 skip=2
SESSION frmtmb.learn forward|reverse|preload files=16/16 pass=474 fail=0 err=0 skip=13
SESSION frmtmb.ode forward|reverse|preload files=11/11 pass=549 fail=0 err=0 skip=1
SESSION frmtmb.sample forward|reverse|preload files=82/87 pass=4050 fail=0 err=0 skip=18
SESSION frmtmb.survival forward|reverse|preload files=16/16 pass=572 fail=0 err=0 skip=0
SESSION frmtmb.spatial forward|reverse|preload files=8/8 pass=478 fail=0 err=0 skip=6
```

(The extension lines are folded: each order printed the same counts.
`files=306/319` counts files with at least one result; a file skipped
whole has none, as in earlier rounds' records.)

## OpenBLAS

`dev/ciharden-openblas.sh 0.3.26 lapack` (build log
`dev/silent-openblas-build.log`: Rblas 74 symbols to OpenBLAS, Rlapack
616), then the changed and nearest test files with that R,
`OPENBLAS_NUM_THREADS=4` (`dev/silent-par2.sh dev/silent-jobs-ob.txt`,
logs `dev/silent-ob-log/`):

```
RESULT test-dharma-cens.R pass=15 fail=0 error=0
RESULT test-drop1-specials.R pass=24 fail=0 error=0
RESULT test-methods-audit.R pass=58 fail=0 error=0
RESULT test-mi-rows-read.R pass=65 fail=0 error=0
RESULT test-newdata-aterms.R pass=57 fail=0 error=0
RESULT test-osa-inference.R pass=34 fail=0 error=0
RESULT test-parse.R pass=57 fail=0 error=0
RESULT test-skew-normal-start.R pass=70 fail=0 error=0
RESULT test-subset-rate.R pass=93 fail=0 error=0
RESULT test-hmm-stationary.R pass=10 fail=0 error=0
RESULT test-hmm.R pass=120 fail=0 error=0
RESULT test-cox-bhaz.R pass=147 fail=0 error=0
RESULT test-cox-bounds.R pass=16 fail=0 error=0
RESULT test-cox.R pass=45 fail=0 error=0
```

0.3.32 was not run: no new assertion sits near a rounding-decided
branch. `dev/ciharden-out/` was deleted afterwards.

## R CMD check

`R CMD check --as-cran`, once per changed package, built and checked in
`dev/silent-check/` (`dev/silent-check.ps1`; pandoc and TinyTeX on
PATH, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, `R_LIBS` = lane library,
rellib-r20, user library; logs `dev/silent-check/check-*.log`; the
`.Rcheck` trees and tarballs were deleted after):

```
frmtmb 0.81.0 (worktree)    Status: 1 NOTE   tests 32m
frmtmb.latent 0.7.3          Status: 1 NOTE   tests 25s
frmtmb.survival 0.2.3        Status: 1 NOTE   tests 18s
```

Each NOTE is the expected "checking HTML version of manual ... NOTE /
Skipping checking math rendering: package 'V8' unavailable".

## Files touched

Core `R/`: `diagnostics.R` (650), `frame.R` (`loglik_nobs()` replaces
`mi_rows_gone()`, 1562), `methods-fit.R` (`logLik()` and
`df.residual()`, 1562), `multiple.R` (pooled anova count, 1562),
`confint.R` (anova count 1562; drop1 scope, `add1.frmtmb_fit()`, 4235),
`frmtmb-package.R` (`importFrom(stats, add1)`), `parse.R` (4675),
`conditional-effects.R` (720), `predict.R` (two refusal texts and two
Rd paragraphs about `dharma_residuals()`, 650; shared area of ppdraw
and names: message and doc hunks only). Generated: `NAMESPACE`,
`man/dharma_residuals.Rd`, `man/drop1.frmtmb_fit.Rd`,
`man/frmtmb_register_aterm.Rd`, `man/residuals.frmtmb_fit.Rd`,
`man/simulate.frmtmb_fit.Rd`. `NEWS.md`.

Core tests: new `test-dharma-cens.R`, `test-drop1-specials.R`; added
blocks in `test-subset-rate.R`, `test-parse.R`,
`test-newdata-aterms.R`; one assertion changed in
`test-osa-inference.R` (the refusal text it matched).

frmtmb.latent: `R/hmm.R`, `NEWS.md`, new
`tests/testthat/test-hmm-stationary.R`.
frmtmb.survival: `R/cox-family.R`, `man/cox.Rd`, `NEWS.md`, new
`tests/testthat/test-cox-bounds.R`.

`dev/`: `test-backlog.md` (seven CLOSED lines, two filed items),
`silent-findings.md`, and the `dev/silent-*` scripts.

No `frm_compat()` row, ledger row or bin-2 manual row is changed. By
construction none flips: `nobs()` is unchanged, brms has no `drop1()`,
`add1()` or `dharma_residuals()`, and brms cannot fit the
interval-censored Cox model (brms-36).

## Bump and floors

Minor bump for core (behavior changes marked BREAKING; a new exported
S3 method `add1.frmtmb_fit`). frmtmb.latent: patch; floor stays
frmtmb (>= 0.80.0), since `frmtmb_family(post = )` and the stationary
escape exist there. frmtmb.survival: minor (BREAKING default knots on
fits with interval or truncation bounds past the data); floor stays
frmtmb (>= 0.75.0). The latent test asserts an escape that core 0.80.0
runs; it was not run against 0.80.0.

## What I did not do

- The one-pass stationary escape (filed).
- A coincident-length addition-term value on new rows (filed).
- No gated brms tier: no changed path reaches a brms fit.
- No OpenBLAS 0.3.32 pass (see above).
- `drop1()` of a `me()` term stays without statistics rather than
  being computed as a nested model with the coefficient held at 0:
  `frm()` has no way to hold one coefficient fixed (no `map`, no
  `constant()` prior), and adding one is beyond an S item.
- No user decision is needed for any item.

## Fix round (review `dev/reviews/2026-11-11-silent.md`, verdict FIX)

Each finding was reproduced on the lane build before it was fixed, with
the reviewer's scripts (`dev/reviews/silent-review-*.R`); my reruns
write `dev/silent-fix-log/`.

### MAJOR 1 (3874): a lower truncation bound moved the default knots

Reproduced (`dev/silent-fix-log/repro-coxtrunc.log`, seed 11, n = 200,
entry = time / 2, times shifted by 2):

```
lane default knots 1.000221 10.92904  internal 2.77058  frmtmb logLik -215.6089451  hand -215.6089455
brms default knots 1.825371 10.92904  hand -189.1194011
```

The reviewer is right: below the lower knot the model has no hazard for
any row, so `H(entry) = 0` there is exact, and only the upper side was
defective. Fix: `cox_read_bounds()` reads only `cens_y2` on code-2 rows
and finite `trunc_ub`; `bhaz_spec()` widens and checks only the upper
knot. After the fix (`dev/silent-fix-log/lane-coxtrunc.log`):

```
lane default knots 1.825371 10.92904  internal 2.77058  frmtmb logLik -189.0598133  hand -189.1194011
```

The knots are brms's again. frmtmb's logLik is 0.06 above the
reviewer's hand optimum on the same knots. The hand fit puts the
baseline weights through a softmax, which cannot reach the simplex
boundary, and a maximum likelihood fit of this baseline often has a
weight at exactly 0 (`?cox`). So the test does not compare optima. It
evaluates the hand-written likelihood (splines2 bases, `H(entry) = 0`
below the lower knot) at the fit's own estimates, and that agrees to
1e-8 relative. The test also asserts that the default knots are
identical to brms's, and that the same knots given explicitly are
accepted and give the identical fit. `trunc(ub = )` past the data now
has its own test (default upper knot widened, given knots refused).
The delayed-entry test passes on base too: base never widened, so it
pins the regression rather than a base defect. The cox Rd (the "brms's
bhaz() addition term" section and `@param Boundary.knots`), the refusal
text and the survival NEWS bullet say "upper" now.

### MAJOR 2 (4235): add1() of a new bar term refitted a different model

Reproduced (`dev/silent-fix-log/repro-add1re.log`): the row `1 | g`
read Df 2, LRT 128.62, the refit of `y ~ (x + 1 | g)`, where the hand
`y ~ x + (1 | g)` gives Df 1, LRT 133.5373. Fix: `add1.frmtmb_fit()`
turns a formula scope into labels itself, in either branch, and wraps
every label that holds `|` and is not already parenthesized, so the
formula scope and step()'s character scope both refit `~ . + (1 | g)`.
The row is named `(1 | g)`. step() builds its update from that name,
so it adds the term in parentheses too. After
(`dev/silent-fix-log/lane-add1re.log`): `(1 | g)  1 649.73 133.54`.
The test checks both spellings against the hand refit. drop1() has no
analogue: its scope leaves out bar terms, and a bar term given as scope
is refused as "not a subset of the term labels" (asserted).

### MINOR 1 (1562): univariate `y | mi()`

Reproduced (`dev/silent-fix-log/repro-mianova-{base,lane}.log`): base
`df.residual` 117 and anova against the complete-case fit refused ("got
120, 116"); lane 113, and anova runs. NEWS now says so. New test in
`test-mi-rows-read.R`: `nobs` attribute and `df.residual()` equal the
complete-case fit's, and `anova(mi fit, y ~ 1 complete-case)` gives the
complete-case pair's Chisq (1e-6 relative). On base it fails with
`df.residual` 117 against 107 and the anova refused ("got 120, 110").

### MINOR 2 (1562): anova() compares a count

Reproduced (`dev/silent-fix-log/repro-nobs-lane.log`, case E: "ran").
Filed under "Filed in the round of 2026-11-11 / From lane silent", not
fixed, as instructed.

### MINOR 3 (650): `method = "traditional"` and `rotation`

Reproduced (`dev/silent-fix-log/repro-trad.log`): PIT KS p 0.86,
traditional 2.32e-25, rotation 2.27e-14 on a correct fit with 94 of 300
rows right-censored. Fix: `dharma_censor_args()` refuses both, with a
reason, when the fit has a censored row; a fit without one passes them
to DHARMa as before (asserted). The Rd section says so; NEWS too.

### MINOR 4 (1562): which count each reader gives

`dev/silent-fix-nobs.R` (log `dev/silent-fix-log/lane-nobs.log`):
`nobs() 80  insight::n_obs() 80  logLik nobs 72  BIC - AIC(k = log(72)) 0`.
New section "The number of observations" in `?logLik_penalized` (the
page of `logLik()`, `AIC()` and `BIC()`) names both counts and their
readers.

### Tests of the fix round

Changed test files, base and lane (`dev/silent-fix-log/tests/`):

```
base test-dharma-cens.R      pass=12 fail=8  error=0
base test-drop1-specials.R   pass=8  fail=20 error=1
base test-mi-rows-read.R     pass=66 fail=2  error=1
base test-cox-bounds.R       pass=13 fail=6  error=0
lane test-dharma-cens.R      pass=20 fail=0  error=0
lane test-drop1-specials.R   pass=31 fail=0  error=0
lane test-mi-rows-read.R     pass=69 fail=0  error=0
lane test-cox-bounds.R       pass=19 fail=0  error=0
```

On base every new block fails except the delayed-entry case (see
MAJOR 1). The add1 block fails there on the row name `1 | g` and its
missing statistics; the behavioral failure (LRT 128.62 on Df 2) is in
the reproducer log.

Core and survival suites, one file per R process, 3 at a time, lane
library first (`dev/release/run-par.sh silent-fix dev/silent-fix.jobs
3`, log `dev/release/silent-fix.log`; every `lib:` line names
`wt-silent-lib`):

```
SILENT-FIX ran 335 of 335
frmtmb files=319 pass=23565 fail=0 err=0 skip=239 warn=0
frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0
```

Survival in one session (`dev/release/run-sessions.sh silent-fixsess 3
frmtmb.survival`, log `dev/release/silent-fixsess.log`): forward,
reverse and preload each `files=16/16 pass=575 fail=0 err=0 skip=0
warn=0`, no LEAK line.

OpenBLAS 0.3.26 with LAPACK (rebuilt; Rblas 74 symbols to OpenBLAS,
Rlapack 616; `dev/silent-ob2-log/`), the fix round's test files and the
cox files: test-dharma-cens 20/0, test-drop1-specials 31/0,
test-mi-rows-read 69/0, test-cox-bounds 19/0, test-cox 45/0,
test-cox-bhaz 147/0 (pass/fail, error 0 each). `dev/ciharden-out/`
deleted again.

R CMD check `--as-cran` of the fix round (`dev/silent-check.ps1`, logs
`dev/silent-check/check-frmtmb.log` and `check-frmtmb.survival.log`;
trees and tarballs deleted after): frmtmb `Status: 1 NOTE` (tests 30m),
frmtmb.survival `Status: 1 NOTE` (tests 20s); each NOTE is the V8
math-rendering one. frmtmb.latent did not change in the fix round, so
its earlier check stands.

Files the fix round touched: `R/confint.R` (`add1.frmtmb_fit()`, drop1
Rd), `R/diagnostics.R` (`dharma_censor_args()`, Rd), `R/methods-fit.R`
(Rd section only), `man/drop1.frmtmb_fit.Rd`,
`man/dharma_residuals.Rd`, `man/logLik_penalized.Rd`, `NEWS.md`,
`tests/testthat/test-dharma-cens.R`,
`tests/testthat/test-drop1-specials.R`,
`tests/testthat/test-mi-rows-read.R`; frmtmb.survival
`R/cox-family.R`, `man/cox.Rd`, `NEWS.md`,
`tests/testthat/test-cox-bounds.R`; `dev/test-backlog.md` (the 3874
CLOSED line, the filed anova item); this file; `dev/silent-fix-nobs.R`.
Bump and floors are unchanged.
