# Lane hollow: three ledger rows (round of 2026-11-09)

Base e0cce906, frmtmb 0.80.0, frmtmb.sample 0.27.0. Brief:
`C:/Users/adf44/source/r/round-briefs/r081-hollow.md`. Status: DONE,
uncommitted. This file is the record; every number below is pasted
from the log named beside it (`dev/hollow-out/`, `dev/release/`,
local).

Summary:
- Row 1 (`brmsfit-methods:394`): `NA` is the right answer on a fit,
  by decision with its reason now in NEWS; the printout's false reason
  on a MAP fit is fixed. Stays divergence (hollow), reason rewritten.
- Row 2 (`brmsfit-methods:469`): draws carry `$data`; the
  leave-one-out functions accept newdata that holds exactly the fitted
  rows (policy confirmed by the user, 2026-11-10). Moves to pass.
- Row 3 (`log_lik:603`, `posterior_predict:440`): decision, confirmed.
- Totals: 634 of 823 to 635 of 823; without decisions 634 of 726 to
  635 of 724. Bin 1 stays 397/49/48.
- Fix round after review `dev/reviews/2026-11-10-hollow.md` (FIX): B1
  (a variable read outside the frame was not compared), m1 (type
  change), m2 (the subset's reason), m3 (Row 1's reason) and m4
  (`make_conditions()` of draws) are fixed; see "Fix round" under Row
  2, and Row 1 (b).
- Bumps: frmtmb patch (a printout, documentation, and
  `make_conditions()` of draws); frmtmb.sample minor (a new element and
  a BREAKING acceptance). Floor: frmtmb.sample needs no new core
  feature, so `frmtmb (>= 0.80.0)` stays.
- `R CMD check --as-cran`, rerun on the fix-round build: frmtmb "Status:
  1 NOTE", frmtmb.sample "Status: 1 NOTE"; each NOTE is the expected
  "Skipping checking math rendering: package 'V8' unavailable" (logs
  `dev/hollow-out/check-core-2.txt`, `check-sample-2.txt`; tests OK,
  17 min and 75 s; the first round's `check-core.txt` and
  `check-sample.txt` gave the same Status). The check trees are
  deleted.
- frmtmb.sample in one session three ways on the fix-round build
  (`dev/release/hollow-fixsess.log`): 3 of 3, each "files=87/87
  pass=4384 fail=0 err=0 skip=1 warn=0"; 38 LEAK lines, 0 without a
  load tag.

## Row 3: custom-family dispatch is a decision (done)

Rows `log_lik:603` and `posterior_predict:440` (bin 2P, groups G3 and
G4) were "cannot transfer, not-ported, cost". The user decided on
2026-11-09 that they are a decision.

Check on the base build (rellib-r19, frmtmb 0.80.0, frmtmb.sample
0.27.0), script `dev/hollow-row3-check.R`, data seed 603, sampler seed
440, log `dev/hollow-out/row3-check.txt`. A custom family carries a
marker in each slot: `lpdf` adds 0.125 to a unit-variance normal log
density, and `sim` returns `mu + 1000`. brms-style functions
`log_lik_mycust(i, prep)` and `posterior_predict_mycust(i, prep)` that
stop when called are defined in the global environment.

```
dim log_lik: 100 30  dim epred: 100 30  dim predict: 100 30  ndraws: 100  nobs: 30
log_lik vs the family's lpdf at each draw's mu, max abs diff: 0
log_lik vs the same density WITHOUT the marker, min abs diff: 0.125
posterior_predict - epred, range: 1000 1000
sim slot calls: 100
```

So `log_lik()` and `posterior_predict()` of draws are the family's own
`lpdf` and `sim` (draws x observations, exact), and neither reaches a
function found by name. The reason holds.

Change: `dev/brmsport-bin2-survey.R` (`bp_custom`) now gives both rows
class `brms-internal` (a decision under `ct_kind` in
`dev/brmsport-bin2-ledger.R`; that map is unchanged) with the reason in
the ledger's words; the comment above `bin2port_g34` and the G3 group
reason no longer list the custom family among the cost rows. I reran
the survey script (user library, brms 2.23.0, log
`dev/hollow-out/survey-regen2.txt`); `git diff` shows exactly the two
survey rows and the one G3 group cell changed. The written record for
section 2.6 is the ledger row's reason, which the standard accepts ("the
brms-suite ledger gives a REASON"); no other file needs it.

Totals, `dev/hollow-totals.R` on the tracked ledgers (no recorded bin-2
run exists in any tree, so the ledger script itself is rerun at
consolidation), log `dev/hollow-out/totals.txt`:

```
base (tracked ledgers):            bins 1 and 2 pass 634 of 823; without decisions 634 of 726
after Row 3:                       bins 1 and 2 pass 634 of 823; without decisions 634 of 724
```

## brms 2.23.0 on rows 1 and 2

Script `dev/hollow-brms-probe.R` (user library only), seed 394, log
`dev/hollow-out/brms-probe.txt`.

Row 1, `hypothesis(fit1, "Intercept = 0", class = "sd", group =
"visit")` on `brmsfit_example1`: the fixture was fitted with
`sample_prior = "yes"` and proper priors (`cauchy(0,2)` on `sd`,
`normal(0,2)` on `b`); it carries 25 prior draws per class column
(`prior_sd_visit` among them) beside 25 posterior draws.

```
       Hypothesis    Estimate   Est.Error     CI.Lower    CI.Upper   Evid.Ratio
1 (Intercept) = 0 3.202570893 1.206439106 0.7883270776 4.906793149 0.1169630956
density_ratio(post, prior, 0) = 0.1169630956
posterior kde at 0 = 0.01615646064
prior kde at 0 = 0.1381329773
```

brms's number is the Savage-Dickey ratio: the kernel density of the 25
posterior draws at 0 over that of the 25 prior draws at 0. A
directional row is the posterior odds (`Intercept > 0`: `Inf`, all 25
draws above 0).

Row 2, fixture 6 (`volume ~ Trt + gp(Age, by = Trt, gr = TRUE)` and
`count ~ Trt + Age`, no rescor): `fit6$data` is a data.frame of 40 rows
and the model's 4 variables (`volume Trt Age count`) with a `terms`
attribute.

```
loo(fit6) estimates:                    elpd_loo -184.0806004720  p_loo 14.1590059834
loo(fit6, newdata = fit6$data):         elpd_loo -184.0810121448  p_loo 14.1592622575
log_lik identical: FALSE  max abs diff: 0.001109178
newdata = first 20 rows: 20 pointwise rows, max abs diff vs loo(fit6) rows 1:20: 0.0003753988
```

brms computes `log_lik(fit, newdata)` and runs PSIS-LOO on it. With the
fitting data that is `loo(fit)` up to brms's re-prediction of the
`gp()` on new rows (4.1e-4 in the elpd); on a subset of the fitting
rows it is those rows' pointwise elpd. On rows the model was not fitted
on the same arithmetic runs, and its result is not a leave-one-out
quantity: a weight `1 / p(y_new | theta)` removes an observation the
posterior never contained.

## Row 1: `brmsfit-methods:394`, Evid.Ratio of a point hypothesis

### What frmtmb reports, and whether it is right

Script `dev/hollow-row1-probe.R` on the base build (rellib-r19), data
seed 394 (120 rows, 12 groups), sampler seed 3940 (4 chains of 2000,
1000 warmup), log `dev/hollow-out/row1-probe-base.txt`.

(a) Maximum-likelihood fit, flat priors. `Evid.Ratio` and `Post.Prob`
are `NA` for `x = 0`, `x > 0` and `sd_g__Intercept = 0`. RIGHT: under a
flat prior the prior has no density, so a Savage-Dickey ratio does not
exist; brms itself gives `NA` without prior draws, and a flat class has
no prior draws (frmtmb.sample's `test-evidence-ratio.R` pins that on
draws). A directional row: brms's ratio is the posterior odds from
draws. A fit has none; its frequentist answer is the one-sided p in
`attr(h, "test")`. Kept `NA`.

(b) MAP fit, `normal(0, 1)` on `b` and `cauchy(0, 2)` (lb 0) on `sd`.
The columns are `NA` too, and the printout said "Evid.Ratio and
Post.Prob are NA, because a maximum-likelihood fit has no posterior",
which is FALSE for a MAP fit. The value is right; the stated reason is
not. The only posterior density a MAP fit has is the normal
approximation at the mode, and at a tested point, which is usually in
the tail, it depends on the scale of the approximation:

```
x = 0, normal on b_x's own scale: 0.00598509
x > 0, normal posterior odds: 17685.3
sd = 0, normal on the sd scale: 0.00751481
sd = 0, normal on the log-sd scale (the optimizer's): density of sd at 0 is 0, ratio 0
```

against (c) on the same model and priors:

```
1    (x) = 0 ... Evid.Ratio 4.94862e-03
1 (Intercept) = 0 (sd_g) ... Evid.Ratio 2.93424e-06
```

For the slope the normal gives 1.21 times the draws' ratio; for the
group sd, 2560 times on the sd scale. Decision: keep `NA` on every fit
and say why. Fixed: the printout (now true for both kinds of fit) and
`?hypothesis` (the reason per case and the route to brms's number).
NEWS gives the reason, which section 2.6 needs.

FIX ROUND (review `dev/reviews/2026-11-10-hollow.md`, item 4): the
first wording ("the normal approximation changes with the scale, so it
gives no single number") was WITHDRAWN. The Savage-Dickey ratio itself
does not depend on the chart (the Jacobian cancels between the
posterior and prior densities); only its normal approximation does, and
frmtmb already reports the sd's Wald interval on the sd scale, so the
chart is chosen. The reason that holds for every parameter class is the
tail error: the ratio needs the posterior density AT the point, which
is usually in the approximation's tail, where it was 2560 times the
draws' ratio for the sd. NEWS, `?hypothesis` and the ledger row
(`dev/brmsport-verdicts-manual-fit.tsv`, `brmsfit-methods:394`) now say
that: a fit does not report a Bayes factor from the normal
approximation, for any parameter.

NOT TAKEN: the ratio for identity-chart rows only (coefficients), where
the normal approximation is the single one and gave 0.00599 against
the draws' 0.00495 (1.21 times) on this model. It would be a number
brms never computes from a fit, its error would still grow into the
tail, and a sd row beside it would stay `NA`. The coordinator kept
`NA` on ML and MAP fits; the alternative is recorded here for the user.

(c) Draws with proper priors and `sample_prior = "yes"`: implemented in
frmtmb.sample before this round (`R/evidence-ratio.R`). On the probe
model, frmtmb's `Evid.Ratio` is IDENTICAL (`identical()` TRUE) to
`brms:::density_ratio()` on the same posterior and prior draws, for
`x = 0` (0.004948622335) and for the sd (2.934235824e-06). The prior
draws of the sd follow the half-Cauchy(0, 2) (KS p 0.4477, 4000
draws). `frm_sample()` of the MAP fit took the fit's priors without
restating them. Nothing was missing here.

### The row

The bin-1 fixture `fit1` is a maximum-likelihood fit; brms's is sampled
with proper priors and prior draws. So the right answer on frmtmb's
fixture is `NA`, and `is.numeric(NA_real_)` holds for that reason, not
because frmtmb computes brms's number. The row stays `divergence`,
class `hollow` (the bin-1 ledger script accepts a held run with a
non-pass verdict only under that class), with the reason rewritten to
meet section 2.6: NEWS ("Hypothesis tests") now gives it. Hand-written
row changed: `dev/brmsport-verdicts-manual-fit.tsv`,
`brmsfit-methods:394`. The fixture's own sd estimate is not used.

Pin: `tests/testthat/test-map.R`, "hypothesis() on a MAP fit does not
call it maximum likelihood". Base (rellib-r19): fail=6 (the three
printout assertions, on the MAP and on the ML fit); lane: pass=16
fail=0 (the file). Logs `dev/hollow-out/t-map-base.txt`,
`t-map-lane.txt`.

## Row 2: `brmsfit-methods:469`, `loo(fit6, newdata = fit6$data)`

### Decision (a): draws carry `$data`

Readers of `$data`: in frmtmb.sample, none (grep of `R/`). In core,
`make_conditions()` reads `x[["data"]]` from any non-fit list that has
it, so `make_conditions(draws, ...)` read the model frame, which fails
on a transformed term (review item 5: "`x` is not a column of the
data. It has: y, poly(x, 2)"; base failed on every model). Fix round:
core's `make_conditions()` now sends draws through their fit, as for a
fit (`R/brms-utilities.R`); test "make_conditions() of draws reads the
fit's data, as of the fit" in `tests/testthat/test-brms-utilities.R`,
error=1 on the unfixed lane build, pass after (logs
`dev/hollow-out/t-util-prefix.txt`, `t-util-fix.txt`). No
extension reads `$d` or another prefix of `data` from a draws object
(grep of every extension's `R/`), so no partial match moves. brms's
suite reads `fit$data` on 10 lines (`fit$data[1:10, ]` as newdata,
`fit1$data$visit`, `model.frame(fit1) == fit1$data`, loo's newdata); on
draws every one read `NULL`. Decision: yes. `frm_sample()` stores `data
= fit[["data"]]`, the fit's model frame, the same object in memory.
Cost: a saved draws object carries the frame once more (serialization
does not deduplicate; the fit already carries it twice, see `?frm`).

### Decision (b): `loo(newdata =)` on the fitted rows follows brms

brms computes `loo()` on `log_lik(newdata =)`. On the fitted rows that
is `loo(fit)`, up to its gp re-prediction. On other rows it runs the
same arithmetic, and the result is not a leave-one-out quantity. The
recorded policy's reason ("a leave-one-out weight is the observation's
own column of log_lik(), so it exists only for a row the model was
fitted on") therefore holds for other rows and argues FOR the fitted
rows. Decision: accept `newdata` that holds exactly the fitted rows,
and give the answer without `newdata`; refuse every other `newdata` as
before, with one sentence added that names the accepted form.

CONFIRMED BY THE USER, 2026-11-10 (relayed by the round's
coordinator): the leave-one-out functions accept `newdata` only when it
holds exactly the fitted rows, and a subset of the fitted rows stays
refused.

"Exactly the fitted rows" (`draws_newdata_is_fitted()`,
`extensions/frmtmb.sample/R/methods-draws.R`), as of the fix round:
the row count of the model frame, and either the frame object itself
(`identical(newdata, fit$data)`: `x$data`, `model.frame()`), or data
that holds every variable the model reads, each equal to its fitted
value (`frame$data_vars`, every column of `frame$raw_vars`, every
variable of the frame's terms, every frame column that is a variable),
with every frame column it also holds equal too. Equality is exact
(`all.equal(tolerance = 0)`, with the same class on both sides);
factors and characters compare by label only when both sides are a
factor or a character. A rebuilt frame is NOT compared:
`poly(x, 2)` rebuilt from its stored coefficients differs from the
fitted basis by a mean relative 2.82e-16 (`dev/hollow-terms-probe.R`),
so an exact test would refuse the fitting data, and a tolerant one
could not tell other data that differs by 2 ulps. A variable with no
stored value refuses (fails closed).

Applied in `log_lik()` (and so `logLik()` and pp_check's `loo_*`
types), `psis()`, `loo()`, `waic()` and `add_criterion()`
(`loo_criterion()`), the `loo_*()` functions and `loo_R2()`
(`loo_draws_opts()`), and `loo_model_weights()` and
`posterior_average()` (`draws_model_weights()`), through
`loo_arg_asks()`.

On the bin-2 stand-in for fixture 6 (built as
`helper-brms-suite-draws.R` builds it, sampler seed 20260917), script
`dev/hollow-row2-probe.R`, logs `dev/hollow-out/row2-base.txt` and
`row2-lane.txt`:

```
base: names(fit6): stanfit draws fit reparam warmup ...   is.null(fit6$data): TRUE
base: loo(fit6, newdata = fit6$data): elpd_loo -183.528521804 identical to loo(fit6): TRUE   (NULL newdata)
base: loo(fit6, newdata = model.frame(fit6)): ERROR loo() cannot honor `newdata` ...
base: loo(fit6, newdata = d), the data fitted: ERROR loo() cannot honor `newdata` ...
lane: names(fit6): stanfit draws fit data reparam ...      is.null(fit6$data): FALSE
lane: loo(fit6, newdata = fit6$data): elpd_loo -183.528521804 identical to loo(fit6): TRUE
lane: loo(fit6, newdata = model.frame(fit6)): ... identical to loo(fit6): TRUE
lane: loo(fit6, newdata = d), the data fitted: ... identical to loo(fit6): TRUE
lane: loo(fit6, newdata = d with one response moved): ERROR loo() cannot honor `newdata` ... Newdata that holds exactly the rows the model was fitted on (x$data) is accepted, and gives the answer without newdata.
```

So the row now holds because frmtmb computes the elpd brms computes on
the fitted rows. It moves from divergence (hollow) to pass: its manual
verdict line is removed from `dev/brmsport-bin2-verdicts-manual.tsv` (a
held run with a manual verdict of a non-held class stops the ledger
script). The policy rows it reaches, `brmsfit-methods:508` and `:544`
(new rows with `allow_new_levels = TRUE`), stay refused; their reason
now reads "refuse newdata other than the fitted rows, as loo() does".

Pin: `extensions/frmtmb.sample/tests/testthat/test-loo-newdata-fitted.R`
(data and sampler seed 469; the transformed-term case, data seed 4691).
Base: pass=21 fail=4 error=2 (`$data` absent, the fitted rows refused,
the new sentence absent, the helper absent); lane: pass=40 fail=0. The
guard's absent cases are in the file: a moved response, `x` times
(1 + 2e-15) in one row, a subset, a reordering, an extra row and a
missing response are each refused by `loo()`, `log_lik()` and
`loo_epred()`; `x + 1e-9` behind `poly()` and `age + 0.5` behind
`offset()` are not taken for the fitted rows. Logs
`dev/hollow-out/t-nf-base.txt`, `t-nf-lane.txt`.

Five existing assertions passed the fitting data to provoke a refusal.
They now pass other rows (a subset, or the data less its first row) and
assert the same refusal: `test-loo.R` (logLik), `test-loo-predict.R`
(loo_model_weights, loo_R2), `test-loglik-autocor.R` (two), and
`test-ppcheck-loo.R`.

### Fix round (review `dev/reviews/2026-11-10-hollow.md`, FIX)

B1 (BLOCKER, fixed). The first helper compared only the frame's columns
when newdata held them all, so a variable read outside the frame was
never compared. `frm(y | cens(cc, y2) ~ x)` holds `y, x, cc` in its
frame and `y2` beside it; the data with `y2[4] + 3` was accepted, and
`log_lik(ds, newdata = d1)` returned the fitted matrix where row 4's
density on the moved interval is -3.558393 against -3.625350
(reviewer's `review-hollow-probe3.log`). The same with `fit$data` plus a
different `x` on `y ~ poly(x, 2)`. Now every path compares every
variable the model reads, `frame$raw_vars` included, and refuses one
that is absent, except for the frame object itself.

m1 (fixed). A numeric fitted column against a factor or character one
is refused; labels compare only when both sides are factor or character,
and other columns need the same class.

m2 (fixed). A subset of the fitted rows stays refused, now with its
true reason in NEWS, the refusal and the docs: rows are matched by
value, and a subset cannot be placed among the fitted rows
unambiguously where rows repeat and a per-row latent value (`me()`,
`mi()`) or an autocorrelation neighbor enters. The row-name match for
models without those terms is filed in `dev/test-backlog.md`, "Filed in
the round of 2026-11-09 / From lane hollow".

Seen to fail on the unfixed lane build (log
`dev/hollow-out/t-nf-prefix.txt`): `test-loo-newdata-fitted.R` pass=46
fail=6. The six: the frame plus a different `x` (poly), the moved `y2`
in the data, the frame plus the moved `y2`, a plain copy of the frame
without `y2`, a factor `k`, a character `k`. After the fix: pass=53
fail=0 (`t-nf-fix.txt`).

The reviewer's case lists rerun on the fixed build (scripts
`dev/review-hollow-probe1.R`, `probe2.R`, `probe3.R`, `probe5.R`
unchanged, logs `dev/hollow-out/rev-probe<k>-fix.txt`), diffed against
the reviewer's logs. Every line is the same except these:

```
== probe1
<    ACCEPTS change to y2          (cens model; now refused)
== probe2
< cens: upper bound of an interval row changed, accepted  TRUE
> cens: upper bound of an interval row changed, accepted  FALSE
< numeric fitted, factor newdata with same labels         TRUE
< numeric fitted, character newdata                       TRUE
> numeric fitted, factor newdata with same labels         FALSE
> numeric fitted, character newdata                       FALSE
< poly: fit$data plus a different x column                TRUE
> poly: fit$data plus a different x column                FALSE
== probe3
< cens changed bound: accepted; identical to log_lik(ds): TRUE
> cens changed bound: refused: log_lik() does not take newdata. ...
< make_conditions(draws, 'x'): make_conditions(): `x` is not a column of the data. It has: y, poly(x, 2)
> make_conditions(draws, 'x'): -1.01354  0.09759  1.20872
== probe5
(no difference)
```

So every case the reviewer saw refused is still refused, and every
true-data case the reviewer saw accepted is still accepted: the data
itself, `fit$data`, a tibble, a data.table, extra columns, `g` as
character, reordered factor levels, integer for double, `na.omit()`
of data whose NA rows the fit dropped, NA in an unused column, a
formula-environment variable given as a column, and in probe1 every
model's data and `fit$data`.

## Ledger rows that move

`dev/hollow-totals.R` (the rules of `dev/brmsport-bin2-ledger.R` on the
tracked ledgers), log `dev/hollow-out/totals.txt`:

```
base (tracked ledgers):            bins 1 and 2 pass 634 of 823; without decisions 634 of 726
                                   bin 1 pass/divergence/cannot transfer 397/49/48
after Row 3:                       bins 1 and 2 pass 634 of 823; without decisions 634 of 724
after Rows 1 and 2 as well:        bins 1 and 2 pass 635 of 823; without decisions 635 of 724
                                   bin 1 pass/divergence/cannot transfer 397/49/48
                                   divergence (all bins) 84; hollow class 1
```

- `log_lik:603`, `posterior_predict:440` (bin 2P): cannot transfer,
  not-ported (cost), to cannot transfer, brms-internal (decision).
- `brmsfit-methods:469` (bin 2T, frmtmb.sample): divergence (hollow)
  to pass.
- `brmsfit-methods:394` (bin 1): stays divergence (hollow), reason
  rewritten. Bin 1 stays 397/49/48.
- `brmsfit-methods:508`, `:544`: verdict unchanged, reason edited.

## Test runs

One file per process, 3 at a time, every suite (core and the eight
extensions), lane library ahead of rellib-r19, `NOT_CRAN=true`,
`FRMTMB_BRMS_FIT_TESTS=true`: `bash dev/release/run-par.sh hollow-full
dev/hollow-out/full-jobs.lst 3`, log `dev/release/hollow-full.log`:
"HOLLOW-FULL ran 510 of 510". Every `lib:` line names `wt-hollow-lib`
for frmtmb (510 of 510).

```
frmtmb files 317 pass 26312 fail 0 err 0 skip 15 warn 0
frmtmb.sample files 87 pass 4371 fail 0 err 0 skip 1 warn 0
frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0
frmtmb.eam files 33 pass 1854 fail 0 err 0 skip 3 warn 0
frmtmb.latent files 12 pass 687 fail 0 err 0 skip 2 warn 0
frmtmb.learn files 16 pass 545 fail 0 err 0 skip 2 warn 0
frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0
frmtmb.spatial files 8 pass 561 fail 0 err 0 skip 0 warn 0
frmtmb.survival files 15 pass 552 fail 0 err 0 skip 0 warn 0
```

The 29 skips are the gated tiers: 15 "set FRMTMB_SCALE_TESTS=true", 13
"set FRMTMB_DRMTMB_FIT_TESTS=true", 1 "set FRMTMB_FUZZ=true". The known
red of the base, `test-se-tier.R:323` and `:327`, did not fail here.

One session per suite, three ways (`bash dev/release/run-sessions.sh
hollow-sess 3`, lane library ahead), log `dev/release/hollow-sess.log`:
"HOLLOW-SESS ran 27 of 27", every SESSION line fail=0 err=0 warn=0:

```
SESSION frmtmb forward files=317/317 pass=26312 fail=0 err=0 skip=15 warn=0
SESSION frmtmb reverse files=317/317 pass=26462 fail=0 err=0 skip=15 warn=0
SESSION frmtmb preload files=317/317 pass=26555 fail=0 err=0 skip=17 warn=0
SESSION frmtmb.sample forward files=87/87 pass=4371 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample reverse files=87/87 pass=4371 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample preload files=87/87 pass=4371 fail=0 err=0 skip=1 warn=0
```

The other seven extensions' 21 sessions give their per-file counts
above. 194 LEAK lines, 0 without a `[with load of ...]` tag.

Fix round, one file per process, 3 at a time, every frmtmb.sample file
and core's `test-map.R`, `test-brms-utilities.R`, `test-refit-data.R`
and `test-brms-suite-brmsfit-helpers.R` (the files that reach
`make_conditions()`): `bash dev/release/run-par.sh hollow-fix
dev/hollow-out/fix-jobs.lst 3`, log `dev/release/hollow-fix.log`,
"HOLLOW-FIX ran 91 of 91", every `lib:` line names `wt-hollow-lib`:

```
frmtmb.sample files 87 pass 4384 fail 0 err 0 skip 1 warn 0
frmtmb files 4 pass 172 fail 0 err 0 skip 0 warn 0
```

The one skip is the scale tier. The earlier whole-repository run above
is of the build before the fix round.

OpenBLAS: not run. No numerical code changed (a printout, a list
element, and an exact equality test on data), and the new assertions
compare two results computed in one process.

## Not done, and why

- A normal-approximation evidence ratio on a MAP fit: rejected, see
  Row 1 (b).
- `newdata` as a SUBSET of the fitted rows: brms gives those rows'
  pointwise elpd. Not accepted: rows are matched by value, and a subset
  cannot be placed among the fitted rows unambiguously where rows
  repeat and a per-row latent value (`me()`, `mi()`) or an
  autocorrelation neighbor enters. It is refused with that reason, and
  the row-name match for models without those terms is filed in
  `dev/test-backlog.md`.
- `pp_mixture(newdata =)` still refuses the fitted rows: it is not a
  leave-one-out function, and its refusal has its own reason.
- The generated port files (`tests/testthat/test-brms-suite-*.R`,
  `dev/brmsport-bin2/`) and the generated ledgers are not edited; the
  consolidation regenerates them.
