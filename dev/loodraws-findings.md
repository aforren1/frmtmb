# Lane loodraws: brms's loo_*() functions on draws

Round opened 2026-10-20, base 4854e425 (frmtmb 0.72.0, frmtmb.sample
0.20.0). Scope: frmtmb.sample only. Core is untouched.

Status: fix round 1 done (section "Fix round 1" at the end). Sections
before it describe the first submission; where fix round 1 changed a
result, the section says SUPERSEDED and the new result is at the end.

## What changed

frmtmb.sample:

- `R/loo-predict.R` (new): `loo_epred()`, `loo_linpred()`,
  `loo_predict()`, `loo_predictive_interval()`, `loo_R2()` and
  `loo_model_weights()`: generics (bare `UseMethod()`, owner formals),
  `frmtmb_draws` methods, `frmtmb_fit` refusals that name
  `frm_sample()`. Shared internals: `loo_draws_opts()` (sorts `...`),
  `draws_psis_of()` (the one `loo::psis()` call, with `loo_r_eff()`),
  `loo_draws_weights()`, `loo_draws_expect()`, `loo_expect_value()`
  (brms's `E_loo_value()`), `loo_R2_boot()` (brms's `.loo_R2()`),
  `draws_model_weights()` (brms's `loo_model_weights.brmsfit()`).
- `R/loo.R`: `psis.frmtmb_draws()` goes through `draws_psis_of()`;
  `loo_matrix()` passes `...` to `log_lik()` (for `draw_ids`);
  `add_criterion(criterion = "loo_R2")` stores `loo_R2(summary =
  FALSE)` instead of refusing; its doc. `brms_resp_index()`, brms's
  `validate_resp()` for `bayes_R2()` and `loo_R2()` alike (one message
  template).
- `R/postfit-draws.R`: `pa_weights()` (posterior_average's
  `"stacking"`/`"pseudobma"`) calls `draws_model_weights()`, the path
  `loo_model_weights()` takes (brief item 5). BREAKING in NEWS: it no
  longer reads a stored `loo` criterion. Without one the weights are
  identical: stacking 0.81370648034831772 / 0.18629351965168228 and
  pseudobma 0.83117464540834785 / 0.16882535459165210 by both routes
  (`dev/loodraws-probe-pa.R`, seeds 1, 2, 5; log
  `dev/loodraws-log/probe-pa.txt`).
- `R/generic-owners.R`: six owner entries; count comment 24 -> 30.
- `R/sample.R`: the laplace-refusal list names the `loo_*()` functions.
- `NAMESPACE`, `man/` (4 new Rd: `sample-loo_predict`, `sample-loo_R2`,
  `sample-loo_model_weights`, `frmtmb_fit-loo_predict`; 3 regenerated).
- `NEWS.md`: a `# frmtmb.sample (development version)` section.
- Tests: `tests/testthat/test-loo-predict.R` (new, 13 blocks);
  `test-generic-collision.R` (six names in `own_generics` and
  `owner_table`; a new block, rstantools attached first, brms absent);
  `test-add-criterion.R` (the `loo_R2` refusal assertion becomes a
  storage assertion).

dev/: `upstream-bugs.md` (brms-30), this file, and the `dev/loodraws-*`
scripts named below.

## Generic ownership (brief item 2)

Measured with `dev/loodraws-probe-generics.R` on the installed
versions (rstantools 2.7.1, loo 2.10.1, brms 2.23.0): the environment
of each exported function.

| name | defined by | brms's export comes from |
|---|---|---|
| loo_R2, loo_linpred, loo_predict, | rstantools | rstantools |
| loo_predictive_interval | rstantools | rstantools |
| loo_epred | rstantools AND brms (two closures) | brms |
| loo_model_weights | loo | loo |

So `loo_epred` has rival owners, ordered `c("brms", "rstantools")`
(brms's own, which is where `loo_epred.brmsfit` is registered), as
`rhat` and `posterior_samples` are. Each `frmtmb_draws` and
`frmtmb_fit` method is registered in this package's table and in
every owner's table (`S3method(brms::loo_epred, )` and
`S3method(rstantools::loo_epred, )` both).

Dispatch measured by `test-generic-collision.R` (child processes, the
probe reproduces UseMethod()'s table lookup): brms attached after,
before, loaded never attached, imported by a package, detach and
reattach; every owner attached in both orders; no owner usable (shadow
packages); and the new block, rstantools attached first with brms
absent, where `get("loo_epred")` is rstantools' generic and every draws
method is reached.

Measured directly with `dev/loodraws-probe-dispatch.R <arm> <order>`
(hand-built draws of a 40-row gaussian fit, seed 1; logs
`dev/loodraws-log/dispatch-<arm>-<order>.txt`), the generic each name
resolves to and what a call on draws returns:

| order | lane | base (rellib-r10) |
|---|---|---|
| brms attached after | loo_epred from brms, four from rstantools, loo_model_weights from loo; all six return (40 x 1, 40 x 2, 1 x 4, weights named ds,ds) | five "no applicable method for ... frmtmb_draws"; loo_model_weights reaches loo's default: "'x' must be a list of matrices or a list of 'psis_loo' objects." |
| brms not loaded | five from frmtmb.sample, all return | "could not find function" for all six |
| rstantools attached first | five from rstantools (loo_epred too), all return | five "no applicable method"; loo_model_weights not found |

## Exactness (evidence item 1)

`dev/loodraws-probe-exact.R` (seed 1 data, 60 rows; samplers seed 11,
2 chains of 300; log `dev/loodraws-log/probe-exact.txt`): on
gaussian, Poisson, random-intercept and `mvbf()` (resp `y2`) draws,
58 checks, 57 TRUE (the one FALSE is the thinned pp_check line,
expected, below), among them:

- `psis(ds)` weights identical to `loo::psis(-ll, r_eff =
  relative_eff(exp(ll), chain_id))` by hand, and to
  `loo(ds, save_psis = TRUE)$psis_object`.
- `loo_epred()` and `loo_linpred()` for mean, var and quantile
  identical to `loo::E_loo()` by hand, and `loo_epred()` identical to
  `brms:::E_loo_value()` on the same matrices (12 of 12).
- `loo_predict(type = "quantile")` identical to `E_loo()` of
  `posterior_predict()` at one seed; `loo_predictive_interval(prob =
  0.9)` identical to `loo_predict(probs = c((1 - 0.9) / 2, 1 - (1 -
  0.9) / 2))`. Note `(1 - 0.9) / 2 != 0.05` in doubles: brms forms the
  probabilities this way too, so the interval is NOT `probs = c(0.05,
  0.95)` bitwise (it differed on 3 of 4 models when the probe used 0.05).
- mv without `resp`: dim 60 x 1 x 2, dimnames `y, y2`; each layer
  identical to `E_loo()` of that response's `posterior_epred()` under
  the joint weights.
- `loo_model_weights()` stacking and pseudobma (seed 5) identical to
  `loo::loo_model_weights()` on the `log_lik()` matrices with the
  `relative_eff()` list; class `stacking_weights`; names deparsed.
- `loo_R2()` against `brms:::.loo_R2()` at the same seed: max relative
  difference 7.98e-15, 7.26e-15, 2.7e-15, 8.45e-15 (36, 33, 12 and 38
  eps of the largest value). The difference is the arithmetic, not the
  bootstrap: brms forms the weighted sums with `rowSums(sweep())`, which
  builds three more S x N matrices; this uses `w %*% cbind(y, y^2, e,
  e^2)`. The test pins the method against a by-hand reference with the
  same arithmetic (`identical()`).
- pp_check's `lw` for `loo_pit_overlay` on every draw identical to
  `weights(psis(ds), log = TRUE)`.

### pp_check() and loo() disagree on thinned draws (pre-existing; SUPERSEDED by fix round 1, decision (a))

On a `draw_ids` subset, `pp_check()`'s weights take `relative_eff()` on
ONE chain (brms's `r_eff_log_lik()`), while `loo()` and `psis()` drop
r_eff (it becomes 1) because the subset has no chain structure: max
|lw difference| 0.0391 on the random-intercept draws above. Both are
documented decisions of earlier lanes (sampfix; loo). The `loo_*()`
functions follow `loo()` and `psis()`, so `loo_epred(ds, ndraws =)`
equals `loo(ds, ndraws =, save_psis = TRUE)`'s weights. NOT changed
here: it would move `loo(ndraws =)`'s numbers or `pp_check()`'s, and
the choice between them (brms uses one chain everywhere) belongs to
consolidation. Needs a decision.

## Parity with brms 2.23.0 (evidence item 2; SUPERSEDED by the same-draws identity, fix round 1)

`dev/loodraws-parity.R` (data seed 1, 60 rows; brms `brm(chains = 4,
iter = 2000, seed = 1)`, frmtmb.sample `frm_sample(chains = 4, iter =
2000, warmup = 1000, seed = 1)`; logs `dev/loodraws-log/parity-*.txt`).
z is the per-row difference over sqrt(MCSE_f^2 + MCSE_b^2), with the
self-normalized importance-sampling MCSE of each side's own weights
divided by loo's r_eff. Rows of one model share the draws, so the z
values are correlated; the mean z is close to one parameter's MC
difference.

```
gaussian loo_epred mean    n=60 mean z=  0.26 sd z= 0.53 max|z|= 1.50 |z|>3: 0
gaussian loo_linpred mean  n=60 mean z=  0.26 sd z= 0.53 max|z|= 1.50 |z|>3: 0
poisson loo_epred mean     n=60 mean z= -0.97 sd z= 0.39 max|z|= 1.46 |z|>3: 0
poisson loo_linpred mean   n=60 mean z= -0.91 sd z= 0.41 max|z|= 1.38 |z|>3: 0
ranef loo_epred mean       n=60 mean z=  0.19 sd z= 1.20 max|z|= 2.83 |z|>3: 0
ranef loo_linpred mean     n=60 mean z=  0.19 sd z= 1.20 max|z|= 2.83 |z|>3: 0
mvbf loo_epred mean        n=60 mean z=  0.16 sd z= 0.41 max|z|= 1.23 |z|>3: 0
mvbf loo_linpred mean      n=60 mean z=  0.16 sd z= 0.41 max|z|= 1.23 |z|>3: 0
```

loo_R2 (seed 7), frmtmb vs brms Estimate, and the difference over
brms's Est.Error (the bootstrap SD, the scale the summary reports):
gaussian 0.1631 / 0.1614 (0.018), poisson 0.0990 / 0.0963 (0.024),
ranef 0.5222 / 0.5281 (-0.079), mvbf resp y2 0.0267 / 0.0272 (-0.008);
mv without resp, rows `R2y R2y2` on both, 0.1626 / 0.1657 and 0.0248 /
0.0254 (brms Est.Error 0.091, 0.068).

Shapes equal brms's on every model: `loo_epred` 60 x 1, quantile
columns `q25 q75`, `loo_predictive_interval` columns `q5 q95`,
`loo_R2` 1 x 4 summary and 4000 x 1 draws with brms's column names,
mv `loo_epred` 60 x 1 x 2 with dimnames `mean` and `y, y2`.

loo_model_weights (gaussian `y ~ x` vs `y ~ x + (1 | g)`): stacking
0.0000 / 1.0000 on both; pseudobma (seed 5) 0.0002 / 0.9998 frmtmb,
0.0001 / 0.9999 brms. Class `stacking_weights` on both.

## Decisions, with the measurement

(Decision 1 and decision 6 are SUPERSEDED by fix round 1: brms's
default values pass, and the generics moved to core.)

1. newdata (brief item 4). Every `loo_*()` refuses `newdata`,
   `re_formula`, `re.form`, the new-level arguments and `nug`, in one
   sentence that names the function called, as `loo()` refuses them.
   NULL values (brms's defaults) pass. brms passes `newdata` to `loo()`
   AND to the prediction, so its weights are computed from the new
   rows' own response (rows 508 and 544 call it with
   `allow_new_levels = TRUE`). `loo_model_weights()` refuses it too.
   `loo_R2()` refuses it, where `bayes_R2()` takes it: bayes_R2's R2 of
   new rows is a prediction score, while a LOO weight needs the row's
   own log-likelihood column.
2. `ndraws`/`draw_ids` choose the draws once for both halves. brms
   draws two separate random subsamples (upstream-bugs brms-30:
   2 of 100 draws shared; the brms result equals subsample-1 weights on
   subsample-2 predictions, max abs diff 0; the same-draw answer
   differs by up to 8.7 percent of its largest value;
   `dev/loodraws-brms-ndraws.R`, log `dev/loodraws-log/brms-ndraws.txt`).
3. `...` is sorted, not swallowed: `ndraws`/`draw_ids` to both halves,
   the prediction function's own arguments to it, `r_eff`/`cores` to
   `loo::psis()`, `moment_match`/`reloo` refused when TRUE (as `loo()`
   refuses them), `k_threshold`, `pointwise`, `moment_match_args`,
   `reloo_args` accepted and unused. Anything else is an error naming
   it. brms 2.23.0 ignores `scale` (row 526): brms's
   `loo_linpred(fit, scale = "response")` is identical to
   `loo_linpred(fit)` (`parity-compare.txt`, "scale ignored ... TRUE").
   Refused here by name. `save_psis` (brms errors on it, a duplicate
   argument) is refused by name too.
4. Multivariate without `resp`: brms's layout, one layer per response,
   joint weights; with `resp`, that response's own `log_lik()`. Built
   in `loo_draws_expect()` by calling the prediction once per response,
   because frmtmb.sample's `posterior_epred()`/`posterior_predict()`
   with no `resp` return the FIRST response only (see defects).
5. `loo_R2()`: brms's Bayesian bootstrap with `rexp()` in brms's order,
   seed restored afterwards, clipping to [-1, 1], per-response
   `log_lik(resp =)`. Ordinal and categorical refused, as `bayes_R2()`
   refuses them (brms refuses categorical and warns on ordinal). A
   stored `add_criterion("loo_R2")` matrix is read when the call passes
   nothing in `...`, `args_epred` or `args_loglik`; brms reads it
   whatever else is passed, which would return draws of another call.
6. ML fits (brief item 6): brms has no ML fit; each of the six refuses
   on a `frmtmb_fit` and names the route, as `log_lik()`, `loo()`,
   `bayes_R2()` do. The methods live in frmtmb.sample, because the
   generics do (the brief keeps core untouched). So with frmtmb alone
   attached there is no `loo_predict()` at all (as for `kfold()` and
   `reloo()` today), and with rstantools attached and frmtmb.sample not
   loaded, `loo_predict(fit)` reaches rstantools' generic with no
   method. Moving the six generics to core, as `bayes_R2()` and
   `log_lik()` are, would close both; it touches core's `R/loo.R`,
   `R/generic-owners.R`, `NAMESPACE` and man/. Not done; flag for
   consolidation.
7. Compat rows: `R/compat.R` has no row for `loo()`, `waic()`,
   `log_lik()`, `bayes_R2()` or `pp_check()`; its post-fit features are
   `fitted`, `predict`, `residuals`, `simulate` and the like
   (feature list read with grep). So none was added for the six, and
   core is unchanged.

## Time and memory at 4000 draws x 1e4 rows (brief item 3)

`dev/loodraws-perf.R` (data seed 1, gaussian `y ~ x`, 1e4 rows; both
samplers 4 chains x 1000 draws, seed 1; log `dev/loodraws-log/perf.txt`;
loo 2.10.1, `mc.cores = 1`). Memory is gc()'s "max used" after a reset,
minus what was live before the call. While it ran, another lane's R
process (`dev/cdfs-dbg3.R`) held 45 GB of the 63 GB, so the clock is
loaded; the control, frmtmb's path run again, spread 54.9 to 63.5 s.

Arm A, the same matrices (frmtmb's `log_lik()` and
`posterior_predict()`, 4000 x 10000; they took 7.3 s and 6.6 s),
three interleaved rounds:

| arm | s (rounds 1, 2, 3) | min s | peak MB |
|---|---|---|---|
| frmtmb: `draws_psis_of()` + `loo_expect_value()` | 69.4, 55.6, 60.6 | 55.6 | 2681 |
| brms: `relative_eff()` + `loo(save_psis = TRUE)` + `E_loo_value()` | 67.3, 67.4, 68.0 | 67.3 | 2681 |
| control (frmtmb again) | 54.9, 63.5, 59.7 | 54.9 | 2681 |

The two results are identical. brms's path runs the whole `loo()`
(pointwise elpd and its diagnostics) to get the psis object; frmtmb's
runs `loo::psis()` alone. The peak is the same, 2.7 GB, set by
`E_loo()`'s per-column work and the weights matrix.

Arm B, end to end, each package on its own fit: frmtmb `loo_predict(ds)`
76.8 and 106.9 s, peak 3188 and 2657 MB; brms `loo_predict(fit)` 112.2
and 111.4 s, peak 2242 and 3117 MB. Min: frmtmb 76.8 s, brms 111.4 s.
One `log_lik()` per call on both sides (brms's runs inside `loo()`).

## Before arm (evidence item 3)

`test-loo-predict.R` on rellib-r10 (log
`dev/loodraws-log/tests-base-loo-predict.txt`): pass=4 err=11 of 13
blocks. Nine blocks fail with "could not find function" (the weak
form). The behavioral failures: `add_criterion(ds, "loo_R2")` refused
with "add_criterion(criterion = "loo_R2") is not available for frmtmb
draws: loo_R2() is not implemented here"; and the mock block fails on
the missing internal. The two blocks that pass on the base are the
invariants this lane builds on (`psis()` == `loo()`'s weights; pp_check
on all draws), kept as guards. The brms-attached before arm is the
dispatch table above: "no applicable method" on rellib-r10.

## Bin-2 rerun, brmsfit-methods block (evidence item 4)

`FRMTMB_PORT_LIB=<lane lib> Rscript dev/brmsport-bin2-run.R
frmtmb.sample dev/brmsport-bin2/frmtmb.sample/test-brms-suite-bin2-methods.R
dev/loodraws-bin2-log/rec-frmtmb.sample-methods.tsv`: RESULT pass=58
fail=0 err=0 skip=0 blocks=13 (record mode). The 17 rows, all
`absent` in the ledger:

| row | held | new class |
|---|---|---|
| 487, 490 (loo_R2) | TRUE | pass |
| 497, 511, 512 (loo_epred) | TRUE | pass |
| 508 (loo_epred newdata) | FALSE, stale on the refusal | divergence: policy (newdata) |
| 519, 522 (loo_linpred) | TRUE | pass |
| 523 (loo_linpred(fit4) dim 3) | FALSE: dim(llp)[3] is NA | divergence: inherited (see defects) |
| 526 (loo_linpred scale) | FALSE, stale on the refusal | divergence: policy (`scale` named) |
| 533, 547, 548 (loo_predict) | TRUE | pass |
| 544 (loo_predict newdata) | FALSE, stale on the refusal | divergence: policy (newdata) |
| 555 (loo_predictive_interval) | TRUE | pass |
| 562, 563 (loo_model_weights) | TRUE | pass |

13 of 17 flip to pass. The verdict files are not edited.

## Defects found, not fixed

1. `posterior_epred()`, `posterior_linpred()` and `posterior_predict()`
   on multivariate draws with `resp = NULL` return the first response's
   S x N matrix; brms returns S x N x nresp (measured: frmtmb 4000 x 60,
   brms 4000 x 60 x 2, `parity-compare.txt`). The `loo_*()` functions
   build brms's layout themselves.
2. `posterior_linpred()` on draws of an ordinal model with a `cs()`
   term returns S x N; brms returns one layer per threshold (row 523:
   brms's `dim(llp)[3]` is 3; frmtmb's is NA). `loo_linpred()`
   inherits it.
3. `posterior_linpred.frmtmb_draws()` never calls `frm_check_dots()`,
   so an unknown argument is ignored in silence
   (`posterior_linpred(ds, junk = 1)` returns the matrix,
   `dev/loodraws-probe-shapes.R`, log `probe-shapes.txt`), while
   `posterior_epred()` refuses it. `loo_linpred()` checks its own
   `...`, so it is not affected.
4. pp_check() vs loo() r_eff on thinned draws, above.

## Not done, and why

- Core generics for the six (decision 6).
- Moment matching and reloo: draws carry no Stan program to refit; the
  refusal is `loo()`'s.
- Unifying the thinned-draws r_eff convention (needs a decision).

## Runs

All with `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
FRMTMB_RELLIB=<lane lib>`: core and frmtmb.sample from the lane
library, the other six extensions from rellib-r10 over the lane's core
(every `lib:` line checked: 253 frmtmb and 62 frmtmb.sample files from
`wt-loodraws-lib`, every extension's core from `wt-loodraws-lib`).

One file per process, all eight suites, `dev/release/run-par.sh
loodraws-all dev/loodraws-jobs-all.txt 6` (log
`dev/release/loodraws-all.log`), 411 of 411 ran:

```
frmtmb          files=253 pass=20654 fail=0 err=0 skip=14 warn=0
frmtmb.coupling files=11  pass=551   fail=0 err=0 skip=5  warn=0
frmtmb.eam      files=32  pass=1844  fail=0 err=0 skip=3  warn=0
frmtmb.latent   files=11  pass=368   fail=0 err=0 skip=2  warn=0
frmtmb.learn    files=15  pass=501   fail=0 err=0 skip=2  warn=0
frmtmb.ode      files=11  pass=549   fail=0 err=0 skip=1  warn=0
frmtmb.sample   files=62  pass=3168  fail=3 err=0 skip=1  warn=0
frmtmb.spline   files=16  pass=594   fail=0 err=0 skip=1  warn=0
```

The 3 frmtmb.sample failures were mine and are fixed:
`test-generic-collision.R` counted the six new `frmtmb_fit` rows
(36 vs 30; the guard now counts `frmtmb_draws` rows), and
`test-message-uniqueness.R` found three repeated templates (the loo
requirement, "Models have different number of observations.",
"Invalid argument 'resp'"; now one helper each, `loo_need_pkg()` and
`brms_resp_index()`, and a message of its own for the row-count
mismatch). After the fix, frmtmb.sample and frmtmb.eam again
(`loodraws-sample`, 94 of 94): frmtmb.sample files=62 pass=3170 fail=0
err=0 skip=1 warn=0; frmtmb.eam files=32 pass=1844 fail=0 err=0 skip=3
warn=0. Every skip is a tier gate (`FRMTMB_DRMTMB_FIT_TESTS`,
`FRMTMB_FUZZ`, `FRMTMB_SCALE_TESTS`). The core and other extension
suites were not rerun after the fix: it touched frmtmb.sample's R/
only (two helpers, one message), and the 253 core files passed on the
build before it.

One session per suite, three orders, `dev/release/run-sessions.sh
loodraws-sess 6` (log `dev/release/loodraws-sess.log`), on the fixed
build, 24 of 24, fail=0 err=0 in every one:

```
frmtmb           pass forward 20654, reverse 20752, preload 20848
frmtmb.coupling  551 / 551 / 551
frmtmb.eam       1844 / 1844 / 1844
frmtmb.latent    368 / 368 / 368
frmtmb.learn     501 / 501 / 501
frmtmb.ode       549 / 549 / 549
frmtmb.sample    3170 / 3170 / 3170
frmtmb.spline    594 / 594 / 594
```

Untagged LEAK lines: rstan's Rcpp module classes in the global
environment (`.__C__Rcpp_stan_fit4model<hash>_`, `.__C__C++Object`,
`.__global__`, `.requireCachedGenerics`) after files that compile or
load a Stan model, in core (bcm-*, brms-*, emmeans, gp-by and others)
and in frmtmb.sample (test-default-priors-fit.R, test-loo.R,
test-sampling-ported.R); and `+package:emmeans` on the search path
after core's test-brms-suite-emmeans.R. None is a file of this lane.
Base arm for frmtmb.sample (main checkout's tests, rellib-r10, forward;
log `dev/release/loodraws-sessbase.log`): files=61/61 pass=3080 fail=0
err=0 skip=1, with the same three untagged LEAK lines. So they are
pre-existing, which contradicts lane-env's "clean in all three orders"
on 4854e425; not verified on base for core.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
`OPENBLAS_NUM_THREADS=4`, logs `dev/loodraws-log/openblas-*.txt`):
test-loo-predict.R pass=86, test-loo.R pass=142, test-postfit-draws.R
pass=68, test-ppcheck-loo.R pass=32, each fail=0 err=0 skip=0 warn=0.
The R copy is deleted.

R CMD check --as-cran of frmtmb.sample (`dev/loodraws-check.ps1`, log
`dev/loodraws-check/check.log`): first run `Status: OK`, but its
tarball held 25 stray `tests/testthat/_problems/` files from the
before-arm run; deleted, and checked again: `Status: OK`, tarball with
0 `_problems` entries.

## Fix round 1 (review dev/reviews/2026-10-20-loodraws.md)

Coordinator decisions: (a) brms's r_eff convention everywhere; (b) the
generics and their fit refusals move to core, `kfold` and `reloo` too.

### What changed, by finding

- Decision (a), finding 2. `loo_r_eff()` (frmtmb.sample `R/loo.R`) is
  brms's `r_eff_log_lik()`: `relative_eff()` on the chains for every
  draw, on ONE chain for a subset, 1 for combined imputations, all 1
  under brms's warning when one is NA. `pp_check()`'s
  `draws_ppc_psis_fit()` (`R/methods-draws.R`) now calls it, so
  `loo()`, `psis()`, the six and `pp_check()` share one function.
  `waic()` takes no r_eff (brms's neither). Size on `loo(ndraws =)`
  (`dev/loodraws-reff-change.R`, data seed 1, 60 rows, 4 x 1000 draws,
  seed 1; log `dev/loodraws-log/reff-change.txt`), elpd_loo new minus
  old:

  ```
  y ~ x           ndraws  500: 0.0002     (se 5.3)
  y ~ x           ndraws 1000: -0.000265  (se 5.22)
  y ~ x           ndraws 3000: -0.000178  (se 5.26)
  y ~ x + (1 | g) ndraws  500: 0.00645    (se 5.03)
  y ~ x + (1 | g) ndraws 1000: 0.00025    (se 4.93)
  y ~ x + (1 | g) ndraws 3000: 0.000815   (se 4.91)
  ```

  one-chain r_eff ranged 0.570 to 1.311. In each case the new estimate
  is identical to `loo::loo(ll, r_eff = <one-chain relative_eff>)`.
  On brms's own draws and the same `draw_ids` (the reviewer's
  construction, `dev/loodraws-same.R`, log `dev/loodraws-log/same.txt`):
  `loo()` elpd brms -127.3531032 frmtmb -127.3531032 (gauss) and
  -123.8635194 both (ri), diff 0; `loo_epred(draw_ids =)` rel diff 0;
  `pp_check()` lw frmtmb vs brms max abs diff 0, and frmtmb's
  `pp_check()` vs its `loo()` 0. BREAKING in frmtmb.sample's NEWS.
- Decision (b). Core `R/loo-generics.R` (new): the generics
  `loo_epred`, `loo_linpred`, `loo_predict`, `loo_predictive_interval`,
  `loo_R2`, `loo_model_weights`, `kfold`, `reloo`, and their
  `frmtmb_fit` methods. The six refuse naming
  `<fn>(frmtmb.sample::frm_sample(fit))`; `kfold()` and `reloo()` say
  that neither package refits the sampler and name `frm_bootstrap()`
  and `AIC()` (naming `frm_sample()` would send the user to another
  refusal). Core `R/generic-owners.R`: eight owner entries
  (`loo_epred = c("brms", "rstantools")`, the two-owner rule core
  already has for `ngrps`); `NAMESPACE`; `man/loo_predict.Rd`; NEWS.
  frmtmb.sample: the generics and fit methods removed
  (`R/loo-predict.R`, `R/loo.R`), the eight names out of
  `sample_generic_owners` (now 22), re-exported from frmtmb
  (`R/reexports.R`). The kfold/reloo gap measured first, on the base
  build: with frmtmb alone "could not find function", with loo or brms
  attached "no applicable method ... frmtmb_fit" (the reviewer's
  dispatch logs and core test block below, base arm). Not moved, the
  same gap: `loo_moment_match()` and `loo_subsample()` (outside this
  instruction; flagged).
- Finding 1. `loo_default_value()` (`R/loo-predict.R`): `NULL`,
  `allow_new_levels`/`allow.new.levels = FALSE`, `sample_new_levels =
  "uncertainty"` count as absent in `loo_draws_opts()`,
  `draws_model_weights()` and `loo()`/`waic()` (`loo_criterion()`).
  `draw_ids` reaches `log_lik()` in `loo()`, `waic()` and
  `loo_model_weights()` as in the other four.
- Finding 3. `loo_resp_name()` takes `ya` (brms) or `y_a` (frmtmb) in
  `loo_epred()`, `loo_linpred()`, `loo_predict()`,
  `loo_predictive_interval()`, `loo_model_weights()` and `psis()`;
  `loo_R2()` already took brms's `ya` and keeps it alone, as
  `bayes_R2()` does: a brms pin (`test-brms-pins.R`, "multivariate
  VarCorr() and bayes_R2() have a row per response") asserts that
  `bayes_R2(ds, resp = "y_a")` errors as brms errors. A first spelling
  that let both functions take `y_a` too broke that pin (fail=1 in
  the full run below) and was reverted. The mv third dimension is
  named in brms's spelling. Untouched: `log_lik()`, `posterior_*()`
  and `loo()` keep frmtmb's `y_a` only, as before.
- Finding 4. `add_criterion_args()` (`R/loo.R`): each criterion gets
  the names it takes; a name no requested criterion takes is refused
  by name ("no requested criterion (loo) takes the argument `seed`").
- Finding 5. `psis.frmtmb_draws()` sorts `...` with `loo_draws_opts()`:
  `draw_ids` honored, `r_eff`/`cores` to `loo::psis()`, unknown names
  refused.
- Finding 6. A gated block in `test-loo-predict.R` (two brms fits,
  `chains = 2, iter = 1000, seed = 1`, data seed 3, 80 rows): brms's
  draws in a `frmtmb_draws` with brms's stanfit; `expect_identical()`
  for log_lik (attributes aside), `loo_epred()`/`loo_linpred()` mean,
  var, quantile, `loo_epred(draw_ids =)`, stacking and pseudobma weights
  and classes; `loo_R2()` within `N * eps * max * 4`; the stacking pair
  asserted non-degenerate. NEWS states the identity, with the
  reviewer's construction's numbers rerun on the lane build
  (`same.txt`): stacking 0.15475477 / 0.84524523 on both, every
  loo_epred/loo_linpred rel diff 0, 53 dimnames comparisons all TRUE,
  loo_R2 draws largest rel diff 1.08e-14.
- Finding 7. `test-loo-predict.R` "posterior_average() weights do not
  read a stored loo"; sample `test-generic-collision.R` "the loo_*()
  names on draws reach their methods in every order" and "the draws
  methods on frmtmb's loo_*() names have every twin"; core
  `test-generic-collision.R` "the loo_*() names refuse on a fit in
  every attach order" and the eight names in `brms_shared`.
- Finding 9. NEWS and `sample-loo_predict` wording: one log_lik per
  call "per response for loo_R2(), per model for loo_model_weights()";
  "subset" for `ndraws`/`draw_ids`, with "`frm_sample(thin =)` draws
  keep their chains"; the parity sentence replaced by the identity;
  the resp spelling and r_eff in NEWS.
- Finding 8: no change (a tier difference, as the review says).

### Before arm of the new tests (rellib-r10)

`dev/release/loodraws-fix1base.log`: sample test-loo-predict pass=6
fail=3 err=16; sample test-generic-collision pass=76 fail=9 err=1;
core test-generic-collision pass=79 fail=9. Behavioral failures among
them: `psis(d, draw_ids = ids)` dim 120 x 40 where 50 x 40 is asked;
`loo(draw_ids =)` refused; `posterior_average()` weights 0.999999908 /
9.2e-08 with a stored loo against 0.999999998 without;
`add_criterion("loo_R2")` refused; the dispatch block "no applicable
method for 'loo_epred' applied to an object of class "frmtmb_draws""
through `rstantools::loo_epred`; the core block `REF_loo_epred: FALSE`
and the other seven (could not find function / no applicable method).

### Runs

All on the lane library (core and frmtmb.sample from it; the six other
extensions from rellib-r10 over the lane's core; `lib:` lines checked).

One file per process, all eight suites, gated
(`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`; `dev/release/run-par.sh
loodraws-fix1all`, 411 of 411; totals `dev/loodraws-log/fix1all-totals.txt`):

```
frmtmb          files=253 pass=20668 fail=1 err=0 skip=14 warn=0
frmtmb.coupling files=11  pass=551   fail=0 err=0 skip=5  warn=0
frmtmb.eam      files=32  pass=1844  fail=0 err=0 skip=3  warn=0
frmtmb.latent   files=11  pass=368   fail=0 err=0 skip=2  warn=0
frmtmb.learn    files=15  pass=501   fail=0 err=0 skip=2  warn=0
frmtmb.ode      files=11  pass=549   fail=0 err=0 skip=1  warn=0
frmtmb.sample   files=62  pass=3234  fail=1 err=0 skip=1  warn=0
frmtmb.spline   files=16  pass=594   fail=0 err=0 skip=1  warn=0
```

The two failures, both fixed after this run: core
`test-arg-refusal.R` "no registered S3 method has a `...` it never
touches" named the eight `frmtmb_fit` refusals (their body was a call
to a helper the guard did not know; `loo_fit_refusal()` and
`loo_refit_refusal()` are now single `frm_stop()` statements and join
`ar_refusers`, and the "refusal helpers really are unconditional" block
checks them); and sample `test-brms-pins.R` (the `bayes_R2(resp =
"y_a")` pin, finding 3 above). Rerun after the fix
(`loodraws-fix1b`, 99 of 99, totals `dev/loodraws-log/fix1b-totals.txt`):
frmtmb.sample files=62 pass=3235 fail=0 err=0 skip=1 warn=0; frmtmb.eam
files=32 pass=1844 fail=0 err=0 skip=3 warn=0; core test-arg-refusal 119,
test-generic-collision 88, test-message-uniqueness 6, test-brms-postfit
78, test-api-spellings 36, each fail=0 err=0 skip=0 warn=0. The other
248 core files and five extensions were not rerun one per process after
that fix (it changed one helper's body and two test files); the
one-session runs below cover them on the final build.

The gated same-draws block was then made robust to the BLAS, after an
OpenBLAS 0.3.26 LAPACK run failed 18 of its `identical()` assertions.
`dev/loodraws-blas-same.R` (logs `dev/loodraws-log/blas-same-ref.txt`,
`blas-same-openblas.txt`), largest difference in eps of the largest
value, reference BLAS / OpenBLAS: log_lik 0 / 0.94 (gauss), 0 / 2.06
(ri); loo_epred mean 0 / 3.11, 0 / 7.93; var 0 / 38.3, 0 / 233; quantile
0 / 0.70, 0 / 4.84; draw_ids 0 / 3.54, 0 / 5.19; loo_R2 draws 49.4 / 43.9,
31.8 / 30.5; stacking weights 0 / 7958 (0.1364898406 / 0.8635101594 on
every side, to 10 digits). brms forms the predictor with `X %*% b`, frmtmb
row by row, so the identity is a property of reference BLAS. The test
now bounds each difference by `N * eps` (N = 80 rows) of the largest
value, of the largest squared mean for a variance, of four times the
largest value for loo_R2, and by `sqrt(eps)` (optim()'s default
reltol) for the weights. With that: test-loo-predict.R pass=155 fail=0
err=0 skip=0 warn=0 under reference BLAS
(`dev/loodraws-log/fix1-ref-loo-predict.txt`) and under OpenBLAS 0.3.26
LAPACK (`openblas-fix1-test-loo-predict.txt`); test-loo.R 142,
test-postfit-draws.R 68, test-ppcheck-loo.R 32, each clean under
OpenBLAS. The OpenBLAS R copy is deleted.

One session per suite, three orders, final build:

```
gated (FRMTMB_BRMS_FIT_TESTS=true), dev/release/loodraws-fix1sessg.log, 24 of 24
frmtmb           20675 / 20773 / 20869 pass (forward/reverse/preload)
frmtmb.coupling  551 x3; eam 1844 x3; latent 368 x3; learn 501 x3
frmtmb.ode       549 x3; sample 3251 x3; spline 594 x3
ungated, dev/release/loodraws-fix1sessu.log, 24 of 24
frmtmb           18273 / 18371 / 18467 pass, files 243/253, skip 191/191/192
frmtmb.coupling  551 x3; eam 1838 x3; latent 368 x3; learn 430 x3
frmtmb.ode       549 x3; sample 2979 x3 (57/62 files, skip 9); spline 594 x3
```

fail=0 err=0 warn=0 in all 48. Ungated: 0 untagged LEAK lines. Gated:
the untagged lines are rstan's Rcpp module classes in the global
environment and `+package:emmeans` (core's test-brms-suite-emmeans.R),
the classes the review's finding 8 lists, now also after
test-loo-predict.R, whose gated block compiles two brms models; no
other untagged line.

R CMD check --as-cran (`dev/loodraws-check.ps1 <pkg>`, logs
`dev/loodraws-check/<pkg>/check.log`): frmtmb `Status: 1 NOTE` (the
V8 math-rendering NOTE on the HTML manual that lane-rules.md names as
expected); frmtmb.sample `Status: OK`, checked again after the last
test edit, `Status: OK`.

Bin-2 brmsfit-methods block on the fixed build (`dev/brmsport-bin2-run.R`
with the lane library; log `dev/loodraws-bin2-log/run-fix1-*`): RESULT
pass=58 fail=0 err=0 skip=0 blocks=13. The 17 rows as before: 13 hold,
508 and 544 (newdata, policy) and 526 (`scale`, policy) stale on the
refusal, 523 (`cs()` linpred layers, inherited) fails. One row outside
the 17 moved: 469, `loo(fit6, cores = 1, newdata = fit6$data)`, now
HOLDS. On a `frmtmb_draws`, `fit6$data` is NULL (there is no `data`
element), so the call passes brms's default `newdata = NULL`, which
`loo()` now accepts; it was refused by name. The row is hollow, not a
pass: consolidation should class it so. 471 moves from stale to a
plain FALSE (its ledger reason, loo 2.10.1's `loo_compare()` returning
a string, is unchanged).

## Fix round 2 (re-check in dev/reviews/2026-10-20-loodraws.md)

1. Gated same-draws tolerances (`test-loo-predict.R`, "on brms's own
   draws, every loo_*() is brms's to rounding"). Each bound is now a
   stated multiple of what the run measures between the two sides: dl,
   the largest log-likelihood difference, and dx, the largest relative
   difference of `posterior_epred()`. A leave-one-out summary sums over
   the S draws (1000; 334 for the subset), each term's weight erring by
   at most 2 dl and its value by dx, plus eps of rounding: mean bound
   S (2 dl + dx + eps) max|x|; variance, two S-term sums of squares,
   2 S (2 dl + 2 dx + eps) max(x^2) / (1 - max sum w^2); quantile,
   loo's linear interpolation in the cumulative weight, S (2 dl + eps)
   L + (dx + eps) max|x| with L the largest gap-over-weight slope
   between neighboring draws; loo_R2, the N-row bootstrap replayed at
   seed 7 to bound var_err's and var_y's errors. The model weights are
   bounded through the pointwise elpd difference dp: stacking through
   its objective (optim()'s reltol bounds the change in the objective,
   so the two optima are compared there, within reltol of its size
   plus 2 N dp), pseudo-BMA by 2 N dp plus eps per bootstrap draw.
   The comments say so. With the 0.20.0 r_eff rule mocked back in
   (`dev/loodraws-oldreff-check.R`, log
   `dev/loodraws-log/oldreff-check.txt`): RESULT pass=150 fail=3; the
   gated block fails on the subset `loo_epred(draw_ids =)` for both
   models, "0.0004985620601 > 0.0000000000002" and "0.0017014655429 >
   0.0000000000003", and the one-chain block fails on `r_eff`.
2. Core `NEWS.md`: the bullet rewrapped to 80 columns (longest line now
   72).
3. `loo_moment_match()` and `loo_subsample()` moved to core with the
   eight: generics in `R/loo-generics.R`, owner entries (`loo`) in
   `R/generic-owners.R` (43 names), `frmtmb_fit` refusals.
   `loo_moment_match()` uses the refit refusal (it moment-matches and
   then refits); `loo_subsample()` a new `loo_unbuilt_refusal()` that
   names the full `loo()` of draws and `AIC()`. Both helpers are single
   `frm_stop()` statements in `ar_refusers`. The generic signature is
   `(x, ...)`: loo 2.10.1's own `loo_moment_match` generic is
   `function(x, ...)` (printed), and the guard "each generic's formals
   are its first owner's" holds it to that; the `loo` argument belongs to
   loo's `loo_moment_match.default(x, loo, ...)`, where `loo_moment_match(fit)`
   used to land and die. A `frmtmb_fit` method in loo's table means the
   default is never reached. frmtmb.sample: the two generics removed,
   re-exported, `sample_generic_owners` now 20. Tests: core
   `brms_shared` gains both names, the fit block checks `REF_mm` and
   `REF_ss` in four attach orders (10 refusals); the sample dispatch
   block checks that `loo::loo_moment_match(ds)` and
   `loo::loo_subsample(ds)` reach the draws refusals in four orders; the
   sample twin check covers both. On rellib-r10
   (`dev/release/loodraws-fix2base.log`): core test-generic-collision
   pass=79 fail=9 with `REF_mm: FALSE` and `REF_ss: FALSE` in every
   order; sample test-generic-collision pass=82 fail=11 err=1, among
   them `REFUSE_mm`/`REFUSE_ss` and the twin check.

### Runs (fix round 2)

Gated one file per process, core and frmtmb.sample
(`dev/release/loodraws-fix2par.log`, 315 of 315; totals
`dev/loodraws-log/fix2par-totals.txt`):

```
frmtmb        files=253 pass=20678 fail=0 err=0 skip=14 warn=0
frmtmb.sample files=62  pass=3257  fail=0 err=0 skip=1  warn=0
```

One session per suite, three orders (`loodraws-fix2sessg.log`,
`loodraws-fix2sessu.log`), fail=0 err=0 warn=0 in all 12:

```
gated    frmtmb 20678 / 20776 / 20872 (skip 14/14/15); frmtmb.sample 3257 x3 (skip 1)
ungated  frmtmb 18276 / 18374 / 18470 (243/253 files, skip 191/191/192);
         frmtmb.sample 2987 x3 (57/62 files, skip 9)
```

Untagged LEAK lines: ungated 0; gated 96, every one rstan's Rcpp
module classes or `+package:emmeans` (none other).

R CMD check --as-cran: frmtmb `Status: 1 NOTE` (V8 math rendering,
expected); frmtmb.sample `Status: OK`.
