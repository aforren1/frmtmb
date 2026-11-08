# Lane loglik, round r080 (2026-11-06): logLik() without the prior

Base 3cd2a0c8 (frmtmb 0.79.1), worktree `frmtmb-wt-loglik`, branch
`wt-r80-loglik`, library `C:/Users/adf44/source/r/wt-loglik-lib` (core and
frmtmb.latent from this tree) in front of `rellib-r18`. Nothing is
committed. This section is the record of round r080. The record of the
earlier lane of the same name (0.75.0, `log_lik()` for residual
covariance terms) follows it unchanged, under "Earlier lane loglik". The
earlier lane's scripts `dev/loglik-*.R` without `r80` in the name are its
own, except the ones this lane added: `dev/loglik-base.R`,
`dev/loglik-names.R`, `dev/loglik-routes.R`, `dev/loglik-run-tests.R`,
`dev/loglik-run-par.sh`, `dev/loglik-run-session.R`,
`dev/loglik-run-sessions.sh` and every `dev/loglik-r80-*` file. Logs are gitignored; every number below is
pasted from the log named beside it.

STATUS: fix round 1 done (see "Fix round 1"); core suite and R CMD check rerun there.

## Decision: the split, with one refusal

The user's split, not the fallback. `logLik()` on an `frm()` fit is the
log-likelihood of the data at the estimates: the maximized value less the
log prior at the estimates. `logLik_penalized()` (new export) is the
maximized value, with attribute `log_prior`. `AIC()`, `BIC()`,
`extractAIC()` and so `step()` refuse a prior with a density.

One route cannot split exactly: `REML = TRUE` with a prior that reads the
fixed effects (class `b`, `Intercept`, an `nlpar`, or a centered
Intercept's offset). REML integrates `beta` with the Laplace
approximation, so that prior sits inside the integral. There `logLik()`
stops (class `frmtmb_loglik_inner_prior`) and `logLik_penalized()` still
answers, with `log_prior = NA`. The rule in code (`fit_log_prior()`,
`R/methods-fit.R`): resolve the fit's prior entries again, collect every
component an entry or its centering offset reads, and refuse if any is an
integrated component (`b`, `miss`, and `beta` under REML). A profiled
`beta` (`frmtmb_control(profile = TRUE)`) is not integrated: TMB's
`profile` maximizes it in the inner problem with no Laplace term
(`TMB::MakeADFun` source, `tmb_sparse_izamd` zeroes the profiled block
of the Hessian before the log-determinant, then adds back
`sum(profile)/2 * log(2 pi)`). No route of `frm()` puts a prior on `b` or
`miss` today; the rule covers it if one appears.

### Brief's table, rebuilt (`dev/loglik-base.R`, seed 4101)

My construction (`n = 200`, 20 groups, `y = 1 + 0.5 x + u_g + e`) is not
the brief's (its logLik for no prior was -297.44; no script with seed
4101 and that number exists in `dev/`), so the levels differ. The pattern
reproduces: `normal(0, 100)` lowers logLik by 5.5241, and
`log(100 * sqrt(2 * pi)) = 5.5242`.

Before (rellib-r18, `dev/loglik-base-before.log`):

    none             logLik -305.7011162851  AIC 619.4022326  b_x 0.526750
    normal(0, 100)   logLik -311.2252388774  AIC 630.4504778  b_x 0.526750
    normal(0, 10)    logLik -308.9240271702  AIC 625.8480543  b_x 0.526724
    normal(0, 1)     logLik -306.7581141894  AIC 621.5162284  b_x 0.524193
    normal(0, 0.1)   logLik -313.5847876904  AIC 635.1695754  b_x 0.349713
    normal(0, 100)   data logLik at MAP via no-prior tape -305.7011162850
    normal(0, 10)    data logLik at MAP via no-prior tape -305.7011163555
    normal(0, 1)     data logLik at MAP via no-prior tape -305.7017863499
    normal(0, 0.1)   data logLik at MAP via no-prior tape -308.8534760685

After (lane build, `dev/loglik-base-after.log`; AIC refusal text cut):

    none             logLik -305.7011162851  AIC 619.4022326  pen -305.7011163
    normal(0, 100)   logLik -305.7011162850  AIC ERROR: AIC() counts each ...  pen -311.2252389
    normal(0, 10)    logLik -305.7011163555  AIC ERROR: ...  pen -308.9240272
    normal(0, 1)     logLik -305.7017863499  AIC ERROR: ...  pen -306.7581142
    normal(0, 0.1)   logLik -308.8534760685  AIC ERROR: ...  pen -313.5847877

Each new `logLik()` equals the no-prior tape at the MAP estimates to the
10 printed decimals, and `logLik_penalized()` is the old `logLik()`.

glmmTMB 1.1.15.2 on the same data with `normal(0, 100)` on `x`
(`dev/loglik-names.R`, `dev/loglik-names.log`): `logLik` -305.7011162849
with no prior and -311.2252388773 with it, and `AIC(m1)` 630.4505 with no
word. So glmmTMB's `logLik()` is the penalized value, frmtmb's old one.

## First question: every route that adds a prior

`dev/loglik-routes.R`, seeds 4101 to 4106. For each route, a MAP fit `fp`
and the same model with no prior `f0` share one parameter layout, so
`-f0$obj$fn(fp$opt$par)` is the data-only objective of the same route at
the MAP estimates. `split = -fp$opt$objective + nlp(fp$estimates)`, with
`nlp` the fit's own `neg_log_prior_fn()` over its resolved entries. Two
routes need their own data value:

- profile: `f0`'s profiled tape maximizes `beta` again (row "tape
  re-max", 1.55e-03 off, which is the wrong comparison), so the data
  value is the plain Laplace tape at `(beta, theta)` of the MAP fit;
- importance: the fit's own corrected-objective closure is found on its
  tape and evaluated with `neg_log_prior = NULL`, at the same draws and
  offsets.

Lane build (`dev/loglik-routes-after.log`; the before log has the same
split and data columns, `dev/loglik-routes-before.log`):

    A laplace gaussian       pen -305.4110095949 log_prior -5.2257905702 split -300.1852190247 data -300.1852190247 rel 0.00e+00 comps [beta,theta]
    A laplace poisson        pen -312.8105778688 log_prior -6.7961136010 split -306.0144642678 data -306.0144642678 rel 0.00e+00 comps [beta,theta]
    B profile (tape re-max)  pen -305.4110095949 log_prior -5.2257905676 split -300.1852190273 data -299.7207296004 rel 1.55e-03 comps [beta,theta]
    B profile                pen -305.4110095949 log_prior -5.2257905676 split -300.1852190273 data -300.1852190273 rel 0.00e+00 comps [beta,theta]
    C quadrature             pen -139.5810059831 log_prior -8.0746023977 split -131.5064035854 data -131.5064035812 rel 3.15e-11 comps [beta,theta]
    D importance             pen -312.7310342313 log_prior -6.7869910242 split -305.9440432072 data -305.9440432072 rel 3.72e-16 comps [beta,theta]
    E REML, sd prior         pen -304.6727739640 log_prior -2.1499311665 split -302.5228427975 data -302.5228427975 rel 0.00e+00 comps [theta]
    E REML, b prior          pen -308.1547693900 log_prior -5.1993163607 split -302.9554530293 data -302.3746957269 rel 1.92e-03 comps [beta,theta]
    F exact gp               pen -17.0086185014  log_prior -2.0048234697 split -15.0037950318  data -15.0037950318  rel 0.00e+00 comps [theta]
    G latent residual        pen -139.6909608505 log_prior -2.4636062049 split -137.2273546456 data -137.2273546455 rel 4.09e-13 comps [theta]
    H me()                   pen -239.6356354382 log_prior -4.2507389064 split -235.3848965318 data -235.3848965318 rel 1.21e-16 comps [beta,logsdme]
    I mo() simo              pen -87.9070014583  log_prior -0.4658718298 split -87.4411296285  data -87.4411296285  rel 0.00e+00 comps [zeta1]
    J nonlinear              pen -5.1818325534   log_prior -1.2011591578 split -3.9806733957   data -3.9806733957   rel 2.12e-15 comps [beta]
    max rel over separable routes: 3.15e-11

`logLik()` on the lane build returns `split` on every separable row (the
log's last column) and refuses row "E REML, b prior". Route I carries the
simplex chart's Jacobian inside the Dirichlet prior; it splits off with
the prior. Route C's 3.15e-11 is the Gauss-Kronrod rule's own rounding on
two tapes (4.2e-9 absolute); every other route is 4.1e-13 or less. That
the split is the data value is an identity wherever the prior reads only
outer parameters (the inner problem and its log-determinant do not see
it); these rows are the numerical check of it on each route.

## What changed

Core:

- `R/methods-fit.R`: `logLik.frmtmb_fit` subtracts the log prior and
  refuses where the prior is inside the integral; new
  `logLik_penalized()` (generic and `frmtmb_fit` method),
  `AIC.frmtmb_fit`, `BIC.frmtmb_fit` (refuse any fit in the call with a
  density prior, then `NextMethod()`), helpers `fit_integrated_comps()`,
  `fit_log_prior()`, `loglik_inner_prior_refusal()`,
  `ic_prior_refusal()`; one Rd, `man/logLik_penalized.Rd`, that documents
  `logLik()`, `logLik_penalized()`, `AIC()`, `BIC()`. `summary()` prints
  `logLik` and `log prior` and no AIC or BIC under a density prior, and
  marks `(MAP)` only for a density (`fit_prior_density()`), not for
  bounds alone.
- `R/sugar.R`: `extractAIC.frmtmb_fit` refuses (so `step()` does).
- `R/fit-end.R`: the discrete check runs under a prior, on the data value.
- `R/allfit.R`, `R/par-template.R`: compare `logLik_penalized()`.
- `R/multiple.R`: `anova.frmtmb_multiple` refuses D3 and likelihood D2
  under a density prior.
- `R/refit-data.R`, `R/confint.R`: refusal text and `?anova.frmtmb_fit`
  no longer say logLik holds the prior.
- `R/fit.R`: `@param prior` and one comment.
- `R/compat.R`: notes of two rows (statuses unchanged):
  `REML x prior` (conditional) and `importance x prior` (works).
- `R/frmtmb-package.R`: `AIC` and `BIC` added to `@importFrom stats`
  (see "A defect R CMD check found").
- `NEWS.md`: two BREAKING bullets and four others under
  `# frmtmb (development version)`.

frmtmb.latent:

- `R/hmm-starts.R`: the four `logLik()` reads become
  `frmtmb::logLik_penalized()`; `@return` says so. `man/hmm_starts.Rd`
  regenerated. `NEWS.md` bullet.

Tests:

- new `tests/testthat/test-loglik-prior.R` (9 blocks; the 9th,
  `frm_allfit()`, was added after the suite runs: 49 passes alone);
- `tests/testthat/test-refit-data.R`: `step()` now stops at
  `extractAIC()`, so the expected message is the AIC refusal's;
- `tests/testthat/test-beta-binomial-density.R`: one `expect_silent()`
  pinned "skipped under a prior", which the brief reverses; removed (it
  also passed a plain data frame as `prior`, which the resolver refuses),
  with a pointer to the new test;
- `tests/testthat/test-map.R`, `test-nl-rtmb-scope.R`,
  `test-prior-compat.R`, `test-brms-priors.R` (gated): each pinned
  "logLik() holds the prior" (see "Tests" below);
- `extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R`: one block.

Name of the accessor. brms 2.23.0 exports `log_posterior()` (draws of
`lp__`) and `log_lik()` (pointwise data part); glmmTMB has no separate
accessor (its internal `prior_nll`; `logLik()` is penalized); lme4 has no
priors and names the penalized criterion `REMLcrit()` / `devcomp()`'s
`pwrss` ("penalized weighted residual sum of squares";
`dev/loglik-names.log`). `log_posterior` is taken in this ecosystem:
frmtmb.sample exports it for draws, with brms's meaning. So
`logLik_penalized()`: the penalized log-likelihood the fit maximized,
next to `logLik()` in the index. It equals `logLik()` with no prior.

## Tests: seen to fail on rellib-r18, pass on the lane build

`dev/loglik-run-tests.R` (`dev/release/run-tests.R` with `LIB` set to
rellib-r18), one file per process; "before" is rellib-r18 alone, "after"
puts the lane library first. From `dev/loglik-fixbase.log` and
`dev/loglik-fix.log`:

    before  RESULT test-map.R pass=3 fail=0 err=1 skip=0 warn=0
    before  RESULT test-nl-rtmb-scope.R pass=102 fail=0 err=1 skip=0 warn=0
    before  RESULT test-prior-compat.R pass=198 fail=0 err=1 skip=0 warn=0
    before  RESULT test-loglik-prior.R pass=6 fail=13 err=4 skip=0 warn=0
    after   RESULT test-map.R pass=8 fail=0 err=0 skip=0 warn=0
    after   RESULT test-nl-rtmb-scope.R pass=105 fail=0 err=0 skip=0 warn=0
    after   RESULT test-prior-compat.R pass=202 fail=0 err=0 skip=0 warn=0
    after   RESULT test-loglik-prior.R pass=42 fail=0 err=0 skip=0 warn=0

Every block of `test-loglik-prior.R` fails on the base build by
behavior, not only by the missing symbol (`dev/loglik-fixbase-out/`):
the closed-form likelihood check is off by 0.0566 relative; the Laplace
split by 0.0176, quadrature by 0.0401, REML with an sd prior by 0.0162;
`logLik()` does not refuse under REML with a b prior; `AIC()`, `BIC()`,
`extractAIC()` do not refuse; `summary()` has AIC and BIC; the discrete
check gives no warning under a prior; `anova()` of `frm_multiple()` fits
with a prior answers. The three edited files fail on the base build by
the missing symbol only: their new `logLik()` lines are the same claim as
`test-loglik-prior.R`'s first block.

The latent test, three arms (`dev/loglik-r80-latent-arms.R`, logs
`dev/loglik-t-out/latent-arm-*.txt`): base (both packages from
rellib-r18) errors on the missing export; mixed (lane core, rellib-r18
frmtmb.latent) fails three expectations, the table's first `logLik`
-245.7 against the penalized -248.3; lane passes 4 of 4. So the
frmtmb.latent change is needed by the core change, and the floor follows.

Tests changed to the new meaning, each seen to fail on the lane build
before the edit (`dev/loglik-core.log`, `dev/loglik-gated.log`):

    RESULT test-map.R pass=6 fail=1 ...   -logLik(map) vs raw + nlp: 2.760495 off
    RESULT test-nl-rtmb-scope.R pass=104 fail=1 ...   tape 1476.3 vs -logLik 1473.6
    RESULT test-prior-compat.R pass=200 fail=1 ...   12.673717 off
    RESULT test-brms-priors.R pass=101 fail=1 ... (gated)   |logLik(fit) - logLik(no prior)| 0.04, wanted > 1

Each pinned "logLik() holds the prior"; each now reads
`logLik_penalized()` for that identity, and `test-map.R` and
`test-prior-compat.R` also assert that `logLik()` is the raw tape at the
MAP. Two of them had an absolute tolerance on the line I rewrote
(`1e-6`); it is now relative (`1e-8` of the value). Gated rerun:
`RESULT test-brms-priors.R pass=102 fail=0 err=0 skip=0 warn=0`
(`dev/loglik-t-out/gated-brms-priors.txt`).

## Audit of every reader

Each reader of `opt$objective` in `R/` (grep `\$objective` plus
`[["objective"]]`, 49 lines at 3cd2a0c8) and of `logLik()` in core and
the eight extensions, and the value it wants.

### `opt$objective`: optimizer internals, want the penalized value, unchanged

These compare or report the objective the optimizer minimizes. Under a
prior that objective holds the log prior, and that is what they must
compare. None changed.

| Where | What it does |
|---|---|
| `R/allfit.R:39` | bobyqa wrapper returns `fval` as `objective` |
| `R/autoscale.R:341`, `:460` | pre-fit verdict; choose between the scaled and plain fits |
| `R/fit.R:1224`, `:1745`, `:1810` | verbose output of the objective and a gain |
| `R/fit.R:1954`, `:2037`, `:2694`, `:2700` | best of several attempts; optimizer wrapper |
| `R/fit.R:3034` to `:3048`, `:3303` to `:3316` | Newton polish gain; recovery starts |
| `R/fit.R:3361` to `:3399` | the larger-sd escape |
| `R/fit.R:3556` to `:3771`, `:3859` to `:4024` | Newton polish; `mo_search()` |
| `R/fit.R:5160` | `gp_ok`: objective against the gp warm start's marginal |
| `R/gp-start.R:343` to `:398` | gp warm-start scan |
| `R/importance.R:1541`, `:1547` | corrected objective at the reported anchor |
| `R/multiple.R:541` (`dev_at_pooled()`) | D3 deviance at the pooled vector; refused under a prior now (below) |

### Readers that changed

| Where | Before | Wants | Now |
|---|---|---|---|
| `logLik.frmtmb_fit` (`R/methods-fit.R`) | `-objective` | data | `-objective - log prior`; refuses where the prior is inside the integral |
| `logLik_penalized()` (new) | none | penalized | `-objective`, attribute `log_prior` |
| `AIC`, `BIC` (stats defaults through `logLik()`) | penalized, silent | refuse | `AIC.frmtmb_fit`, `BIC.frmtmb_fit` refuse any fit in the call with a density prior |
| `extractAIC.frmtmb_fit` (`R/sugar.R`), so `step()` | penalized | refuse | refuses |
| `summary.frmtmb_fit` | `logLik`, `AIC`, `BIC`; `(MAP)` for any prior | data | `logLik`, `log prior`, AIC/BIC `NA` and not printed; `(MAP)` only for a density |
| `discrete_loglik_warning()` (`R/fit-end.R`) | skipped under any prior | data | runs on the data value; skipped under REML as before |
| `frm_allfit()` (`R/allfit.R`) | `logLik()` | penalized (each optimizer maximizes it) | `logLik_penalized()`, with a line saying so |
| `nl_restart_round()`, `nl_fit_off_saddle()` (`R/par-template.R:767`, `:840`) | `logLik()` | penalized (restarts judged on what they maximize) | `logLik_penalized()` |
| `anova.frmtmb_multiple` D2 (likelihood), D3 (`R/multiple.R:506`) | `2 * objective`, silent under a prior | refuse | refuses, as `anova.frmtmb_fit` does |
| `hmm_starts()` (frmtmb.latent `R/hmm-starts.R:399`, `:431`, `:442`, `:537`) | `logLik()` | penalized (modes ranked on what each refit maximized; `grad_rel` is the objective's gradient) | `frmtmb::logLik_penalized()` |

Without the last change, the lane's core with rellib-r18's frmtmb.latent
ranks the refits on the data value: the test fails three ways
(`dev/loglik-r80-latent-arms.R mixed`, below).

### Readers that need no change

| Where | Why |
|---|---|
| `anova.frmtmb_fit`, `drop1()` (`R/confint.R`) | already refuse a density prior (`fit_prior_density()`); message reworded |
| `deviance.frmtmb_fit` (`R/sugar.R`) | `-2 * logLik()`: now the data deviance, as lme4's |
| `nl_local_message()` (`R/par-template.R:944`) | reached only with no prior: under any prior `use_curvature = FALSE`, so `look$local` is empty and the message never runs |
| `nl_saddle_warning()`, restart messages | print the penalized value they compared, labeled "log-likelihood" as before; not reworded (not done, below) |
| `pooled_order()`, `pooled_lrt_parts()` df reads (`R/multiple.R:382`, `:504`) | read `attr(logLik(), "df")` only |
| `loo()`, `waic()`, `LOO()`, `log_lik()` messages (`R/loo.R`) | point to `AIC()`, which now explains itself under a prior |
| `R/brms-postfit.R:38`, `R/autocor-seam.R:233` examples | fits with no prior |
| frmtmb.learn trace identity `sum(log p) == logLik(fit)` | now holds under a prior too |
| frmtmb.survival `rp-check.R`, `royston-parmar.R` | comments and messages about a floored likelihood, no prior |
| frmtmb.sample | builds its own objectives; reads neither `logLik()` nor `opt$objective` of a fit |
| frmtmb.spatial, .coupling, .eam, .ode | no reader |
| equivalence survey L1 and L4 (`dev/equiv-lib.R:132`, `:402`) | under a prior they read `-f$opt$objective`, the penalized value, which this change does not touch; without one `logLik()`, unchanged. No edit needed |
| `frm_compat()` rows | notes of `REML x prior` and `importance x prior` updated (statuses unchanged) |

## Equivalence survey: the 21 sampled pairs that carry a prior

L1 and L4 read `-f$opt$objective` under a prior and `logLik()` without one
(`dev/equiv-lib.R:132`, `:402`), so the split cannot move them; the rerun
checks that by construction. `dev/loglik-r80-equiv-ids.R` lists the pairs
whose spec has a prior (21 of 172); `dev/loglik-r80-equiv-run.sh` reruns
them with `FRMTMB_EQUIV_LIB` = lane lib; rellib-r18;
`dev/loglik-r80-equiv-compare.R` compares with `results.tsv`'s
`verdict_measured`, `level` and `l1_const` (`dev/loglik-r80-equiv-compare.log`):

       id    rec_verdict    now_verdict rec_level now_level       rec_l1
     p001     equivalent     equivalent        L1        L1 -2.84217e-14
     p002         defect         defect        L1        L1  3.99897e-01
     p006     equivalent     equivalent     L1+L2     L1+L2  0.00000e+00
     p008     equivalent     equivalent        L1        L1  5.68434e-14
     p015     equivalent     equivalent        L1        L1  0.00000e+00
     p023 cannot measure cannot measure      <NA>      <NA>  6.70659e-01
     p035     equivalent     equivalent     L1+L2     L1+L2 -3.41061e-13
     p036     equivalent     equivalent     L1+L2     L1+L2 -5.68434e-14
     p037     equivalent     equivalent  L1+L2+L4  L1+L2+L4 -1.13687e-13
     p038     equivalent     equivalent        L1        L1 -1.13687e-13
     p043     equivalent     equivalent        L1        L1 -6.67200e-01
     p050     equivalent     equivalent        L1        L1  1.13687e-13
     p058     equivalent     equivalent        L1        L1 -1.13687e-12
     p071 cannot measure cannot measure      <NA>      <NA>           NA
     p089 cannot measure cannot measure      <NA>      <NA>           NA
     p095     equivalent     equivalent        L1        L1 -3.41061e-13
     p103     equivalent     equivalent        L1        L1  3.28420e+00
     p134     equivalent     equivalent        L1        L1  5.19146e-01
     p140     equivalent     equivalent        L1        L1  0.00000e+00
     p144     equivalent     equivalent        L1        L1  0.00000e+00
     p145     equivalent     equivalent        L1        L1 -1.07461e-01
           now_l1 same
     -2.84217e-14 TRUE
      3.99897e-01 TRUE
      0.00000e+00 TRUE
      5.68434e-14 TRUE
      0.00000e+00 TRUE
      6.70659e-01 TRUE
     -3.41061e-13 TRUE
     -5.68434e-14 TRUE
     -1.13687e-13 TRUE
     -1.13687e-13 TRUE
     -6.67200e-01 TRUE
      1.13687e-13 TRUE
     -1.13687e-12 TRUE
               NA TRUE
               NA TRUE
     -3.41061e-13 TRUE
      3.28420e+00 TRUE
      5.19146e-01 TRUE
      0.00000e+00 TRUE
      0.00000e+00 TRUE
     -1.07461e-01 TRUE
    pairs: 21  verdict and level unchanged: 21  changed: 0 

## Suites

One file per process (`dev/loglik-run-par.sh`, 3 jobs, `NOT_CRAN=true`,
lane library first). Every log's `lib:` line names the lane core;
frmtmb.latent is the lane's, the other seven extensions rellib-r18's
(their code is unchanged). Totals generated from the logs:

    core (dev/loglik-core.log): files=309 pass=22549 fail=3 err=0 skip=239 warn=0
    frmtmb.coupling (dev/loglik-ext.log): files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam (dev/loglik-ext.log): files=33 pass=1848 fail=0 err=0 skip=5 warn=0
    frmtmb.latent (dev/loglik-ext.log): files=12 pass=687 fail=0 err=0 skip=2 warn=0
    frmtmb.learn (dev/loglik-ext.log): files=16 pass=474 fail=0 err=0 skip=13 warn=0
    frmtmb.ode (dev/loglik-ext.log): files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample (dev/loglik-ext.log): files=85 pass=3962 fail=0 err=0 skip=18 warn=0
    frmtmb.spatial (dev/loglik-ext.log): files=8 pass=478 fail=0 err=0 skip=6 warn=0
    frmtmb.survival (dev/loglik-ext.log): files=15 pass=552 fail=0 err=0 skip=0 warn=0
    gated subset, FRMTMB_BRMS_FIT_TESTS=true (dev/loglik-gated.log): files=38 pass=4468 fail=1 err=0 skip=0 warn=0

The three core failures and the one gated failure are the four tests
that pinned the old meaning ("Tests" above); each was edited and passes
on rerun (`dev/loglik-fix.log`, `dev/loglik-final.log`,
`dev/loglik-t-out/gated-brms-priors.txt`):

    RESULT test-compat-extensions.R pass=2 fail=0 err=0 skip=0 warn=0
    RESULT test-compat-register.R pass=103 fail=0 err=0 skip=0 warn=0
    RESULT test-compat.R pass=759 fail=0 err=0 skip=0 warn=0
    RESULT test-prior-compat.R pass=202 fail=0 err=0 skip=0 warn=0
    RESULT test-loglik-prior.R pass=49 fail=0 err=0 skip=0 warn=0
    RESULT test-map.R pass=8 fail=0 err=0 skip=0 warn=0
    RESULT test-brms-priors.R pass=102 fail=0 err=0 skip=0 warn=0

The gated subset is the 38 core files that mention a prior and one of
logLik, AIC, BIC, summary, brms_lp_check, extractAIC, step, allfit,
anova or drop1 and are gated (`dev/loglik-r80-gated.jobs`); skip 0.

One session per suite, three orders (`dev/loglik-run-sessions.sh`,
`dev/loglik-sessions.log`); every LEAK line is tagged
`[with load of ...]`, none is this lane's:

    SESSION frmtmb forward files=296/309 pass=22554 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb reverse files=296/309 pass=22704 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb preload files=296/309 pass=22797 fail=0 err=0 skip=241 warn=0
    SESSION frmtmb.latent forward files=12/12 pass=687 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent reverse files=12/12 pass=687 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent preload files=12/12 pass=687 fail=0 err=0 skip=2 warn=0
    SESSIONS ran 6 of 6

files=296/309: the 13 others skip at file level without the gate
(test-brms-suite-*.R and three *-brms.R files): they hold gated tests only.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
`dev/loglik-openblas.log`), the four files with exactness checks:

    RESULT test-map.R pass=8 fail=0 err=0 skip=0 warn=0
    RESULT test-nl-rtmb-scope.R pass=105 fail=0 err=0 skip=0 warn=0
    RESULT test-prior-compat.R pass=202 fail=0 err=0 skip=0 warn=0
    RESULT test-loglik-prior.R pass=42 fail=0 err=0 skip=0 warn=0

## A defect R CMD check found, fixed

The first `R CMD check --as-cran` of frmtmb.latent
(`dev/loglik-r80-check-latent.log`, first run) gave 3 WARNINGs and a NOTE
from one cause: frmtmb's namespace did not load with only base attached,
`object 'BIC' not found whilst loading namespace 'frmtmb'`. The new
`S3method(AIC, frmtmb_fit)` and `S3method(BIC, frmtmb_fit)` need the
generics visible in the namespace, and frmtmb imported neither from
stats. No test saw it, because every test session attaches stats. Fixed
by adding `AIC BIC` to the `@importFrom stats` line of
`R/frmtmb-package.R` (so `NAMESPACE` imports them). After the fix, with
`R_DEFAULT_PACKAGES=NULL`, `loadNamespace("frmtmb")` and
`loadNamespace("frmtmb.latent")` both load (`dev/loglik-r80-install3.log`
and the command in the lane's report).

## Cost of the split

Superseded by "Fix round 1", m8: the entries are now kept on the fit.

`dev/loglik-r80-timing.R` (seed 4101, one process, arms interleaved,
blocks over 1.2 s, minimum of 5 rounds; `dev/loglik-r80-timing.log`):

    control_noprior  min 134.3 us  median 180.7 us
    prior            min 2392.6 us  median 4013.7 us
    penalized_prior  min 2685.5 us  median 3681.6 us
    one MAP fit: 0.050 s

`logLik()` under a prior resolves the prior again on every call
(`resolve_prior_input()`, as `fit_prior_density()` already did for
`anova()`), about 2.3 ms more per call, or 1/20 of one fit of this
model. `summary()` pays it three times. Not cached: a cache keyed on the
fit's estimates would have to be invalidated by every path that replaces
`fit$opt`, and the cost is small against any fit.

## Bump and floor

- frmtmb: minor. Two BREAKING changes (`logLik()` under a prior; AIC,
  BIC, extractAIC refuse) and a new export.
- frmtmb.latent: minor, `hmm_starts()` ranks on a different value under a
  prior. Floor: frmtmb at the release that exports `logLik_penalized()`
  (the version after 0.79.1). I did not edit either DESCRIPTION: the
  round sets versions at consolidation. Roxygen 8.1.1 rewrote
  `Config/roxygen2/version` to 8.1.1 in both DESCRIPTIONs; I set it back
  to 8.1.0, so neither DESCRIPTION is in the diff.
- The other seven extensions: no change, no floor. Their suites pass on
  the lane core (above).

## Not done, and why

Fix round 1 closed the first three items and the restart label below
(M2, m6, m7); the gated-tier and seed items stand.

- The restart messages of a nonlinear fit (`R/par-template.R`,
  "log-likelihood X") print the penalized value under a prior and still
  call it "log-likelihood". The number is the one that was compared, as
  before this change; rewording five messages was outside the split.
- `AIC(lm_fit, map_fit)` (a non-frmtmb fit first) dispatches to
  `stats:::AIC.default`, which calls `logLik()` on each and gives the
  data value with no refusal; so does `AIC(logLik(map_fit))`, and any
  package that builds AIC from `logLik()` (performance, MuMIn). Closing
  those would need a marked `logLik` class, which `print`, `all.equal`
  and other packages read; I judged it not worth it.
- No `logLik_penalized()` method for `frm_multiple()` fits; the generic
  gives "no applicable method". `logLik()` on them already refuses.
- The change to the nonlinear restarts (`logLik_penalized()` in
  `nl_restart_round()`) has no test that fails without it: on the base
  build `logLik()` was the penalized value, so the old and new code agree
  there, and under the lane's `logLik()` the old code would rank restarts
  on the data value. Same for `frm_allfit()`: its new test is an identity
  check, not a seen-to-fail test.
- Gated tier: the 38-file subset above, not the whole gated tier, and no
  extension's gated tests.
- The brief's own seed-4101 numbers (-297.44 ...) were not reproduced:
  no script in `dev/` has them; mine are from my own construction.

## Process notes

- Once, four R processes ran for about a minute: a single-file rerun of
  `test-map.R` while three suite processes ran.
- `dev/loglik-r80-equiv-run.sh` first piped the id list into the loop,
  so `wait` returned at once and the driver printed "done" while the last
  three pairs ran; they finished and wrote their results (21 of 21).
  Fixed in the script.
- The first core `R CMD build` copied `dev/ciharden-out` while I deleted
  it, and failed ("copying to build directory failed"). I stopped that
  task with TaskStop; its `R CMD build` children lived on, and
  `Stop-Process` on them was denied by the permission system, so I let
  them end on their own (they did) and reran the checks one after the
  other.
- `dev/ciharden-out/` is deleted.

## R CMD check --as-cran (final code)

`dev/loglik-r80-rcheck.sh`, built and checked in `dev/loglik-check/`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
`R_LIBS` = lane lib; rellib-r18; user library. Run after the
`importFrom` fix, so its test step ran the whole core suite on the final
code (`Running 'testthat.R' [41m] OK`):

    frmtmb         Status: 1 NOTE   (dev/loglik-r80-check-frmtmb.log)
    frmtmb.latent  Status: 1 NOTE   (dev/loglik-r80-check-latent.log)

Both NOTEs are "Skipping checking math rendering: package 'V8'
unavailable", the expected one. The check trees are deleted; the logs
stay.

## Fix round 1 (review `dev/reviews/2026-11-06-loglik.md`, verdict FIX)

Every item, with what was done and the evidence. "Pre-fix" is the
lane's code before this round, as the reviewer built it in
`rev-loglik-lib` (`REV_MUT` unset); "mutant" is the same build with one
`REV_MUT` switch.

- **M1, pkgdown index.** `logLik_penalized` added under "Inference and
  diagnostics" in `_pkgdown.yml`, after `anova.frmtmb_fit`.
  `dev/loglik-r80-pkgdown.R` (`dev/loglik-r80-pkgdown.log`): pkgdown
  2.2.1, `check_pkgdown: OK` ("No problems found"); the review measured
  "1 topic missing from index" before.
- **M2, criteria that dispatch elsewhere.** Under a prior with a density
  `logLik()` now has `df = NA` and a `df_why` that says why
  (`prior_df_why`, `R/methods-fit.R`); the value is still the data
  log-likelihood. Every internal reader of the count now calls
  `loglik_df()`: `pooled_order()` and `pooled_lrt_parts()`
  (`R/multiple.R`), `anova.frmtmb_fit` and `drop1()` (`R/confint.R`),
  `extractAIC()` (`R/sugar.R`). `dev/loglik-r80-dfna-probe.R`, seed 201
  (`dev/loglik-r80-dfna-probe.log`, lines shortened):

      print(logLik(fp))            'log Lik.' -192.5685 (df=NA)
      AIC(fp)                      ERROR <frmtmb_prior_ic> ...
      AIC(lm, fp)                     df      AIC | fl  3 392.2962 | fp NA       NA
      BIC(lm, fp)                     df      BIC | fl  3 400.6586 | fp NA       NA
      AIC(logLik(fp))              [1] NA
      BIC(logLik(fp))              [1] NA
      anova(fp, fp2)               ERROR <frmtmb_error> anova() compares log-likelihoods, and model 1 has a prior ...
      lmtest::lrtest(fp, fp2)      #Df blank; LogLik -192.57, -192.52; Chisq 0.096; Pr(>Chisq) blank
      print(frm_allfit(fp))        original fit: logLik_penalized -196.7122664, code 0 ...
      anova(m1, m0) D1             statistic 0.07733, df1 1, p 0.7809 (the Wald rule still runs)
      anova(m1, m0) D3             ERROR <frmtmb_error> anova() compares log-likelihoods, and a pooled fit has a prior ...
      logLik_penalized(m1)         ERROR <frmtmb_error> logLik_penalized() has no pooled answer ...
      df.residual(fp)              [1] 116
      summary(fp) Method line      " Method: ML (MAP)   logLik: -192.569   log prior: -4.14372 "
      frmtmb.sample functions whose body calls logLik(: none

  Nothing errors unhelpfully: `lmtest::lrtest()` prints the statistic
  with no df and no p-value. The review measured `AIC(lm, fp)` and
  `AIC(logLik(fp))` at 301.5597279 for `fp` before. NEWS and
  `?logLik_penalized` now say which calls stop and which get `NA`.
  Pinned by the block "criteria that dispatch elsewhere read NA, never a
  number" of `test-loglik-prior.R` (`AIC(lm, fp)`, `BIC(lm, fp)`,
  `AIC(logLik(fp))`, and the no-prior control).
- **M3, nlstart merge hazard.** Not in this tree. Consolidation must
  change `nl_sign_probe()` in `frmtmb-wt-nlstart/R/par-template.R`
  (its lines 90 and 91, `ll0 <- as.numeric(stats::logLik(fit))` and
  `ll1 <- ... stats::logLik(r)`) to `logLik_penalized()`, because the
  probe selects a refit on the penalized objective (its lines 62 to 68),
  and label its `frmtmb_nl_multistart` message "penalized
  log-likelihood" under a prior (`fit_penalized(fit)`), as
  `nl_fit_off_saddle()` now does.
- **m1, negative control for the discrete check.** The reviewer's case A
  (50 zero counts, `normal(-5, 0.01)` on the Intercept: penalized +3.35,
  data -0.34) is in the block "the discrete check runs under a prior" as
  `expect_no_warning()`. Run against the review build
  (`dev/loglik-t-out/mut-*.txt`): with `REV_MUT=discpen` that block
  fails, "Actually got a <frmtmb_discrete_loglik> with message: The
  log-likelihood is 3.34934, above 0"; with `REV_MUT` unset it passes.
  The family type is now checked first in `discrete_loglik_warning()`.
- **m2, Jacobian rows.** New block "the split holds through every
  chart's Jacobian": `mo()` with a `simo` Dirichlet, `ar` with a prior on
  class `ar`, exact `gp()` with `sdgp` and `lscale` priors, and
  `(1 + x | g)` with `lkj(4)` and an `sd` prior. Each `logLik()` equals
  the no-prior tape at the MAP to `sqrt(eps)` relative, with a nonzero
  log prior.
- **m3, uniform priors.** `prior_entry_constant()`: an entry is a
  constant when its density is `uniform()` and it is read with no change
  of variables (scale not `sd`, `natural` or `ordthres`; a `trans` chart
  is not of kind `uniform`). `fit_prior_density()` ignores such entries,
  so AIC, BIC, `anova()`, `drop1()` and the summary treat the fit as
  bounded ML. Where the reasoning fails the refusal stays: a uniform
  density on an `sd`, or on any parameter read through a log-Jacobian,
  is not constant on the optimizer's scale and moves the estimate. New
  block "a uniform prior on a coefficient is bounds, not a density":
  `uniform(-5, 5)` on `b` gives the no-prior slope to 1e-6, a finite df,
  `AIC()` equal to the no-prior one to 1e-8, `log_prior` = log(1/10) and
  method "ML"; `uniform(0, 10)` on `sd` is refused, with `df = NA`.
- **m4, meaning of `log_prior`.** Kept: the value the optimizer adds,
  with every log-Jacobian and no truncation normalizer. New section "The
  log prior term" in `?logLik_penalized` names the charts and works the
  difference from brms's `lprior` for `normal(0, 1)` on an `sd`:
  `+ log(s)` here, `+ log(2)` in brms (the review's `brmscode.log` and
  `discrete-after.log`). NEWS says it in one sentence.
- **m5, summary under REML with a b prior.** Prints
  `logLik_penalized: <value>` and "logLik: none, as the prior on the
  fixed effects is inside the REML integral". Pinned in the REML block.
- **m6, labels.** `nl_fit_off_saddle()` and `nl_saddle_warning()` say
  "penalized log-likelihood" when the maximized value differs from
  `logLik()` (`fit_penalized()`, which also covers a uniform prior's
  constant). `print.frmtmb_allfit()` names the column, the original
  line and the spread `logLik_penalized` in that case. Pinned by the
  block "the nonlinear restart message names the penalized value" (the
  reviewer's seed 501 construction): the message says "penalized
  log-likelihood", and its last number is `logLik_penalized()` of the
  fit to 1e-5 relative, away from `logLik()`.
- **m7.** `logLik_penalized.frmtmb_multiple` refuses through
  `multiple_refuse()`, as `logLik()` does. Pinned.
- **m8, cost.** The resolved prior entries stay on the fit
  (`fit$resolved_prior`; not `prior_*`, which `$prior` would
  partial-match), so `fit_log_prior()` is one closure call;
  `fit_prior_entries()` resolves again only for an object without the
  slot. `dev/loglik-r80-cost.R` (the review's script with three arms,
  seed 401, 7 blocks, run base, prefix, after, prefix, after;
  `dev/loglik-r80-cost.log`); `big_prior` is the review's 24-row prior.
  `proc.time()` ticks at 10 ms, so these are block means, min / median:

      arm      big_prior logLik()   big_prior summary()   big_prior print(fit)
      base     0.00 / 0.00 ms       2.00 / 4.00 ms        4.00 / 6.00 ms
      prefix   3.00 / 5.50 ms       18.00 / 22.00 ms      16.00 / 24.00 ms
      after    0.00 / 0.00 ms       2.00 / 6.00 ms        4.00 / 6.00 ms
      prefix   5.50 / 6.00 ms       18.00 / 22.00 ms      18.00 / 22.00 ms
      after    0.00 / 0.00 ms       2.00 / 4.00 ms        2.00 / 6.00 ms

  The fix-round build is back at the base build's cost.
- **m9, floor.** For consolidation: frmtmb.latent's DESCRIPTION floor
  must rise from `frmtmb (>= 0.79.0)` to the frmtmb release that exports
  `logLik_penalized()`; `hmm_starts()` calls it.

Tests of this round seen to fail on the pre-fix build
(`dev/loglik-t-out/mut-none.txt`: `RESULT test-loglik-prior.R pass=68
fail=4 err=3`): the REML summary block (2 fails), the dispatch block (df
not NA, then no `df_why`), the uniform block (`AIC()` refused), the
`frm_multiple()` block ("no applicable method"), the restart-message
block (the message says "log-likelihood").

Fix-round runs, one file per process, lane library first (`dev/loglik-fix2.log`,
`dev/loglik-core2.log`, `dev/loglik-ext2.log`; every `lib:` line names the lane core):

    RESULT test-loglik-prior.R pass=88 fail=0 err=0 skip=0 warn=3
    RESULT test-map.R pass=8 fail=0 err=0 skip=0 warn=0
    RESULT test-nl-rtmb-scope.R pass=105 fail=0 err=0 skip=0 warn=0
    RESULT test-prior-compat.R pass=202 fail=0 err=0 skip=0 warn=0
    RESULT test-refit-data.R pass=122 fail=0 err=0 skip=0 warn=0
    RESULT test-beta-binomial-density.R pass=93 fail=0 err=0 skip=0 warn=0
    RESULT test-arg-refusal.R pass=122 fail=0 err=0 skip=0 warn=0
    (test-loglik-prior.R rerun with explicit bounds on the uniform rows:
    RESULT test-loglik-prior.R pass=88 fail=0 err=0 skip=0 warn=0)
    core suite: files=309 pass=22600 fail=0 err=0 skip=239 warn=0
    frmtmb.coupling: files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam: files=33 pass=1848 fail=0 err=0 skip=5 warn=0
    frmtmb.latent: files=12 pass=687 fail=0 err=0 skip=2 warn=0
    frmtmb.learn: files=16 pass=474 fail=0 err=0 skip=13 warn=0
    frmtmb.ode: files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample: files=85 pass=3962 fail=0 err=0 skip=18 warn=0
    frmtmb.spatial: files=8 pass=478 fail=0 err=0 skip=6 warn=0
    frmtmb.survival: files=15 pass=552 fail=0 err=0 skip=0 warn=0

R CMD check --as-cran of core on the fix-round code
(`dev/loglik-r80-check-frmtmb2.log`): `Running 'testthat.R' [37m] OK`,
`Status: 1 NOTE` (V8 math rendering). frmtmb.latent is unchanged in
this round; its suite passes on the new core (above).

## After the short check (SHIP), two text items

- NEWS and `?logLik_penalized` now name exactly which calls give `NA`
  (`AIC(lm_fit, map_fit)`, `BIC(lm_fit, map_fit)`, `AIC(logLik(map_fit))`),
  which refuse (`AIC()`, `BIC()`, `extractAIC()`, `step()` dispatched on
  the fit), and that a package counting df itself is not covered
  (`insight::get_loglikelihood()`, so `performance`). Text only.
- `?anova.frmtmb_multiple` and the NEWS bullet on pooling now say that
  under a prior the Wald rules (D1, D2 with `use = "wald"`) are Wald tests
  at the MAP, consistent with `summary()`'s Wald intervals, with nominal
  size not guaranteed. Text only. Rds rendered with `Rd2txt` and read.

# Earlier lane loglik

# Lane loglik: log_lik() for residual covariance terms, and the 13 bin-2 sar() and fcor() rows

Base 70476dc3 (frmtmb 0.75.0), worktree `frmtmb-wt-loglik`, branch
`wt-loglik`, library `C:/Users/adf44/source/r/wt-loglik-lib` (core,
frmtmb.sample and frmtmb.spatial from this tree) in front of
`rellib-r13`. Nothing is committed. Logs are under `dev/loglik-*.txt`,
`dev/loglik-log/` and `dev/release/loglik-*` (all gitignored); every
number below is pasted from them, with its script and seed.

The brief named `C:\Users\adf44\source\r\frmtmb-wt-sarfcor\dev\sarfcor-brmsrows.R`.
That worktree was removed during the round; the tracked copy
`dev/sarfcor-brmsrows.R` in this tree is the one meant (the
coordinator's correction).

## Part A. brms's leave-one-out conditional density

### A1. What brms 2.23.0 does (run, not recalled)

brms's `R/log_lik.R` (`dev/brms-suite/brms/R/log_lik.R`, lines 246 to
412 and 1078 to 1100): for `ar()`, `ma()`, `arma()` (`cov = TRUE`),
`cosy()` and `unstr()` one call per time-series GROUP
(`log_lik_gaussian_time()`, `begin_tg:end_tg`); for `fcor()` and both
SAR types one call over all rows (`stopifnot(i == 1)`). Per draw, a
dense `solve()` of the group's covariance `C`, then
`g = C^-1 e`, `cbar = diag(C^-1)`, `yloo = y - g / cbar`,
`sdloo = 1 / sqrt(cbar)`; for student `dfloo = nu + k - 1` and
`sdloo^2 = (nu + e'g - g^2 / cbar) / (nu + k - 1) / cbar`
(`student_t_cov_factor()`, Shah et al. 2014 eq. 6), `k` the group size.
For the lag SAR, `e = y - (I - rho W)^-1 mu`. `get_cov_matrix_ac()`
reads lags by position within the group (only `unstr()` carries
`Jtime`), and `sigma ~ x` enters as `outer(sigma, sigma) * cor`.

brms was run on every case (`dev/loglik-brms-probe.R`, seed 1, log
`dev/loglik-brms-probe.txt`; grouped student cases again in
`dev/loglik-brms-probe2.R`, `dev/loglik-brms-probe2.txt`): a brms fit
whose draws are two chosen parameter vectors (Fixed_param, one chain
each; `dev/loglik-brms-lib.R`), and `log_lik()` against the conditional
computed directly from the joint covariance by the partition formula,
with no precision matrix. 24 rows in 5 groups of sizes 6, 5, 4, 6, 3,
rows shuffled:

    ar gaussian                draws 2  rows 24  max|brms - dense| / max|dense| = 3.981e-16
    ma gaussian                draws 2  rows 24  max|brms - dense| / max|dense| = 2.876e-16
    arma gaussian              draws 2  rows 24  max|brms - dense| / max|dense| = 3.597e-16
    cosy gaussian              draws 2  rows 24  max|brms - dense| / max|dense| = 2.774e-16
    fcor gaussian              draws 2  rows 24  max|brms - dense| / max|dense| = 6.253e-16
    fcor student               draws 2  rows 24  max|brms - dense| / max|dense| = 6.570e-16
    lag gaussian               draws 2  rows 24  max|brms - dense| / max|dense| = 4.111e-16
    lag student                draws 2  rows 24  max|brms - dense| / max|dense| = 4.396e-16
    err gaussian               draws 2  rows 24  max|brms - dense| / max|dense| = 3.378e-16
    err student                draws 2  rows 24  max|brms - dense| / max|dense| = 3.768e-16
    ar sigma~x gaussian        draws 2  rows 24  max|brms - dense| / max|dense| = 5.196e-16
    unstr gap gaussian         draws 2  rows 22  max|brms - dense| / max|dense| = 2.670e-16
    ar gap gaussian: brms vs by-position 2.955e-16, vs by-level 6.752e-02
    ar gaussian draw 1: sum of conditionals -25.055181, joint -26.097025

and, with the reference taken group by group (probe2):

    ar student                 draws 2  rows 24  max|brms - dense| / max|dense| = 5.532e-16
    ma student                 draws 2  rows 24  max|brms - dense| / max|dense| = 4.569e-16
    arma student               draws 2  rows 24  max|brms - dense| / max|dense| = 5.160e-16
    cosy student               draws 2  rows 24  max|brms - dense| / max|dense| = 4.352e-16
    ar sigma~x student         draws 2  rows 24  max|brms - dense| / max|dense| = 5.991e-16
    unstr gap student          draws 2  rows 22  max|brms - dense| / max|dense| = 5.879e-16

So brms's `log_lik()` is the exact conditional for every term and both
families; no brms defect here. The first probe's grouped student rows
differed by 0.10 to 0.17: its reference was one multivariate t over all
24 rows, a different model from brms's one t per group (each group
draws its own chi-square), so those numbers say nothing about brms and
probe2 replaces them. brms's `log_lik()` works on the error SAR, where
its `posterior_epred()` fails (brms-39, already filed). The ARMA(1,1)
covariance brms uses has gamma0 = (1 + 2 phi theta + theta^2) / (1 -
phi^2) at unit innovations, which the probe's reference assumed and
which agreed.

### A2. The route, and why it is exact and cheap

`R/autocor-loglik.R`, `autocor_row_loglik(fit, dpv, resp)`, exported
in the sampling API (`R/sampling-api.R`). Everything is read off the
precision `Q` of the response; neither `Q` nor the covariance of more
than one group is formed:

- the `ar()` family (`ar`, `ma`, `arma` with `cov = TRUE`, `cosy`,
  `unstr`): `Q = D^-1 R^-1 D^-1` per group, `D = diag(sigma)`. One
  `chol2inv(chol(R_p))` per PATTERN of present time levels (the
  patterns the objective already uses), shared by every group of the
  pattern, applied to the `k x G` matrix of standardized residuals in
  one product. Unequal group sizes, missing time points and `sigma ~ x`
  come from the patterns and the per-row sigma;
- a registered term (`fcor()`, `sar()`, any `frmtmb_register_autocor()`
  term): the whitening map is affine, `u = K (y - m)`, so `Q = K'K`,
  `g = K'u`, `cbar = colSums(K^2)`, `e'Qe = u'u`. `K` comes from a new
  documented slot, `frmtmb_register_autocor(jacobian =)`; `fcor()`
  gives `L^-1 D^-1`, frmtmb.spatial's `sar()` gives the sparse
  `(I - rho W) / sigma` (both types: they differ in the mean only).
  Without the slot, `K = B^-1` with `B` from `unwhiten()` of the
  identity: exact, dense;
- the density is written in the residual `g / cbar` itself (which is
  `y - yloo`), so no `y` is subtracted from a nearly equal `yloo`.

The review's estimate (`dev/reviews/2026-10-24-sarfcor.md`, finding 12)
was this route; brms's grouping and its per-group t were the parts to
verify, and A1 verified them.

Cost, measured (`dev/loglik-cost.R`, seed 5, log `dev/loglik-cost.txt`;
arms interleaved in one process, each block repeated past 1.2 s, min of
5 rounds, a control arm that repeats the first; run alone on the
machine's other load):

    sar n = 1600  nnz(W) = 6240
    agreement: sar vs brms 9.107e-16, sarfb vs sar 6.960e-15 (max abs / max)
    sar, per draw, min over 5 rounds:
      sar        1.924 ms per draw (reps 1024)  alloc    1809120 bytes
      sarfb   1730.000 ms per draw (reps    1)  alloc  207899968 bytes
      brms    4120.000 ms per draw (reps    1)  alloc  102560792 bytes
      ctrl       1.895 ms per draw (reps  512)  alloc    1809120 bytes
      brms / sar = 2141.6, sarfb / sar = 899.2, ctrl / sar = 0.98
    ar agreement with brms's algorithm: 6.193e-16
    ar(cov = TRUE), 400 groups x 10, per draw:
      ar         0.635 ms per draw (reps 2048)  alloc     729536 bytes
      brms      17.813 ms per draw (reps   64)  alloc   15536048 bytes
      ctrl       0.605 ms per draw (reps 2048)  alloc     729536 bytes
      brms / ar = 28.1, ctrl / ar = 0.95
    fcor agreement with brms's algorithm: 1.843e-15
    fcor(), 600 rows, per draw:
      fcor       2.637 ms per draw (reps  512)  alloc    5844960 bytes
      brms      97.500 ms per draw (reps   16)  alloc   14482512 bytes
      ctrl       2.695 ms per draw (reps  512)  alloc    5844960 bytes
      brms / fcor = 37.0, ctrl / fcor = 1.02

`sar` is the term with its jacobian slot, `sarfb` the same map
registered without it (the fallback), `brms` brms's algorithm
transcribed (`log_lik_gaussian_lagsar()` and `log_lik_gaussian_time()`,
dense). Allocation is `Rprofmem()` bytes per call. The control ratios
(0.95 to 1.02) bound the clock's noise; the ratios that matter are 28
to 2142. The fallback is 900 times the slot on the 1600-row lattice,
which is why `sar()` fills the slot; a term registered without one
still gets exact values. `fcor()` allocates two `n x n` matrices per
draw (`K` and `K * K`); a precomputed `colSums(L^-1 ^2)` would save one
of them, not done because the term is dense by nature and is 37 times
brms already.

### A3. Agreement with brms at the same draws

`dev/loglik-brms-agree.R` (seed 4, log `dev/loglik-brms-agree.txt`):
frmtmb.sample's `log_lik()` on hand-set draws (no sampler) against
brms's `log_lik()` on a brms fit whose three chains sit at the same
three parameter vectors. 43 rows in 8 groups of sizes 7, 5, 6, 4, 7,
3, 6, 5, rows shuffled; the `cosy()` and `unstr()` rows on the same
data with two time points removed (41 rows):

    ar gaussian            3 x 43  max|frmtmb - brms| / max|brms| = 2.033e-15
    ar student             3 x 43  max|frmtmb - brms| / max|brms| = 8.183e-16
    ma gaussian            3 x 43  max|frmtmb - brms| / max|brms| = 6.853e-16
    ma student             3 x 43  max|frmtmb - brms| / max|brms| = 7.070e-16
    arma gaussian          3 x 43  max|frmtmb - brms| / max|brms| = 1.185e-15
    arma student           3 x 43  max|frmtmb - brms| / max|brms| = 7.988e-16
    ar sigma~x gaussian    3 x 43  max|frmtmb - brms| / max|brms| = 6.163e-16
    ar sigma~x student     3 x 43  max|frmtmb - brms| / max|brms| = 4.288e-16
    cosy gaps gaussian     3 x 41  max|frmtmb - brms| / max|brms| = 6.024e-16
    cosy gaps student      3 x 41  max|frmtmb - brms| / max|brms| = 4.040e-16
    unstr gaps gaussian    3 x 41  max|frmtmb - brms| / max|brms| = 9.479e-16
    unstr gaps student     3 x 41  max|frmtmb - brms| / max|brms| = 5.220e-16
    fcor(V) gaussian       3 x 43  max|frmtmb - brms| / max|brms| = 1.210e-15
    fcor(V) student        3 x 43  max|frmtmb - brms| / max|brms| = 9.476e-16
    sar(W) gaussian        3 x 43  max|frmtmb - brms| / max|brms| = 6.347e-16
    sar(W) student         3 x 43  max|frmtmb - brms| / max|brms| = 5.366e-16
    sar(W, type = 'error') gaussian 3 x 43  max|frmtmb - brms| / max|brms| = 5.928e-16
    sar(W, type = 'error') student 3 x 43  max|frmtmb - brms| / max|brms| = 4.822e-16

The translation between the two parameterizations is the one the
log-density tier uses: brms's `ar()`, `ma()` and `arma()` sigma is the
innovation sd, frmtmb's the marginal one (`sigma_brms = sigma /
sqrt(gamma0)`, and for `sigma ~ x` the intercept shifts by
`log sqrt(1 - ar^2)`); brms's `Lcortime` is the Cholesky factor of
frmtmb's `unstr()` matrix; `lagsar`/`errorsar` are rho. The ar() rows
use data with no interior gap: with one, brms (lags by position) and
frmtmb (lags by time level, `?frmtmb-autocor`) fit different models,
and the probe's "ar gap" line shows brms's values are the by-position
conditional and differ from the by-level one by 6.8e-2. `cosy()` has
one correlation for every pair and brms's `unstr()` carries the time
index, so those agree with gaps.

The same comparisons are tests (A5). The ML-estimate check
(`dev/loglik-check.R`, seed 3, `dev/loglik-check.txt`) compares
`autocor_row_loglik()` at each fit's estimates with the dense reference
built from the fit's own covariance, for 10 formulas x 2 families on
ragged groups with an interior gap: 1.3e-16 to 8.7e-16, except
`unstr()`, whose estimates on 7 levels over 8 groups are near singular
(false convergence; the reference itself is NaN), so the tests set
`unstr()`'s correlations by hand. The same log shows the row sums are
not the joint log density, e.g. `fcor` gaussian row sum -64.1559,
joint -62.6940.

### A4. The seam slot, and core's log_lik() on an ML fit

`frmtmb_register_autocor(jacobian = function(sigma, theta, aux))`:
documented in a new section "The Jacobian" of its help page. When the
frame is built a given slot is checked once on numbers against the map
(`K (y - m) = whiten(y)$u`, relative 1e-6), and a wrong one is refused
by the term's name. frmtmb.spatial fills it through the exported
function alone: no `frmtmb:::` (grep of `extensions/frmtmb.spatial/R`
for `frmtmb:::` finds none).

Core's `log_lik()` on a `frmtmb_fit` keeps its refusal: it refuses for
EVERY model ("the pointwise likelihood at every posterior draw ... this
is a maximum-likelihood fit"), and making it answer for the covariance
terms alone would make the one family of models where a fit has a
pointwise density the exception. The computation at one parameter
vector is core's exported `autocor_row_loglik()`, which on an ML fit is
the value at the estimates (the A3 check runs it that way). The
coordinator decides whether the brief meant more.

### A5. Tests

New:

- `tests/testthat/test-autocor-loglik.R` (core): the dense conditional
  for ar, ar(p = 2), ma, arma, cosy, ar + (1 | g), ar with sigma ~ x
  (gaussian and student), unstr at hand-set correlations, fcor with
  and without sigma ~ x and student; the row sum against the joint
  density (it differs by more than 1 percent at ar = 0.8, and equals it
  at ar = 0, where each conditional is the row's own density; the joint
  is checked against `build_objective()` too); the jacobian slot and the
  fallback against fcor(); a wrong jacobian refused, with fcor()'s own
  slot as the absent case; the refusal on a response without a
  covariance term.
- `extensions/frmtmb.sample/tests/testthat/test-loglik-autocor.R`:
  `log_lik()` of hand-set draws against the dense conditional for ar,
  ma, arma (with the correlation from `stats::ARMAacf()`, not from
  frmtmb), cosy and unstr with gaps, fcor with sigma ~ x, gaussian and
  student, each column identical to core's function at the draw; a
  registered term without a jacobian slot; `loo()`, `waic()`,
  `loo_epred()`, `loo_predict()` on a covariance block; and, gated,
  brms's `log_lik()` at the same draws for every row of A3 but sar.
- `extensions/frmtmb.spatial/tests/testthat/test-sar.R`, two new
  blocks: the dense conditional for both types, both families, and the
  row sum against the joint; gated, brms's `log_lik()` at three
  parameter vectors per type and family.

Changed (the refusals removed, each pin now an agreement test):

- `extensions/frmtmb.sample/tests/testthat/test-loo.R`, "log_lik()
  refuses likelihoods with no per-observation column": the
  `ar(cov = TRUE)` arm asserted the refusal; it now checks the columns
  against the dense conditional of each subject at four draws and that
  `loo()` is finite. The mixture arm is unchanged.
- the same file, "loo() runs on a cov = FALSE fit and cov = TRUE still
  refuses", renamed "... and on a cov = TRUE one": the three refusal
  messages it pinned (ar, unstr, and that unstr's message offers no
  cov = FALSE form) are replaced by: each cov = TRUE column is core's
  function at the draw, it differs from the cov = FALSE column at
  ar = 0.5 and equals it at ar = 0, `loo()` is finite, and unstr gives
  finite columns.

Seen to fail on rellib-r13 (logs `dev/loglik-log/base-*.txt`):

    RESULT frmtmb test-autocor-loglik.R pass=0 fail=0 err=7 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=0 fail=0 err=6 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=132 fail=0 err=2 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=93 fail=0 err=2 skip=0 warn=0

frmtmb.sample's errors are the behavioral refusal ("log_lik() needs a
likelihood that factors into one term per observation, and the
response 'y' does not: it carries ar(t, g, cov = TRUE) ..."), for every
block, the registered term without a jacobian among them. Core's and
frmtmb.spatial's are the weak form, "could not find function
autocor_row_loglik"; the behavioral failure for `sar()` is the same
frmtmb.sample refusal, which the registered lag term of
test-loglik-autocor.R shows.

On the lane build (one file per process, gated):

    RESULT frmtmb test-autocor-loglik.R pass=47 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=93 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=146 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=109 fail=0 err=0 skip=0 warn=0

## Part B. The 13 bin-2 sar() and fcor() rows

Adapter (`dev/brmsport-bin2-adapters.R`, section "Residual correlation
terms"; spec `dev/brmsport-bin2-adapters-spec.R`): brms's
`log_lik_<family>_<term>(1, prep)` and
`posterior_predict_<family>_<term>(1, prep)` for `lagsar`, `errorsar`
and `fcor` become `brms_adapt_log_lik()` and
`brms_adapt_posterior_predict()`, which build a frmtmb model of the
prep's shape (response `Y`, one coefficient per row set to brms's `mu`
of the draw, the term over brms's matrix, brms's sigma, nu and rho per
draw) and return frmtmb's `autocor_row_loglik()` or `sim_draw()`.
ADAPTATION, recorded in the spec header: the SAR blocks'
`Msar = diag(10)` becomes `brms_adapt_sar_path(10)`, the
row-standardized path graph over the same 10 rows, because frmtmb
refuses `diag(10)` by design (`I - rho I` makes rho and the scale one
parameter). The substitution applies to the setup line, so brms's own
function in the value checks reads the same matrix. Blocks added to
`bin2_adapted`: "log_lik for SAR models runs without errors", "log_lik
for FCOR models runs without errors", "log_lik for FCOR student-t models
works as expected", "loglik for SAR models runs without errors" (brms's
label in tests.posterior_predict.R), "posterior_predict for FCOR models
runs without errors".

`Rscript dev/brmsport-bin2-gen.R`: "wrote 6 files; 238 assertion ids;
13 without a verdict; 149 adapted; 71 value checks". Only
`dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R` (+199 lines)
and `-posterior_predict.R` (+90) changed. The 13 without a verdict are
the 13 rows; the ledger and `dev/brmsport-bin2-verdicts.tsv` were not
touched, so the rows are "unclassified" until consolidation.

Recorded on the lane build (`FRMTMB_PORT_LIB` = the lane library,
`R_LIBS_SITE` = rellib-r13 for the extensions the other blocks load;
logs `dev/loglik-log/run-frmtmb-*.txt`, records
`dev/loglik-log/rec-frmtmb-*.tsv`):

    RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=97 fail=0 err=0 skip=0 blocks=15
    RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=61 fail=0 err=0 skip=0 blocks=14
    assert frmtmb log_lik:150 unclassified TRUE
    value  frmtmb log_lik:150  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 4.44e-16, max rel 1.44e-16
    assert frmtmb log_lik:152 unclassified TRUE
    value  frmtmb log_lik:152  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 3.12e-16
    assert frmtmb log_lik:157 unclassified TRUE
    value  frmtmb log_lik:157  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 4.44e-16, max rel 1.44e-16
    assert frmtmb log_lik:159 unclassified TRUE
    value  frmtmb log_lik:159  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 3.06e-16
    assert frmtmb log_lik:174 unclassified TRUE
    value  frmtmb log_lik:174  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 1.97e-16
    assert frmtmb log_lik:176 unclassified TRUE
    value  frmtmb log_lik:176  HOLDS: 30 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 1.62e-16
    assert frmtmb log_lik:230 unclassified TRUE
    assert frmtmb posterior_predict:111 unclassified TRUE
    assert frmtmb posterior_predict:113 unclassified TRUE
    assert frmtmb posterior_predict:118 unclassified TRUE
    assert frmtmb posterior_predict:120 unclassified TRUE
    assert frmtmb posterior_predict:133 unclassified TRUE
    assert frmtmb posterior_predict:135 unclassified TRUE

(the fifth column is `held`). Every one of the 13 holds, so each would
be a pass at consolidation; the 6 value checks hold against brms's own
`log_lik_*()` on the same prep at 3.1e-16 relative or better, and
log_lik:230 is itself a value assertion (brms's naive multivariate-t
conditional, `expect_equal(ll, ll_student)`) and holds. The 4 SAR
posterior_predict rows and the 4 SAR log_lik rows pass through the
path-graph adaptation; the 2 fcor posterior_predict rows and the 3
fcor log_lik rows pass as written.

The base arm (rellib-r13, `dev/loglik-log/base-rec-frmtmb-*.tsv`): the
7 log_lik rows do not hold ("'autocor_row_loglik' is not an exported
object"), the 6 posterior_predict rows hold, as the sarfcor replay
found at 0.75.0. So the port moves 13 rows from gap to pass, 7 of them
because of part A.

Not ported, and now portable: the two "log_lik for ARMA models runs
without errors" rows (log_lik:130, 133; survey: "brms's row is the
leave-one-out conditional density ... which frmtmb does not compute")
read brms's prep of per-draw covariance matrices, which an adapter
would have to map onto a time-series fit; outside the brief's 13.

## Open items and defects

- No new brms defect. brms-39 (`posterior_epred()` on an error SAR)
  stands; brms's `log_lik()` on the error SAR works.
- Gaps under `ar()`, `ma()`, `arma()`: frmtmb's conditional is the
  conditional of its own (by-level) model, which differs from brms's
  by-position one exactly as the two likelihoods differ; documented in
  `?frmtmb-autocor`, "The pointwise log-density".
- `fcor()` allocates two `n x n` matrices per draw (A2).

## Suites

Every suite the change can reach: core and all eight extensions, the
lane's core, frmtmb.sample and frmtmb.spatial in front of rellib-r13
(the `lib:` lines show `wt-loglik-lib` for those three and rellib-r13
for the other six).

One file per process, gated (`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`,
`dev/release/run-par.sh loglik-gated`, 8 at once; per-file logs
`dev/release/loglik-gated-files/`): "LOGLIK-GATED ran 452 of 452",
totals from the per-file RESULT lines:

    frmtmb           files 280 pass 22243 fail 0 err 0 skip 15 warn 0
    frmtmb.coupling  files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam       files 32 pass 1844 fail 0 err 0 skip 3 warn 0
    frmtmb.latent    files 12 pass 681 fail 0 err 0 skip 2 warn 0
    frmtmb.learn     files 15 pass 501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode       files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample    files 68 pass 3423 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial   files 8 pass 544 fail 0 err 0 skip 0 warn 0
    frmtmb.survival  files 15 pass 552 fail 0 err 0 skip 0 warn 0

The 29 skips are in other tiers' gates: test-drmtmb-agreement.R (13,
"set FRMTMB_DRMTMB_FIT_TESTS=true"), each package's test-scale.R,
test-fuzz.R and test-curve-scale.R. None is on a path this lane
changed.

One session per suite, three orders (`NOT_CRAN=true`,
`dev/release/run-sessions.sh loglik-sess 8`): "LOGLIK-SESS ran 27 of
27", fail 0 and err 0 in every session, 46 LEAK lines, every one tagged
"[with load of ...]", none untagged:

    SESSION frmtmb forward files=270/280 pass=19807 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=270/280 pass=19957 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=270/280 pass=20050 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample forward files=63/68 pass=3125 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3125 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3125 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=461 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=461 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=461 fail=0 err=0 skip=6 warn=0

frmtmb.sample's 5 files with no RESULT in a session are its five
`test-brms-suite-*.R` files, which skip as a whole outside the suite
tier; core's 10 are the same kind.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
"Rblas: 74 to OpenBLAS, 2 to Rblasref", "Rlapack: 616 to OpenBLAS, 10
to Rlapackref"; `OPENBLAS_NUM_THREADS=4`; gated; logs
`dev/loglik-log/ob-*.txt`):

    RESULT frmtmb test-autocor-loglik.R pass=47 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-fcor.R pass=57 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=93 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=146 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=109 fail=0 err=0 skip=0 warn=0

0.3.32 was not run: the new code has no branch that rounding decides
(no threshold on a computed value), and every assertion is a ratio to
the reference's own size with margins of 1e5 or more over the observed
differences. The OpenBLAS copy (`dev/ciharden-out/`) is deleted.

## R CMD check --as-cran

Once per changed package, built and checked in `dev/loglik-check/`
(`dev/loglik-rcheck.sh`; `R_LIBS` the lane library, rellib-r13, the
user library; pandoc and TinyTeX on PATH;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, `NOT_CRAN=true`; logs
`dev/loglik-log/check-*.txt`):

    frmtmb          Status: 1 NOTE
    frmtmb.sample   Status: 1 NOTE   (tests 344 s, OK)
    frmtmb.spatial  Status: 1 NOTE   (tests 29 s, OK)

Each NOTE is the expected "Skipping checking math rendering: package
'V8' unavailable". `dev/loglik-check/` and the bin-2 runner's staging
directory are deleted.

## Bumps and floors

- frmtmb: minor (a new export in the sampling API,
  `autocor_row_loglik()`, and a new argument of an exported seam,
  `frmtmb_register_autocor(jacobian =)`).
- frmtmb.sample: minor (`log_lik()` and the leave-one-out functions on
  residual covariance terms); floor the new core version, for
  `autocor_row_loglik()`.
- frmtmb.spatial: minor (`sar()` reaches `log_lik()`); floor the new
  core version, for `jacobian =` (an older core refuses the unknown
  argument when the package loads).

## Overlaps with the other lanes

`R/autocor-seam.R` (acprior touches `frmtmb_autocor_pars()` there; my
edits are the registration's `jacobian` argument, the map check and
`fcor_jacobian()`), `R/autocor.R` (doc sections only),
`R/sampling-api.R` (one alias and one export, each on its own line),
`NAMESPACE`, `NEWS.md`, `man/frmtmb-autocor.Rd`,
`man/frmtmb-sampling-api.Rd`, `man/frmtmb_register_autocor.Rd`,
frmtmb.sample (`R/loo.R`, `NEWS.md`, `man/sample-log_lik.Rd`,
`man/sample-loo.Rd`, `tests/testthat/test-loo.R`), frmtmb.spatial
(`R/sar.R`: the registration call gains `jacobian =` on its own line;
`NEWS.md`, `man/sar.Rd`, `tests/testthat/test-sar.R`), and
`dev/brmsport-bin2-adapters*.R` with the two regenerated bin-2 files.

## Fix round 1

Review `dev/reviews/2026-10-27-loglik.md`, verdict FIX. Numbered as the
review.

### 1. The Jacobian check at a generic point

Seen first on the build under review: the reviewer's
`dev/review-loglik-seam.R` (seed 81), run unedited on the lane library
before the fix (`dev/loglik-log/rv-seam-before.txt`):

    rvnosig sigma ~ x                              ACCEPTED, rel err vs dense = 1.047e-01
    rvnosig constant sigma 1.7                     ACCEPTED, rel err vs dense = 3.040e-01
    rvrows sigma ~ x                               ACCEPTED, rel err vs dense = 4.357e-01
    rvrows constant sigma 1.7                      ACCEPTED, rel err vs dense = 2.544e-16
    rvnorho rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 2.322e-01
    rvtrans rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 8.844e-02
    rvright rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 1.391e-15

The check is now `autocor_seam_check_jacobian()` (R/autocor-seam.R).
It requires `K (y - m) = whiten(y)$u` for two directions of `y - m`
(`sin(i)` and `cos(1.7 i)`), at two points. The first is the start
with sigma = 1. The second is a generic point: sigma 1.37, or
`exp(0.3 sin(i))` on row `i` under `sigma_varies = TRUE`, and every
parameter moved from the start by `0.4 k`, with k its position (0.04,
and then 0, when the map cannot be evaluated there). The same script
after the fix (`dev/loglik-log/rv-seam-after.txt`; each refusal reads
"The '<name>' term's map is wrong: jacobian() must give the n x n
matrix K with K (y - m) = whiten(y)$u, m the mean unwhiten(0), at every
sigma and parameter value; ..."):

    rvnosig sigma ~ x                              REFUSED
    rvnosig constant sigma 1.7                     REFUSED
    rvrows sigma ~ x                               REFUSED
    rvrows constant sigma 1.7                      REFUSED
    rvnorho rho = tanh(0.8)                        REFUSED
    rvtrans rho = tanh(0.8)                        REFUSED
    rvright rho = tanh(0.8)                        ACCEPTED, rel err vs dense = 1.391e-15

"rvrows constant sigma 1.7" is refused too: the term is registered
with `sigma_varies = TRUE`, and the slot is wrong for that whatever
sigma the model has.

Why one generic point is enough: a slot passes only if its error
matrix `K_slot - K` annihilates both directions at that point. The
slips that matter are functions of sigma, of the parameters, or of the
orientation of `K`:

- sigma: caught by a sigma that is not 1 and differs by row, which
  also tells rows from columns;
- the parameters: caught by moving every parameter off the start;
- the orientation: a transpose differs on two generic directions
  unless `K` is symmetric, and then `K' = K` and the values are right.

A slot that is right at both points and wrong elsewhere (a branch on
the parameter value) is still possible. No finite set of points
excludes it, and such a slot is not a function of the map. The help
section "The Jacobian" states the points. The help example now says
that `Linv / sigma` holds for one sigma, and that a sigma per row
divides the columns.

New test block in `tests/testthat/test-autocor-loglik.R`, "a jacobian
slot wrong away from the start point is refused". It refuses the four
slips by name: no division by sigma; rows instead of columns under
`sigma_varies = TRUE`; rho ignored; the transpose on an asymmetric
row-standardized W. As the absent cases, it accepts the right slot of
each map and checks it against the dense conditional; fcor()'s and
sar()'s own slots pass in every fit.

### 2. The 8 SAR rows: class "adapted-input"

`dev/brmsport-bin2-adapters-spec.R` gains `bin2_adapted_input`: the 8
ids (log_lik:150, 152, 157, 159; posterior_predict:111, 113, 118, 120),
each with the reason "pass on an adapted input: brms's prep uses Msar =
diag(10), and the port replaces it by the row-standardized path graph
over the same 10 rows (brms_adapt_sar_path()). frmtmb refuses brms's
literal input by design: I - rho I makes rho and the scale one
parameter, so sar() needs a weights matrix with a negative and a
positive eigenvalue. ...".

`dev/brmsport-bin2-ledger.R` gives a held row on that list verdict
"pass", class "adapted-input" and that reason. The reason goes into
the verdicts it writes, so the generated tier file carries it. The
script stops on an entry whose row does not hold, and the summary
states how many passes are adapted. For a dry run, the script now
takes `FRMTMB_BIN2_REC_DIR` and `FRMTMB_BIN2_LEDGER_OUT`.

Dry run: all six bin-2 files were recorded on the lane build
(`dev/loglik-bin2-record.sh`, records in `dev/loglik-bin2-log/`), and
the ledger was written into `dev/loglik-ledger-out/` (log
`ledger-run.txt`). `git status` shows no tracked verdict, ledger or
summary file changed. The 13 rows:

    id	outcome	class	kind	value_check
    log_lik:150	pass	adapted-input	shape	holds
    log_lik:152	pass	adapted-input	shape	holds
    log_lik:157	pass	adapted-input	shape	holds
    log_lik:159	pass	adapted-input	shape	holds
    log_lik:174	pass		shape	holds
    log_lik:176	pass		shape	holds
    log_lik:230	pass		value
    posterior_predict:111	pass	adapted-input	shape
    posterior_predict:113	pass	adapted-input	shape
    posterior_predict:118	pass	adapted-input	shape
    posterior_predict:120	pass	adapted-input	shape
    posterior_predict:133	pass		shape
    posterior_predict:135	pass		shape

The dry-run summary:

    | 2P | pass | - | - | 162 |
    | 2P | pass | adapted-input | - | 8 |
    8 of the bin-2 passes hold on an adapted input (class adapted-input: the port replaced brms's input, which frmtmb refuses by design); 220 are literal.
    | pass | 392 | 44 | 184 | 620 |

The tracked summary has `| pass | 392 | 44 | 171 | 607 |`. The gain is
13 rows: 5 literal passes (the fcor rows log_lik:174, 176, 230 and
posterior_predict:133, 135) and 8 adapted. Every other row of the dry
run matches the tracked verdicts; the ledger stops on any difference.

### 3. The dense fallback warns above 1000 rows

Measured (`dev/loglik-fallback-cost.R`, seed 7,
`dev/loglik-fallback-cost.txt`; sar()'s lag map registered without its
slot, over a path graph; minimum of 5 calls, peak memory of one call):

    n =   500  call 0.030 s  peak Mb 21  4000 draws 0.03 h  silent
    n =  1000  call 0.500 s  peak Mb 48  4000 draws 0.56 h  silent
    n =  1500  call 0.720 s  peak Mb 105  4000 draws 0.80 h  warns
    n =  2000  call 1.670 s  peak Mb 155  4000 draws 1.86 h  warns
    n =  3000  call 4.950 s  peak Mb 345  4000 draws 5.50 h  warns

The reviewer's map measured 4.3 s and 306 Mb at 2000; its whiten and
unwhiten are dense too.

I chose a warning, not a refusal: the values are exact, and 345 Mb per
draw at 3000 rows is not a memory hazard on a machine that fits such a
model. Above `autocor_loo_dense_max` = 1000 rows, core raises
`frmtmb_autocor_dense_fallback`, which names the term and `jacobian =`.
frmtmb.sample's `log_lik()` passes on the first one per call and
muffles the rest. Tests:

- core, "a term without a jacobian slot warns above 1000 rows": 1001
  rows warn, 1000 do not;
- frmtmb.sample, "log_lik() says once per call that a term takes the
  dense route": 3 draws at 1001 rows give exactly 1 warning.

### 4. A test that sar() uses its slot

test-sar.R, "sar() reads its precision off its sparse jacobian slot,
not the dense fallback": `sar_jacobian()` returns a `sparseMatrix`; a
1024-row sar() raises no fallback warning; and the same map registered
without the slot warns. Seen to fail on a mutant
(`dev/loglik-mutant-sarslot.R`, `dev/loglik-mutant-sarslot.txt`): sar()
registered again without the slot, with test-sar.R run against it:

    slot after the mutant registration: TRUE (TRUE = no slot)
    RESULT pass=91 fail=1 err=0
    FAILING:
      sar() reads its precision off its sparse jacobian slot, not the
               dense fallback

### 5. The gap difference in `?frmtmb.sample::log_lik`

One paragraph added to "Residual covariance terms". Under `ar()`,
`ma()` and `arma()`, a group that skips a time point inside its range
differs from brms: frmtmb reads lags by time level and brms by
position. On one such group in the reviewer's case the difference is
17 percent (`dev/review-loglik-log/brms.txt`, "rel = 1.721e-01"),
against 7.8e-16 on the groups without a gap. `cosy()` and `unstr()`
agree.

### 6. The cross-reference

sample-loo's link now points to `[sample-log_lik]`.

### 7. The loo_*() functions against the loo package

test-loglik-autocor.R, "loo(), waic() and the loo_*() functions on a
covariance block are the loo package's on the same log_lik() matrix".
For `ar(cov = TRUE)` and `cosy()`, with 40 hand-set draws each:

- `loo()` equals `loo::loo(ll, r_eff = relative_eff(exp(ll), one
  chain))`, in the estimates and the pointwise elpd;
- `waic()` equals `loo::waic(ll)`;
- `loo_epred()` equals `loo::E_loo(posterior_epred(ds), psis(-ll))`,
  and also the weighted mean `colSums(w * posterior_epred(ds))` written
  out;
- `loo_predict()` equals `E_loo()` of `posterior_predict()` drawn from
  the same random state.

There is no brms reference for these. brms's values would need brms's
draws through brms's `loo_*()`, which read the same `log_lik()` (A3
agrees with brms to 2.0e-15) through the same loo package.

### 8. test-loo.R

The `ar(cov = TRUE)` arm is now its own block, "log_lik() of an ar(cov
= TRUE) block is the conditional of each row given the rest of its
subject", at the end of the section "log_lik() is the conditional
density, exactly". The refusal block keeps the mixture arm.

### 9. theta_map

`autocor_seam_loo_parts()` applies `autocor_seam_theta_on()` once, for
`whiten()` and for the slot. The fallback gets the raw parameters,
because `autocor_seam_B()` applies the map itself.

### Not done

- Under `ar()`, `ma()` or `arma()`, `frm()` gives no message when a
  group skips a time level inside its range. `autocor_warn_gaps()`
  warns only on non-consecutive GLOBAL levels, so a ported brms script
  fits the by-level model where brms fits the by-position one. The
  defect predates this lane; filed here, not fixed.
- brms accepts `se()` and `mi()` beside `ar(cov = TRUE)`, where frmtmb
  refuses them (review, "What I checked"). This also predates the lane.

### The reruns

The first full gated run of fix round 1 (`dev/release/loglik-fr1-gated`,
"LOGLIK-FR1-GATED ran 452 of 452") found 2 failures in core's
test-message-uniqueness.R: "duplicated error templates: The '' term's
map  for". The new `autocor_seam_check_jacobian()` defined a second copy
of the map check's refusal. It now receives the map check's own `bad()`.
Core was reinstalled, and its 280 files and its three sessions were run
again. The extension results below come from the first run; the fix
changes no code that they reach other than core's message wiring.

One file per process, gated (`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`;
core from `dev/release/loglik-fr1b-gated`, "LOGLIK-FR1B-GATED ran 280
of 280"; the extensions from `loglik-fr1-gated`):

    frmtmb           files 280 pass 22257 fail 0 err 0 skip 15 warn 0
    frmtmb.coupling  files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam       files 32 pass 1844 fail 0 err 0 skip 3 warn 0
    frmtmb.latent    files 12 pass 681 fail 0 err 0 skip 2 warn 0
    frmtmb.learn     files 15 pass 501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode       files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample    files 68 pass 3434 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial   files 8 pass 547 fail 0 err 0 skip 0 warn 0
    frmtmb.survival  files 15 pass 552 fail 0 err 0 skip 0 warn 0

The skips are the other tiers' gates, as before.

One session per suite, three orders, `NOT_CRAN=true` (core:
`loglik-fr1b-sess`; sample and spatial: `loglik-fr1-sess`). No LEAK
line is untagged; 20 and 37 LEAK lines are all "[with load of ...]":

    SESSION frmtmb forward files=270/280 pass=19821 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=270/280 pass=19971 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=270/280 pass=20064 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.sample forward files=63/68 pass=3136 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3136 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3136 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=464 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=464 fail=0 err=0 skip=6 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=464 fail=0 err=0 skip=6 warn=0

OpenBLAS 0.3.26 with LAPACK ("Rblas: 74 to OpenBLAS, 2 to Rblasref",
"Rlapack: 616 to OpenBLAS, 10 to Rlapackref"; `OPENBLAS_NUM_THREADS=4`;
gated; logs `dev/loglik-log/ob-fr1-*.txt`; the copy is deleted again):

    RESULT frmtmb test-autocor-loglik.R pass=61 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-fcor.R pass=57 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=104 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loo.R pass=146 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=112 fail=0 err=0 skip=0 warn=0

The two bin-2 files on the final build (`dev/loglik-log/fr1-run-*.txt`,
records `fr1-rec-*.tsv`):

    RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=97 fail=0 err=0 skip=0 blocks=15
    RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=61 fail=0 err=0 skip=0 blocks=14

All 13 rows hold (held = TRUE), and the 6 value checks hold at 3.1e-16
relative or better, as in Part B.

R CMD check --as-cran (`dev/loglik-log/fr1-check-*.txt`): frmtmb,
frmtmb.sample and frmtmb.spatial each "Status: 1 NOTE". The NOTE is
"Skipping checking math rendering: package 'V8' unavailable".
`dev/loglik-check/` is deleted.

## Fix round 2

The reviewer's "Fix round 1 check" (`dev/reviews/2026-10-27-loglik.md`)
says SHIP, with three gaps. All three are fixed.

### 1. The generic point steps parameters both ways

The gap: the generic point moved the parameters upward only (0.4, then
0.04, then 0). The reviewer's `rvoneside` has rho = theta on (-1, 0],
starts at 0, its map stops for theta > 0, and its slot ignores rho. It
was accepted, and `log_lik()` was 3.465e-01 off at rho = -0.5
(`dev/review-loglik-log/fr1-seam2.txt`).

New test blocks in `tests/testthat/test-autocor-loglik.R`:

- "the generic point also steps a parameter down from its start": the
  reviewer's construction with the slot that ignores rho, which must be
  refused. The absent case is the same map with the right slot, which
  is accepted and equals the dense conditional at rho = -0.5.
- "a slot whose map cannot move off its start is checked there only,
  and says so": a map defined at its start alone, with a right slot.

Seen to fail on the build before this round
(`dev/loglik-log/fr2-before-autocor-loglik.txt`):

    RESULT frmtmb test-autocor-loglik.R pass=62 fail=2 err=1 skip=0 warn=0
    ── 1. Failure ('test-autocor-loglik.R:487:3'): the generic point also steps a pa
    Expected `ll_unfit(y ~ x + oneside(W), data = d, data2 = list(W = W))` to throw a error with class <frmtmb_error>.
    ── 2. Failure ('test-autocor-loglik.R:543:3'): a slot whose map cannot move off
    Expected `w` to inherit from "frmtmb_autocor_jacobian_unchecked".

The first failure shows the slip accepted; the absent case passed on
that build too.

The fix, in `autocor_seam_check_jacobian()`:

- For a term with parameters, the steps are 0.4, -0.4, 0.04 and -0.04
  times each parameter's position, in that order.
- If no moved point can be evaluated, the start is checked with the
  generic sigma, and a warning of class
  `frmtmb_autocor_jacobian_unchecked` names the term. It says that the
  slot was checked only at the start of its parameters, that a slot
  which ignores a parameter would pass, and that `log_lik()` reads its
  precision off that slot alone.
- It warns and does not refuse, because the slot may be right: refusing
  a correct model is worse than a check that says what it could not do.

The help section "The Jacobian" and the function's own doc state the
steps and the warning.

The reviewer's script `dev/review-loglik-fr1-seam2.R` (seed 101),
rerun unedited on the fixed build (`dev/loglik-log/fr2-rv-seam2.txt`):
every case got its expected outcome, and `rvoneside` is now refused.

    rvoneside (ignores rho) asymmetric W         expected REFUSED  got REFUSED  The 'rvoneside' term's map is wrong: jacobia
    rvok asymmetric W                            expected ACCEPTED got ACCEPTED rel err vs dense = 1.459e-16
    rvsym symmetric W                            expected ACCEPTED got ACCEPTED rel err vs dense = 2.927e-16
    rvtsym symmetric W                           expected ACCEPTED got ACCEPTED rel err vs dense = 2.927e-16

The other slips (`rvsign`, `rvcache`, `rvfirst`, `rvsym` and `rvtsym` on
an asymmetric W, `rvfrozen0`, `rvfrozen8` and `rvfrozen14`) stay
refused.

### 2. The gap paragraph

`?frmtmb.sample::log_lik`, "Residual covariance terms", now says what
was measured. The case has one such group, 22 rows and 3 draws. The
largest difference over the matrix is 17 percent of brms's largest
value, and the rows of the groups without a gap agree to 8e-16
(`dev/review-loglik-log/brms.txt`, "rel = 1.721e-01" and "rel on
groups with no interior gap = 7.844e-16").

### 3. The sar-uses-its-slot test reads the threshold

test-sar.R builds the lattice as `k = ceiling(sqrt(lim + 1))`, where
`lim` is `frmtmb:::autocor_loo_dense_max`. It asserts `n > lim` before
it asserts that no warning is raised, so a change to the threshold
cannot make the block pass on a case below it.

### Reruns

The changed files, one per process, in both tiers (logs
`dev/loglik-log/fr2-{plain,gated}-*.txt`). The plain tier's skips are
the brms-gated blocks.

    plain  RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    plain  RESULT frmtmb test-autocor-loglik.R pass=66 fail=0 err=0 skip=0 warn=0
    plain  RESULT frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
    plain  RESULT frmtmb.sample test-loglik-autocor.R pass=76 fail=0 err=0 skip=2 warn=0
    plain  RESULT frmtmb.spatial test-sar.R pass=93 fail=0 err=0 skip=4 warn=0
    gated  RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb test-autocor-loglik.R pass=66 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb.sample test-loglik-autocor.R pass=104 fail=0 err=0 skip=0 warn=0
    gated  RESULT frmtmb.spatial test-sar.R pass=113 fail=0 err=0 skip=0 warn=0

Core one-session runs (`NOT_CRAN=true`, `dev/release/loglik-fr2-sess`,
"LOGLIK-FR2-SESS ran 3 of 3"). There are 20 LEAK lines, all of them
tagged:

    SESSION frmtmb forward files=270/280 pass=19826 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=270/280 pass=19976 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=270/280 pass=20069 fail=0 err=0 skip=206 warn=0

OpenBLAS 0.3.26 with LAPACK ("Rblas: 74 to OpenBLAS", "Rlapack: 616 to
OpenBLAS"; gated; `dev/loglik-log/ob-fr2-*.txt`; the copy is deleted):

    RESULT frmtmb test-autocor-loglik.R pass=66 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-autocor-seam.R pass=72 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb test-message-uniqueness.R pass=6 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.sample test-loglik-autocor.R pass=104 fail=0 err=0 skip=0 warn=0
    RESULT frmtmb.spatial test-sar.R pass=113 fail=0 err=0 skip=0 warn=0

R CMD check --as-cran (`dev/loglik-log/fr2-check-*.txt`): frmtmb.sample
and frmtmb.spatial "Status: 1 NOTE"; frmtmb "Status: 1 NOTE". Each NOTE
is the V8 math-rendering note. `dev/loglik-check/` is deleted.
