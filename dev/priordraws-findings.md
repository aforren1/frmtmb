# Lane priordraws: sample_prior, prior_draws(), stored warmup

Round opened 2026-10-27, base 435e7ce6 (frmtmb 0.76.0, frmtmb.sample
0.23.0). Worktree `frmtmb-wt-priordraws`, library
`C:/Users/adf44/source/r/wt-priordraws-lib`, base build `rellib-r14`.
Logs are in `dev/priordraws-log/` (gitignored); every number below
names its script, and the counts are pasted from
`dev/priordraws-summarize.R`'s output (section 9).

STATUS: DONE, fix round 1 included (sections 13 and 14).

## 1. What brms 2.23.0 does (run, not recalled)

Scripts: `dev/priordraws-brms-stancode.R` (stancode only, log
`brms-stancode.txt`), `dev/priordraws-brms-names.R` (brms fits, sets 1
to 3, logs `brms-names-[123].txt`), `dev/priordraws-brms-only.R` and
`dev/priordraws-brms-only2.R` (logs `brms-only.txt`, `brms-only2.txt`),
`dev/priordraws-brms-examples.R` (the stored example fits, log
`brms-examples.txt`). brms source read from `dev/brms-suite/brms/R`
(`stan_rngprior()`, `rename_prior()`, `reorder_pars()`,
`prior_draws.brmsfit()`, `ignore_prior()`, `evidence_ratio()`,
`density_ratio()`, `stan_prior()`'s `prior_only` check,
`.as_draws_list()`, `as.mcmc.brmsfit()`, `update.brmsfit()`).

- `sample_prior = "yes"`: one `prior_` draw per `lprior` statement of
  the Stan program, in generated quantities, with the statement's `_rng`
  and rejection sampling for a bounded parameter. A class-wide prior on
  a vector is one draw (`prior_b`); a vector with any
  coefficient-specific prior has one statement, and one draw, per
  coefficient (`prior_b_x`, `prior_b_z`; case B). A flat prior has no
  statement and no draw (default `b`, `bs`, `bsp`; cases A, D, M). A
  separate Stan vector gets its own draw from the class prior
  (`prior_bsp` beside `prior_b`; D2). Names after `rename_pars()`:
  `prior_sd_g`, `prior_sd_g__Intercept` (per coefficient),
  `prior_cor_g` (one element `[1, 2]` of an `lkj_corr_rng()` draw),
  `prior_simo_mom1[k]` (a Dirichlet draw), `prior_sds_sx`,
  `prior_Intercept_sigma`, `prior_b_sigma`, `prior_sigma_y`,
  `prior_b_y`, `prior_rescor`, `prior_theta[k]` (mixture weights,
  Dirichlet), `prior_b_a` (nlpar), `prior_nu` (truncated at 1),
  `prior_hu` (on [0, 1]), `prior_phi`, `prior_ar` (unbounded under
  `cov = FALSE`). A GEV's `xi` is brms's transformed `tmp_xi`, which
  `stan_rngprior()` excludes. `gp()`: `prior_sdgp` and
  `prior_lscale__1_gpx` (case GP).
- Order: `reorder_pars()` puts the class `prior` after the parameters,
  and within it the parameters whose Stan name (before renaming) ends in
  `_Intercept` first. Only the location's `prior_Intercept` moves
  (`prior_sd_g__Intercept` was `prior_sd_1__1`; case B), and the rest
  keep the Stan program's order: per response, nlpars, then each dpar
  in the family's order with b, Intercept, simo, bsp, bcs, bs, sds and
  ac, then mixture weights, rescor, group-level sd and cor.
- `sample_prior = "only"`: the model block's likelihood sits inside
  `if (!prior_only)`; `lprior`, the group-level `std_normal(z)` and
  `me()`'s measurement density `normal_lpdf(Xn | Xme, noise)` stay
  outside. `mi()`'s response densities are inside, so the missing
  values are flat. No `prior_` draws in generated quantities, so
  `prior_draws()` of an "only" fit is a 0 x 0 data frame. A parameter
  with no prior statement is refused at `stancode()`: "Sampling from
  priors is not possible as some parameters have no proper priors.
  Error occurred for parameter 'b'." ('sd_1' for a flat sd, 'ar' for a
  flat ar). A vector with a coefficient-specific prior on one
  coefficient and none on the others is NOT refused (it has a prior
  statement), and the flat coefficients are sampled from nothing.
- `prior_draws(x, variable)`: each name takes the prior column whose
  name is the longest match at its start (`b` for `b_x`), in random
  order (`sample(ndraws)`); `b_Intercept` of a model with a centered
  intercept warns "Sampling from the prior of an overall intercept is
  not possible by default..." and is left out, but only after the
  prefix match found something (`max(matches) == -1 ||
  ignore_prior()`).
- `hypothesis()`: `prior_draws(x, variable = parsH)` for every
  hypothesis, used when it has a column for each name; a point
  hypothesis's `Evid.Ratio` is `density_ratio(post, prior, point = 0)`,
  NA with no message without prior draws. brms gives ratios for
  `sigma = 1` and `sd_g__Intercept = 1` (case "yes fit").
- Warmup: `save_warmup` is rstan's (default TRUE), passed through
  `brm(...)`. `as_draws*(inc_warmup = TRUE)` adds the saved warmup;
  with `save_warmup = FALSE`, `warmup2` is 0 and `inc_warmup = TRUE`
  returns the post-warmup draws (100 of 100 rows; "nowarm" cases).
  `as.mcmc(inc_warmup = TRUE)` refuses with `combine_chains = TRUE` and
  with `thin >= 2`. `rhat()` and `neff_ratio()` pass `...` to
  `as_draws_array()`. `nuts_params()` and `log_posterior()` pass `...`
  to bayesplot's stanfit methods, which take `inc_warmup`.
- `update(fit, sample_prior = "yes")` is supported (recompiles).
- The stored example fits 1 and 2 have `sample_prior = "yes"` with
  priors `normal(0,2)` on b, `normal(0,3)` on b of sigma and
  `cauchy(0,2)` on sd (lb 0); and `normal(2,2)` and `normal(0,3)` on
  the b of nlpars a and b.

## 2. What changed

Core (`R/objective.R`, `R/sampling-api.R`, `R/confint.R`):

- `build_objective()` reads `frame$prior_only`: with TRUE it leaves out
  every response density (and the `mi()` measurement term) and keeps
  the group-level and `me()` densities, brms's `if (!prior_only)`. It
  is the documented objective seam (`?frmtmb-sampling-api`), useful on
  its own to anyone sampling a model's prior. Only a copy of a frame
  carries it.
- `plot()` of a hypothesis overlays the density of the prior draws
  (second of `colors`) unless `ignore_prior = TRUE`, as brms does.

frmtmb.sample:

- `R/prior-draws.R` (new): `check_sample_prior()`, the plan of brms's
  prior statements (`prior_plan()`, read off `validate_prior()`'s table
  of the effective prior, `brms_coef_table()`, the blocks and the
  draws' names), direct draws (`prior_draws_make()`, `rprior_real()`:
  inverse-CDF truncation, LKJ marginal `2 Beta(a, a) - 1`,
  `a = eta - 1 + M / 2`, Dirichlet by gammas), the `"only"` refusals,
  `prior_draws()`.
- `frm_sample(sample_prior = )`; the warmup kept as `x$warmup`; the
  call's arguments as `x$sample_args`; `x$sample_prior`,
  `x$prior_info`.
- `draws_raw_array()`, `draws_as_array()`, `as.mcmc()`, `rhat()`,
  `neff_ratio()`: `inc_warmup`.
- `R/evidence-ratio.R` rewritten: brms's `evidence_ratio()` from prior
  draws; the analytic denominator and its refusals are gone.
- `update.frmtmb_draws()` samples again.
- Generic owners: `prior_draws` (brms). brms's deprecated `prior_samples()`
  is not defined here: brms's calls `UseMethod("prior_draws")`, which
  reaches the method registered in brms's table. A generic of our own
  could not follow it: the owner bindings need a bare `UseMethod()`
  generic with a `brmsfit` method in brms's table, and brms has no
  `prior_samples.brmsfit` (test-generic-collision.R's BRMSCONTROL
  counted 21 of 22 when it was defined).
- The ported suite's draws fixtures 1 and 2 sample with brms's example
  priors and `sample_prior = "yes"`, as brms's stored example fits were
  made (tests/testthat/helper-brms-suite-draws.R).

## 3. The two routes to prior draws, measured

Script `dev/priordraws-route.R`, log `route.txt` (model B, cauchy(0, 1)
on one sd among its priors; seed 21, 4 chains of 2000, 4000 draws):

    only run: 4000 draws, 3.5 s elapsed, 0 divergent
    posterior run without prior draws 2.8 s, with 2.5 s
    direct draws of 4000 rows: 0.0090 s each (20 repetitions)
      b_x              run: ess_bulk/n 1.31 ess_tail/n 0.68 KS(thin 10) p=0.238 | direct: ess/n 1 by construction, KS p=0.052
      b_z              run: ess_bulk/n 1.54 ess_tail/n 0.72 KS(thin 10) p=0.886 | direct: ess/n 1 by construction, KS p=0.747
      sigma            run: ess_bulk/n 1.14 ess_tail/n 0.60 KS(thin 10) p=0.713 | direct: ess/n 1 by construction, KS p=0.523
      sd_g__x          run: ess_bulk/n 1.31 ess_tail/n 0.52 KS(thin 10) p=0.439 | direct: ess/n 1 by construction, KS p=0.334
      sd_g__Intercept  run: ess_bulk/n 1.31 ess_tail/n 0.46 KS(thin 10) p=0.637 | direct: ess/n 1 by construction, KS p=0.555
      sd_g__x upper quantiles, exact: 63.7 636.6 | run: 105.2 806 | direct: 68.6 324.9

Decision, on these numbers and on brms's semantics:

- `"yes"`: direct draws. brms's `prior_` columns are independent
  generated-quantities draws, not a second chain, so direct draws are
  the same object. They cost 0.009 s for 4000 rows against 3.5 s for a
  run; a run's tail effective size is 0.46 to 0.72 per draw, the direct
  draws' is 1 by construction (the class-level `prior_b` and
  `prior_cor_g` statements are not even parameters a run could sample).
  The posterior run itself took 2.8 s without and 2.5 s with the prior
  draws (one run each, the difference is noise).
- `"only"`: a sampler run on core's `prior_only` objective. brms samples
  every parameter of the model from its prior, the group-level
  coefficients included; direct draws could only do that with a
  generator per covariance structure (car, spde, gp, ar1, ...), and the
  objective already has every one of their densities. Section 4 has
  its agreement with brms and with the exact prior.

## 4. Evidence against brms

Script `dev/priordraws-evidence.R` (parts 1 to 4, logs
`evidence-[1234].txt`) and `dev/priordraws-evidence2.R` (parts 5 and 6,
logs `evidence-[56].txt`); data seed 1 (part 4: seed 4), sampler seeds
in each script. Names: `dev/priordraws-names-frm.R`, log
`names-frm.txt`. Counts generated by `dev/priordraws-summarize.R`
(section 9).

- Names: 17 of 17 models give brms's `prior_draws()` names in brms's
  order (cases A, B, B2, B3, C2, H, D, D2, I, O, R, P, N, M, E, K, S of
  section 1), and 2 of 2 more fitted on both sides by
  `dev/priordraws-names2.R` (log `names2.txt`): smooths in two
  predictors (`Intercept sds_sx sds_sz Intercept_sigma sds_sigma_sx`)
  and `ar` priors on both responses of a multivariate model. The first
  run of names2 found a defect, since fixed: `s(x)` in the location and
  `s(x)` in sigma share a label, and the plan keyed smooth terms by label
  alone, so the second got no column. The ported fixture 1 gives the
  14 names brms's suite expects (`dev/priordraws-fixture.R`, log
  `fixture.txt`), and an `ar()` model `Intercept ar sigma`, brms's
  stancode's.
- Part 1 (frmtmb's draws against the exact truncated prior, one-sample
  KS, 2000 draws a column): 39 columns of 9 models, p from 0.125 to
  0.986, none below 0.05.
- Part 2 (against brms's own draws, two-sample KS, 4000 a side): 14
  columns of models B and D2, p from 0.033 to 0.962, 2 below 0.05
  (B `cor_g` and D2 `sigma`, both 0.033). Part 5 reads brms's own
  draws against the exact prior: 14 columns at two seeds, p from 0.120
  to 0.913; frmtmb's same two columns against the exact prior in part 1
  have p 0.347 and 0.230. Two of 14 below 0.05 has probability about
  0.15 under the null, so this is read as chance.
- Part 3 (`"only"`, `y ~ x + (1 | g)`, `normal(0, 1)` on b, 4 chains
  of 4000 a side): 24 z-scores (mean, 5, 50 and 95 percent quantiles of
  `b_Intercept`, `b_x`, `sigma`, `sd_g__Intercept`,
  `r_g[1,Intercept]` and of `yrep[, 1]`), max |z| 3.00 (the 95 percent
  quantile of `sd_g__Intercept`: brms 7.749, frmtmb 8.376, exact
  7.956). frmtmb's Rhat max 1.0026, 0 divergences.
- Part 6 (`"only"`, two more seeds, each side against the exact prior
  quantiles with `posterior::mcse_quantile()`): frmtmb 18 z-scores, max
  |z| 1.29; brms 18 z-scores, max |z| 2.43. So part 3's 3.00 is the two
  sides' noise added, not a disagreement of frmtmb with the prior.
- Part 4 (Evid.Ratio, `y ~ x + (1 | g)` with `normal(0, 1)` on b and
  `student_t(3, 0, 1)` on sd and sigma, 4 chains of 6000 a side, error
  of each side from 8 blocks of each density estimate):

        x = 0                  brms ER   2.6281 (se 0.0731) | frmtmb ER   2.7762 (se 0.1032) | z   1.17
        x = 0.2                brms ER   8.1709 (se 0.2109) | frmtmb ER   7.8632 (se 0.1662) | z  -1.15
        sigma = 1              brms ER   1.9518 (se 0.0709) | frmtmb ER   1.9589 (se 0.0517) | z   0.08
        sd_g__Intercept = 0.5  brms ER   2.8601 (se 0.0606) | frmtmb ER   2.9952 (se 0.0557) | z   1.64

- Identity: draws made without `sample_prior` are bitwise those of
  rellib-r14 under fixed seeds (`dev/priordraws-identity.R`, run by
  `dev/priordraws-identity.cmd` pinned to logical core 0 with `start
  /affinity 1`, 3 runs per arm): 4 models (fit route with correlated
  group terms, formula route with `mo()` and student, `ar(cov = TRUE)`
  through the residual-correlation chart, `laplace = TRUE`), one md5
  each over all 6 runs.


## 5. Seen to fail on rellib-r14

- `tests/testthat/test-prior-only-objective.R` (core): lane 8 pass;
  rellib-r14 5 failures of 8 (lines 24, 30, 38, 55, 64: the objective
  still holds the response density). Logs `t-core-prioronly-*.txt`.
- The three new frmtmb.sample files on rellib-r14 fail in the weak form
  only, because the tests pass `sample_prior`, which 0.23.0 forwards to
  rstan ("passing unknown arguments: sample_prior"): test-prior-draws
  0 pass, 1 fail, 14 errors; test-warmup-draws 4 errors;
  test-evidence-ratio 7 pass, 8 errors. Logs `t-test-*-base.txt`.
- The behavioral form, built with arguments 0.23.0 takes
  (`dev/priordraws-seenfail.R`, logs `seenfail-base.txt` and
  `seenfail-lane.txt`, seed 7, 2 chains of 1200 with 400 warmup):

  | case | brms | rellib-r14 | lane |
  | --- | --- | --- | --- |
  | A `dim(as_draws_array(inc_warmup = TRUE))` | 1200 x 2 | refused | 1200 x 2 x 17 |
  | B `as.mcmc(inc_warmup = TRUE)` rows | 1200 | refused | 1200 |
  | C `Evid.Ratio` of `x = 0`, no prior draws | NA | 4.024 | NA |
  | D `update(ds, sample_prior = "yes")` | prior_ columns | refused | `prior_Intercept prior_b prior_sigma prior_sd_g` |
  | E `brms::prior_draws(ds)`, no prior draws | 0 x 0 | 1600 x 0 | 0 x 0 |
  | F `Evid.Ratio` of `sigma = 1`, prior draws | a number | NA + warning | 12.906 |

- The ported rows 794, 797 and 798 fail on rellib-r14 by their ledger
  verdict ("cannot transfer": no prior draws) and hold on the lane
  (section 7).

## 6. Tests

New files:

- core `tests/testthat/test-prior-only-objective.R` (8 assertions): the
  `prior_only` objective is the group-level density exactly, reads no
  coefficient, leaves the fit's frame alone, and differs from the full
  objective by the response log-likelihood; me()'s two densities stay.
- frmtmb.sample `test-prior-draws.R` (57 assertions): `sample_prior`
  values; brms's names and order (model B, mo()+student, smooths in two
  predictors); each column against its exact prior by one-sample KS
  (with a control that the half-t check rejects a full t); truncation
  of a bounded `b`; Dirichlet and `nu`; `prior_draws(variable = )`'s
  longest match, random order, `pars` warning and the intercept rule;
  no prior draws without "yes"; prior draws leave the posterior draws
  bitwise as they were; the seed fixes them; `update(sample_prior =
  "yes")` reproduces the direct call bitwise; the `"only"` refusals
  (brms's words; a partly flat vector; `mi()`) and their absent case (a
  model with every class proper samples); `"only"` gives `b_x` its
  `normal(0, 1)` prior and a prior predictive of the right shape; a
  simplex of one weight is drawn as 1; a flat `sd` has no prior column
  and `"only"` refuses it naming `sd_1`.
- frmtmb.sample `test-warmup-draws.R` (33 assertions): `inc_warmup` on
  every accessor, the warmup's rows are the sampler's, brms's
  `as.mcmc()` refusal, `save_warmup = FALSE`, and every other method
  identical to the same object stripped of its warmup.

Changed files: `test-evidence-ratio.R` (rewritten for brms's ratio from
prior draws; adds the plot overlay), `test-draws-methods.R` and
`test-draws-spellings.R` (`inc_warmup` is answered now),
`test-loo.R` (`update()` of fake draws), `test-generic-collision.R`
(`prior_draws` in the owner table), `helper-brms-suite-draws.R`
(fixtures 1 and 2 sample with brms's example priors and prior draws).

Runs, all on the lane's library (`lib:` lines checked), one file per
process through `dev/release/run-par.sh` with `FRMTMB_RELLIB`, tallied by
`dev/priordraws-tally.R` (logs `tally-*.txt`):

- ungated, all nine suites, 458 files (`priordraws-ungated.log`):
  28510 pass, 5 fail, 0 error, 250 skip. The 5 failures were this lane's
  and are fixed: `fam$dpars` read with `$` (core and sample
  test-bracket-access), the old `inc_warmup` refusal
  (test-draws-spellings), and one refusal text written twice
  (test-message-uniqueness, 2). Every skip is a gate.
- frmtmb.sample's 72 files and core's two touched files again after the
  fixes (`priordraws-sample2.log`): 3405 pass, 0 fail, 0 error, 14 skip.
- the eight files touched by the last fix (`priordraws-final.log`):
  361 pass, 0 fail.
- gated (`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS`,
  `FRMTMB_FUZZ`), the 49 files that skipped ungated
  (`priordraws-gated.log`): 3888 pass, 0 fail, 0 error, 15 skip (the
  `test-scale.R` files and `test-curve-scale.R`, gated by their own
  switch); and the 15 files that skip whole ungated, the ported brms
  suite (`priordraws-gated2.log`): 631 pass, 3 fail, 0 error, 0 skip.
  The 3 failures are the ledger rows that now hold (section 7).
- OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
  lapack`): `test-prior-only-objective.R` 8 pass (log
  `t-ob-prioronly.txt`), `test-prior-draws.R` 53 pass and
  `test-warmup-draws.R` 33 pass (logs `t-ob-test-*.txt`).
- One session per suite, three orders (`dev/release/run-sessions.sh
  priordraws-sessions 6`, ungated, `FRMTMB_RELLIB` = the lane's library;
  log `dev/release/priordraws-sessions.log`): 27 of 27 sessions, every
  one fail 0 and error 0. frmtmb forward 272/282 files, 20074 pass,
  204 skip; reverse 20224 pass; preload 20317 pass, 206 skip (the
  file-level skips run no test, so the files count 272). coupling 551,
  eam 1838, latent 681, learn 430, ode 549, sample 3362 (67/72 files,
  14 skip), survival 552, spatial 478 pass in each order. Every LEAK
  line is tagged "[with load of ...]": extensions a test loads, and
  `rxTempDir` and `GDAL_DATA` from rxode2's and sf's loads.

## 7. Ported ledger

`extensions/frmtmb.sample/tests/testthat/test-brms-suite-methods.R`,
gated, on the lane (log
`dev/release/priordraws-gated2-files/frmtmb.sample--test-brms-suite-methods.txt`):
58 pass, 3 fail, every failure a row that now holds against its
recorded verdict:

- `brmsfit-methods:794` (prior_draws names): cannot transfer -> pass.
- `brmsfit-methods:797` (prior_draws(fit1, "b_Trt1") dimnames):
  cannot transfer -> pass.
- `brmsfit-methods:798` (`sort(prs1$b) == sort(prs2$b_Trt)`): cannot
  transfer -> pass.

Rows that keep their verdict, with a reason the consolidation should
rewrite:

- `brmsfit-methods:107` stays "cannot transfer", but no longer for the
  recorded reason (no warmup): it reads `fit1$fit@sim$iter`, which in
  brms is the stanfit and here the frmtmb fit. `dim(as.mcmc(fit1,
  inc_warmup = TRUE)[[1]])` is 75 x the variable count here, brms's
  shape (fixture: 1 chain of 75 with 50 warmup).
- `brmsfit-methods:802` stays "cannot transfer": `prior_draws()` on
  `fit1$fit`, brms's stanfit, which here is the frmtmb fit and carries
  no draws.
- `brmsfit-methods:595` (`nsamples(incl_warmup = TRUE)`) passed before
  and passes.

`dev/brms-api-diff.md` (c): the `inc_warmup` item is marked closed.

## 8. Not done, and defects found

Not done, with the reason:

- `gp()` prior draws. brms 2.23.0 names them `prior_sdgp` (whatever the
  term) and `prior_lscale__1_gpx` (case GP), which `rename_prior()`
  produces by accident of its regular expressions; a second gp term
  would give brms duplicate names. Not emulated; no column is drawn.
- (Done in fix round 1: coefficient-specific ordinal threshold priors
  are drawn as brms's `prior_Intercept__k`.)
- Two group-level terms on one grouping factor with different ids
  (`(1 | g) + (0 + x | g)`): brms makes one statement per id and names
  both `prior_sd_g`, a duplicate; this package makes one statement per
  grouping factor.
- `"only"` refuses a vector of coefficients some of which are flat,
  where brms samples the flat ones from nothing (section 1). Deliberate:
  the result has no proper prior.
- `"only"` refuses missing values of an `mi()` term, where brms samples
  them flat. Deliberate, for the same reason. (Fix round 1: a `mi(sdy)`
  response with no missing value is sampled, as brms does.)
- `"only"` refuses `me()` models: brms needs proper `meanme` and `sdme`
  priors, and frmtmb has no prior slot for those classes (fix round 1).
- `prior_samples()` is not defined in this package (section 2).
- `hypothesis(scope = "ranef" | "coef")` gives no prior draws, as in
  brms (its `hypothesis_coef()` reads a data frame with no `prior_`
  columns).
- No runtime guard against an older core: frmtmb.sample with frmtmb
  0.76.0 would ignore `prior_only` and sample the POSTERIOR under
  `"only"`. The DESCRIPTION floor has to rise to the core that carries
  the seam (section 10); the floor is the guard.
- (Fix round 1: the warmup is no longer stored twice; it stays in the
  stanfit and is converted on demand.)

Defects found, not fixed (outside this lane):

- `frm_sample()`'s default-prior disclosure prints `b (flat), as brms
  leaves slopes` when the call itself gives class `"b"` a prior
  (`dev/priordraws-disclosure.R` on rellib-r14, log `disclosure.txt`).
  The prior applies; the message is wrong. Pre-existing.
- `extensions/frmtmb.sample/R/methods-draws.R`: roxygen warns "@noRd
  must not be followed by any text" at the `draws_pp_structured()`
  block (already in test-backlog.md, 0.76.0 consolidation). Left alone
  to keep this lane's diff to what it needs.
- The ported helper `brms_shim_update()` drops `sample_prior` with the
  other sampling arguments (`brms_sampling_args`), so a ported
  `update(fit, sample_prior = )` row would not reach it. No such row is
  ported today.

## 9. Generated counts

`dev/priordraws-summarize.R`, output `dev/priordraws-log/summary.txt`:

    ## names (names-frm.txt)
    SUMMARY same 17 of 17
    ## part 1: 39 columns of 9 models, KS p from 0.125 to 0.986; below 0.05: 0
    ## part 2: 14 columns, two-sample KS p from 0.033 to 0.962; below 0.05: 2
    ## part 5: brms's own draws vs exact, 14 columns, p from 0.120 to 0.913
    ## part 3: 24 z-scores, max |z| 3.00
    ## part 4: 4 hypotheses, z from -1.15 to 1.64
    ## part 6: brms 18 z-scores against the exact prior, max |z| 2.43
    ## part 6: frmtmb 18 z-scores against the exact prior, max |z| 1.29
    ## identity re_fit: 6 runs, distinct md5 1 (base 1, lane 1)
    ## identity mo_student: 6 runs, distinct md5 1 (base 1, lane 1)
    ## identity ar_cov: 6 runs, distinct md5 1 (base 1, lane 1)
    ## identity laplace: 6 runs, distinct md5 1 (base 1, lane 1)

Which build: parts 1 to 6 and the identity runs ran before four later
changes to `R/prior-draws.R` (no `prior_samples()`, `fam[["dpars"]]`,
smooth terms keyed by predictor, `ar` rows deduplicated per response)
and the shared refusal helper. None of them is on the path of a draws
object made without `sample_prior`, and the names were run again on
the final build: 17 of 17 and 2 of 2. Part 1 run again on the final
build (log `evidence-1-final.txt`) prints its 39 KS lines identical to
the first run's.

## 10. Bump, floor, overlaps

- frmtmb: MINOR. A new documented field of the objective seam
  (`frame$prior_only`) and the prior overlay of `plot()` of a
  hypothesis. Nothing else in core changes behavior.
- frmtmb.sample: MINOR, with BREAKING items (the evidence ratio,
  `inc_warmup`, `update()`). Floor: the core release that carries
  `prior_only`. Without it, `sample_prior = "only"` would sample the
  posterior (section 8), so the floor is not optional.
- Files other lanes of this round are likely to touch too: core
  `NEWS.md`, `R/confint.R` (hypothesis), `R/objective.R`
  (build_objective), frmtmb.sample `R/sample.R` (frm_sample's body and
  roxygen), `R/methods-draws.R` (hypothesis, as.mcmc, rhat, update),
  `R/draws-brms.R`, `NAMESPACE`, `NEWS.md`,
  `tests/testthat/test-generic-collision.R` (the owner lists; the
  new entry, `prior_draws`, is on its own line in each).

## 11. R CMD check --as-cran

`dev/priordraws-check.ps1` (modeled on `dev/release/run-check.ps1`:
built with vignettes, `NOT_CRAN=true`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
`R_LIBS` = the lane's library, rellib-r14, the user library), in
`dev/priordraws-check/<package>`, logs `dev/priordraws-check/*.log`:

- frmtmb.sample: `Status: 1 NOTE`, the HTML manual's "Skipping checking
  math rendering: package 'V8' unavailable", which is this machine's.
  Tests 387 s, OK.
- frmtmb: `Status: 1 NOTE`, the same V8 NOTE. Tests 22 min, OK. The
  first attempt ran while the frmtmb.sample check wrote into
  `dev/priordraws-check/`, and `R CMD build` of core, which copies the
  whole tree before it applies `.Rbuildignore`, failed on a locked file
  there ("problem copying ... Permission denied"); it was run again
  alone.
- frmtmb.sample was checked again after the last two edits of
  `R/prior-draws.R` (section 12): `Status: 1 NOTE` (V8), tests 395 s,
  OK.

## 12. Overlap with the reslot and mointerop lanes

Read from their worktrees (read-only) at the coordinator's request.

- reslot writes brms's (flat) as an effective prior row with
  `prior = ""`: a `"flat"` density (the map's Jacobian), an LKJ made
  flat on the Cholesky factor (`dist$flat = TRUE`), the autocor
  defaults on their interval (a uniform, `flat_interval`), and a
  bounds-only user row that replaces a default. `validate_prior()`'s
  table writes no density for a row with `prior = ""`
  (`fill_prior_table()` sets `prior` only from a non-empty string and
  vectorizes only from non-flat rows), so each such slot reads `(flat)`
  there. The plan reads the table, so `"yes"` draws nothing for it, as
  brms draws nothing for a flat parameter, and `"only"` refuses with
  brms's message naming brms's Stan parameter (`b`, `sd_1`, `ar`,
  `L_1` for a correlation, `sigma`, ...), as brms does
  (`dev/priordraws-brms-only.R` measured `b`, `ar` and `sd_1`; `L_1`
  is read from `stan_prior()`, not run). The backstop on internal
  parameters now skips a `"flat"` or flat-LKJ entry, so it cannot count
  one as a proper density (edit after the coordinator's note). Measured
  on this tree with a bounds-only `set_prior("", class = "sd", lb = 0)`:
  no `prior_sd_g` column, and `"only"` refused with "Error occurred for
  parameter 'sd_1'" (test-prior-draws.R, 57 pass). Not run on reslot's
  tree; the consolidation should run test-prior-draws.R on the merge.
- mointerop: a two-value `mo()` has `simo_mox1[1]` fixed at 1; brms
  declares `simplex[1]`, and `dirichlet_rng()` of length one is 1, so
  `prior_simo_mox1[1]` is 1 in every draw. The plan counts the draws'
  `simo_mox1[` columns (1), and `prior_draws_make()` returns 1 for a
  simplex of one weight (edit after the note; a unit test in
  test-prior-draws.R). A one-value `mo()` is refused by mointerop's
  `frm_sample()`, since brms cannot sample it either; with no `simo`
  column the plan has no statement for it. This lane does not touch
  `posterior_predict()` of a `set_rescor(TRUE)` fit.

## 13. For consolidation (review 2026-10-27, findings 1 to 3)

Recorded here, not worked around in this lane
(`dev/reviews/2026-10-27-priordraws.md`):

1. Default `ar`/`ma` under `cov = TRUE` and `cosy` get prior draws on
   this lane alone (the 0.76.0 effective prior holds their interval
   default as a proper uniform), and `"only"` samples them; brms has no
   statement for them and refuses `"only"` naming `ar`, `ma` or `cosy`.
   reslot's `(flat)` rows fix it (the reviewer ran the combination:
   `yes: Intercept b sigma` and brms's three refusals). Consolidation
   must merge reslot with or before this lane and ADD a test to
   test-prior-draws.R: under default priors, no `prior_ar`,
   `prior_ma`, `prior_cosy`, and `"only"` refused naming `ar`, `ma`,
   `cosy`. It fails on this lane alone and passes on the merge.
2. `sigma` under first-order `cov = TRUE` is frmtmb's marginal sd while
   its prior (and `prior_sigma`) is on brms's innovation sd, so
   `"only"` `sigma` and the `Evid.Ratio` of `sigma` disagree with brms
   (O3: z -7.95 on the mean). The sigmainn lane changes the semantics.
   Consolidation must rerun O3 (`dev/review-priordraws-only.R brms 2`
   and `frm 2`) after sigmainn and add a test that `"only"` `sigma`
   under `ar(cov = TRUE)` matches its prior; if sigmainn does not land,
   document that `prior_sigma` is the innovation sd's.
3. The Depends floor: frmtmb.sample must require the core release that
   carries `prior_only`; on 0.76.0 core `"only"` samples the posterior
   (review: `mean(b_x)` 0.448). Raise the floor in DESCRIPTION and run
   test-prior-draws.R in the floor script (`dev/rel073-floor.sh`); it
   fails 2 assertions on the stale core, so the floor job catches it.

Also from the review's list: run test-prior-draws.R on the merge of
priordraws, reslot, sigmainn and mointerop; fit `y ~ mo(m2)` with a
two-level `m2` under `"yes"` on the merge and check
`prior_simo_mom21[1]` is 1 in every draw; ledger 794, 797, 798 to
pass, 107 and 802 keep "cannot transfer" with the new reasons
(section 7).

## 14. Fix round 1

Review 2026-10-27, findings 4 to 11.

- 4, `mi(sdy)`. brms 2.23.0 keeps `normal_lpdf(Y[Jme] | Yl[Jme],
  noise[Jme])` outside `if (!prior_only)` (`dev/review-priordraws-mi.R`,
  log `dev/review-priordraws-log/mi.txt`). Core `build_objective()`
  now keeps the `mi(sdy)` measurement term under `prior_only`
  (`R/objective.R`), and `?frmtmb-sampling-api` and NEWS say so.
  `frm_sample("only")` refuses only a model with missing values of an
  `mi()` term (latent rows without a measurement), naming the count
  and the variable ("5 missing values of xm have no density without
  the likelihood ..."). Tests: core test-prior-only-objective.R (the
  term equals `-sum(dnorm(y, Yl, sdy))` and reads no coefficient);
  test-prior-draws.R (`y | mi(sdy) ~ x` samples under `"only"`, `b_x`
  is its prior and a latent value is `normal(y[1], 0.3)`, each within
  4 Monte Carlo errors; the missing-value refusal's wording).
- 5, `me()`. brms 2.23.0 (`dev/priordraws-brms-me.R`, log
  `brms-me.txt`): `default_prior()` has `meanme` and `sdme` flat
  (`sdme` lb 0); `"yes"` draws `prior_meanme_1` and `prior_sdme_1`
  only when they have proper priors (brms does not rename them), after
  the predictor statements; `"only"` with the defaults is refused
  "Error occurred for parameter 'meanme_1'"; with `normal(0, 2)` on
  meanme and `exponential(1)` on sdme it samples, keeping
  `normal_lpdf(Xn | Xme, noise)` and `std_normal(zme)` outside the
  block. I did NOT add the classes to core's resolver: the change is
  not self-contained (set_prior's class lists, the prior table rows,
  `resolve_priorlist()`, the brmsprior translation and its refusal
  text, frm_sample()'s defaults, MAP fits) and sits in `R/priors.R`,
  which reslot rewrites. Instead `"only"` on a model with `me()` terms
  is refused with brms's words for 'meanme_1' followed by "frmtmb has
  no prior slot yet for brms's classes "meanme" and "sdme" ...".
  `"yes"` draws no meanme or sdme, brms's answer under its (flat)
  defaults. `?frm_sample` (Prior draws) and NEWS say this; the core
  seam still keeps both `me()` densities (true, and what a later slot
  needs). To file: brms's `meanme` and `sdme` as prior classes in core.
- 6, ordinal thresholds. brms 2.23.0 with `normal(0, 2)` on
  `Intercept` `coef = 1` (`brms-me.txt`): one statement per threshold,
  `prior_b`, `prior_Intercept__1` (normal(0, 2)), `__2` and `__3` (the
  class's student_t(3, 0, 2.5)); a flat class leaves only `__1`. The
  plan now does that (K from the draws' `b_Intercept[k]`), and these
  columns do not move to the front, since brms's reorder reads the
  generated-quantity name `prior_Intercept__1`. Test: names
  `b Intercept__1 Intercept__2 Intercept__3` on brms's data, and the
  first two against their exact priors (KS).
- 7, warmup stored once. `frm_sample()` leaves the warmup in the
  stanfit alone, with a small recipe (`warmup_recipe`: laplace, mo,
  ncp blocks, chart bounds) and the prior draws of the warmup rows
  (`warmup_prior`); `draws_warmup()` converts it on demand through the
  same `sample_finish()` the draws are made with. A GEV xi chart reads
  each draw through the objective's tape, which the draws do not keep,
  so there (and for several imputations, each with its own fit) the
  warmup is converted at sampling time as before. Measured on the
  reviewer's 400-group case (`dev/priordraws-size.R`, data seed 1,
  sampler seed 2, 2 chains of 1000; logs `size-before.txt`,
  `size-after.txt`): `object.size()` 14.6 MB before, 11.5 MB after
  (stanfit 6.6 MB, draws 3.1 MB, warmup copy 3.1 MB to none);
  `saveRDS()` 18.8 MB to 16.0 MB; `save_warmup = FALSE` 8.4 MB both.
  The `inc_warmup = TRUE` array (1000 x 2 x 406) is `identical()`
  before and after; it takes 0.19 s against 0.07 s. Tests
  (test-warmup-draws.R): no stored copy; the recipe applied to the
  post-warmup iterations gives the draws back bitwise; every earlier
  `inc_warmup` assertion holds.
- 8, mismatch. `inc_warmup = TRUE` with a warmup whose variables or
  rows do not fit the draws is refused ("the saved warmup does not fit
  these draws (...)"). Tests construct both: draws with a column
  dropped, a warmup with a row dropped.
- 9, `update(ds, 5)` says "the first unnamed argument is the new
  formula, and this one is of class numeric. Name the arguments of
  frm_sample() to change ...". Test added.
- 10, tolerances. test-prior-draws.R: `|sd(b_x) - 1|` within 4
  `posterior::mcse_sd()`. test-evidence-ratio.R: the KDE at 5 of
  normal(5, 1) draws against its exact expectation
  `dnorm(0, 0, sqrt(1 + h^2))`, within 4 first-order standard errors
  `sqrt(f R(K) / (n h))`.
- 11, line endings. The working copy is CRLF (`core.autocrlf=true`, no
  `.gitattributes`; `git ls-files --eol` shows `w/crlf` on the other
  touched files). `R/evidence-ratio.R` and `test-evidence-ratio.R` were
  LF and are CRLF now. No other file converted.

Seen to fail on the lane's build before this round (the reviewer's
read-only copy `wt-priordraws-rv-lib/floor/frmtmb.sample`, on the
current core; `dev/priordraws-runtest-prev.R`, logs `t-prev-*.txt`):
test-prior-draws.R 56 pass, 1 fail, 5 errors (the mi refusal's reason,
`mi(sdy)` refused, the me message, the ordinal names (`b` only), the
update message); test-warmup-draws.R 27 pass, 3 fail, 3 errors (the
stored copy, the recipe, the mismatch returned silently).

Tests after the round (lane library, one file per process,
`dev/priordraws-fr1.jobs`, log `dev/release/priordraws-fr1.log`): 19
files, 1196 pass, 0 fail, 0 error, 2 skip (test-loo's gated two).
Core: test-prior-only-objective 10, test-bracket-access 36,
test-message-uniqueness 6, test-me 70, test-mi-bsp-names 11,
test-brms-parity-defects 172, test-edgecases 64, test-case-studies 31.
frmtmb.sample: test-prior-draws 65, test-warmup-draws 39,
test-evidence-ratio 41, test-draws-methods 174, test-draws-spellings
97, test-loo 138, test-generic-collision 94, test-message-uniqueness 6,
test-bracket-access 1, test-parallel-chains 4, test-sample-direct 137.

One session per suite, three orders, core and frmtmb.sample
(`dev/release/run-sessions.sh priordraws-fr1-sessions 6 frmtmb
frmtmb.sample`, log `dev/release/priordraws-fr1-sessions.log`): 6 of 6
sessions, fail 0, error 0 in each. frmtmb forward 20076, reverse
20226, preload 20319 pass (272/282 files); frmtmb.sample 3380 pass in
each order (67/72 files, 14 skip). Every LEAK line is tagged "[with
load of ...]".

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`, logs `t-ob2-*.txt`; the R copy deleted afterwards):
test-prior-only-objective 10, test-prior-draws 65, test-evidence-ratio
41, test-warmup-draws 39, all pass.

R CMD check --as-cran, one package at a time, after every edit of the
round (`dev/priordraws-check.ps1`, logs `dev/priordraws-check/*.log`):
frmtmb.sample `Status: 1 NOTE` (tests 358 s); frmtmb `Status: 1 NOTE`
(tests 18 min). The NOTE is the HTML manual's "Skipping checking math
rendering: package 'V8' unavailable" in both.
