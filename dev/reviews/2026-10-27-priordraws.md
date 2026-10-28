# Review: lane priordraws (2026-10-27)

Reviewer of `frmtmb-wt-priordraws` (branch `wt-priordraws` at 435e7ce6,
uncommitted). Lane build `wt-priordraws-lib`, base `rellib-r14`. Every
script below is `dev/review-priordraws-*.R` in this tree; logs are in
`dev/review-priordraws-log/` (gitignored). The scratch combination with
reslot was built in the session scratchpad and installed into
`C:/Users/adf44/source/r/wt-priordraws-rv-lib` (the only install this
review made); `wt-priordraws-rv-lib/floor/` holds a plain copy of the
lane's installed frmtmb.sample for the floor arm.

## Verdict: FIX

No blocker. The core seam, the `"yes"` draws, the evidence ratio, the
warmup storage and `update()` hold up under every check below. Three
items need the consolidation (findings 1 to 3). Findings 4 to 6 are
small fixes in this lane. Findings 7 to 11 are minor or nits.

## What was verified (passes)

- Tests, lane library, one file per process, `NOT_CRAN=true`,
  `FRMTMB_BRMS_FIT_TESTS=true` (`dev/review-priordraws-runall.sh`,
  logs `t-*.txt`, every `lib:` line the lane's):

  | file | pass | fail | error | skip | escaped warnings |
  | --- | --- | --- | --- | --- | --- |
  | core test-prior-only-objective.R | 8 | 0 | 0 | 0 | 0 |
  | core test-message-uniqueness.R | 6 | 0 | 0 | 0 | 0 |
  | core test-bracket-access.R | 36 | 0 | 0 | 0 | 0 |
  | sample test-prior-draws.R (run 1, 2, 3) | 57, 57, 57 | 0 | 0 | 0 | 0 |
  | sample test-warmup-draws.R | 33 | 0 | 0 | 0 | 0 |
  | sample test-evidence-ratio.R | 41 | 0 | 0 | 0 | 0 |
  | sample test-draws-methods.R | 174 | 0 | 0 | 0 | 0 |
  | sample test-draws-spellings.R | 97 | 0 | 0 | 0 | 0 |
  | sample test-loo.R | 146 | 0 | 0 | 0 | 0 |
  | sample test-generic-collision.R | 94 | 0 | 0 | 0 | 0 |
  | sample test-bracket-access.R | 1 | 0 | 0 | 0 | 0 |
  | sample test-message-uniqueness.R | 6 | 0 | 0 | 0 | 0 |
  | sample test-brms-suite-methods.R (gated) | 58 | 3 | 0 | 0 | 0 |

  The 3 failures are the ledger rows brmsfit-methods:794, 797 and 798
  reporting "now HOLDS but is recorded as 'cannot transfer'", as the
  lane says. Core test-prior-only-objective.R on rellib-r14: 3 pass,
  5 fail (seen to fail, confirmed).
- Flakiness: test-prior-draws.R three fresh processes, 57/57 each. Its
  KS and z checks run on fixed seeds (`with_sample_seed()` restores the
  session's RNG), so the prior-draw checks are deterministic; the only
  sampler-dependent bounds are `|mean| / mcse < 4` and
  `|sd - 1| < 0.1` on 2 x 1000 draws (finding 10).
- Names on 11 models the lane did not test, brms 2.23.0 fits against
  frm_sample(), `sample_prior = "yes"`, data seed 1, sampler seed 5, 2
  chains of 3000 with 1000 warmup (`review-priordraws-yes.R`, compare
  with `review-priordraws-compare.R`, log `compare.txt`): 10 of 11
  identical in name and order. T1 multivariate without rescor and
  resp-specific priors (`b_y Intercept_y sigma_y b_z_x Intercept_z
  sigma_z sd_g`), T2 `sigma ~ x` with an `Intercept_sigma` prior and a
  group-specific student_t sd, T3 two correlated terms with
  `lkj(1.5)` class-wide and `lkj(3)` on a 3-term group, T4 student_t
  `b` with `lb = 0` plus a coef prior, T5 nonlinear with an nlpar RE
  and `lb = 0` on an nlpar, T6 negbinomial `shape` gamma and logistic
  `b`, T7 zero-inflated `zi` beta, T9 `ma`, T10 `ar(cov = TRUE)`, T11
  `cosy`. T8 differs (finding 6).
- Distributions: two-sample KS, 4000 a side, 42 matched columns: p from
  0.022 to 0.996, 1 below 0.05 (T6 `sd_g`, 0.022; 2.1 expected by
  chance), 0 below 0.01. Truncation (student_t `b` at 0, sd, sigma,
  nlpar at 0, `ar`/`ma` on their declared interval), LKJ marginals for
  dimension 2 and 3, and the centered `prior_Intercept` (default
  `student_t(3, median, 2.5)`) all agree.
- `"only"` against brms, 4 chains of 3000 (1000 warmup), seed 9
  (`review-priordraws-only.R`): O1 weights + offset + `(1 | g)`, O2
  censoring, O3 `ar(cov = TRUE)`, O4 nonlinear, O6 truncation. 120
  z-scores (mean, 5, 50, 95 percent quantile, MCSE from `posterior`),
  of parameters, `r_g[1,Intercept]` and three prior predictive columns.
  |z| > 3: four. Three are O3 `sigma` (finding 2). The fourth is O1
  `b_Intercept` q95, z -3.46; against the exact prior
  (`review-priordraws-o1o6.R`) brms is the one off (z -3.29, frmtmb
  1.52). O6: brms's `posterior_predict()` returns 59 non-finite draws
  of 24000 for `trunc(lb = -3)`, frmtmb 0; on brms's finite draws the
  three yrep columns have KS p 0.572, 0.738, 0.571. All Rhat < 1.004.
- Refusals under `"only"`, brms `stancode()` against frm_sample() on 17
  models (`review-priordraws-onlyrefuse.R`, log `onlyrefuse.txt`;
  mixture redone with brms's families in `acdefault-*.txt`): same
  verdict and same parameter name in 14 (`bs`, `b_sigma`, `ar` with
  `cov = FALSE`; s(x), cumulative, cor, zip, bernoulli, student RE,
  mo, no intercept, `0 + Intercept`, offset only and mixture all
  sample on both). Different in `ar(cov = TRUE)` and `cosy` (finding
  1).
- Evidence ratio, same draws (`review-priordraws-er.R`, log `er.txt`):
  frmtmb's `hypothesis()` against brms's own `evidence_ratio()` and
  `find_vars()`, reading prior draws for every hypothesis in brms's
  RNG sequence, bitwise identical for `x = 0`, `x > 0`, `x < 0.3`,
  `x - z = 0`, `x + z = 0.5`, `z = x`, `class = "sd"` (`g__Intercept =
  0.5`, `g__Intercept > 1`) and `class = ""` (`sigma = 1`,
  `sigma^2 = 1`). The prior of `x - z` has sd 1.392 (two independent
  normal(0, 1) columns: 1.414). `Intercept = 0` gives NA with brms's
  warning; draws without prior draws give NA with no message.
- No leak of `prior_only` (`review-priordraws-probe.R`, log
  `probe.txt`, P1): `only$fit$frame$prior_only` is NULL, also after
  saveRDS/readRDS; `log_lik()` of an `"only"` fit equals an independent
  normal density (max relative difference 0); `frm()` afterwards and
  `build_objective()` of the fit's frame are unaffected.
- Prior columns do not change any method (P2): summary, fixef, ranef,
  posterior_predict, log_lik, loo, bayes_R2, conditional_effects
  identical between `"yes"` and `"no"` under the same seed;
  `mcmc_plot()` defaults leave `prior_` out; `posterior_summary()` and
  `rhat()` include the 4 `prior_` rows, as brms's do.
- Warmup isolation (P3): 20 methods identical with and without
  `x$warmup` (ranef, log_lik, posterior_predict, posterior_epred, loo,
  bayes_R2, pp_check data, conditional_effects, hypothesis,
  nuts_params, log_posterior, ndraws, variables, print, as_draws_df,
  fitted, predict, coef, waic, loo_R2). Thinning (P4): `thin = 3,
  warmup = 400, iter = 1000` gives 134 + 200 rows per chain, rstan's
  own `warmup2` and brms's shape; `as.mcmc(inc_warmup = TRUE)` refuses
  with brms's words; `nsamples(incl_warmup = TRUE)` 668.
- `update()` (P6): `update(no, sample_prior = "yes")` and the fit
  route's equivalent reproduce the direct call bitwise; `seed`,
  `chains`, `iter` and `sample_prior = "only"`/`"no"` rerun as asked;
  formula, `formula.`, `newdata` are refused with a clear message.
- Reslot overlap: lane + reslot's diff (dev/ excluded; three rejected
  hunks, all roxygen or NEWS, left out) installed into
  `wt-priordraws-rv-lib`: test-prior-draws.R 57/0, test-evidence-ratio.R
  41/0, test-warmup-draws.R 33/0, core test-prior-only-objective.R 8/0
  (logs `t-combo-*.txt`). Reslot's `(flat)` rows break nothing, and fix
  finding 1.
- NEWS: BREAKING items listed (evidence ratio, `inc_warmup`,
  `update()`), no lane names, no `dev/` paths, no em dashes; new code
  and tests within 80 columns.

## Findings

### 1. major: on the lane alone, default `ar`/`ma` (`cov = TRUE`) and `cosy` get prior draws and "only" samples them; brms has neither

Evidence: `review-priordraws-acdefault.R`, logs `acdefault-lane.txt`
and `acdefault-rv.txt`, `y ~ x + ar(time = t, gr = h, cov = TRUE)`,
`cosy(...)`, `ma(..., cov = TRUE)` with `normal(0, 1)` on b and the
default autocor prior. brms 2.23.0's `stancode(sample_prior = "yes")`
has `prior_b`, `prior_Intercept`, `prior_sigma` only (no statement for
the bounded flat `ar`, `ma`, `cosy`), and `"only"` is refused
"Error occurred for parameter 'ar'" (`'cosy'`). Lane build:

    yes: Intercept b ar sigma     (ar uniform on [-1, 1])
    yes: Intercept b cosy sigma
    yes: Intercept b ma sigma
    only: sampled 200             (all three)

Lane + reslot: `yes: Intercept b sigma` and the three refusals in
brms's words. The plan reads the effective prior, which on 0.76.0
holds the interval default as a proper uniform. Also seen in
`onlyrefuse.txt` (`ar_cov_default`, `cosy_default`).

Fix: consolidation merges reslot with or before this lane. Add a test
to test-prior-draws.R for the three classes (no `prior_ar`, `prior_ma`,
`prior_cosy` under default priors; `"only"` refused naming `ar`, `ma`,
`cosy`); it fails on this lane alone and passes on the combination,
which pins the dependency.

### 2. major (overlap with sigmainn): `sigma` under first-order `cov = TRUE` is not brms's, so `"only"` and `Evid.Ratio` of `sigma` disagree with brms

Evidence: O3, `y ~ x + ar(time = t, gr = h, cov = TRUE)`,
`exponential(1)` on sigma, 4 x 2000 draws a side: brms `sigma` mean
1.007, q95 2.975; frmtmb mean 1.232, q95 3.633; z -7.95 (mean), -4.48
(q50), -5.99 (q95). The prior predictive agrees (yrep z within 1.5)
and `ar[1]` agrees (z 0.04), so the prior is on the innovation sd and
the reported `sigma` column is the marginal sd. Under `"yes"`,
`prior_sigma` is drawn from the innovation-sd prior while the
posterior `sigma` column is the marginal sd, so
`hypothesis(x, "sigma = 1", class = "")` divides densities of two
different quantities. Pre-existing semantics (0.76.0's
`autocor_innov()`), which the sigmainn lane changes.

Fix: at consolidation, after sigmainn, run O3 again
(`review-priordraws-only.R brms 2` and `frm 2`) and add a test that
`"only"` sigma under `ar(cov = TRUE)` matches its prior. If sigmainn
does not land, document that `prior_sigma` is the innovation sd's.

### 3. major (consolidation): the Depends floor is the only guard against a stale core

Evidence: `review-priordraws-floor.R`, log `floor.txt`: the lane's
frmtmb.sample on rellib-r14's frmtmb 0.76.0, `y ~ x + (1 | g)`,
`normal(0, 1)` on b, `"only"`: `mean(b_x)` 0.448, `sd(b_x)` 0.138, no
warning (the posterior, not the prior). test-prior-draws.R on that
stack: 55 pass, 2 fail (lines 228 and 229), so the floor job catches
it. DESCRIPTION still reads `frmtmb (>= 0.76.0)`.

Fix: raise the floor to the core release that carries `prior_only`
and run test-prior-draws.R in the floor script. No runtime guard is
needed: once the floor is raised, `library()`/`loadNamespace()` refuse
an older core, so a stale core in a user library cannot load the new
extension.

### 4. minor: the `mi()` measurement term is dropped, where brms keeps it, and the `mi()` refusal gives a wrong reason

Evidence: `review-priordraws-mi.R`, log `mi.txt`. brms 2.23.0
`stancode(bf(y | mi(sdy) ~ x), sample_prior = "only")`:

    if (!prior_only) {
      target += normal_id_glm_lpdf(Yl | Xc, Intercept, b, sigma);
    }
    target += lprior;
    target += normal_lpdf(Y[Jme] | Yl[Jme], noise[Jme]);

The measurement term is in `model_prior`, outside the block, and the
model has no missing values, so brms's prior is proper. Core's
`build_objective()` guards that term with `!prior_only`
(R/objective.R, `mm_$se`), and `?frmtmb-sampling-api` and NEWS say the
seam is brms's block. frm_sample() refuses this model with "without
the likelihood the missing values have no density at all", which is
false here.

Fix: keep the `mm_$se` term under `prior_only` (drop the
`&& !prior_only`), and refuse only an `mi()` model that has missing
values without a measurement density. Add a core test for the term.

### 5. minor: no `me()` model can reach `"only"`, but the docs say it does

Evidence: O5 (`only-frm-3.txt`): `bf(y ~ me(xo, sdx))` with proper
`b`, `Intercept`, `sigma` priors is refused "Error occurred for
parameter 'meanme_mexo'", and `set_prior(..., class = "meanme")` is
refused by core ("frmtmb ... has no prior slot for it"), so the user
cannot act on the refusal. brms samples it once `meanme` has a prior.
`?frm_sample` (Prior draws) and frmtmb.sample's NEWS describe
`"only"` sampling `me()` values; core keeps their densities for a path
no user can reach.

Fix: state the limit in `?frm_sample` and NEWS, and make the refusal
say that frmtmb has no prior slot for `meanme`, `sdme`. A prior slot
is out of this lane's scope; file it.

### 6. minor: coefficient-specific threshold priors get no draws; brms's names are now measured

Evidence: T8, `cumulative`, `normal(0, 2)` on `Intercept`, `coef = 1`:
brms `b Intercept__1 Intercept__2 Intercept__3 Intercept__4`; lane
`b`. The lane left these out as "not measured".

Fix: draw `prior_Intercept__k` in brms's position (after `b`), each
from its own row or the class prior, or list the gap in NEWS with
brms's names.

### 7. minor: the warmup is stored twice

Evidence: P7 in `probe.txt`, `y ~ x + (1 | g)` with 400 groups, 2
chains of 1000: `save_warmup = TRUE` draws 3.1 MB, warmup 3.1 MB,
stanfit 6.6 MB, total 14.6 MB; `save_warmup = FALSE` total 8.4 MB. The
stanfit already holds the warmup; `x$warmup` is a second copy in draws
coordinates, and `frm_sample()` extracts the whole array a second time
to get it. brms reads the warmup from the stanfit on demand.

Fix (optional): convert the warmup from `x$stanfit` when
`inc_warmup = TRUE` is asked for, or document the cost in `?frm_sample`.

### 8. minor: `draws_raw_array(inc_warmup = TRUE)` drops the warmup silently on a mismatch

`R/draws-brms.R`: when `x$warmup` has other columns or a row count not
divisible by the chains, the post-warmup array is returned with no
message. No current path triggers it, but a later method that edits
`x$draws` would make `inc_warmup = TRUE` quietly mean FALSE.

Fix: refuse in that case; return the draws alone only when
`x$warmup` is NULL.

### 9. nit: `update(ds, 5)` says "a new formula or new data"

`probe.txt` P6. An unnamed non-formula argument binds to `formula.`.
Fix: refuse unnamed non-formula arguments with the "named arguments"
message.

### 10. nit: two fixed relative tolerances in tests

test-prior-draws.R `abs(sd(bx) - 1) < 0.1` and test-evidence-ratio.R
`abs(kde / dnorm(0) - 1) < 0.1`. Both are about 4 standard errors at
their draw counts and will not flake, but neither is stated against a
measured error. Fix: compare with `posterior::mcse_sd()` (and a stated
kernel bias bound), or say in a comment why 0.1 is safe.

### 11. nit: `evidence-ratio.R` and `test-evidence-ratio.R` are LF in a CRLF working copy

Git warns on `git diff`. Harmless; normalize at consolidation.

## What the consolidation must run

- test-prior-draws.R on the merge of priordraws, reslot, sigmainn and
  mointerop, plus the new test of finding 1.
- O3 of `review-priordraws-only.R` against brms after sigmainn
  (finding 2).
- mointerop's two-value `mo()`: the unit test `prior_draws_make()`
  with a simplex of one weight gives 1; also fit `y ~ mo(m2)` with a
  two-level `m2` under `"yes"` on the merge and check
  `prior_simo_mom21[1]` is 1 in every draw, as brms's.
- The floor script with the raised floor (finding 3).
- Ledger: 794, 797, 798 to pass; 107 and 802 keep "cannot transfer"
  with the lane's new reasons.

## Short check of fix round 1 (2026-10-28)

### Verdict: SHIP

Findings 4 to 11 are fixed and verified. Section 13 of the findings
lists what the consolidation must do, completely (one addition below,
optional). Two nits, neither blocks.

### The library is the tree

`dev/review-priordraws-libmatch.R`: every function in the tree's `R/`,
sourced and compared by `deparse()` (sources removed) with the
installed namespace in `wt-priordraws-lib`: frmtmb 2104 functions, 0
differ, 0 missing; frmtmb.sample 306, 0 differ, 0 missing.

### Tests, lane library, one file per process

`dev/review-priordraws-runall2.sh`, logs `fr1-*.txt`, every `lib:` line
the lane's, 0 escaped warnings:

| file | pass | fail | error | skip |
| --- | --- | --- | --- | --- |
| sample test-prior-draws.R, runs 1, 2, 3 | 65, 65, 65 | 0 | 0 | 0 |
| sample test-evidence-ratio.R, runs 1, 2, 3 | 41, 41, 41 | 0 | 0 | 0 |
| sample test-warmup-draws.R | 39 | 0 | 0 | 0 |
| sample test-draws-methods.R | 174 | 0 | 0 | 0 |
| sample test-message-uniqueness.R | 6 | 0 | 0 | 0 |
| core test-prior-only-objective.R | 10 | 0 | 0 | 0 |
| core test-me.R | 70 | 0 | 0 | 0 |
| core test-message-uniqueness.R | 6 | 0 | 0 | 0 |

### Item by item

- 4, `mi(sdy)`: verified. `dev/review-priordraws-fr1.R` (brms and frm
  sides, data seed 1, `y | mi(sdy) ~ x`, `sdy = 0.4`, 4 chains of 3000
  with 1000 warmup, seed 9) and `dev/review-priordraws-fr1cmp.R`:
  `b_Intercept`, `b_x`, `sigma` and three prior predictive columns, 24
  z-scores (mean, 5, 50, 95 percent quantiles), max |z| 2.39. The
  first latent value has mean 0.204 against `y[1]` = 0.209 and sd
  0.395 against 0.4. With 2 missing `y`, brms's `stancode()` accepts
  the model (it samples them flat), and frm_sample() refuses with
  "2 missing values of y have no density without the likelihood",
  which is true; the same for a plain `mi()`. Without missing values
  there is no refusal.
- 5, `me()`: verified. brms 2.23.0 `default_prior()` has `meanme` flat
  and `sdme` flat with lb 0; `"only"` with the defaults is refused
  "Error occurred for parameter 'meanme_1'"; with `normal(0, 2)` and
  `exponential(1)` it is accepted, and `"yes"` then draws
  `prior_meanme_1` and `prior_sdme_1`. frm_sample's refusal gives
  brms's words for `meanme_1` and then says frmtmb has no prior slot
  for `meanme` and `sdme`, which is true and tells the user why no
  prior helps. `"yes"` gives `Intercept bsp sigma`, brms's names under
  its flat defaults. `?frm_sample` says the same. Making the two
  classes a follow-up is right: the change touches `set_prior()`'s
  class lists, the prior table, the resolver, the brmsprior
  translation and frm_sample's defaults, all in `R/priors.R`, which
  reslot rewrites in this round.
- 6, ordinal thresholds: verified. T8 (`dev/review-priordraws-yes.R
  frm 3` again, compare with `dev/review-priordraws-compare.R`):
  names and order identical to brms, `b Intercept__1 Intercept__2
  Intercept__3 Intercept__4`. Two-sample KS against brms: p 0.241,
  0.200, 0.033, 0.962, 0.355. Against the exact priors (one-sample,
  frmtmb): `__1` normal(0, 2) p 0.286; `__2` to `__4` student_t(3, 0,
  2.5) p 0.061, 0.916, 0.981 (brms's own: 0.811; 0.197, 0.835, 0.072).
  The 0.033 is chance.
- 7, warmup stored once: verified. `dev/review-priordraws-warmcmp.R`
  run by `dev/review-priordraws-warmcmp.cmd` pinned to logical core 0,
  old arm the previous build (`wt-priordraws-rv-lib/floor`, which
  stored the warmup) and new arm the lane library, data seed 1,
  sampler seed 13, 2 chains of 300 with 200 warmup; compared by
  `dev/review-priordraws-warmcmp2.R`. Six models the lane did not use:
  correlated random effects with prior draws, `mo()` + student (simplex
  and `nu`), `ar(cov = TRUE)` (the bounded chart), GEV (fit route, the
  xi chart), `thin = 2`, and two imputations. In all six the draws,
  the whole `as_draws_array(inc_warmup = TRUE)` array and
  `rhat(inc_warmup = TRUE)` are `identical()` old to new, and
  `as.mcmc(inc_warmup = TRUE)` and `as.data.frame(inc_warmup = TRUE)`
  have the same shape (and the `thin` refusal). The new objects store
  no warmup except GEV and the imputations, which convert at sampling
  time, as the lane says. `object.size()` falls by 18 to 70 KB on
  these small models (re_cor 1723 to 1653 KB); GEV is unchanged.
- 8, mismatch: verified. Draws with a column dropped give "inc_warmup =
  TRUE: the saved warmup does not fit these draws (7 variables and 300
  rows against 6 variables over 2 chains) ..."; `inc_warmup = FALSE`
  on the same object works.
- 9, `update(ds, 5)`: "the first unnamed argument is the new formula,
  and this one is of class numeric. Name the arguments ...";
  `update(ds, y ~ x)` keeps the new-model refusal.
- 10, tolerances: both are now stated against a measured error
  (`posterior::mcse_sd()`, and the kernel estimate's exact expectation
  and first-order standard error, deterministic under `set.seed(3)`).
  At 4 errors the false-failure rate is about 6e-5 per check; three
  runs of each file passed.
- 11, line endings: `R/evidence-ratio.R` and `test-evidence-ratio.R`
  are CRLF now.

### Section 13

Complete and precise against this review's findings 1 to 3 and its
consolidation list: the reslot merge order and the test that pins it,
the O3 rerun after sigmainn with its fallback, the floor and the floor
script, the four-lane test-prior-draws.R run, the two-value `mo()`
check, and the ledger rows. Optional addition: run
`dev/review-priordraws-acdefault.R rv` on the merge (it prints the
names and refusals of finding 1 directly).

### New findings

12. nit: `extensions/frmtmb.sample/R/sample.R` is LF in the working
    copy now (`git ls-files --eol`: `i/lf w/lf`; the other touched
    files are `w/crlf`), and git warns on `git diff`. The section 14
    note "No other file converted" is not quite right. The diff is not
    inflated (277 insertions, 24 deletions). Normalize at
    consolidation.
13. Outside this lane, to file: `frm_sample(bf(ye ~ x), family =
    gen_extreme_value(), ...)` from a formula fails on rellib-r14 too
    with "frm_linpred() needs a fitted model" (`dev/review-priordraws-gevbase.R`);
    the fit route works.
