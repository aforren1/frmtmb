# Lane aterms2: brms's subset(), index(), mi(x, idx = ), rate() and cat()

Worktree `C:\Users\adf44\source\r\frmtmb-wt-aterms2`, branch
`wt-aterms2`, base `1f40800d` (frmtmb 0.65.0, frmtmb.sample 0.13.0).
Private library `C:/Users/adf44/source/r/wt-aterms2-lib`; "before" is
`C:/Users/adf44/source/r/rellib-r3`. brms 2.23.0 and rstan 2.32.7 from
the user library. Every script below is in `dev/`, and its log beside it
as `dev/aterms2-log-*.txt`.

## 1. What brms 2.23.0 does, read from its source

Read from the installed namespace (`dev/aterms2-brms-read*.R`, logs
`dev/aterms2-log-brms-read*.txt`) and from `stancode()` / `standata()`
(`dev/aterms2-brms-probe.R`, `dev/aterms2-brms-probe2.R`), and from
fits (`dev/aterms2-brms-fitprobe.R`, seed 11, and
`dev/aterms2-brms-fitprobe2.R`, seed 12).

- `cat(x)`: `resp_cat()` stores `thres = "x - 1"` with no group, and
  `data_response()` warns "Addition argument 'cat' is deprecated. Use
  'thres' instead." `cat()` beside `thres()` warns and `cat()` wins
  (`terms_ad()`: `x$thres <- x$cat`). Valid for the ordinal families.
- `rate(denom)`: families poisson, negbinomial, negbinomial2, geometric
  (their `$ad`). `denom` must be numeric and positive ("Rate
  denomiators should be positive."). Stan code: `poisson_log_lpmf(Y |
  mu + log_denom)` under the log link without cens/trunc,
  `poisson_lpmf(Y | mu .* denom)` otherwise; negbinomial
  `neg_binomial_2_log_lpmf(Y | mu + log_denom, shape .* denom)`;
  geometric the same with shape `1 .* denom`; negbinomial2
  `inv(sigma) .* denom`. R side: `multiply_dpar_rate_denom()` in
  `posterior_epred_*`, `posterior_predict_*` and `log_lik_*`, which
  multiplies mu and the shape. `posterior_linpred()` has no `log(denom)`.
  `fitted()` equals `fitted(dpar = "mu") * denom` (ratio 1 to 1 in the
  fit probe). Newdata without `denom` is refused ("can neither be found
  in 'data' nor in 'data2'"); `conditional_effects()` holds it at its
  mean.
- `subset(s)`: every family. `subset_data()` keeps the rows where
  `as.logical(s)` is TRUE; refuses NA ("Subset variables may not contain
  NAs."), a wrong length, and an empty result. `na_omit()` drops a row
  for an NA only if some response that uses the variable includes the
  row. Group levels come from the WHOLE data (`N_1 = 5` with a level
  present only outside the subset). `me()` with subset is refused.
  After the fit, `fitted()`, `predict()`, `log_lik()` and `loo()` refuse
  a missing or multiple `resp` ("Argument 'resp' must be a single
  variable name for models using addition argument 'subset'."), return
  `N_resp` rows with it, and on `newdata` return the rows where the
  subset is TRUE; `nobs()` is the whole data (40).
- `index(id)` and `mi(x, idx = ref)`: `idxl_<resp>_<x>_<k> =
  match(ref, index(x) on x's rows)`, and the predictor reads
  `Yl_x[idxl[n]]`. Refusals: "Response 'x' needs to have an 'index'
  addition term", "Could not match all indices in response 'x'.",
  "Index of response 'x' contains duplicated values.", "mi() terms of
  subsetted variables require the 'idx' argument to be specified." and
  "mi() terms in subsetted formulas require the 'idx' argument to be
  specified." The coefficient is `rename("mi(x,idx=g1)")`, `mixidxEQg1`.
  On newdata the indices are matched within newdata.

## 2. What changed

The design follows brms: one model frame holds every row, which keeps
brms's NA rule, and each response then takes its own rows of it before
anything is built, so its designs, smooth bases, grouping levels and
addition-term values come from its own rows.

- `R/subset.R` (new): the subset and index evaluation with brms's
  checks, the NA rule (`subset_na_rows()`), the per-response rows,
  the refusals (`subset_check_spec()`), `mi_idx_rows()` (brms's
  `idxl`), `mi_idx_newdata()`, and `subset_resp_check()` /
  `subset_newdata()` for the post-fit methods.
- `R/parse.R` (addition-term hunks only): `core_aterms` gains `rate`,
  `subset`, `index`, `cat`; `row_aterms`; the `cat()` branch and the
  cat/thres double-count refusal; `mi_pred_args()` and its use in the
  two predictor `mi()` branches.
- `R/frame.R`: `rhs_resp` beside `rhs_comb` in the frame formula
  loop; the NA branch uses `subset_na_rows()` under subset; the
  univariate filter; `mf`/`n` set to a response's rows at the top of
  the per-response loop and of Phase 1 and restored after each; the
  `cens_y2` rows; the `rate()` positivity check; `idx_expr`/`idxl` and
  brms's label on `mi_info`; Z's row count from its component; the
  `subset_rows` element of the frame; the |ID| level message.
- `R/families.R`: `rate_mu()`, `rate_mean()`, `rate_shape()`; poisson,
  negbinomial and geometric read `rate` in lpdf, lcdf, mean, variance,
  deviance, truncated mean and simulator; `check_accepted_aterms()`
  ignores `row_aterms`.
- `R/objective.R`: `idxl` in `lp_eta_fixed()`; the zero-column eta and
  the `cs()` matrix take the design's rows, not the frame's.
- `R/predict.R`: `idxl` in `patch_mo_cols()`; `eval_dpars()` and
  `lp_eta_design()` take the design's rows; `mean_is_mu()` accepts
  `rate_mean()` and `has_rate()` routes a rate response through the
  mean; `aterms_for_newdata()` needs `rate` and skips the row terms;
  the resp check and newdata filter in `frm_linpred()` and `fitted()`;
  `mi_idx_newdata()` in `pred_design()`.
- `R/predict-brms.R`: `predict()`'s resp check and newdata filter; the
  exposure follows newdata in `predict_simulate()`.
- `R/conditional-effects.R`: `row_aterms` skipped; a subset variable is
  held at TRUE on the grid, as brms's `prepare_conditions()`; `has_rate`
  in the two mean-display tests.
- `R/sandwich.R`: `vcov_cluster()` refuses a multivariate subset model.
- `R/compat.R`: features `rate()`, `subset()`, `index()` with rules;
  the row terms are exempt from the derived allow-list refusals.
- `R/sampling-api.R`: `subset_resp_check()` and `subset_newdata()`
  exported and documented there.
- `R/fit.R`: `?frm` sections "Exposures, rate()" and "Response subsets,
  subset() and index()", and the `cat()` paragraph.
- frmtmb.sample: `R/loo.R` (`log_lik()`'s resp check, one column per
  row of the asked response), `R/methods-draws.R` (resp checks in
  `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`,
  `predictive_error()`, newdata filter, and the exposure in
  `posterior_predict(newdata =)`).
- Tests: `tests/testthat/test-subset-rate.R` (new, ungated),
  rows 24, 25 and 25b of `test-brms-likelihood.R` (gated),
  `extensions/frmtmb.sample/tests/testthat/test-subset-rate-draws.R`
  (new), and brms's names `denom` and `idxl_*` in
  `brms_standata_view()` of `tests/testthat/helper-brms-suite.R`, with
  the same hunk in the sample package's generated copy.
- Docs: `NEWS.md`, the sample package's `NEWS.md`,
  `vignettes/brms-migration.Rmd`, `man/frmtmb-sampling-api.Rd`,
  `man/frm.Rd` (roxygenised).

## 3. Validation

### 3.1 brms's own compiled program, check A and C

`dev/aterms2-lpcheck.R` (seeds 24, 25, 26), log
`dev/aterms2-log-lpcheck.txt`; the same shapes are rows 24, 25 and 25b
of `test-brms-likelihood.R`. `measured_const` is brms's `log_prob`
minus frmtmb's log density at frmtmb's estimates, `max_grad` brms's
gradient there (the inner block for check C):

| shape | measured_const | max_grad |
|---|---|---|
| poisson, `rate(time)` | 1.13687e-12 | 0.000504 |
| poisson("sqrt"), `rate(time)` | -4.54747e-13 | 3.79e-08 |
| negbinomial, `rate(time)` | 5.68434e-13 | 1.21e-05 |
| negbinomial, `rate(time)`, `shape ~ x` | -3.41061e-13 | 2.29e-06 |
| geometric, `rate(time)` | 1.13687e-13 | 3.33e-06 |
| gaussian + poisson, each `subset()`, NA outside y2's rows | 0 | 0.000783 |
| `mi(x, idx = g1)` with x `mi() + index(g2) + subset(s)`, 3 latent values (check C) | 0 | 6.94e-16 |

Row 25 adds `(1 | g)` to the subsetted gaussian response as a check C
shape. The whole gated file: 464 pass, 0 fail, 0 skip
(`dev/aterms2-log-t-brms-likelihood.txt`).

### 3.2 Identities without Stan (`test-subset-rate.R`)

Each is an IDENTITY, checked as a residual:

- `rate()` on poisson against `offset(log(time))`: `obj$fn` and
  `obj$gr` `identical()` at a perturbed point (seed 1).
- Negative binomial and geometric `rate()`: logLik against
  `dnbinom(mu * d, size = shape * d)` and `size = d` written out,
  within `1e3 * eps` relative; the pearson residual against
  `mu d (1 + mu / shape)`.
- A subset model against the sum of the separate fits on their rows,
  at a perturbed point of the joint vector (seed 2), within `64 * eps`.
- `mi(x, idx = )` with x observed against `lm`-style `y ~ x[idx] + w`
  plus `x ~ 1` on x's rows, at a perturbed point (seed 5); and
  `mi(x, idx = g2)` with `index(g2)` naming each row against plain
  `mi(x)`: `obj$fn` `identical()` (seed 6).
- The univariate filter against the fit on `d[d$s1, ]`: `identical()`.

### 3.3 Nothing else moved

`dev/aterms2-bitwise.R`: `fn`, `gr`, `logLik`, `fitted()` and pearson
residuals in hexadecimal at a fixed perturbed point, base against lane,
for ten models the changes pass through without the new terms (poisson
with `(1 | g)`, poisson sqrt, poisson `cens()`, negbinomial with
`shape ~ z`, negbinomial identity, geometric, a two-response model with
`(1 | p | g)`, `mi()`, `sratio()` with `cs()`, a threshold-only
cumulative model with `(1 | g)`). `diff` of the 44 lines is empty
(`dev/aterms2-log-bitwise-base.txt`, `-lane.txt`).

### 3.4 Each new path fails without itself

`dev/aterms2-mutants.R`, log `dev/aterms2-log-mutants.txt`: a mutant
replaces one internal function for the session.

| path | lane | mutant |
|---|---|---|
| `has_rate()`: predict(newdata) at exposure 1000, draws / mu | 1000.075 | 1.009387 |
| `subset_na_rows()`: n_obs of the NA-outside-the-subset model | 120 | 90 |
| `frame_rows()`: the subset model | fits | "missing value where TRUE/FALSE needed" |
| `mi_idx_rows()`: logLik minus the separate fits | 2.862635e-10 | -49.03075 |

### 3.5 Before and after

- `test-subset-rate.R` on rellib-r3 (after punch round 1): 17 of 17
  blocks error, 13 on "Addition term `rate()` / `subset()` / `cat()`
  is not supported" (6 rate, 6 subset, 1 cat) and 4 on "mi() in a
  predictor takes one variable name"
  (`dev/aterms2-log-t-subset-rate-base.txt`). This is
  the refusal the base build makes; the behavioural form is in 3.4.
- `test-subset-rate-draws.R` with the lane's core and the BASE
  frmtmb.sample (`dev/aterms2-mixlib` held a copy of rellib-r3's
  frmtmb.sample, `cp -r rellib-r3/frmtmb.sample dev/aterms2-mixlib/`,
  first on the path through `dev/aterms2-testfile-mix.R` with
  `ATERMS2_MIX=1`; the copy was deleted afterwards so that no installed
  package sits in the tree): 6 failures and 303 warnings
  (`dev/aterms2-log-t-subset-rate-draws-mix.txt`). `log_lik()` of the
  subset model returns 60 recycled columns with "longer object length
  is not a multiple" instead of refusing, `log_lik(resp = "y1")` has
  60 columns instead of 30, and `posterior_predict(newdata =)` at
  exposure 1e4 draws around mu instead of mu * 1e4 (ratio 0.00010).
  With the lane's frmtmb.sample: 17 pass.

## 4. brms's ported suite

`test-brms-suite-standata.R`, gated, one process each arm
(`dev/aterms2-log-g-suite-standata-{base,lane}.txt`). Base: 87 pass.
Lane: 77 pass and 10 fail, every one of the ten "now HOLDS but is
recorded as 'cannot transfer'", and nothing else moved. The rows this
change flips, for lane `defects` to re-verdict in
`dev/brmsport-verdicts*.tsv`:

`standata:171` (cat), `standata:266`, `standata:267`, `standata:268`,
`standata:269` (subset), `standata:620`, `standata:634`,
`standata:641`, `standata:649` (mi idx, index, subset refusals),
`standata:1047` (rate).

`standata:1047` holds because `brms_standata_view()` now gives the
frame's `rate` values brms's name `denom`. `standata:620` is a
different case, corrected in punch round 1 after the review
(`dev/reviews/2026-09-29-aterms2.md`, claim 7): without the helper
edit it ALSO holds, vacuously, because `sdata$idxl_y_x_1` is NULL
there, `all(NULL %in% 9:5)` is TRUE, and `brms_hollow_null()` does
not inspect `expect_true(all(...))`. With the edit the value read is
`8 7 8 8 5 8 9 7 7 6`, identical to brms's own `idxl_y_x_1`, so the
lane's verdict is real. The harness hole predates this lane; the
coordinator is telling lane `defects`. The other eight rows needed
no harness change.

## 5. Test runs and checks

One test file per R process through `dev/aterms2-testfile.R` (the
package attached, as `dev/release/run-tests.R` has it), driven by
`dev/aterms2-suite.ps1`; one log per file in
`dev/aterms2-suite-<tag>/` and the RESULT lines in
`dev/aterms2-suite-<tag>-summary.txt`.

- Ungated, every test file of frmtmb and frmtmb.sample (219 files,
  `NOT_CRAN=true`, the final build): 219 of 219 ran; pass 14293, fail
  0, error 0, warn 0, skip 159
  (`dev/aterms2-suite-ungated-summary.txt`). Every skip is in a file
  of the gated tiers (brms, bcm, drmtmb, fuzz, rl, sampler) except one
  in `test-scale.R`. An earlier run of the same list, made with a
  runner that did not attach the package and before the `cat()` /
  `thres()` refusal was narrowed, failed in `test-conditions.R`,
  `test-data2.R`, `test-id-kron.R` (functions not found without the
  package attached) and `test-thres.R` (`thres(4) + thres(gr = g)`
  reached the new cat message instead of "Duplicated addition term");
  both causes were fixed before this run.
- Gated (`FRMTMB_BRMS_FIT_TESTS=true`), 24 files: the 10 core and 5
  sample `test-brms-suite-*.R`, `test-brms-agreement.R`,
  `test-brms-likelihood.R`, `test-brms-methods.R`,
  `test-brms-priors.R`, `test-brms-port.R`, `test-subset-rate.R`, and
  sample's `test-loo.R`, `test-sampling-ported.R`,
  `test-subset-rate-draws.R`: 24 of 24 ran, 0 skip in any file, and
  the only failures are the ten flipped rows of section 4
  (`dev/aterms2-suite-gated-summary.txt`).
- `test-brms-likelihood.R` on rellib-r3: 440 pass, 3 errors, the three
  new rows, on "`rate()` / `subset()` is not supported" and "mi() in a
  predictor takes one variable name"
  (`dev/aterms2-log-t-brms-likelihood-base.txt`).
- `R CMD check --as-cran` (`dev/aterms2-check.ps1`, logs
  `dev/aterms2-check/{core,sample}/check.log`): frmtmb `Status: 1
  NOTE`, the V8 math-rendering note on the HTML manual; frmtmb.sample
  `Status: OK`. Rerun after punch round 1 on the final tree, with the
  same two results (the logs are the rerun's).

## 6. Decisions, and what was NOT done

- **Three intentional divergences, one cause.** brms keeps every
  level of the WHOLE data in a subsetted response; frmtmb builds a
  response's design on its own rows and keeps the levels those rows
  carry, so it creates no parameter that the data cannot place. The
  review measured all three (`dev/aterms2-rev-log-02a.txt`):
  1. Grouping levels. A level with no rows adds nothing to the Laplace
     marginal likelihood, so the fits agree (3.1 row 25 uses data where
     the sets are equal; the review's `gr(cov = A)` case agrees to
     2.7e-12); `ranef()` lists fewer levels, and `|ID|`-linked terms
     whose rows carry different levels are refused with a message
     naming subset().
  2. Factor predictor levels. A level that occurs only outside the
     subset is an all-zero column `fc` in brms's `X`, a coefficient
     only its prior places; frmtmb drops the column with its
     rank-deficiency message, and the fit equals the fit on `d[s, ]`.
  3. Ordinal categories. A category that occurs only outside the
     subset gives brms 3 thresholds where frmtmb fits 2; the extra one
     has no data above it, the case `?frm` "Ordinal thresholds" says
     needs a prior to be held.
  `?frm` "Response subsets" and the migration vignette now name all
  three (punch round 1; before, they said one difference remained).
  Global levels would need each component's Z, each X and each
  threshold count built from a level set computed outside the
  response's rows; not done, by choice.
- `subset()` with `rescor`, `me()` (as brms), `na.exclude`, and in a
  multivariate model with an autocorrelation term or a structured
  family are refused by name, and so is `vcov_cluster()` on a
  multivariate subset model.
- `conditional_effects()` of a model with `mi(x, idx = )` stops with
  "Could not match all indices" on its grid, whose idx and index values
  are reference values. Not changed.
- `cat()` on a non-ordinal family is refused with `thres()`'s message.
  `cat(x = 4)` is accepted, as brms accepts it (punch round 1).
- brms's `negbinomial2` does not exist in frmtmb, so its `rate()` does
  not either.
- `rate()` with `cens()` or `trunc()` on poisson runs through the CDF
  at `mu * d`, which is brms's code, but no brms comparison was run
  (compat row "untested").
- The multivariate default of `conditional_effects()` covers the first
  response only; that predates this lane.

## 7. Defects found, not fixed

Found by the review on BOTH builds, so not this lane's; recorded,
not fixed (`dev/aterms2-rev-log-06b-ceoffset.txt`,
`dev/aterms2-rev-log-12-preexisting.txt`):

- `conditional_effects()` on `y ~ x + offset(log(time))` fails with
  "non-numeric argument to mathematical function".
- `emmeans()` on the same model fails with "undefined columns
  selected".
- An addition term given an expression, `weights(wt * 2)`, fails
  with "invalid model formula in ExtractVars"; so does
  `rate(time * 2)`, which brms accepts.

## 8. Punch round 1

The review (`dev/reviews/2026-09-29-aterms2.md`) found the lane
MERGEABLE with four minors and three notes. What was done:

1. **`rate()` exposure on newdata.** `aterms_for_newdata()`
   (`R/predict.R`) repeats the frame's positivity check, with brms's
   sentence; NA stays allowed and predicts NA. It covers `fitted()`,
   `predict()` and frmtmb.sample's `posterior_predict()`, which all read
   newdata's exposure through it. New block "newdata's exposure must be
   positive" in `test-subset-rate.R`, with the absent case (1e-3
   predicts `mu * 1e-3`). On rellib-r3 the block errors on "`rate()` is
   not supported"; the behavioural form is `dev/aterms2-mutant-newrate.R`
   (seed 24, log `dev/aterms2-log-mutant-newrate.txt`): the lane refuses
   time -2 and 0, and `aterms_for_newdata()` without the check returns
   an expected count of -3.130614 at -2 and 0 at 0.
2. **`?frm` subset text.** Rewritten to name the three differences
   (section 6); same in `vignettes/brms-migration.Rmd` and the compat
   row.
3. **`nobs()` as brms's.** `nobs(fit)` is the data's rows before any
   `subset()` cut, for a univariate subset model too (the frame keeps
   `nobs_data`); `nobs(fit, resp = "y1")` returns that response's rows
   in a multivariate model and the data's rows in a univariate one,
   and refuses a name that is no response. `resp` sits AFTER the dots,
   so `nobs(fit, 7)` is still refused as an unnamed argument
   (`test-arg-refusal.R`), and brms's positional `nobs(fit, "y1")` is
   refused the same way; named, it works. The test row that expected
   the old refusal (`nobs(fit, resp = "y")`, "same nobs") is removed
   from `test-arg-refusal.R`. The likelihood keeps its own count:
   `logLik()`'s `nobs` attribute, which `BIC()` reads, and
   `df.residual()` stay the rows the likelihood sums over, and
   `simulate()` draws one value per fitted row (it read `nobs()` and
   now reads the frame). frmtmb.sample's `nobs()` passes `resp` to the
   fit's method. Tests in `test-subset-rate.R` and
   `test-subset-rate-draws.R`.
4. **Section 4, row 620.** Corrected above.

Notes: a character subset variable is refused, as brms refuses it,
with a factor of the same values accepted (test added); `cat(x = 4)`
is accepted and `cat(k = 4)` refused (test added); the `predict()`
behind the two statistical thresholds of the rate block is seeded
(seed 11). The three failures found on both builds are in section 7.

Reruns (one process each, `dev/aterms2-suite-punch1*-summary.txt`):
`test-subset-rate.R` 85 pass, `test-arg-refusal.R` 120, `test-compat.R`
680, sample `test-subset-rate-draws.R` 19, all 0 fail, 0 error, 0 warn,
0 skip; also `test-methods-audit.R` 58, `test-methods.R` 64,
`test-mv-gaps.R` 57, `test-multivariate.R` 23, `test-structure.R` 135,
`test-simulate-ergonomics.R` 50, `test-predict-newdata.R` 12, sample
`test-draws-methods.R` 148; gated `test-brms-likelihood.R` 464,
`test-brms-suite-standata.R` 77 pass and the same 10 "now HOLDS", sample
`test-brms-suite-methods.R` 61 (this second group on the build before
`resp` moved after the dots, which only `test-arg-refusal.R` reaches).
On rellib-r3, `test-subset-rate.R` has
17 of 17 blocks erroring (section 3.5).

## 9. Punch round 2

The re-check (`dev/reviews/2026-09-29-aterms2.md`, "Re-check after
punch round 1") found one blocker, from round 1: frmtmb.eam's
`test-family.R`, "a term a family does not accept reads refused, not
untested", failed 5 times on the lane (wiener, gddm, lba, rdm,
wiener_gng) and passed 51 of 51 on rellib-r3. Its premise, that every
aterm outside a family's `accepts_aterms` reads "refused", became
stale when core's compat table began to read `subset()` and `index()`
"works" for every family. The cell is right: the reviewer measured
`subset()` on `wiener()` and `lba(3)` equal to the fit on `d[s, ]`.

- Fix, in `extensions/frmtmb.eam/tests/testthat/test-family.R` only:
  the row terms are exempt from the refused check, read from core's
  own `frmtmb:::row_aterms` rather than a literal, and the test now also
  asserts that the family's table holds exactly `length(row_aterms)`
  such cells and that each reads "works", so the exemption cannot pass
  on an empty set. The test therefore needs a frmtmb that has
  `row_aterms` (this change); on rellib-r3's frmtmb it errors on the
  missing object.
- The same pattern in the other extensions: every test that reads
  `frm_compat_features()` or `kind == "aterm"`
  (`frmtmb.coupling/test-surface.R`, `frmtmb.learn/test-surface.R`,
  `frmtmb.ode/test-compat.R`, `frmtmb.sample/test-compat-preflight.R`,
  `frmtmb.eam/test-wiener-cdf.R`) was read; none asserts refusal over
  every unaccepted aterm, and the full run below confirms no other
  failure.
- Optional notes, both done: `?frm` now says the log-likelihood's count
  is the rows of the fitted frame, which in a multivariate model
  includes a row that no response's subset uses; and the fit's
  character-subset refusal applies on `newdata` too
  (`subset_newdata()`, tested in `test-subset-rate.R`).

**Process note.** The lane rule asks for every suite that can reach a
changed path. Round 1 and punch round 1 ran only core and
frmtmb.sample; the change to core's compat table reached frmtmb.eam's
suite, which was not run until the re-check. Every extension's
ungated suite now runs against this build.

Full ungated run, one file per process, `NOT_CRAN=true`, lane frmtmb
and frmtmb.sample from `wt-aterms2-lib`, the other extensions from
rellib-r3 (their `R/` is unchanged in this worktree), tests from the
worktree (`dev/aterms2-suite-list-full.txt`,
`dev/aterms2-suite-full/`, per-file lines in
`dev/aterms2-suite-full-bypkg.txt`, counts generated into
`dev/aterms2-suite-full-counts.txt`): 311 of 311 files ran.

```
frmtmb files=183 pass=12360 fail=0 error=0 skip=155 warn=0
frmtmb.coupling files=11 pass=542 fail=0 error=0 skip=5 warn=0
frmtmb.eam files=29 pass=1743 fail=0 error=0 skip=3 warn=0
frmtmb.latent files=10 pass=359 fail=0 error=0 skip=2 warn=0
frmtmb.learn files=15 pass=429 fail=0 error=0 skip=13 warn=0
frmtmb.ode files=11 pass=547 fail=0 error=0 skip=1 warn=0
frmtmb.sample files=37 pass=1965 fail=0 error=0 skip=4 warn=0
frmtmb.spline files=15 pass=553 fail=0 error=0 skip=1 warn=0
```

frmtmb.eam `test-family.R`: 61 pass, 0 fail. Core `test-subset-rate.R`:
86 pass.

## 10. Versions

A minor bump for frmtmb (new addition terms and two new sampling-API
exports). frmtmb.sample needs the frmtmb that exports
`subset_resp_check()` and `subset_newdata()`, so its floor moves to the
next frmtmb release; a minor bump for it too.
frmtmb.eam's `test-family.R` now reads `frmtmb:::row_aterms`, so its
test suite needs this frmtmb; its package code is unchanged, so no
version bump or floor change for its code alone.

`R CMD check --as-cran` after punch round 2 (`dev/aterms2-check.ps1`,
`R_LIBS` = the lane library, then rellib-r3): frmtmb `Status: 1 NOTE`
and frmtmb.eam `Status: 1 NOTE`, both the V8 math-rendering note on the
HTML manual (`dev/aterms2-check/{core,eam}/check.log`). frmtmb.sample
did not change in this round; its punch round 1 check stands
(`Status: OK`).
