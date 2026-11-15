# Lane parity2: brms surface parity (round of 2026-11-13)

Base ec8178d3, frmtmb 0.82.0, frmtmb.sample 0.29.0. Brief:
`C:/Users/adf44/source/r/round-briefs/r083-parity2.md`. Worktree
`frmtmb-wt-parity2`, library `wt-parity2-lib` (core and frmtmb.sample
from the worktree; the other seven extensions read from rellib-r21).
Every change is uncommitted. Logs are local, under `dev/parity2-out/`;
every number below is pasted from the log named beside it.

STATUS: fix round of review `dev/reviews/2026-11-13-parity2.md` (FIX)
done; see "Fix round" at the end, which supersedes the earlier
sections where they differ (1497 above all: a fit keeps "gaussian").
Bumps: frmtmb minor (two new sampling-API exports, BREAKING behavior);
frmtmb.sample minor (BREAKING behavior). Floor: frmtmb.sample's
DESCRIPTION still says `frmtmb (>= 0.82.0)`; it must rise at release to
the frmtmb release that carries `ce_snl()` and `latres_err()`, which
the consolidation sets with the version. On 0.82.0's core, a CE at
`re_formula = NULL` on draws stops with "unused argument (ce_snl(dots))"
and latent draws with 'could not find function "latres_err"' (review,
`floor.txt`).

## Scripts

- `dev/parity2-data.R`: the data both arms read (seed 1, n = 150, the
  triage's check2.R data plus columns; seed 11 for the random-intercept
  data, seed 3 for drop_unused_levels, seed 5 for multinomial/ordinal).
- `dev/parity2-brms-probe.R` (3922, 3919), `dev/parity2-brms-probe2.R`
  (everything after 3919), `dev/parity2-brms-probe3.R` (ordinal draw
  margin under ndraws, err's place): brms 2.23.0, `fixed_param` fits
  at init 0 for names and layouts; one sampled fit (2 chains of 1000,
  seed 1) for 1497 and 5331. Logs `brms-probe.txt`, `brms-probe2.txt`,
  `brms-probe3.txt`.
- `dev/parity2-frm-probe.R base|lane`: frmtmb on the same items. Logs
  `frm-base.txt` (rellib-r21) and `frm-lane*.txt`.
- `dev/parity2-run-test.R`, `dev/parity2-runjobs.sh`: one test file per
  R process, at most 3 at once (copies of sampfront's).
- `dev/parity2-install.ps1`: roxygenise and install into the lane
  library. `dev/parity2-roxy.R`: roxygenise with full output.
- `dev/parity2-lpinspect.R`, `dev/parity2-center.R`, `dev/parity2-rr.R`:
  what a fit's linear predictors hold, the intercept's `center` record,
  and the traceback of 4337.
- `dev/parity2-backlog-edit.ps1`: the CLOSED lines and the filed items.

## Item 3922: conditional_effects() default effects

brms's rule (`get_all_effects.brmsterms()`, read in
`dev/brms-suite/brms/R/conditional_effects.R:524`): for every dpar,
then every nlpar, each fixed-effect term label and each `cs()` term is
one effect of the variables it reads (`get_var_combs()`); each `sp`
(mo/me/mi), `sm` and `gp` term gives its variables joined by `*`, so
alone and in pairs; a nonlinear predictor gives its covariates alone and
in pairs; reserved `intercept` is removed; sets of more than two
variables are dropped. frmtmb used the selected dpar's variables plus
the leading pair of each interaction.

| formula | brms 2.23.0 | frmtmb base | frmtmb lane |
|---|---|---|---|
| `y ~ exp(x):z` | x:z | x z x:z | x:z |
| `y ~ t2(x, z)` | x z x:z | x z | x z x:z |
| `y ~ mo(k):z` | k z k:z | k z | k z k:z |
| `yo ~ cs(x) + z`, sratio | z x | z:cats__ | z:cats__ x:cats__ |
| `bf(y ~ x, sigma ~ exp(z))` | x z | x | x z |
| `bf(y ~ x) + bf(y2 ~ z)` | y.y_x y2.y2_x y.y_z y2.y2_z | same | same |
| `bf(y ~ mi(x2) + z) + bf(x2 \| mi() ~ w)` | y.y_z x2.x2_z y.y_x2 x2.x2_x2 y.y_w x2.x2_w | same | same |
| `y ~ x * z * w` | x z w x:z x:w z:w | same | same |
| `y ~ a * exp(b * x)`, `a ~ 1 + w`, `b ~ z` | x w z | same | same |
| `bf(y ~ x + (1 + z \| g), sigma ~ w)` | x w | x | x w |
| `y ~ poly(x, 2) + z + x:w` | x z x:w | x z w x:w | x z x:w |

Sources: `brms-probe.txt` (both `get_all_effects()` and the names of a
`fixed_param` fit's `conditional_effects()`), `frm-base.txt`,
`frm-lane1.txt`. The cs() row differs only in the key: frmtmb's ordinal
default is the per-category display (a recorded decision, ledger row
brmsfit-methods:217), whose key is brms's `x:cats__`; the effects are
brms's two, in brms's order. The mv and mi rows already matched; the
triage compared `get_all_effects()` with the output names and so
counted them as differences.

Change: `ce_star()`, `ce_lp_effects()`, `ce_default_effects()` in
`R/conditional-effects.R` replace `ce_lp_pairs()`/`ce_plot_pairs()`;
`ce_grids_build()` takes them when `effects` is NULL, over every
response (brms reads the whole model, so `resp = "y2"` on the mv model
gives `y2.y2_x y2.y2_z`). Kept deliberately: a factor-smooth's own
grouping factor stays out (frmtmb drops that term at the population
level; filed LOW), and mo/me/mi terms come in that order rather than
written order (filed LOW).

Test `tests/testthat/test-ce-default-effects.R`: base pass=3 fail=8;
lane pass=11 fail=0 (`t1/`).

## Item 1497: new levels follow sample_new_levels = "uncertainty"

brms 2.23.0, `y ~ x + (1 | g)` on the seed-11 data, sampled fit, band
at x = -1, 0, 1 with `re_formula = NULL` (`brms-probe2.txt`):

| rule | width at -1, 0, 1 |
|---|---|
| brms default | 1.971 1.948 2.031 |
| brms "uncertainty" | 1.971 1.948 2.031 |
| brms "gaussian" | 3.633 3.527 3.559 |
| frmtmb base, boot 200, seed 2 (gaussian) | 2.040 2.163 2.279 |
| frmtmb lane, boot 200, seed 2 (uncertainty) | 1.883 1.938 2.064 |

brms's rule (`get_new_rdraws()`, prepare_predictions.R:1024): per draw
and new level, the effects of one observed level picked at random
(within the row's `by` level); "old_levels" picks once per call.

Change: `ce_plan_draw()` (R/ce-levels.R) reads the plan's rule; the plan
carries it (`ce_level_plan(..., sample_new_levels)`) with an
environment for the "old_levels" pick; `ce_plan_key()` includes it, so
boot draws of one rule are refused under the other; `ce_snl()` reads
and checks the argument; `ce_dots()` lets it through. The fit method
and frmtmb.sample's draws method pass it. Blocks whose levels are their
structure: "uncertainty" lends a `gr_cov`/`gr_prec` level as brms does,
and `car`/`spde` stay at 0 (brms does not draw a new level there). The
Wald band is unchanged: it adds the gaussian variance of a new level.

Tests: `tests/testthat/test-ce-new-levels-snl.R` (every draw is one of
the six fitted modes; old_levels one mode; gaussian none) base pass=1
fail=2 err=2, lane pass=10 fail=0. frmtmb.sample
`test-draws-brms-surface.R`, first block: per draw the shift of the
new-level curve from the population curve is one of that draw's six
`r_g` values, to 64 ulps; under "gaussian" none is.

## Item 3919: s(exp(x))

brms 2.23.0: fits; `bs_sexpx_1`, `sds_sexpx_1`, `s_sexpx_1[1..8]`; CE
default `x` (`brms-probe.txt`). frmtmb base: "The model uses `exp(x)`,
which is not a column of `data`". Cause: `dpar_frame_rhs()` (R/frame.R)
put a smooth's term into the frame formula with `as.name()`. Change:
`str2lang()`. Lane: variables `b_Intercept bs_sexpx_1 sds_sexpx_1
sigma`, CE `x` (`frm-lane1.txt`). Test (in
`test-brms-surface-mixed.R`): the logLik equals that of `s(ex)` with
`ex = exp(x)` to 1e-6 relative, fitted values on new data too.

## Item 4160: attr(x, "min") of a mo() variable

brms 2.23.0 (`brms-probe2.txt`): 1..4 with `attr(k, "min") <- 0` fits
Jmo = 4, Xmo_1 = 1 2 3 4; new data 1..4 are coded 0 1 2 3 and 0 is
refused ("Invalid values in variable 'k': '0'"), because
`validate_newdata()` overwrites the attribute with the observed minimum
(brms-40, dev/upstream-bugs.md, not followed). An attribute above the
data gives codes -1 0 1 2 in brms. frmtmb base: min 1, D 3.

Change (R/frame.R, the mo() loop): an integer variable's `attr(, "min")`
is read from the data at the fit; one whole number no larger than the
smallest value, else refused by name. The frame keeps it, so new data
are coded from it (predict.R's `mo_min()` is unchanged). Lane: min 0,
D 4, codes 1..4; new data 0..4 answer. The declared category holds no
row, so its simplex cell and the intercept trade off, and the SE check
names them; that is brms's model, which the user asked for.

## Item 3373: REML df

lme4 1.x and frmtmb on `Reaction ~ Days + (Days | Subject)` REML
(`frm-base.txt`): lme4 AIC 1755.628, df 6; frmtmb base AIC 1751.628,
df 4, logLik -871.8141 on both. Change: `n_outer_est()` counts beta
under REML, `mo_const_df()` subtracts every constant mo() coefficient,
`frame_n_outer()` no longer drops beta. Lane (`frm-lane2.txt`): AIC
1755.628, df 6, BIC 1774.786, equal to lme4's on every digit printed.

## Item 1560: a later lf() on a dpar bf() sets

brms 2.23.0: `bf(y ~ x, sigma ~ z) + lf(sigma ~ w)` gives `sigma ~ w`
(`brms-probe2.txt`), with "Replacing initial definitions of parameters"
(brmsformula.R:671-676). frmtmb base: "lf() sets 'sigma', which the bf()
it is added to already sets". Change: `bform_replace_pars()` in R/bf.R
drops the earlier formula, fixed value, equation or nonlinear body and
gives brms's message, for `lf()` and `nlf()`. The other half of the
backlog item (`update(fit, bf())` pooling) is not done; filed again.

## Item 5337: drop_unused_levels = FALSE

brms 2.23.0 `make_standata(y ~ f + (1 | g), drop_unused_levels =
FALSE)`: N_1 7, K 3, X columns Intercept fb fc with column sums 60 30 0;
with `(1 + f | g)` N_1 7, M_1 3; default N_1 6, K 2; draws carry
`b_fc` and `r_g[7,Intercept]` (`brms-probe2.txt`). frmtmb base: fixef
Intercept fb, ranef rows 1..6.

Change: `reformulas::mkReTrms(drop.unused.levels = drop_unused_levels)`
in R/frame.R. Lane: ranef rows 1..7, level 7 at 0, the logLik unchanged
(test). Not changed: the zero column `fc` stays dropped as an aliased
column. A coefficient of a column of zeros has no maximum likelihood
estimate, and brms draws it from its prior, an improper one under its
flat default. Filed as a decision, recommendation keep.

## Item 5315: b_Intercept under 0 + Intercept

brms 2.23.0 `variables()` (`brms-probe2.txt`):
- `bf(y ~ 0 + Intercept + x, sigma ~ z)`: b_sigma_Intercept b_Intercept
  b_x b_sigma_z
- `y ~ 0 + Intercept + x`: b_Intercept b_x sigma
- `bf(y ~ x, sigma ~ 0 + Intercept + z)`: b_Intercept b_x
  b_sigma_Intercept b_sigma_z

frmtmb base gave `b_Intercept b_sigma_Intercept b_x b_sigma_z` on all
three (`parity2-center.R` output). Change: `brms_variable_sort()`
(R/brms-names.R) treats the intercept of a linear predictor with
`center = FALSE` (0 + Intercept or bf(center = FALSE)) as an ordinary
column, as it already did for a nonlinear parameter's. Lane: brms's
order on all three (test).

## Item 5323: insight::get_predicted() on multinomial()

brms 2.23.0: `insight::get_predicted()` fails with "'length = 2' in
coercion to 'logical(1)'" (inside insight; `brms-probe2.txt`), with and
without `predict = "expectation"`. frmtmb base: class get_predicted,
length 180 (60 rows x 3 categories flattened). Change: refused by name
for every family whose expectation is a matrix (`ce_matrix_mean()`:
multinomial and the simplex families), before any `predict` type.

## Item 5325: ordinal posterior_epred() draw margin

brms 2.23.0 (`brms-probe2.txt`): cumulative fit, dim 20 60 3, dimnames
`"1".."20"`, NULL, `"1" "2" "3"`; on newdata the same; categorical:
NULL, NULL, categories. Change: frmtmb.sample's `posterior_epred()`
names the first margin `"1".."S"` for an ordinal family. Three
assertions of `test-sampling-ported.R` that pinned NULL now pin brms's.

## Item 5331: loo_*() on a subset of the fitted rows

brms 2.23.0, sampled fit, `newdata = dg[1:10, ]` (`brms-probe2.txt`):
`loo_epred` and `loo_linpred` 0 difference from the fitted rows 1:10;
`loo_predict` 0.0428 and `loo_predictive_interval` 0.128 (they
simulate); `loo_R2(newdata = )` answers (R2 0.613). Change
(frmtmb.sample `loo-predict.R`): `loo_draws_opts(subset = TRUE)` lets a
subset through, `loo_draws_expect()` cuts the prediction to the rows and
takes the PSIS weights of `log_lik(newdata = )`, `loo_R2()` cuts y, the
prediction and the log-likelihood. Test: identical to the full answer's
rows on a reordered subset with a repeated row (`c(5, 3, 20, 3)`), for
loo_epred, loo_linpred, loo_predict and loo_predictive_interval at one
seed; moved rows still refused.

## Item 5300: err[k] of a latent residual

brms 2.23.0 (`brms-probe2.txt`), poisson `ar(time = t, gr = g, cov =
FALSE)` and `cov = TRUE`, 48 rows: `err[1..48]` saved, `zerr` not saved
(count 0) on both. Change: core `latres_err()` (exported in the
sampling API: `sderr * zerr` under cov = FALSE, `latres_expand()`
otherwise); frmtmb.sample `draws_latres_cols()` adds `err[k]` as a
derived block after the sampled columns, and drops it again when a draw
is read back. Test: on the cov = FALSE model `err[k] = sderr * zerr[k]`
to 64 ulps, err is in `variables()`, and not in the internal matrix.
Not done (larger than names, filed again): the corrected estimate and
logLik, `get_joint_cov()` on a new row, M7.

## Item 2980: frm_multiple() + int_conditions

brms 2.23.0 `brm_multiple()` + `conditional_effects(effects = "x:z",
int_conditions = list(z = c(-1, 1)))`: answers, z = -1 1
(`brms-probe2.txt`). frmtmb base: "formal argument "int_conditions"
matched by multiple actual arguments". Cause: the j > 1 imputations got
the caller's `int_conditions` in `...` and the first imputation's grid
as `int_conditions =`. Change (R/multiple.R): a wrapper drops the
caller's two from `...`.

## Item 4337: default_prior() on rr()

frmtmb only (brms has no rr()). Cause: `prior_table()` built a data
frame of class "sd" rows from `block_sd_prefix()`, which has 0 rows for
an rr() block (its loadings carry its scale), beside a length-1 group.
Change (R/priors.R): skip a block with no sd. Lane: the table lists the
b rows and no sd row for g (test).

## Later additions

- 1560 also for `nlf()`: brms 2.23.0 replaces on all three orderings
  (`bf() + nlf()`, `nlf() + nlf()`, `nlf() + lf()`), each with
  "Replacing initial definitions of parameters 'sigma'"
  (`dev/parity2-brms-nlf.R`, `brms-nlf.txt`). `test-nlf.R` pinned the
  three refusals; it now pins brms's replacement.
- 5325 under subsetting (`brms-probe3.txt`): `ndraws = 5` gives
  `"1".."5"` and `draw_ids = c(3, 7, 12)` gives `"1" "2" "3"`, so the
  names count the returned draws, as `seq_len()` here does.
- 5300's place (`brms-probe3.txt`): brms lists `... sderr Intercept
  lprior lp__ err[1] ... err[48]`, err at positions 8 to 55 of 55, last:
  `reorder_pars()` does not list class err, and an unlisted class sorts
  after lp. `brms_class_levels` gains "err" after "lp".
- 3922 reached three existing assertions that pinned the old defaults:
  `test-brms-methods.R` ("conditional_effects(dpar =) enumerates ...",
  which recorded brms's `x z` as a divergence; now equal, with the flat
  z panel compared value for value), `test-nlf.R` (a nonlinear sigma's
  covariate beside mu's: `x z`, brms's rule). 1497 reached
  `test-ce-levels.R`, frmtmb.sample `test-postfit-draws.R` (two
  tests) and `test-laplace-draws.R`, whose assertions compare with the
  Wald band or a brms run under the gaussian law; they now ask for
  `sample_new_levels = "gaussian"`, and the laplace test also pins that
  "uncertainty" is refused on laplace draws (it reads the observed
  levels' effects, which such draws do not hold). 5331 reached
  `test-loo-predict.R` and `test-loo-newdata-subset.R`, which pinned the
  refusal of a subset; they now pin the refusal of moved rows.
- The rounding check of the one new warning assertion (the attr-min fit
  names its flat parameters): `dev/parity2-ulp.R`, y * (1 + 2e-15 k),
  k = 0..24, data seed 9: "SE warning on 25 of 25" on reference BLAS
  (`ulp.txt`) and on OpenBLAS 0.3.26 with LAPACK (`ulp-ob.txt`), no
  other warning on any k.

## Files touched

Core: `R/conditional-effects.R` (default effects, `ce_snl()`,
`ce_dots()`, docs), `R/ce-levels.R` (`ce_level_plan()`,
`ce_plan_draw()`, `ce_plan_key()`), `R/frame.R` (smooth terms in the
frame, `attr(, "min")`, `mkReTrms(drop.unused.levels)`), `R/bf.R`
(`bform_replace_pars()`, `plus_bf()`), `R/parse.R` (a comment),
`R/methods-fit.R` (REML df), `R/brms-names.R` (`brms_variable_sort()`,
`brms_class_levels`), `R/priors.R` (one `next` in `prior_table()`),
`R/multiple.R`, `R/insight.R`, `R/autocor-latent.R` (`latres_err()`),
`R/autocor.R` (docs), `R/sampling-api.R` (exports `ce_snl`,
`latres_err` on their own lines, docs), `NAMESPACE`, `man/` (4 Rd),
`NEWS.md`, `dev/test-backlog.md`. Tests: new
`test-ce-default-effects.R`, `test-ce-new-levels-snl.R`,
`test-brms-surface-mixed.R`; changed `test-brms-methods.R`,
`test-ce-levels.R`, `test-nlf.R`, `test-portability.R`.

frmtmb.sample: `R/conditional-effects-draws.R`, `R/draws-brms.R`
(`draws_latres_cols()`), `R/loo-predict.R`, `R/methods-draws.R`
(ordinal draw margin, the refusal note), `man/` (2 Rd), `NEWS.md`.
Tests: new `test-draws-brms-surface.R`; changed `test-laplace-draws.R`,
`test-loo-newdata-subset.R`, `test-loo-predict.R`,
`test-postfit-draws.R`, `test-sampling-ported.R`.

Shared areas touched outside the lane's own: `R/priors.R` (priors2's;
one line in the sd-row loop, no theta line), `R/sampling-api.R` (two
export lines), `R/brms-names.R` (two lines). No compat rule and no
ledger file changed; the ledger row `standata:1133` (the zero column of
an unused factor level) stays true, since that column is still
dropped.

## Runs

All on the lane library (core and frmtmb.sample from the worktree,
the other seven extensions from rellib-r21; each log's `lib:` line).

- Every file of all nine suites, one file per process, 3 at once,
  `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true` (`all1/`,
  `all1-results.txt`): 533 of 533 files with a RESULT line, pass 37485,
  skip 29 (all gated by other variables: drmTMB agreement 13, the scale
  tier 15, the fuzz tier 1). Seven files failed there; each was a test
  pinning a behavior this lane changed (see "Later additions") and each
  was updated and rerun (`t3/`, `t4/`): test-brms-methods 1015 pass,
  test-ce-levels 101, test-nlf 83, test-portability 97, frmtmb.sample
  test-laplace-draws 119, test-loo-newdata-subset 53, test-loo-predict
  153, test-postfit-draws 89, test-sampling-ported 215, all fail 0,
  err 0, warn 0. The drmTMB agreement file under
  `FRMTMB_DRMTMB_FIT_TESTS=true` (REML fits): lane pass=131, base
  pass=131.
- New test files, lane against base (`t1/`, `t2/`, `t5/`):
  test-ce-default-effects lane 11/0 base pass=3 fail=8;
  test-ce-new-levels-snl lane 10/0 base pass=1 fail=2 err=2;
  test-brms-surface-mixed lane 31/0 base pass=1 fail=11 err=6;
  frmtmb.sample test-draws-brms-surface lane 21/0 base pass=4 fail=4
  err=3. Every base failure is behavioral except where an argument did
  not exist (`sample_new_levels`, `ce_level_plan()`'s fifth argument).
- Name and latent-residual files after the err order change (`t4/`):
  17 files, all fail 0.
- One session per suite, three orders (`sessions.log`): frmtmb forward
  pass=24446, reverse 24611, preload 24714; frmtmb.sample 4260 on each;
  fail 0, err 0 on all six; every LEAK line is tagged "with load of".
- OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
  lapack`, `t-ob/`, `t-ob3/`): the ten changed or new files of both
  packages, all fail 0; test-ce-default-effects first showed one
  escaped warning (nlminb code 8 at a measured maximum on the nonlinear
  fit), now allowed by its text; rerun warn=0 on both BLAS. The OpenBLAS
  copy was deleted afterwards.
- Equivalence pairs that reach a changed path (latent residuals, REML),
  `dev/parity2-equiv-out/`: p028, p044, p112, p132, p154, p167
  equivalent L1, p128 equivalent L4, p156 refusal correct; each equal
  to the record's verdict in `dev/equiv-data/results.tsv`.
- R CMD check --as-cran: frmtmb "Status: 2 NOTEs" (examples timing,
  frm_hazard_reads 5.02 s elapsed, a load measurement; V8 for the HTML
  math), tests OK in 65 min; frmtmb.sample "Status: 1 NOTE" (V8), tests
  OK. The core check ran before `test-brms-surface-parity2.R` was
  renamed `test-brms-surface-mixed.R` (content unchanged).
- R process cap: twice a fourth R process ran next to three test
  processes (the first brms probe 3, which I stopped and reran later;
  and `dev/parity2-brms-nlf.R`, about 5 seconds).

## Not done, and why

- 5300's corrected estimate, `get_joint_cov()` on a new row and M7: not
  names, and each a model change of its own; filed again.
- 1560's `update(fit, bf())` pooling: filed again.
- 5337's zero fixed-effect column: a decision, filed with a
  recommendation.
- The ordinal per-category default of `conditional_effects()` is kept
  (an earlier decision); only its effects follow brms.

## Fix round (review dev/reviews/2026-11-13-parity2.md, FIX)

### F1 (1497): a fit's boot band keeps "gaussian"

The review: under "uncertainty" a bootstrap replicate lends the refit's
conditional MODE of a level, where brms lends a posterior DRAW of the
level's effect; modes shrink, so the band is too narrow where shrinkage
is material. Its data: 30 groups of 2 rows, sd_g 1, sigma 2, data seed
100 (ML sd_g 0.941). Widths at x = -1, 0, 1 (review's table): brms
"uncertainty" 4.096 3.963 3.922; fit boot 400 lending modes 2.903 2.815
2.787; fit boot 400 gaussian 4.036 4.344 4.560; Wald 3.967 3.840 3.940.

Option (a), as the coordinator chose: `ce_snl_fit()`
(R/conditional-effects.R) gives the fit's display the default
"gaussian", as `predict()` and `fitted()` have it, and refuses
"uncertainty" and "old_levels" on `band = "boot"` by name (a refit holds
each level's mode, so a lent level carries none of its conditional
variance; use "gaussian" or the draws of frm_sample()). "old_levels" is
refused too, since it lends a mode in the same way. Draws keep brms's
"uncertainty" default (frmtmb.sample passes `ce_snl(dots)`, default
"uncertainty"). Option (b) is filed MEDIUM in the round's backlog
section with the review's widths.

After (`dev/parity2-snl-width.R`, `snl-width.txt`, lane library): the
default boot band at 400 refits, seed 2, is 4.036 4.344 4.560, the
review's gaussian figures to the digit; at 200 refits 4.031 3.690
4.264, ratio to Wald 1.016 0.961 1.082. New test (second block of
`test-ce-new-levels-snl.R`): min(boot 200 / Wald) > 0.8 on that data.
Seen to fail on the lane build before the fix: "Expected `min(wb/ww)` >
0.8. Actual comparison: 0.734" (`t6-prefix/`). The threshold: the width
of a 95 percent band from 200 refits has a Monte Carlo sd of about 7
percent (binomial quantile error with the normal density at 1.96), so
0.8 is three of those below 1 and above the lending band's 0.73.

### F2: a band that draws nothing refuses a rule it would ignore

`ce_snl_fit()`: on `band = "wald"` or `"profile"` any value but
"gaussian" is refused, "this band draws no new level; it adds the
group-level term's variance for one, which is the "gaussian" law"
(`snl-width.txt`: "uncertainty" and "old_levels" on wald and on boot
each refused). Before: `identical()` to the call without it (review).
Test: the same block asserts both refusals per value and that
`sample_new_levels = "gaussian"` on Wald is `identical()` to the
default call (the guard's absent case).

### F3 (3922): a default pair is ordered numeric first

brms 2.23.0 (`dev/parity2-brms-probe4.R`, `brms-probe4.txt`): `y ~ f *
x` gives `f x x:f`; `effects = "f:x"` gives `f:x`. Change:
`ce_numeric_first()` in `ce_grids_build()`'s default branch swaps a
pair whose first variable is a factor, character, logical or grouping
variable and whose second is not (brms's `prepare_cond_data(reorder =
TRUE)`); named effects keep the caller's order. Before: `f x f:x`
(review, both arms). Test: `y ~ f * x` gives `f x x:f` and
`effects = "f:x"` stays `f:x` (test-ce-default-effects.R).

### LOW items

- `` s(`my x`) ``: brms 2.23.0 stops parsing it, "<text>:1:4: unexpected
  symbol 1: my x" (`brms-probe4.txt`). The lane build had fitted it with
  the names ``bs_s`myx`_1``. Now refused by name in `dpar_frame_rhs()`
  ("not a syntactic name ... brms 2.23.0 cannot parse"); test in
  `test-brms-surface-mixed.R`.
- `test-ce-default-effects.R` allows "A Newton step from the stop" on
  one fit alone (`allowed =` per call): `bf(y ~ x + (1 + z | g), sigma ~
  w)`, the fit that warns under OpenBLAS 0.3.32 (`dev/parity2-ob-warn.R`
  under the 0.3.32 copy: 1 warning on that fit, 0 on the eight others).
  The lane had wrongly put the allowance on the nonlinear fit; the
  first OpenBLAS 0.3.26 warning was this fit's too (same 1.49e-06).
- The laplace-draws refusal of a CE whose plan lends a level now says
  `sample_new_levels = "gaussian"` draws the new level from each draw's
  covariance and needs no level's effect (`draws_call_state$snl_lends`,
  set by the CE draws method); `test-laplace-draws.R` asserts the
  sentence.
- NEWS: "on eleven formulas" dropped; the 1497 bullet rewritten: the fit
  bullet is not BREAKING (its default band is unchanged), the draws
  bullet in frmtmb.sample is.
- `dev/upstream-bugs.md`: brms-52, `cumulative()` with `center = FALSE`
  does not compile (stanc, "Ill-typed arguments supplied to function
  "ordered_logistic_glm_lpmf"", the fourth argument the integer 0),
  reproduced in `brms-probe4.txt`. Not posted. No other lane of the
  round had a brms-52 when filed.

### Fix-round runs (lane library, after the fixes)

- Changed and new files plus every core and frmtmb.sample test file
  that names `conditional_effects`, a `loo_` function or
  `sample_new_levels`: 87 files, one per process, gated
  (`FRMTMB_BRMS_FIT_TESTS=true`, `t7/`): pass 10191, skip 1 (the scale
  tier); one failure, `test-effects.R`, whose three-way assertion named
  `Trt:w` and now names brms's `w:Trt` (F3); rerun pass=54 fail=0
  (`t9/`). The direct files (`t8/`, `t10/`): test-ce-new-levels-snl 14
  pass, test-ce-default-effects 13, test-brms-surface-mixed 32, all
  fail 0, warn 0.
- Base arms (rellib-r21, `t12/`): test-ce-default-effects pass=4
  fail=9 (`f * x` gave `f x f:x`), test-ce-new-levels-snl pass=2 fail=2
  err=2 (the width block passes on 0.82.0, whose fit default was already
  gaussian; it failed on the lane build before this round, 0.734),
  test-brms-surface-mixed pass=1 fail=11 err=6.
- OpenBLAS 0.3.32 with LAPACK (`dev/ciharden-openblas.sh 0.3.32
  lapack`, `t-ob332/`, `t10-ob332/`): test-ce-new-levels-snl 14,
  test-ce-default-effects 13 (warn 0 after the allowance moved to the
  fit that warns), test-brms-surface-mixed 32, test-effects 54,
  test-ce-levels 101, frmtmb.sample test-draws-brms-surface 21,
  test-laplace-draws 120, test-postfit-draws 89; fail 0 on all. The copy
  was deleted afterwards.
- One session forward (`sess2-frmtmb.txt`, `sess2-sample.txt`): frmtmb
  files=320/333 pass=24453 fail=0 err=0 warn=10, the ten warnings all
  in test-latres-sparse.R, "cannot create dir ... No space left on
  device" (the shared disk was full for a moment; 144 GB free after);
  that file rerun alone pass=67 warn=0 (`t11/`). frmtmb.sample
  files=87/92 pass=4261 fail=0 err=0 warn=0. Every LEAK line is tagged
  "with load of".
- R CMD check --as-cran: frmtmb "Status: 1 NOTE" (V8 for the HTML
  math; tests OK in 34 min); frmtmb.sample "Status: 1 NOTE" (V8; tests
  OK in 143 s). `dev/parity2-check/` deleted afterwards.

### Files touched in the fix round

`R/conditional-effects.R` (`ce_snl()` default argument, `ce_snl_fit()`,
`ce_numeric_first()`, docs), `R/ce-levels.R` (docs), `R/frame.R` (the
backtick refusal), `NEWS.md`, `dev/test-backlog.md` (the 1497 CLOSED
line, the MEDIUM item), `dev/upstream-bugs.md` (brms-52), `man/`
(roxygenised); tests `test-ce-new-levels-snl.R`,
`test-ce-default-effects.R`, `test-brms-surface-mixed.R`,
`test-effects.R`. frmtmb.sample: `R/methods-draws.R`,
`R/conditional-effects-draws.R`, `tests/testthat/test-laplace-draws.R`.
Scripts: `dev/parity2-brms-probe4.R`, `dev/parity2-snl-width.R`,
`dev/parity2-ob-warn.R`.

## Fix round 2 (short check dev/reviews/2026-11-13-parity2-short.md, FIX)

### F3a: grouping variables only of the terms re_formula keeps

The short check: brms 2.23.0 takes its group variables from
`update_re_terms(formula, re_formula)`, so at the default `re_formula =
NA` there are none; `ce_numeric_first()` counted every grouping
variable as a factor. `y ~ kk * x + (1 | kk)`, integer kk: brms `kk x
kk:x` at NA and `kk x x:kk` at NULL; the lane gave `x:kk` at both.

Change: `ce_kept_group_vars(x, re_formula)` (the grouping variables of
`ce_kept_blocks()`, an `mm()` term's member columns included);
`ce_grids_build()` takes `re_formula` (a new last argument, default NA)
and passes them to `ce_numeric_first()`; the fit method and
frmtmb.sample's draws method pass their `re_formula`. Roxygen of
`ce_numeric_first()` and the sampling-API page state the rule as above.

Test (`test-ce-default-effects.R`, data seed 21, n 200, kk in 1..6):
names at the default and at `re_formula = NULL`. Seen to fail on the
lane build before this round (`t13-prefix/`): "Expected
`names(conditional_effects(fk, resolution = 5))` to be identical to
`c("kk", "x", "kk:x")` ... actual "kk" "x" "x:kk"" (pass=14 fail=1).

### LOW and filed

- The width test's comment now says what the short check measured: over
  30 boot seeds the minimum ratio has mean 0.985, sd 0.086, so 0.8 is
  about 2.1 sd below its mean, 0 of 30 at or below; the fixed seed (2,
  0.961) carries it, and 20 ulp perturbations keep 0.961 on both BLAS.
- Filed in the round's backlog section: MEDIUM (pre-existing), the boot
  band at an OBSERVED level omits the level's conditional variance (55
  to 61 percent of Wald, identical on 0.82.0); LOW, "old_levels" on the
  Wald band by `fitted()`'s route, or a refusal naming `conditions =`.

### Runs (lane library, 2 R processes, gated)

`t14/`: test-ce-default-effects pass=15, test-ce-new-levels-snl 14,
test-effects 54, test-ce-levels 101, frmtmb.sample test-postfit-draws
89, test-draws-brms-surface 21; fail 0, err 0, warn 0 on all.

Files touched: `R/conditional-effects.R` (`ce_kept_group_vars()`,
`ce_grids_build(re_formula =)`, the call site, roxygen),
`R/sampling-api.R` (docs), `man/` (roxygenised),
frmtmb.sample `R/conditional-effects-draws.R` (passes `re_formula`),
`tests/testthat/test-ce-default-effects.R`,
`tests/testthat/test-ce-new-levels-snl.R` (comment),
`dev/test-backlog.md` (two items).
