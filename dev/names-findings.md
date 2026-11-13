# Lane names, round of 2026-11-11: brms names and layouts

Base e8abed0a (frmtmb 0.81.0, `rellib-r20`). Worktree
`frmtmb-wt-names`, library `wt-names-lib`. Brief:
`round-briefs/r082-names.md`. brms 2.23.0 from the user library.

STATUS: done (2026-11-11).

## Instruments

- `dev/names-cases.R`: the twelve cases (data seed 20261111, n = 120),
  one definition for both arms.
- `dev/names-brms.R`: brms 2.23.0 on each case (`chains = 1, iter = 60,
  seed = 1`; names, dims and orders only). Logs
  `dev/names-out/brms-*.txt`.
- `dev/names-brms-code.R`: brms's Stan declaration order, no compile.
  Log `dev/names-out/brms-code.txt`.
- `dev/names-frm.R <base|lane> [sample]`: frmtmb on the same cases, ML
  fit and `frm_sample(chains = 1, iter = 60, seed = 1)`. Logs
  `dev/names-out/frm-{base,lane}[-sample].txt`.

## What brms does (from the logs above)

- `variables()` order is `brms:::reorder_pars()`: by class (`b`, `bs`,
  `bsp`, ..., `sd`, `cor`, `sds`, `sdgp`, `lscale`, the dpars, ...,
  `simo`, `r`, `s`, `zgp`, ..., `meanme`, `sdme`, `corme`), a class `b`
  intercept first, then the Stan declaration order.
- Stan declares the predictors in the FAMILY's dpar order (mu, sigma,
  alpha), but numbers the group ids in the order the dpar formulas are
  WRITTEN. `brms-code.txt`, skew and skew2: with `alpha ~` written
  before `sigma ~`, `r_2_alpha_1, r_3_sigma_1`; written the other way,
  `r_2_sigma_1, r_3_alpha_1`; `b_Intercept, b_sigma_Intercept,
  b_alpha_Intercept` both times.
- `r_` columns: block by block in id order, coefficient-major within a
  block (`r_g[g01,Intercept] ... r_g[g12,Intercept], r_g[g01,ffb]`).
- `me()` hyperparameters: group by group (`meanme_1` is the `gr = g`
  group): `meanme_mexg meanme_mewg meanme_mez`.
- `ranef()` and `coef()`: one entry per grouping factor, whatever
  predictors its terms are on. `coef()` holds every population-level
  coefficient of the model, the group's own first.
- `ranef()` of a model without group-level terms (smooth only, mo only,
  gp only): "The model does not contain group-level effects."
- One-response `posterior_predict()`, `predictive_error()`,
  `predict(summary = FALSE)`, `log_lik()`: dimnames `list(NULL, NULL)`;
  `posterior_epred()` and `posterior_linpred()`: `NULL`.
- `summary()`: `Formula:` then every pform, indented 9 spaces
  (`print.brmsformula()`: pforms then pfix).
- gp draws: `zgp_gpxffa[1..60]`; smooth draws `s_sx_1[1..8]`; mo
  `simo_mok1[1..3]`; latent residuals `err[1..120]` (and no `zerr`).

## The brms vs frmtmb table (generated)

`dev/names-table.R` over the saved runs (brms `dev/names-out/brms-*.rds`,
frmtmb `dev/names-out/frm-{base,lane}[-sample].rds`); output
`dev/names-out/table.md`, pasted verbatim below. "variables() order"
compares the order of the names both sides carry; "ranef()"/"coef()"
compare the entry names and each entry's coefficient columns. Cases: cat
(3158, 3155), skew (3190), dist (3196, 3199), distpop (3196), twoblk
(3205), smooth (3843), mo (3846), gpby (1881, 1884), me (4525), latres
(5087), mv (3848), dblbar (2897).

| case | arm | variables() order | ranef() | coef() |
|---|---|---|---|---|
| cat | base | DIFFERS (6 common) | DIFFERS: brms g [mub_Intercept,muc_Intercept] frmtmb g,g [mub_Intercept; muc_Intercept] | DIFFERS: brms g [mub_Intercept,muc_Intercept,mub_x,muc_x] frmtmb mub,muc [g; g] |
| cat | lane | same (6 common) | same | same |
| cat | base draws | DIFFERS (31 common) | same | same |
| cat | lane draws | same (31 common) | same | same |
| skew | base | DIFFERS (6 common) | DIFFERS: brms h,g,ga [Intercept; alpha_Intercept; sigma_Intercept] frmtmb h,ga,g [Intercept; sigma_Intercept; alpha_Intercept] | DIFFERS: brms h,g,ga [Intercept,sigma_Intercept,alpha_Intercept; alpha_Intercept,Intercept,sigma_Intercept; sigma_Intercept,Intercept,alpha_Intercept] frmtmb mu,sigma,alpha [h; ga; g] |
| skew | lane | same (6 common) | same | same |
| skew | base draws | DIFFERS (33 common) | DIFFERS: brms h,g,ga [Intercept; alpha_Intercept; sigma_Intercept] frmtmb h,ga,g [Intercept; sigma_Intercept; alpha_Intercept] | DIFFERS: brms h,g,ga [Intercept,sigma_Intercept,alpha_Intercept; alpha_Intercept,Intercept,sigma_Intercept; sigma_Intercept,Intercept,alpha_Intercept] frmtmb h,ga,g [Intercept,sigma_Intercept,alpha_Intercept; sigma_Intercept,Intercept,alpha_Intercept; alpha_Intercept,Intercept,sigma_Intercept] |
| skew | lane draws | same (33 common) | same | same |
| dist | base | DIFFERS (6 common) | DIFFERS: brms g [Intercept,sigma_Intercept] frmtmb g,g [Intercept; sigma_Intercept] | DIFFERS: brms g [Intercept,sigma_Intercept,x,sigma_z] frmtmb mu,sigma [g; g] |
| dist | lane | same (6 common) | same | same |
| dist | base draws | DIFFERS (31 common) | same | same |
| dist | lane draws | same (31 common) | same | same |
| distpop | base | DIFFERS (3 common) | brms stops, frmtmb answers | brms stops, frmtmb answers |
| distpop | lane | same (3 common) | both stop | brms stops, frmtmb answers |
| distpop | base draws | DIFFERS (4 common) | both stop | both stop |
| distpop | lane draws | same (4 common) | both stop | both stop |
| twoblk | base | same (5 common) | DIFFERS: brms g [Intercept,x] frmtmb g,g [Intercept; x] | same |
| twoblk | lane | same (5 common) | same | same |
| twoblk | base draws | DIFFERS (30 common) | same | same |
| twoblk | lane draws | same (30 common) | same | same |
| smooth | base | same (4 common) | brms stops, frmtmb answers | brms stops, frmtmb answers |
| smooth | lane | same (4 common) | both stop | brms stops, frmtmb answers |
| smooth | base draws | DIFFERS (5 common) | both stop | both stop |
| smooth | lane draws | same (13 common) | both stop | both stop |
| mo | base | same (3 common) | brms stops, frmtmb answers | brms stops, frmtmb answers |
| mo | lane | same (6 common) | both stop | brms stops, frmtmb answers |
| mo | base draws | same (7 common) | both stop | both stop |
| mo | lane draws | same (7 common) | both stop | both stop |
| gpby | base | same (6 common) | brms stops, frmtmb answers | brms stops, frmtmb answers |
| gpby | lane | same (6 common) | brms stops, frmtmb answers | brms stops, frmtmb answers |
| gpby | base draws | DIFFERS (7 common) | both stop | both stop |
| gpby | lane draws | same (127 common) | both stop | both stop |
| me | base | DIFFERS (12 common) | brms stops, frmtmb answers | brms stops, frmtmb answers |
| me | lane | same (12 common) | both stop | brms stops, frmtmb answers |
| me | base draws | DIFFERS (13 common) | both stop | both stop |
| me | lane draws | same (13 common) | both stop | both stop |
| latres | base | same (4 common) | brms stops, frmtmb answers | brms stops, frmtmb answers |
| latres | lane | same (4 common) | both stop | brms stops, frmtmb answers |
| latres | base draws | same (5 common) | both stop | both stop |
| latres | lane draws | same (5 common) | both stop | both stop |
| mv | base | DIFFERS (8 common) | DIFFERS: brms g [y_Intercept,y2_Intercept] frmtmb g,g [y_Intercept; y2_Intercept] | DIFFERS: brms g [y_Intercept,y2_Intercept,y_x,y2_z] frmtmb y_mu,y2_mu [g; g] |
| mv | lane | same (8 common) | same | same |
| mv | base draws | DIFFERS (33 common) | same | same |
| mv | lane draws | same (33 common) | same | same |
| dblbar | base | same (5 common) | same | DIFFERS: brms g [Intercept,ffb,x] frmtmb g [Intercept,x,ffb] |
| dblbar | lane | same (5 common) | same | same |
| dblbar | base draws | DIFFERS (30 common) | same | same |
| dblbar | lane draws | same (30 common) | same | same |

Every remaining lane difference is a documented frmtmb divergence (filed
under "Filed in the round of 2026-11-11"): `coef()` of a fit without
group-level terms returns the coefficient vectors, and an ML gp() fit
answers `ranef()`/`coef()` with the field. "brms stops, frmtmb answers"
on the base ML rows for `ranef()` is 3843.

## What changed, item by item

Backlog line numbers are the item's first line at e8abed0a.

### 1. variables() and draws order: 3158, 3190, 3196, 3848, 4525, 1884

- Core `brms_variable_sort(fit, names)` (R/brms-names.R) is brms's
  `reorder_pars()` on frmtmb's brms names: class rank
  (`brms_class_levels`, dpars ranked in the families' order), a class
  `b` intercept first, then `r_` keys (block in brms id order,
  coefficient, level) and the `me()` group order, else the input order.
- `brms_block_order()` / `brms_block_keys()`: brms's group-id order,
  response by response, then the dpar/nlpar formulas in WRITTEN order
  (`fit$bform`), `mu` first; ties keep the frame order.
- `hyp_env_vals()` walks the blocks in that order, so `sd_` and `cor_`
  come out in brms's order for `variables()` AND for the draws' derived
  `sd_`/`cor_` columns. Smooth `sds_` keep the frame order.
- `me_hyper_values()` lists the calls group by group (`me$groups[[i]]$K`).
- `variables.frmtmb_fit()` returns the sorted names. Draws:
  `variables.frmtmb_draws()`, `draws_raw_array()` (behind every
  `as_draws_*()`, `as.matrix()`, `as.array()`, `as.data.frame()`,
  `posterior_summary()`, `mcmc_plot()`, `pairs()`) and `as.mcmc(combine_chains
  = TRUE)` permute through `draws_brms_order()`. `x$draws` keeps the
  template order, which every positional reader uses.
- 1884 (gp `sdgp` order in the draws) falls out of the class sort.
- 3196's population order and 3158's categorical order fall out of the
  intercept-first key.
- The sort is exported to frmtmb.sample through the sampling API
  (`R/sampling-api.R`, rawNamespace and alias).

### 2. simo_ in an ML fit's variables(): 3846 (= 4072, 4373)

`hyp_vals_only()` and `hyp_par_cov()` carry the `zeta<j>` components
(`mo_extra_comps()`), and `hyp_env_vals()` puts the weights
(`mo_put_simplex()`), so `variables()` lists `simo_mok1[1..3]` and
`hypothesis()` reads them with the delta method. The draws already had
them; `draws_cov_cols()` skips them (they are in its `have`).

### 3. gp draws: 1881, 1884

`brms_basis_labels()`: an exact `gp()` block's `b` columns are
`zgp_<prefix>_<term>[i]` (brms's `rename_gp()` base, `sdgp_` with `zgp_`).
The values are z in the block's farthest-first order, so they are not
elementwise brms's `zgp` (brms factors in first-appearance order); the
docs say so. `gp(k = )` keeps `b[i]`: its weights carry the spectral
density (`covstruct_core$hsgp$nll` is `dnorm(b, 0, hsgp_sds())`), brms's
`zgp` is standardized.
Also: a smooth's basis weights are now `s_<label>_<k>[i]` (brms's
`s_sx_1[i]`), which 3843 records; frmtmb's smooth `b` is the coefficient
on the same `smooth2random(type = 2)` basis brms uses (homdiag density,
`dnorm(b, 0, sds)`).

### 4. ranef()/coef(): 3155, 3199, 3205, 3843

- `ranef()`: blocks sorted into brms id order, then
  `ranef_merge_groups()` puts every block on one factor into one matrix
  (attributes `term` = the block labels, `term_ncol` = their column
  counts). Blocks whose levels differ, or whose columns would share a
  name, stay separate. `ranef_pick()` returns one block's columns for a
  block label; `print()` and `as.data.frame()` read the per-column
  labels.
- smooth blocks are no longer in `ranef()`, and a model with no
  group-level block stops with brms's message (3843).
- `coef()`: `coef_merge_groups()` builds one data frame per group with
  every `fixef()` coefficient (the group's own columns first, then fixef
  order); the per-dpar nesting is gone.
- NOT done: brms's `ranef()`/`coef()` of a brmsfit are `levels x 4 x
  coefs` summary arrays. The ML methods keep the lme4 `levels x coefs`
  shape (and `condVar`), as documented; the draws methods already return
  brms's arrays.

### 5. Others

- 4573: `pp_dimnames()` in core `predict(summary = FALSE)`; frmtmb.sample
  `posterior_predict()` (one response; so also `predict(summary =
  FALSE)` and `predictive_error()`), and `log_lik()`. `posterior_epred()`
  and `posterior_linpred()` stay `NULL`, as brms's are.
- 5087: the states of `ar(cov = FALSE)` (latres kind `condar`) are named
  `zerr[k]` (`zerr_<resp>[k]` in a multivariate model): one standard
  normal per row in row order, `latres_lpdf()` `-0.5 * sum(b * b)` and
  `latres_expand()` `s * (b + lag_sum)`, which is brms's `zerr`. Other
  latent kinds (ar1 under cov = TRUE, ma, arma, cosy, unstr) keep
  `latent_state[k]`: their states are not brms's `zerr`. NOT done:
  brms's derived `err[k]` columns (a content change, not a name), the
  corrected estimate, `get_joint_cov()` on a new row and M7, which the
  bundle also lists.
- 2897 (= 3203): `confint_varcorr()$block` shows a `||` term as written
  (`block_display_label()`, from the components' `written_text`).
  `term_label` stays the internal key. NOT changed: `VarCorr()` and
  `as.data.frame(VarCorr())`, which key by brms group, and diagnose()'s
  singular table.
- 1865: `insight_param_scale()` maps every component's thresholds when
  several read one block (order = "mu"), with a rectangular delta-method
  map. The sum-to-zero `Intercept` row is the decision of 1832, not
  done.
- 1543: `summary_formula_lines()` and `cat_formula_lines()`: brms's
  `print.brmsformula()` lines (formula, pforms, nlforms, pfix), indent 9.
  Also the `frm_multiple()` summary.
- 1550 (user decision, coordinator 2026-11-11): `data_name_of()`; a
  value is `<class: rows x cols>` (`tbl_df` shown as `tibble`), no
  deparse. `update(newdata = )` goes through the same line in `frm()`.
- 1279 (user decision): CLOSED line only, no code change.

### 6. 3164 (multinomial/categorical K, unused levels, epred margin)

brms 2.23.0 (`dev/names-brms-3164.R`, `dev/names-out/brms-3164.txt`,
seed 3164, n = 90): a factor with levels `a b c unused` gives `ncat = 3`
and `b_mub_Intercept b_muc_Intercept b_mub_x b_muc_x`; a named
`cbind(y1, y2, y3)` multinomial declares `muy2`, `muy3`;
`posterior_epred()` of a categorical fit has dimnames `list(NULL, NULL,
c("a", "b", "c"))`; `posterior_predict()` `list(NULL, NULL)`.

frmtmb (`dev/names-frm-3164.R`, `dev/names-out/frm-3164-{base,lane}.txt`):

| case | base (rellib-r20) | lane |
|---|---|---|
| categorical, unused level | `b_mub_* b_muc_* b_muunused_*` (6 coefs) | `b_mub_Intercept b_muc_Intercept b_mub_x b_muc_x` |
| categorical(K = 3), factor a b c | `b_mu2_* b_mu3_*` | as categorical() |
| multinomial(K = 3), named cbind | `b_mu2_* b_mu3_*` | `b_muy2_* b_muy3_*` |
| draws posterior_epred dimnames | `NULL, "1".."90", a b c` | `NULL, NULL, a b c` |
| draws posterior_predict dimnames | `NULL` | `list(NULL, NULL)` |

Changes: `categorical_y_levels()` drops unused factor levels;
`categorical(K =)`, `multinomial(K =)`, `dirichlet_multinomial(K =)`
carry a `defer` slot (`k_defer()`) that takes the data's names when the
data have `K` categories, and otherwise leaves the K family; the draws'
polytomous `posterior_epred()` leaves the observation margin `NULL`
(ordinal too: the gated test `test-sampling-ported.R` "ordinal
posterior_epred has brms's draws x obs x category shape" now asserts
`identical(dimnames(ep), dimnames(bep))` against brms; it read "brms
... leaves the observation margin unnamed" before). NOT done:
`insight::get_predicted()` on `multinomial()` (filed).

Tests that pinned the old K naming and were changed:
`test-refcat.R` (`fk` now `muy2 muy3`), `test-famgaps.R` (one-hot
multinomial(K = 3) with columns `a b c` is `mub muc`; a factor with one
used level is now "fewer than two categories"),
`test-simulate-density.R` (two `K =` tests now use unnamed columns so
the K family's `mu2`, `mu3` still apply).

## Tests

New files: `tests/testthat/test-brms-layout.R` (core, no extension),
`extensions/frmtmb.sample/tests/testthat/test-brms-draws-order.R`.

### Seen failing on rellib-r20 (base), passing on the lane

`dev/release/names-newbase.log` (base, one file per process) against
`dev/release/names-rerun3.log` (lane):

- base: `RESULT test-brms-layout.R pass=9 fail=27 err=4`; lane
  `pass=57 fail=0 err=0`. Behavioral failures on base, one or more per
  item: the variables() order (dist, distpop, skew, mv, categorical,
  me), `names(ranef())` `h ga g` (skew), `ranef()` refusing `$g` on a
  distributional, categorical and mv fit, no simo_ in variables(),
  `ranef()` answering on a smooth-only and a no-RE fit, `dimnames(p)`
  NULL, the Formula line without `sigma ~ z + (1 | g)`, the empty
  `Data:` line under `do.call()` and `update()` by value, `1 + f | g`
  in confint_varcorr(), `tau_raw` in insight's parameters, the unused
  level kept and `mu2 mu3` under K.
- base: `RESULT test-brms-draws-order.R pass=11 fail=15 err=0 warn=1`;
  lane `pass=28 fail=0 err=0`: the draws' order, `b[i]` for the smooth
  and gp, `latent_state` for ar(cov = FALSE), NULL dimnames.

### Existing tests changed (they pinned the old layout)

core: test-brms-names.R (mv order), test-brms-parity-defects.R
(data_name of a frame by value), test-methods.R (ranef blocks, two
tests rewritten, one added for the still-separate case with
`ranef_merge_groups()` on hand-built matrices), test-group-labels.R
(`coef(fit)[["g1:g2"]]`), test-latres.R (ranef() of a latres-only fit
stops), test-review-v29.R and test-sugar.R (merged ranef), test-refcat.R,
test-famgaps.R, test-simulate-density.R (K naming, unused level).
frmtmb.sample: test-latres.R (`zerr`), test-brms-output.R (the brms shim
is built from the draws in variables() order, as brms stores them),
test-draws-methods.R and test-sample-direct.R (accessors in variables()
order, `x$draws` setequal), test-rescor-draws.R (`unname()` both sides),
test-sampling-ported.R (epred observation margin NULL; the gated brms
comparison asserts `identical(dimnames(ep), dimnames(bep))`).
frmtmb.ode: test-ode.R (`ranef(fit)$id[, "lka_Intercept"]`; the lka and
lke blocks on id are one matrix now).

### Every suite, one file per process (lane library)

`dev/release/names-all.log`, NOT_CRAN=true, 3 at a time, every `lib:`
line `frmtmb 0.81.0 C:/Users/adf44/source/r/wt-names-lib` (512 of 512):

| package | files | pass | fail | err | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 318 | 23550 | 0 | 0 | 239 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 687 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 548 | 0 | 1 | 1 | 0 |
| frmtmb.sample | 88 | 4079 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 556 | 0 | 0 | 0 | 0 |

The one ode error was test-ode.R reading `ranef(fit)[[1L]]` as one
block; fixed in the test, rerun below.

### Gated tier (FRMTMB_BRMS_FIT_TESTS=true), lane library

`dev/release/names-gated.log`: every gated file of every package (35)
plus frmtmb.ode's test-ode.R after its fix, 36 of 36. Results other than
fail=0 err=0:

- test-brms-methods.R fail=1 err=1: two assertions pinned the old
  layout (frmtmb's coef() with the group's own two columns, and
  `ranef()` of a no-RE fit as an empty list). brms's own output in the
  same test is now what frmtmb gives; updated to compare against it.
  Rerun `dev/release/names-gated2.log`: pass=1013 fail=0 err=0.
- test-sampling-ported.R fail=1: my first strengthening asserted
  `identical(dimnames(ep), dimnames(bep))`, and brms names an ordinal
  epred's DRAW margin `"1".."S"` (filed). Narrowed to the observation
  margin; rerun `dev/release/names-gated3.log`: pass=215 fail=0 err=0
  skip=0. So brms's ordinal observation margin is `NULL`, as the lane
  now gives.
- test-drmtmb-agreement.R skip=13: its own gate,
  FRMTMB_DRMTMB_FIT_TESTS, not set (it reads drmTMB, not names).
- test-ode.R: pass=72 fail=0 err=0.

### Equivalence pairs remeasured

No reader in `dev/equiv-*.R` changed; the harness reads `ranef()`
through helper-brms.R's `brms_ranef_block()`, which pools the blocks of
a group itself, so the merged `ranef()` reaches it. Two pairs with
group-level terms on several predictors, remeasured with
`FRMTMB_EQUIV_LIB="<lane lib>;rellib-r20"` into `dev/names-equiv-out`:
p139 (multinomial x us) equivalent L1 "every level agrees", the record's
verdict; p129 (cratio x mvbf) equivalent L4 "every level agrees", the
record's verdict.

### Not run

- OpenBLAS (`dev/ciharden-openblas.sh`): the change is names, orders and
  layouts; the only arithmetic touched is insight's delta-method map for
  order = "mu" (a rectangular `A V A'`, same formula as before) and no
  exactness test changed.
- One-session runs (`run-sessions.sh`): see the end.

## R CMD check --as-cran (once per changed package)

`bash dev/names-check.sh <pkg>`, built and checked in
`dev/names-check/<pkg>/`, R_LIBS = lane lib, rellib-r20, user library,
pandoc and TinyTeX on PATH, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, no
`--no-manual`:

- frmtmb: `Status: 1 NOTE` (HTML manual: "Skipping checking math
  rendering: package 'V8' unavailable", the expected one).
- frmtmb.sample: `Status: 1 NOTE` (the same V8 note).

frmtmb.ode changed only a test file; no check run for it (its one-file
run of test-ode.R is pass=72 fail=0 err=0 above).

## Bump and floor

Minor for frmtmb (BREAKING layout changes: ranef(), coef(), variables()
order, VarCorr() group order, K families' names, categorical()'s unused
level) and minor for frmtmb.sample (BREAKING draws order and names). The
new frmtmb.sample code calls core's `brms_variable_sort()` (sampling
API), so frmtmb.sample's floor is the frmtmb version this round
releases. frmtmb.ode needs no bump (test only).

## One session per suite, three ways

`bash dev/release/run-sessions.sh names-sess 3 frmtmb frmtmb.sample`
(`dev/release/names-sess.log`), lane library, NOT_CRAN=true:

    SESSION frmtmb forward files=305/318 pass=23550 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb reverse files=305/318 pass=23715 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb preload files=305/318 pass=23818 fail=0 err=0 skip=241 warn=0
    SESSION frmtmb.sample forward files=83/88 pass=4079 fail=0 err=0 skip=18 warn=0
    SESSION frmtmb.sample reverse files=83/88 pass=4079 fail=0 err=0 skip=18 warn=0
    SESSION frmtmb.sample preload files=83/88 pass=4079 fail=0 err=0 skip=18 warn=0
    NAMES-SESS ran 6 of 6

Every LEAK line carries `[with load of ...]` (0 untagged). The seven
other extensions' sessions were NOT run (time); their one-file runs are
in the table above, and the lane adds no registry entry and no test
that registers anything.

## Files touched

R: R/brms-names.R, R/confint.R, R/methods-fit.R, R/me.R, R/insight.R,
R/families.R, R/fit.R, R/multiple.R, R/predict-brms.R, R/autocor.R,
R/gp-by.R, R/sampling-api.R; NAMESPACE and man/ by roxygen
(coef.frmtmb_fit.Rd, ranef.Rd, frmtmb-sampling-api.Rd, frmtmb-gp.Rd,
frmtmb-families.Rd, frmtmb-autocor.Rd); NEWS.md.
frmtmb.sample: R/draws-brms.R, R/methods-draws.R, R/loo.R, NEWS.md.
Tests: listed under Tests. dev/test-backlog.md (CLOSED lines, the 1279
decision line, the new section). dev/names-* scripts and this file.

Outside the lane's obvious area: R/predict-brms.R (ppdraw's; one hunk,
`pp_dimnames()` on the one-response `summary = FALSE` return),
frmtmb.sample R/methods-draws.R and R/loo.R (sampfront's; hunks:
`variables()`, `draws_brms_order()`, the summary's row order, the
as.mcmc combine path, one-response pp dimnames, the polytomous epred
margin, log_lik's matrix dimnames), R/fit.R (one hunk: data_name),
R/families.R (K defer, unused level), R/insight.R. R/priors.R not
touched (no printed prior name changed).

## Fix round (review dev/reviews/2026-11-11-names.md, verdict SHIP)

The coordinator asked for fixes to minors 1 to 8 before consolidation;
minor 9 (line endings) is left to consolidation. The user's three open
questions (2-D ML ranef(), coef()/gp ranef() where brms stops, the
`x$draws` order) are unchanged.

1. `pars =` regex order: `draws_accessor_args()` (draws-brms.R) and
   `draws_select_variables()` (methods-draws.R) match against
   `variables(x)`, as brms's `extract_pars(pars, variables(x))`.
   `as.matrix()`, `as.array()`, `as.data.frame()`, `as.mcmc()`,
   `posterior_interval()`, `mcmc_plot()`, `pairs()` and `plot()` with
   `pars =` now select in variables() order. Test added to
   test-brms-draws-order.R ("variables() and the draws objects list
   brms's order": `as.matrix`, `as.array`, `as.data.frame`,
   `posterior_interval` and `as.mcmc` with `pars = "^b_"`). The failing
   case is the review's reproducer (`names-review-draws.log`, every
   shape: "as.matrix(pars = regex) in variables order FAIL").
2. `||` shown one way: `frm_multiple()`'s `pooled_varcorr$grp` goes
   through `block_display_label()` (R/multiple.R), so it and the
   `frm_multiple()` summary read `1 + f || g`, as `confint_varcorr()`
   does. Test added (test-brms-layout.R, the confint_varcorr test:
   pooled grp, the summary print with no `1 + f | g`, and a
   `gr(g, cor = FALSE)` block shown as written). Failing case: the
   review's `names-review-misc-lane.log`.
3. NEWS: BREAKING on the confint_varcorr()/pooled block labels
   (`||` and `gr(cor = FALSE)`), the insight change, and frmtmb.sample's
   epred margin; new bullets for `ranef()` dropping a smooth beside
   group terms and for the draws `summary()` row order; the order bullet
   names `posterior_interval()` and `pars =` and says to select
   `x$draws` by name.
4. dev/brmsport-bin2-verdicts-2t.R: the reasons of brmsfit-methods:136
   and :140 describe one data frame per group with every fixef()
   coefficient. Verdicts and classes unchanged; no count is stated that
   was not measured.
5. vignettes/brms-migration.Rmd: the gp draws are `zgp_<term>[i]`, not
   elementwise brms's. Also frmtmb.sample's vignettes/sampling.Rmd said
   "the same `b[i]` columns" of non-centered draws; it now names the
   group-level columns.
6. The 1279 CLOSED line uses `(lane names): by user decision, ...`.
7. test-brms-layout.R: the simplex sum is compared with 1 (relative);
   the other new tolerances are plain relative `tolerance =` values,
   no longer scaled by a target.
8. ?categorical: with data, K does not fix the names; formulas and
   priors take `mub`, `muc`; `mu2` is refused for factor data; numeric
   codes 2, 3, 4 give `mu3`, `mu4`.
9. ?frm_sample `@return`: `x$draws` is in an internal order, not
   variables()'s; select its columns by name.

Runs (lane library, NOT_CRAN=true, one file per process,
`dev/release/names-fix.log`, 139 of 139): frmtmb.sample's whole suite,
files=88 pass=4087 fail=0 err=0 skip=18 warn=0; core
test-brms-layout.R plus every core file that names confint_varcorr,
frm_multiple or categorical(, files=51 pass=7101 fail=0 err=0 skip=71
warn=0. test-brms-layout.R pass=61, test-brms-draws-order.R pass=36,
both fail=0 err=0.

R CMD check --as-cran after the fix round (`bash dev/names-check.sh`,
same settings as before): frmtmb `Status: 1 NOTE`, frmtmb.sample
`Status: 1 NOTE`; both NOTEs are the HTML manual's "Skipping checking
math rendering: package 'V8' unavailable".

Files touched in the fix round: frmtmb.sample R/draws-brms.R,
R/methods-draws.R, R/sample.R (docs), man/frm_sample.Rd, NEWS.md,
vignettes/sampling.Rmd, tests test-brms-draws-order.R; core
R/multiple.R, R/families.R (docs), man/frmtmb-families.Rd, NEWS.md,
vignettes/brms-migration.Rmd, tests test-brms-layout.R;
dev/brmsport-bin2-verdicts-2t.R, dev/test-backlog.md, this file.
