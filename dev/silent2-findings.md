# Lane silent2: small silent defects

Round opened 2026-11-13, base ec8178d3 (frmtmb 0.82.0), branch
`wt-r83-silent2`, uncommitted. Before: `rellib-r21`. Lane library:
`C:/Users/adf44/source/r/wt-silent2-lib`. Scripts `dev/silent2-*.R`; logs
in `dev/silent2-out/` (local, gitignored). Every number below is pasted
from the log named beside it.

Status: DONE (2026-11-13), then a fix round (section at the end).
Everything is uncommitted.

## Backlog numbers

The brief's items 5 and 6 carry the numbers 5347 and 5353, but at
ec8178d3 those lines start other items. By description:

- "prediction design under the session's contrasts" is backlog 5353;
- "laplace draws, `conditional_effects()` NA band" is backlog 5362;
- 5347 is "`log_lik(newdata =)` of draws on a subset computes every fitted
  column", which the environment file gives to lane parity2 ("loo
  subset"). This lane did not touch it.

The descriptions decide; the CLOSED lines go under 5353 and 5362.

## Reproducers before and after

`dev/silent2-repro.R base|lane`, logs `dev/silent2-out/repro-base.txt`
and `repro-lane-final.txt`. Base (rellib-r21, 0.82.0):

```text
5235 CE predict, resolution 120: ran
5235 CE mean, resolution 120: ran
5235 predict(newdata = 120 rows): ran
  rows 1:2 ratio 6.59135 (x nearly equal)
5235 in-sample fitted(): ran
5246 anova(rows 1..72, rows 9..80): ran
5246 control anova(rows 1..72, rows 1..72): ran
5246 anova(y ~ x, y ~ z), NA in different rows: ran [1 warnings: anova(): the fixed effects of y ~ x and y ~ z are not nested ...]
5246 anova(y ~ x on d, y ~ x + z on rev(y)): ran
5246 control anova(y ~ x, y ~ x + z) on one frame: ran
5353 predict after options() reset: ERROR: Random-effect design for `1 + f | g` does not match the fitted model (columns: (Intercept), fb, fc vs (Intercept), f1, f2)
5362 CE predict, laplace, sigma ~ (1 | h): ran [6 warnings: NAs produced]
  NA cells 100 of 100
5362 posterior_predict on the same draws: ERROR: posterior_predict() needs the random effects, ...
5362 control CE predict, laplace, re_formula = NA: ran
5362 control CE predict, full draws, re_formula = NULL: ran
5362 control CE predict, laplace, mu ~ (1 | h): ERROR: conditional_effects() needs the random effects, ...
print order: variables() sdgp_gpxfa sdgp_gpxfb lscale_gpxfa lscale_gpxfb
print order: print()      sdgp_gpxfa lscale_gpxfa sdgp_gpxfb lscale_gpxfb
```

Lane, the same lines on the final build (`repro-lane-final.txt`):

```text
5235 CE predict, resolution 120: ERROR: Addition term vreal(rep_len(c(1, 3), 120)) could not be evaluated on newdata: it reads no column of newdata, so its values are not newdata's rows. Write it as a column of the data
5235 CE mean, resolution 120: ERROR: conditional_effects() could not evaluate vreal(rep_len(c(1, 3), 120)) on the effect grid: it reads no column of the data, so its values are not the grid's rows. Write it as a column of the data
5235 predict(newdata = 120 rows): ERROR: (as the first line)
5235 in-sample fitted(): ran
5246 anova(rows 1..72, rows 9..80): ERROR: anova() needs fits of the same data: model 2's likelihood reads other response values than model 1's, ...
5246 control anova(rows 1..72, rows 1..72): ran
5246 anova(y ~ x, y ~ z), NA in different rows: ERROR: anova() needs fits of the same data: ...
5246 anova(y ~ x on d, y ~ x + z on rev(y)): ERROR: anova() needs fits of the same data: ...
5246 control anova(y ~ x, y ~ x + z) on one frame: ran
5353 predict after options() reset: ran
  max |diff| against the prediction under the fit's contrasts 0
5362 CE predict, laplace, sigma ~ (1 | h): ERROR: conditional_effects() needs the random effects, and these draws come from frm_sample(laplace = TRUE), ...
5362 control CE predict, laplace, re_formula = NA: ran
5362 control CE predict, full draws, re_formula = NULL: ran
print order: variables() sdgp_gpxfa sdgp_gpxfb lscale_gpxfa lscale_gpxfb
print order: print()      sdgp_gpxfa sdgp_gpxfb lscale_gpxfa lscale_gpxfb
```

The CE line through `method = "predict"` reaches `aterms_for_newdata()`
through `ce_sim_one()`'s caller, so it carries that function's words.

## Item 5235: an addition term of length n by coincidence

`R/predict.R`: new `aterm_rowless(ex, v, newdata)`: a value with more than
one row whose expression names no column of newdata. In
`aterms_for_newdata()` such a value is dropped, so a term the family needs
is refused, with the reason "it reads no column of newdata, so its values
are not newdata's rows. Write it as a column of the data"; any other term
gets the omission warning with "(it reads no column of newdata)".
`R/conditional-effects.R`, `ce_aterms()`: the same test on the grid, with
"it reads no column of the data, so its values are not the grid's rows".

A value of length 1 (`trials(10)`, `trunc(lb = 0)`) is unaffected. The
fitted rows (`newdata = NULL`) read the frame's stored values and are
unaffected: test "a payload whose length matches the new rows by chance
is refused" asserts `fitted(fe) == mu * s` on them. A payload written as
a column follows newdata (the absent case, same test).

## Item 5246: `anova()` compares rows, not their count

`R/frame.R`: new `loglik_rows_key(fit)`: per response (sorted by name),
the response values each density reads (`frame$y`, less `mi_map` unread
rows, through `subset_rows`), ordered by the data's row names (numeric
when every name is a number, else radix order). `R/confint.R`,
`anova.frmtmb_fit()`: after the count check, the keys must be identical,
else "anova() needs fits of the same data: model k's likelihood reads
other response values than model 1's, over the same number of
observations (other rows, or the rows of another data set). Fit the
models to the same rows".

Why values in row-name order and not the row names alone: row names alone
would refuse a fit to `d[complete.cases(d), ]` of a tibble (names restart
at 1) against a fit to `d` that dropped the same rows itself; values
alone (unordered) would pass a bernoulli 1..72 against 9..80 whenever the
counts of ones agree. The test asserts both absent cases: the same rows in
reverse order and under names reset to 1 compare (`anova(a, b, c0)`, and
its chisq equals `2 (logLik(b) - logLik(c0))`).

Not changed: `multiple.R`'s pooled `anova()` (counts per imputation),
which compares fits of the m imputed data sets by design.

## Item 5353: group-level design under the session's contrasts

`R/frame.R`: new `re_term_contrasts(lhs, data, env)`, the contrasts
`model.matrix()` gives the factors of a bar's left side at the fit (a
factor's own `contrasts` attribute, else `options("contrasts")`), stored
on each plain component (`contrasts`, carried through the block builder)
and on an `mm()` term's spec (`mms$contrasts`). New
`re_contrasts_for(contr, mf)` keeps the entries `model.matrix()` can apply
to a model frame. Readers: `pred_design()` (`R/predict.R`, the backlog's
line), `mm_member_designs(use_model_frame = TRUE)` (the `mm()` route) and
`re_term_cnms()` (`R/re-formula.R`, a `re_formula` naming part of a term).

Test "a group-level design is rebuilt under the fit's contrasts"
(test-predict-newdata.R): a fit under `contr.sum`, predicted after the
option is reset, equals its prediction under `contr.sum`, for the plain
term, a partial `re_formula` and an `mm()` term; a fit under the default
predicts the same under either option (the absent case).

## Item 5362: laplace draws, CE predictive band

`extensions/frmtmb.sample/R/conditional-effects-draws.R`: on laplace
draws with `method = "posterior_predict"` the watch reads the expected
values AND every predicted dpar's link-scale predictor on the grid
(`frm_linpred(type = "link")`, one leaf per dpar, cells in the expected
values' order so the design skip applies cell by cell), as
`posterior_predict()`'s watch reads its dpars. Through a level plan the
dpars go through `ce_plan_eval(one = )` and share the new-level cache.

Test "a predictive band reads the group-level terms of every dpar"
(test-laplace-draws.R): gaussian and student `sigma ~ (1 | h)` and
negbinomial `shape ~ (1 | h)` refuse with the re_formula = NA hint; at
re_formula = NA the laplace band equals the full draws' band bit for bit;
the full draws' band is finite. The existing case "(0 + x | g) at x = 0"
(design skip) still passes, so the skip alignment holds.

## Unfiled: `print()` of draws in stored order (filed and closed)

`extensions/frmtmb.sample/R/sample.R`, `print.frmtmb_draws()`: the rows in
`variables()` order (`draws_brms_order()`). Test "print() lists the
parameters in variables() order" (test-brms-draws-order.R), on draws
built from an ML fit of `gp(x, by = f, k = 5)` (no sampler).

## Items 5207, 5167, 5214: ordinal-mixture detection

### The measurement

`dev/silent2-ordmeasure.R <arm> <set>` fits the controls with the
default `frm()` call (jittered starts on): the 40 of
`dev/ordstart-starts.R` (cum_cum and cum_sratio, seeds 1 to 20, n = 300)
and the review's 145 (`dev/reviews/ordstart-review-sweep.R`, sets sweep
104, steep 38, wide 3). Per fit it records the degenerate rule's numbers,
the largest threshold standard error per component from `fixef()`
(natural units), the smallest gap between ordered thresholds, a shift
test (the smallest over the parameters of the larger log-likelihood
change of a shift of +-10) and a per-threshold sharpening (each
threshold's latent distances doubled alone). Summary
`dev/silent2-ordsum.R base`, log `dev/silent2-out/ordsum-base.txt`:

```text
fits 185 by set: forty 40, steep 38, sweep 104, wide 3
degenerate warning on 0.82.0: 21
passed (no degenerate warning, code 0): 163
SE cut 50: passed fits flagged 20 of 163
SE cut 100: passed fits flagged 20 of 163
SE cut 500: passed fits flagged 20 of 163
SE cut 5000: passed fits flagged 9 of 163
largest coefficient SE on passed fits not flagged at 100: cs_sr_acat_seed3 se_t 16.43 se_b 43.56; theta_formula_seed4 se_t 27.21 se_b 36.77
```

The 20 flagged at 50, 100 and 500 alike (cum_cum 16; cum_sratio 5, 6,
14, 16, 18; sr_cr_cloglog 5, 6; acat_acat_K3 1, 5; dpar_formula 1 to 4;
cr_cr_probit_K6 1, 4; acat_cumcll_K7 3, 6; cc_logit_K7 6; steep_n300_s3
K6 seed 2) have a largest threshold SE of 640 (cum_cum 16) to 3.8e33
(the probit's underflow). The other 143 have at most 27.21
(theta_formula seed 4). Nothing lies between 27.21 and 640, so any cut
from 28 to 639 gives the same verdicts on these fits. 100 sits in the
gap, 3.7 times the clean group's top and 6.4 times below the flagged
group's bottom, and a threshold with a standard error of 100 latent
units is placed by no reading of the data. ordstart's suggested 5e3
flags 9 of the 20 and misses 11, among them cum_cum 16 (640), cum_sratio
5 and 6 (4620, 4345) and one of its own six flat points (dpar_formula 4,
4443).

ordstart's six flat points and the 5207 fit on the base arm
(`ordsum-base.txt`): cum_sratio 14, 18, dpar_formula 4, acat_cumcll_K7 6
and cum_sratio 16 have threshold SEs 5053, 5200, 4443, 8193 and 1.75e7
and are flagged; cs_sr_acat 1 has a component-2 SE that fixef() reports
as NaN (flagged on the lane from the Hessian's finite value); cs_sr_acat
4's flat direction is its cs() coefficient (325.5), which the threshold
test does not read (filed). Both cs_sr_acat fits get the standard-error
check's "not available for 4 of 15 parameters" on the base already.

Candidates rejected on these numbers:

- the shift test: the flat thresholds give 7e-11 to 0.0156, but fits with
  a gap of about 1e-8 between two thresholds give 4.6e-4 to 1.1e-3 too
  (the log-gap parameter moves 10 either way at no cost), so it does not
  separate them (`ord-forty-base.txt`, column flat10);
- the per-threshold sharpening: sound fits reach -0.034 to -0.07
  (cum_sratio 4 and 10, cum_cum 4, standard errors 3.9 and below),
  because doubling the distances to a threshold the rows straddle moves
  little; one flat point gains (+0.751, cum_sratio 18)
  (`ord-forty-base.txt`, column tsharp).

### The rule built

`R/families.R`:

- new `mixture_ord_thres_se(fit, resp)`: per component, the largest
  standard error of a threshold the data place (`unident` left out), on
  the latent scale (times the component's median `disc`) and per unit of
  mean case weight, from the covariance of the fit-end standard-error
  check (`fit$cache$sdr$cov.fixed`, else `se_fit_analysis()$V`, which
  check_convergence() and se_check() share and cache) through the
  Jacobian of `ord_threshold_values()` by central differences. `NA` when
  there is no such covariance (a random-effects fit whose check waits,
  `check_se = "ignore"`) or the block is not held whole.
- `mixture_ord_se_cut <- 100`.
- `mixture_ord_degenerate_check()`: a component whose threshold standard
  error is above the cut, and not already a run-off threshold, is
  reported: "a threshold the data do not place: its standard error is X
  on the latent scale, so the likelihood is flat along it".
- `mixture_ord_degeneracy()`: "collapsed" is two ordered thresholds
  within 8 ulps (`diff <= 8 eps max(1, |tau|)`), not the same double
  (backlog 5167's first half). The message now says "the same number up
  to rounding, so the category between them has no probability in it".
- The start ranking (`mixture_ord_start_class()`) is unchanged in rule;
  it reads the new "collapsed" through `mixture_ord_degeneracy()`.

Why weights: the sharpening cut is per unit of mean weight already
(`test-ordinal-mixture.R`, "constant weights do not change the
degenerate verdict": weights 1/500 on a sound fit with largest SE 11.8
would read 264 unnormalized).

### Lane against base on the controls

`dev/silent2-ordcmp.R`, logs `dev/silent2-out/ordcmp-*.txt`:

```text
sets forty sweep steep wide ; fits 185 of base 185 lane 185
logLik moved by > 1e-6: 0
code changed: 0
degenerate warning: base 21 lane 42 ; gained 21 ; lost 0
gained with every component's SE finite and at most 100: 1   (cs_sr_acat_seed1, component 2's fixef() SE NaN)
passed by the lane (code 0, no warning): 142 ; their largest threshold SE on the base arm: 27.21
```

Per set (`ordcmp-forty.txt` and so on): forty gained 6 (cum_cum 16,
cum_sratio 5, 6, 14, 16, 18), sweep 14, steep 1 (steep_n300_s3_K6
seed 2), wide 0. The 8-ulp collapse changed no control fit
(`dev/silent2-ordcollapsed.R`, `ordcollapsed.txt`: "fits 185 ;
collapsed on base 13 ; on lane 13 ; newly collapsed 0"); its case is
the 1e-15 twin, which the test builds directly.

The lane arm of forty and sweep ran on the build before the weight
normalization was added (build 3); steep and wide on build 4. No control
fit has case weights, so the two builds give the same verdicts on them.

### Every rounding

`dev/silent2-ulps.R lane`, log `dev/silent2-out/ulps-lane.txt`: the test's
three fits on x scaled by 1 + 2e-15 k, k = 0..20:

```text
SUMMARY seed16 pass 21 of 21; seed2 sratio pass 21 of 21; seed2 cumulative pass 21 of 21
```

Every k ends at logLik -379.935017 with SE 1.753e+07 (seed 16) and 1.41,
1.06 (seed 2).

### Not done in this group

- 5214 (no jittered starts for ordinal mixtures with random effects):
  the default is the user's (ordstart Q1 and F6, the release record) and
  how often the answer moves was not measured; nothing here changes it.
  Left open.
- The threshold-SE test needs the fit-end covariance. A random-effects
  ordinal mixture whose standard-error check waits for the first SE use
  (`se_check_at_fit()` false) is not tested at fit end. Filed.
- A flat direction on a `cs()` coefficient (ordstart's cs_sr_acat seed 4,
  coefficient 325.5): the threshold-SE test does not read coefficients.
  On base that fit already gets "Standard errors are not available for 4
  of 15 parameters" (`ord-sweep-base.tsv`), so it is not silent. Filed
  as a note.

## Item 5524: floors outside a support for a reader

Done in core: the ordinal mixture's gap floor. `R/families.R`: new
`ord_gap_exact(...)` (`exact_support()` and no advector argument);
`ord_log_interior()` and `ord_gap_floor()` skip the floor for it, and
`ord_touch_nan()` returns 0, so a category between tied thresholds reads
`-Inf` and a crossed `cs()` pair reads NaN (the log of a negative
number). thres.R's grouped path calls the same two functions. Test "an
exact reader reads a closed category as -Inf" (test-ordinal-mixture.R):
with both components' second threshold tied, `row_lpdf(exact = TRUE)`
gives `-Inf` on every category-2 row, the default gives a finite value,
every other row is identical, and `fit$obj$fn()` at that point is finite
(the tape keeps the floor).

Not done, and refiled as a narrower item:

- `log1m_exp_ad()`: its floor binds where `l` rounds to 0, which is an
  empty window (exact `-Inf`) or a window narrower than the tails
  resolve (exact value finite); the doubles do not tell them apart
  (`dev/support-findings.md`, "Why the rows marked not changed stay").
  A switch would trade one wrong value for another.
- The U rows of the survey (lba and rdm densities floored at 1e-300,
  rdm's lcdf, gddm's grid floor, survival's caps and `sv_squeeze()`):
  the exact value is positive and below what the formula resolves, so
  `-Inf` would be as wrong as the floor; they need log-space forms.
- The S rows in the extensions (lba/rdm race at a decision time <= 0,
  wiener_gng go rows at y <= ndt, wiener `allow_unreachable`, the
  contaminant outside `crange`, survival's `sv_floor_pos`): each can read
  `-Inf` for an exact reader, but `exact_support()` is not exported, so
  each extension needs it added to the sampling API first, and each site
  needs its own reader test. That is four files in two extensions plus a
  core export, more than this item's size; refiled.

## Files touched

Package code:

- `R/predict.R`: `aterm_rowless()` (new), `aterms_for_newdata()` (item
  5235); `pred_design()`'s group-level `model.matrix()` (item 5353).
- `R/conditional-effects.R`: `ce_aterms()` (item 5235).
- `R/frame.R`: `loglik_rows_key()` (new, item 5246);
  `re_term_contrasts()`, `re_contrasts_for()` (new), the plain
  component's `contrasts`, `mms$contrasts`, the block builder's copy of
  `contrasts`, `mm_member_designs()` (item 5353).
- `R/confint.R`: `anova.frmtmb_fit()` (item 5246).
- `R/re-formula.R`: `re_term_cnms()` (item 5353).
- `R/families.R`: `ord_log_interior()`, `ord_touch_nan()`,
  `ord_gap_floor()`, `ord_gap_exact()` (new) (item 5524);
  `mixture_ord_degenerate_check()`, `mixture_ord_thres_se()` (new),
  `mixture_ord_se_cut` (new), `mixture_ord_degeneracy()`'s collapse test
  (items 5207, 5167).
- `extensions/frmtmb.sample/R/conditional-effects-draws.R` (item 5362).
- `extensions/frmtmb.sample/R/sample.R`: `print.frmtmb_draws()` (print
  order).

Tests: `tests/testthat/test-newdata-aterms.R`, `test-subset-rate.R`,
`test-predict-newdata.R`, `test-ordinal-mixture.R` (four new tests and
the `thres(4)` test's allowed warning); frmtmb.sample
`tests/testthat/test-laplace-draws.R`, `test-brms-draws-order.R`.

Records: `NEWS.md`, `extensions/frmtmb.sample/NEWS.md`,
`dev/test-backlog.md` (seven CLOSED lines, the section "Filed in the
round of 2026-11-13" with "### From lane silent2", four items),
`dev/silent2-*` (scripts and this file). No man page changed (every new
function is `@noRd`; roxygenise rewrote nothing).

Shared areas named by the environment file: `R/predict.R` (one hunk in
`pred_design()`, one in `aterms_for_newdata()`, one new function before
it), `R/confint.R` (one hunk in `anova.frmtmb_fit()`), frmtmb.sample
`R/` (CE draws, print). `R/fit.R`, `R/se-check.R`, `R/priors.R` and the
`theta` lines are not touched; `mixture_ord_thres_se()` reads
`se_fit_analysis()` without changing it.

## Tests

### Seen to fail on rellib-r21, then pass on the lane

`dev/silent2-runtest.R base|lane <package> <file>`, one process per file,
NOT_CRAN=true; logs `dev/silent2-out/t-base-*.txt`, `t-lane-*.txt`:

```text
test-newdata-aterms.R          base pass=59 fail=4 err=0          lane pass=63 fail=0 err=0 warn=0
test-subset-rate.R             base pass=97 fail=2 err=0 warn=1   lane pass=99 fail=0 err=0 warn=0
test-predict-newdata.R         base pass=15 fail=0 err=1          lane pass=20 fail=0 err=0 warn=0
test-ordinal-mixture.R         base pass=110 fail=3 err=0         lane pass=114 fail=0 err=0 warn=0
frmtmb.sample test-laplace-draws.R     base pass=124 fail=6 warn=24   lane pass=130 fail=0 err=0 warn=0
frmtmb.sample test-brms-draws-order.R  base pass=37 fail=2            lane pass=39 fail=0 err=0 warn=0
```

The base failures are behavioral: the 120-row newdata and grid are read
(4 failed refusals); `anova()` runs on both pairs (2); the prediction
after the option reset errors "Random-effect design ... does not match"
(the 1 error); the flat-threshold fit is silent while its SE is 1.75e7,
the 1e-15 twin is not collapsed, and the exact reader reads -690.8 at a
tie (3); the band is all NA under "NAs produced" (6 failures, 24
warnings); `print()` lists the gp columns pairwise (2).

The lane's first full run of test-ordinal-mixture.R let one warning
escape: the `thres(4)` test's `mixture(cumulative(), sratio())` fit now
gets the degenerate warning for component 2's third threshold, at -16.7
with a standard error of 4846 (`dev/silent2-probe-thres4.R`; fixef():
`mu2_Intercept[3] -16.676 4.846e+03`, near 12.6, cover 1.3e-12). That is
a correct report, so the test now allows it by its words.

### Every suite, one file per process

`bash dev/silent2-full.sh silent2-full dev/silent2-out/jobs-all.txt`
(run-par.sh, 3 at a time, NOT_CRAN and FRMTMB_BRMS_FIT_TESTS true,
FRMTMB_RELLIB the lane library); summary `dev/silent2-fullsum.R`, log
`dev/silent2-out/fullsum.txt`:

```text
frmtmb files=330 pass=27475 fail=0 err=0 skip=15 warn=0 lanelib=330
frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0 lanelib=11
frmtmb.eam files=33 pass=1857 fail=0 err=0 skip=3 warn=0 lanelib=33
frmtmb.latent files=13 pass=697 fail=0 err=0 skip=2 warn=0 lanelib=13
frmtmb.learn files=16 pass=545 fail=0 err=0 skip=2 warn=0 lanelib=16
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0 lanelib=11
frmtmb.sample files=91 pass=4618 fail=0 err=0 skip=1 warn=0 lanelib=91
frmtmb.spatial files=8 pass=612 fail=0 err=0 skip=0 warn=0 lanelib=8
frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0 lanelib=16
no RESULT line: 0
```

`SILENT2-FULL ran 529 of 529`. Every log's `lib:` line names
`wt-silent2-lib` for frmtmb (and for frmtmb.sample in its own files);
the other extensions load from rellib-r21, unchanged. The 29 skips are
the known ones: test-drmtmb-agreement.R 13, six extensions'
test-scale.R files 14, test-fuzz.R 1, test-curve-scale.R 1. The installed core then took a
comment-only change (the roxygen of `mixture_ord_thres_se()` and
`mixture_ord_se_cut`) and was reinstalled before the sessions and the
check.

### OpenBLAS 0.3.26 with LAPACK

`bash dev/ciharden-openblas.sh 0.3.26 lapack` (log
`dev/silent2-out/openblas-build.txt`), OPENBLAS_NUM_THREADS=4, logs
`dev/silent2-out/ob326-*.txt`:

```text
RESULT test-ordinal-mixture.R pass=114 fail=0 err=0 skip=0 warn=0
RESULT test-mixture-starts.R pass=95 fail=0 err=0 skip=0 warn=0
RESULT test-newdata-aterms.R pass=63 fail=0 err=0 skip=0 warn=0
RESULT test-subset-rate.R pass=99 fail=0 err=0 skip=0 warn=0
RESULT test-predict-newdata.R pass=20 fail=0 err=0 skip=0 warn=0
RESULT test-laplace-draws.R pass=130 fail=0 err=0 skip=0 warn=0 (frmtmb.sample)
RESULT test-brms-draws-order.R pass=39 fail=0 err=0 skip=0 warn=0 (frmtmb.sample)
```

`dev/ciharden-out/` deleted after. 0.3.32 was not run: the SE cut sits
between 27 and 640 on the controls, not near a rounding-decided branch,
and the 8-ulp collapse test involves no BLAS.

### Every suite in one session, three orders

`run-session.R` through `dev/silent2-sessions.sh` (run-sessions.sh, two
at a time), then, after the core forward and reverse sessions were under
way, the driver loop was stopped and the other 25 sessions were run by
`dev/silent2-session1.sh` / `dev/silent2-sessions-list.sh` (the same
runner, environment and output files), three at a time. NOT_CRAN and
FRMTMB_BRMS_FIT_TESTS true, the lane library first. Logs
`dev/release/silent2-sess-sessions/<package>--<order>.txt`:

```text
SESSION frmtmb forward files=330/330 pass=27475 fail=0 err=0 skip=15 warn=0
SESSION frmtmb preload files=330/330 pass=27743 fail=0 err=0 skip=17 warn=0
SESSION frmtmb reverse files=330/330 pass=27640 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=33/33 pass=1857 fail=0 err=0 skip=3 warn=0
SESSION frmtmb.eam preload files=33/33 pass=1857 fail=0 err=0 skip=3 warn=0
SESSION frmtmb.eam reverse files=33/33 pass=1857 fail=0 err=0 skip=3 warn=0
SESSION frmtmb.latent forward files=13/13 pass=697 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=13/13 pass=697 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=13/13 pass=697 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=16/16 pass=545 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn preload files=16/16 pass=545 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn reverse files=16/16 pass=545 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward files=91/91 pass=4618 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample preload files=91/91 pass=4618 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample reverse files=91/91 pass=4618 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=612 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial preload files=8/8 pass=612 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial reverse files=8/8 pass=612 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival forward files=16/16 pass=575 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=16/16 pass=575 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=16/16 pass=575 fail=0 err=0 skip=0 warn=0
```

27 of 27 sessions; 194 LEAK lines, every one tagged `[with load of
...]` (an extension's or a Stan model's load), 0 untagged; 0 `DETAIL
FAIL|ERROR|WARN` lines. Each session's `lib:` lines name
`wt-silent2-lib` for frmtmb and frmtmb.sample, rellib-r21 for the rest.

### R CMD check --as-cran

`bash dev/silent2-check.sh <package>` (built and checked in
`dev/silent2-check/<package>/`, R_LIBS the lane library, rellib-r21,
the user library; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; pandoc and
TinyTeX on PATH). Logs kept as `dev/silent2-out/check-frmtmb-00check.txt`
and `check-sample-00check.txt`; the check trees are deleted.

```text
frmtmb         * checking tests ... [68m] OK     Status: 2 NOTEs
frmtmb.sample  * checking tests ... [292s] OK    Status: 1 NOTE
```

The NOTEs: "Skipping checking math rendering: package 'V8' unavailable"
(both, expected), and for core "Examples with CPU (user + system) or
elapsed time > 5s": `frm_hazard_reads` 7.68 s and `dharma_residuals`
6.89 s elapsed (3.25 and 1.07 s user). This lane changed neither
example, and the machine ran fourteen lanes' suites at the time
(lane-rules: "An examples-timing NOTE on this box measures LOAD").
Examples with `--run-donttest` OK (144 s); vignettes rebuilt OK (18 min).

## NEWS and version

NEWS bullets: `NEWS.md` under "# frmtmb (development version)" (six
bullets, four marked BREAKING), `extensions/frmtmb.sample/NEWS.md` under
"# frmtmb.sample (development version)" (two bullets).

- frmtmb: minor bump. `anova()` and new data now refuse what they read
  wrongly, and the ordinal-mixture check warns on 21 of the 185 control
  fits that passed before. No export changes.
- frmtmb.sample: patch bump. It calls only what frmtmb 0.82.0 exports
  (`frm_linpred()`, `ce_plan_eval()`), so its floor stays frmtmb 0.82.0.
- Other extensions: no change.

## For the user

1. The threshold standard-error cut (100 on the latent scale) turns 21
   of the 185 control fits from silent to the degenerate warning; none
   of them moves. My recommendation is to keep it as a warning, as the
   other degenerate cases are: a threshold with a standard error of 640
   or more is not an estimate. The alternative is a message.
2. 5214 (no jittered starts for ordinal mixtures with random effects)
   stays as the earlier rounds left it, the default being yours.
3. The rest of 5524 needs `exact_support()` exported through the
   sampling API before frmtmb.eam and frmtmb.survival can read it; a
   follow-up lane can take the five S sites (filed).

## Fix round (review dev/reviews/2026-11-13-silent2.md, FIX)

Where this section and the ones above
differ, this section is the record. Logs `dev/silent2-out/fx-*`,
`rev-*-after.txt`, `probe-fix-*.txt`, `recost.txt`, `ulps2-lane.txt`.

### F1 (5353, a transformed factor in a bar)

`re_term_contrasts()` now reads each variable of the left side by its
deparsed name on the frame (`factor(cyl)`), and evaluates it on the frame
only when the frame lacks that column; the names are the ones the
prediction's model frame gives, which `re_contrasts_for()` matches.

Reviewer's scripts on the fixed build (`rev-probe-after.txt`,
`rev-contr-after.txt`, `rev-core-after.txt`) against the reviewer's lane
logs:

```text
RE-only (1 + factor(cyl) | g) sum -> helmert, max |diff|   before 0.4826564   after 0
(1+factor(cyl)|g) sum -> helmert (same names)               before 0.538       after 0
(1+factor(cyl)|g) sum -> treatment                          before ERROR       after 0
```

Every other contrasts case of `silent2-review-contr.R` and
`silent2-review-core.R` stays 0: `(1+f|g)` to helmert and treatment,
ordered `(1+o|g)`, character column, `(1+f:x|g)`, `(1+f||g)`,
`gr(g, by = )`, `sigma ~ (1+f|g)`, a data attribute, `re_formula` all
and part; `mm(g, g2)` 2.2e-16 and 4.4e-16 (in-sample against newdata,
rounding), reset 0. (`silent2-review-core.R` stops at its own
`gr(g, by = ch)` case, whose `ch` varies within `g`, as in the
reviewer's log.) New test "a transformed factor in a group-level term
keeps the fit's contrasts" (test-predict-newdata.R): both formulas, fit
under contr.sum, predicted under contr.helmert and contr.treatment.

### F2 (5246, the row key)

`loglik_rows_key()` stores a factor or ordered response by its codes and
a logical one as 0/1, every value as a double (character as character),
and drops the response names when each fit has one response.
`silent2-review-core.R`, lane, before (reviewer's log) and after:

```text
anova y ~ 1 vs I(y) ~ x                  before ERROR   after ran
anova cumulative int vs ordered factor y before ERROR   after ran
anova bernoulli int vs logical y         before ERROR   after ran
```

and `y` against a copy `y2` (new test). Every other row of that script
still "ran". The lane's refusals still refuse (`repro-lane-fix.txt`):
rows 1..72 against 9..80, NA rows of different predictors, `rev(y)`. The
refusal text is unchanged: it is now true wherever it fires. The
reviewer's `cens(c1)` against `cens(c2)` still compares; filed LOW.

### F3 (5207, a lost threshold's message)

`mixture_ord_thres_se()` returns attribute `lost`: whether the SE check
lost (by `outer_par_map()` name) one of the component's raw threshold
parameters. The flag rule is unchanged (SE above 100); the message then
reads "a threshold the data do not place: it has no standard error",
with no number and no "flat". A warning that names an unplaced
threshold carries class `frmtmb_ord_unplaced`, and `fit_end_checks()`
(R/fit-end.R) does not mark the SE check explained for it, so the SE
check's list stands beside it. cs_sr_acat seed 1 (`probe-fix-before.txt`,
`probe-fix-after.txt`):

```text
before: component 2 (a threshold the data do not place: its standard error is 334 on the latent scale, so the likelihood is flat along it) ...   [no SE warning]
after:  component 2 (a threshold the data do not place: it has no standard error) ...
        Standard errors are not available for 4 of 15 parameters. tau_raw2_2, tau_raw2_3, bcs4_2, bcs4_3: ...
```

The same un-suppression made `test-ordinal-mixture.R`'s "the degenerate
check says only what it found" (three components) raise the SE list
(6 of 14 parameters); that test now allows it, by its words.

### F4 (the 5167 CLOSED line)

Reworded: closed in part; the 8-ulp band closes the 1e-15 twin only, the
SE test flags 5 of the six flat moves (not cs_sr_acat seed 4). The
remainder (twins 1e-12 and 1e-8 apart) is filed LOW under "From lane
silent2".

### F5 (5235, a constant value)

New `aterm_scalar()`: a value that reads no column of newdata (or the
grid) and has one distinct element is that element, before the length
check, in `aterms_for_newdata()` and `ce_aterms()`. The refusal of a
varying value now says "its values are the fitted rows' and not
newdata's, even when newdata is the fitted data. Leave newdata out to
predict the fitted rows, or write the term as a column of the data".
`silent2-review-core.R`:

```text
5235 constant rep(2,120) newdata 120   before ERROR   after ran
5235 constant rep(2,120) CE res 120    before ERROR   after ran
5235 constant rep(2,120) CE default    before ran     after ran
5235 fitted(newdata = data)            ERROR both, now with the message above
```

The CE `method = "predict"` route still reaches `aterms_for_newdata()`
and so says "newdata" for a grid (noted by the reviewer; not changed).

### 5214

1. The degenerate warning of a fit with random effects, where no starts
   ran and `starts` was not set, now says "A fit with random effects
   runs no jittered starts by default: run them with
   frmtmb_control(starts = 10)". Test "a degenerate fit with random
   effects names the starts to run" (cum_sratio seed 9 with `(1 | g)`,
   20 groups; and, absent, the no-RE seed 16 fit does not say it).
2. Not shipped. `dev/silent2-recost.R` (reviewer's design: n = 1000, 50
   groups, `y ~ x + z + (1 | g)`, `mixture(cumulative(), cumulative())`,
   minimum of 3 rounds), `recost.txt`:

```text
RECOST seed 1 np 13 evals 248 deferred FALSE fit 4.80s covariance 0.00s ratio (fit + cov) / fit 1.000
RECOST seed 2 np 13 evals 121 deferred TRUE fit 2.30s covariance 0.72s ratio (fit + cov) / fit 1.313
RECOST seed 3 np 13 evals 142 deferred TRUE fit 2.40s covariance 0.86s ratio (fit + cov) / fit 1.358
```

   Above the 1.2 bound on both deferred fits, so the test still reads
   only a covariance the fit built anyway (seed 1 builds it: 4 x 13 = 52
   gradients against a quarter of 248 evaluations). The filed item
   carries the numbers.

### Nits

NEWS: the `log_lik()` bullet is BREAKING; the control-fit counts left
NEWS (they are above, in this file). NEWS also gained the F1, F2, F3, F5
and 5214 facts.

### Merge note for consolidation (famend, R/families.R)

famend's tree changes `mixture_ord_degenerate_check()` at the same lines
(base numbering). Resolve each conflict to keep both:

1. Doc list (base 7049, 7064-7066): "flagged in five cases"; the bullets
   in the order step, run-off, silent2's rewritten collapsed bullet
   ("the same double up to rounding (8 ulps) ..."), silent2's standard
   error bullet, then famend's "A `cs()` effect has closed a row's own
   category" bullet.
2. The flags (base 7088): keep silent2's `se`, `unplaced` and `nose`
   lines and famend's `shut <- (dg$closed < 1e-8) %in% TRUE`, then
   `bad <- step | runoff | unplaced | dg$collapsed %in% TRUE | shut`.
3. The messages (base 7105-7107): silent2's `unplaced` branch (with the
   `nose` alternative) after `runoff`, silent2's collapsed text, then
   famend's `if (shut[i])` branch.

The `frm_warning()` call that follows is changed by silent2 alone (the
random-effects starts text and `class = "frmtmb_ord_unplaced"`) and
merges. After the merge run famend's csmix test and test-ordinal-mixture.R:
a fit that is both `shut` and `unplaced` now also keeps the SE check's
list (F3), which a test allowing only the degenerate warning would see.
### Fix-round tests

New tests seen failing on the build before the fix (the installed lane
build of the first round), then passing (`fx-pre-*.txt`, `fx-post-*.txt`):

```text
test-newdata-aterms.R    before pass=63 fail=0 err=1    after pass=67 fail=0 err=0 warn=0
test-predict-newdata.R   before pass=21 fail=1 err=1    after pass=26 fail=0 err=0 warn=0
test-subset-rate.R       before pass=100 fail=4 err=0   after pass=104 fail=0 err=0 warn=0
test-ordinal-mixture.R   before pass=118 fail=6 err=0   after pass=124 fail=0 err=0 warn=0
frmtmb.sample test-laplace-draws.R      after pass=130 fail=0 err=0 warn=0
frmtmb.sample test-brms-draws-order.R   after pass=39 fail=0 err=0 warn=0
```

Before: the constant payload errors (F5); helmert silently off and
treatment refused (F1); four anova pairs refused (F2); the F3 test's
five assertions and the 5214 test's message fail. The first post-fix
run of test-ordinal-mixture.R had warn=1: the three-component test met
the SE list F3 now keeps; that test allows it (above).

Every rounding (`dev/silent2-ulps2.R`, `ulps2-lane.txt`), x scaled by
1 + 2e-15 k, k = 0..20:

```text
SUMMARY cs1 pass 21 of 21; re9 pass 21 of 21
```

(cs_sr_acat seed 1 ends at -404.991560 with its threshold SE lost on
every k.)

OpenBLAS 0.3.32 with LAPACK, `OPENBLAS_CORETYPE=Haswell`,
`OPENBLAS_NUM_THREADS=4` (`ob332-*.txt`; `dev/ciharden-out/` deleted):

```text
RESULT test-ordinal-mixture.R pass=124 fail=0 err=0 skip=0 warn=0
RESULT test-predict-newdata.R pass=26 fail=0 err=0 skip=0 warn=0
```

Full core suite, one file per process (`dev/release/silent2-fix-core.log`,
`fix-core-sum.txt`):

```text
frmtmb files=330 pass=27500 fail=0 err=0 skip=15 warn=0 lanelib=330
SILENT2-FIX-CORE ran 330 of 330
```
Sessions, forward order (`dev/silent2-session1.sh`, logs
`dev/release/silent2-sess-sessions/<package>--forward.txt`):

```text
SESSION frmtmb forward files=330/330 pass=27500 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.sample forward files=91/91 pass=4618 fail=0 err=0 skip=1 warn=0
```

0 untagged LEAK lines in either (48 and 15 tagged). The first
frmtmb.sample run read pass=4615 warn=2 with an ERROR in
test-brms-draws-order.R, "cannot create dir ... No space left on
device": drive C ran out while fourteen lanes ran (4.25 GB free after,
144 GB later). It was rerun once the disk had room; the line above is
the rerun. The first core check died the same way in `R CMD build`
("No space left on device") and was rerun.

`R CMD check --as-cran` of core on the fix build
(`dev/silent2-out/check-frmtmb-fix-00check.txt`; check tree deleted):

```text
* checking examples ... [88s] OK
* checking examples with --run-donttest ... [100s] OK
* checking tests ... [36m] OK
* checking HTML version of manual ... [34s] NOTE   (V8 unavailable)
Status: 1 NOTE
```

frmtmb.sample's code did not change in the fix round; its check of the
first round stands (Status: 1 NOTE).

Status: FIX ROUND DONE (2026-11-13).
## Fix round 2 (short check dev/reviews/2026-11-13-silent2-short.md, FIX on F2)

The frame stores a factor response as its codes (`frame$y` numeric;
`y_levels` only for ordinal and categorical; `dev/silent2-probe-ylev.R`),
so the key of the first fix round compared codes. Now
`loglik_rows_key()` reads a factor, character or logical response's raw
values off the fit's data (new `loglik_rows_raw()`, by the response
expression) and keys them by their labels, a logical as "0"/"1"; a
numeric response keeps its double key. New `loglik_rows_same()`
compares two keys element by element, a numeric key against a label key
through `as.character()` of the number; `anova()` uses it.

`dev/reviews/silent2-short-f2.R lane`, before (`f2-lane.txt` of the
short check) and after (`dev/silent2-out/short-f2-after.txt`):

```text
                                                      before    after
y ~ 1 vs I(y) ~ x                                     ran       ran
y ~ 1 vs y2 ~ x (copy)                                ran       ran
cumulative int vs ordered factor                      ran       ran
bernoulli int vs logical                              ran       ran
rows 1..72 vs 9..80                                   REFUSED   REFUSED
NA rows of different predictors                       REFUSED   REFUSED
y vs rev(y)                                           REFUSED   REFUSED
y vs y + 1e-12                                        REFUSED   REFUSED
categorical fa (a..d) vs fb (w..z), equal codes       ran       REFUSED
cumulative oa (1,2,3,4) vs ob (1,2,3,5), equal codes  ran       REFUSED
cumulative oa (1..4) vs oc (2..5), equal codes        ran       REFUSED
bernoulli factor no/yes vs lo/hi, equal codes         ran       REFUSED
```

The short check's two side notes (false refusals) now compare: the same
labels with levels reversed, and integers 2..5 against their ordered
factor ("ran" in the after log; the first warns only that the models
are not nested).

Test (test-subset-rate.R, "anova() compares one response by its values,
whatever its form"): new refusals `fa` against `fb` and `oa` (1..4)
against `oc` (2..5), and the integers 2..5 against `oc` compare. On the
build of fix round 1 (`fx2-pre-test-subset-rate.R.txt`): pass=104
fail=3 err=0 warn=1 (the two refusals did not happen, the 2..5 pair was
refused, and the `fa`/`fb` anova ran with its not-nested warning). On
the fix build (`fx2-post-*.txt`):

```text
RESULT test-subset-rate.R pass=107 fail=0 err=0 skip=0 warn=0
RESULT test-predict-newdata.R pass=26 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture.R pass=124 fail=0 err=0 skip=0 warn=0
RESULT test-newdata-aterms.R pass=67 fail=0 err=0 skip=0 warn=0
```

Nits: the contrasts NEWS bullet is BREAKING; `R/families.R` says "four
cases" (the famend merge makes it five, as the merge note says). The
anova NEWS bullet says a factor response compares by its labels.

Files: `R/frame.R` (`loglik_rows_key()`, `loglik_rows_raw()` and
`loglik_rows_same()`), `R/confint.R` (one line in `anova.frmtmb_fit()`),
`R/families.R` (one doc word), `tests/testthat/test-subset-rate.R`,
`NEWS.md`. Only the key, its comparison, a comment and NEWS changed, so
`R CMD check` was not rerun (the brief's condition); the fix-round-1
core check stands (Status: 1 NOTE).

Status: FIX ROUND 2 DONE (2026-11-13).