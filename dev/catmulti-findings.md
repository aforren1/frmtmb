# Lane catmulti: refcat, the bare multinomial(), and the simplex families in mvbf()

Lane catmulti, round opened 2026-10-14, base 47cf659d (frmtmb 0.71.0,
frmtmb.sample 0.19.0). Worktree `frmtmb-wt-catmulti`, branch
`wt-catmulti`, library `C:/Users/adf44/source/r/wt-catmulti-lib` ahead
of `rellib-r9`. Nothing is committed. Every script named is in `dev/`;
its log is in `dev/catmulti-out/` (local, gitignored). Numbers below
are pasted from those logs.

Scope: the categorical, multinomial and simplex parity gaps filed at
0.71.0 (`dev/test-backlog.md`, lane famsimplex's section;
`dev/famsimplex-findings.md`, sections 5 and 8), items 1 to 5 of the
brief, the stretch items of 6, and one item the coordinator added: the
cumulative tied-threshold NaN that lane bin2port found.

## 1. What brms 2.23.0 does, as run

`dev/catmulti-brms-probe.R` (stancode, standata, default_prior, family
objects; log `brms-probe.txt`), `dev/catmulti-brms-fits.R` (six fits,
2 chains x 2000; log `brms-fits.txt`), `dev/catmulti-brms-fits2.R` (a
mixture, a two-category categorical, multinomial newdata; log
`brms-fits2.txt`), `dev/catmulti-brms-nobs.R` (log `brms-nobs.txt`).
Data generators: `dev/catmulti-data.R`.

- `categorical(link = "logit", refcat = NULL)` and
  `multinomial(link = "logit", refcat = NULL)`. `refcat` is
  `as_one_character(allow_na = TRUE)`: `refcat = 2` is the category
  `"2"`, `c("a", "b")` is refused ("Cannot coerce 'refcat' to a single
  character value."), an unknown name is refused ("The reference
  response category must be one of 'a', 'b', 'c'.").
- dpar names: `make_stan_names(paste0("mu", predcats))`. Factor levels
  a, b, c: `mub muc` (default), `mua muc` (`refcat = "b"`),
  `mua mub muc` (`refcat = NA`). Levels `p.q, r_s, t` with
  `refcat = "r_s"`: `mupq mut`.
- multinomial categories are the response columns: `cbind(y1, y2, y3)`
  and a matrix column `Yc` give `muy2 muy3`; an unnamed matrix gives
  `mu2 mu3`; `refcat = "y2"` gives `muy1 muy3`; `NA` gives all three.
  Stan: the reference's predictor is the literal 0 in `mu[n]`.
  `standata()` without `trials()`: "Specifying 'trials' is required
  for this model." One column: "At least 2 response categories are
  required."
- `default_prior()`: `student_t(3, 0, 2.5)` Intercept and flat `b` per
  predicted category, named by the reference (rows `mua`, `muc` under
  `refcat = "b"`; three dpars under `NA`).
- `refcat = NA` fits in brms only through its priors, and with the
  default priors (flat `b`) it does not converge: on 300 rows the
  slopes came out near 519 with sd 696, max R-hat 1.709, 1422
  transitions at the maximum treedepth (`brms-fits.txt`, "cat-refna").
- Post-fit shapes (`brms-fits.txt`):
  - categorical `refcat = "b"`: `fitted()` 300 x 4 x 3, layers
    `P(Y = a), P(Y = b), P(Y = c)`; `posterior_linpred()` 2000 x 300 x 2,
    every dimname NULL; `fitted(scale = "linear")` 300 x 4 x 2, layers
    `eta1, eta2`; under `refcat = NA` both have 3 layers.
  - multinomial: `fitted()` 200 x 4 x 3, layers `P(Y = y1)` ...,
    expected COUNTS (probability times trials); `predict()` 200 x 4 x 3
    layers `y1, y2, y3`; `posterior_predict()` S x N x 3 named
    `y1, y2, y3`; `fitted(newdata)` without the trials column is
    refused ("The following variables can neither be found ... 'N'").
  - `posterior_linpred(transform = TRUE)` with no dpar on all of these:
    "Invalid argument 'dpar'. Valid distributional parameters are:
    'mua', 'muc'".
  - mixture of two gaussians: `posterior_linpred()` 1000 x 300 x 2 and
    `fitted(scale = "linear")` 300 x 4 x 2, both with NULL layer names
    (`brms-fits2.txt`). A two-category categorical stays S x N and
    n x 4.
- mvbf: `dirichlet()` + gaussian fits (fixef below). `logistic_normal()`
  + gaussian DOES NOT COMPILE in brms 2.23.0 (brms-C, section 8).
- `nobs()` of a univariate `subset()` model is 300 of 300 rows in
  brms, and its summary says "Number of observations: 300";
  `log_lik()` has 299 columns (`brms-nobs.txt`).

## 2. What changed

### Core

- `categorical(link, refcat, levels, K)` (brms's order: `refcat`
  second). `fam_categorical_impl(cats, refcat, link, levels)` keeps
  the 1..K response codes and stores the predicted positions
  (`cat_pred`) and `refcat`; `cat_logit_lpdf(pred = )` places each
  predictor at its category, and under `NA` the softmax denominator has
  no implicit zero. The reference-1 arithmetic is unchanged term for
  term (the default `pred` is `2..K`). `cat_probs()` and the simulator
  read `cat_pred`. dpar names are brms's `make_stan_names()`
  (`simplex_dpar_names()`), a BREAKING change for levels with dots or
  underscores.
- `multinomial(link = "logit", refcat = NULL, K = NULL)`. With no `K`
  a deferred family reads the categories off the response columns
  (`multinomial_categories()`), as `categorical()` reads levels;
  `K =` still builds a positional family away from the data (the
  frmtmb.latent `hmm()` component path uses it). `multinomial(3)` is
  refused with a pointer to `K = 3`. `fam_multinomial_impl()` has the
  same `refcat` handling, a `post$mean_fn` of expected counts, and
  simulated counts carry the column names.
- `fitted()` on multinomial: brms's n x 4 x K through the existing
  matrix-mean route (finite-difference delta method, `fit_fd_se()`);
  the Wald bounds are held in [0, trials] (`prob_clamp_quantiles(upper
  = )`).
- `fitted(scale = "linear")` with no dpar stacks every location
  predictor (`mu_stack_dpars()`: categorical, multinomial, simplex,
  mixture) into n x 4 x P, layers `eta1 ...` on a polytomous family
  and unnamed on a mixture.
- `logistic_normal()` in `mvbf()`: `check_mv_extra_family()` admits it.
  `lncor_blocks()` gives each block a reporting prefix, `lncor` in a
  univariate model and `lncor_<resp>` in a multivariate one; summary
  rows `lncor_Y(b,c)`, `variables()`/`hypothesis()` names
  `lncor_Y__b__c`, a `default_prior()` row per response with `resp`,
  and `set_prior(class = "lncor", resp = )` resolves to that response
  (`R/priors.R`, two hunks).
- `brmsfamily()` doc and `?frmtmb-families` "Categorical (nominal)
  responses" rewritten; compat rows `fitted x multinomial` and
  `mvbf x group:simplex` added (`R/compat.R`).
- Stretch: `confint(parm = "lncor__b__c")` and
  `confint(parm = "rescor__y1__y2")` on a 2 x 2 matrix resolve to the
  one chart coordinate (`par_alias_index()`), with the existing
  "internal scale" message; the refusal for a wider block now advises
  `hypothesis(..., class = NULL)` (delta method) or `method = "boot"`
  instead of `method = 'profile'`, which fails there. `drop1()` on a
  family without `mu` reads its scope from the first location
  predictor that takes the main formula.
- Coordinator item: `ord_cumulative_logpmf()` adds `(i1 + iK)` to `a`
  so an end row's masked interior pair is strictly ordered (section 6).

### frmtmb.sample

- `posterior_linpred()` with no dpar stacks the same set
  (`draws_mu_stack()`) into a draws x obs x P array, every dimname
  NULL; `transform = TRUE` there is refused in brms's words plus a
  frmtmb sentence. `?posterior_epred` "Categorical outcomes" updated.
- `draws_lncor_cols()` and the default `lkj(1)` read each response's
  own chart name from `frame$extra_map` and report
  `lncor_<resp>__b__c`.

### Added after the first suite run

- `conditional_effects()` on `multinomial()` takes the per-category
  display, as the simplex families do (`ce_matrix_mean()` in
  `R/conditional-effects.R`); `method = "predict"` there is refused as
  for a simplex family.
- `dharma_residuals()` on `multinomial()` is refused by name
  (`R/diagnostics.R`): with a mean, it had reached the rank code and
  failed with "(subscript) logical subscript too long".

## 3. Evidence

### Fit agreement with brms fits on the same data

frmtmb's ML estimate against the brms posterior (2 chains x 1000
post-warmup draws, brms default priors), `z = (ML - mean) / sd` and the
ratio of standard errors; `dev/catmulti-shapes.R lane`, log
`shapes-lane.txt`, data `dev/catmulti-data.R` seeds 11 and 12:

```
categorical(refcat = "b"), n = 300 (seed 11)
               frmtmb    brms       z se_ratio
mua_Intercept -0.3429 -0.3443  0.0096   0.9919
muc_Intercept -0.7649 -0.7607 -0.0230   0.9578
mua_x         -0.9922 -1.0106  0.1017   0.9819
muc_x         -1.2884 -1.3018  0.0666   1.0139
multinomial(), n = 200, N = 10 (seed 12)
muy2_Intercept  0.2587  0.2585  0.0043   0.9998
muy3_Intercept -0.3852 -0.3841 -0.0163   1.0148
muy2_x          0.7638  0.7683 -0.0682   1.0177
muy3_x         -0.5815 -0.5796 -0.0249   0.9939
multinomial(refcat = "y2")
muy1_Intercept -0.2587 -0.2603  0.0282   1.0044
muy3_Intercept -0.6440 -0.6480  0.0622   1.0066
muy1_x         -0.7638 -0.7636 -0.0041   0.9908
muy3_x         -1.3454 -1.3488  0.0446   1.0086
mvbf(dirichlet Y, gaussian w), n = 300 (seed 13; dev/catmulti-mvln.R)
mub_Y_Intercept  0.3490  0.3501 -0.0209   0.9913
muc_Y_Intercept -0.3812 -0.3804 -0.0151   1.0224
w_Intercept      1.0389  1.0386  0.0067   0.9961
mub_Y_x          0.5876  0.5884 -0.0157   1.0100
muc_Y_x         -0.4469 -0.4481  0.0236   1.0298
w_x              0.5536  0.5541 -0.0097   0.9808
```

Names are identical in every case. brms max R-hat: 1.0048, 1.0009,
1.0044, 1.0039 for those four fits. The first rows of `fitted()` agree
to the third decimal (categorical `P(Y = a)` row 1: 0.38985, brms
0.39007; multinomial `P(Y = y1)` row 1: 3.3033, brms 3.3100).

### brms's Stan program at frmtmb's estimates

`tests/testthat/test-refcat-brms.R` (gated), reported by
`dev/catmulti-lpreport.R`, log `lpreport.txt`: brms's log density with
flat priors minus frmtmb's log-likelihood (`const`), and the largest
gradient component of brms's log density at frmtmb's maximum:

```
LPCHECK const 2.273736754e-13 grad 0.000106 ours -272.8212648   categorical(refcat = "b")
LPCHECK const -1.136868377e-12 grad 0.00016 ours -620.4326575   multinomial()
LPCHECK const -7.958078641e-13 grad 0.00017 ours -620.4326575   multinomial(refcat = "y2")
LPCHECK const -1.47792889e-12 grad 0.000171 ours -610.7012958   multinomial(refcat = NA)
LPCHECK const -1.136868377e-13 grad 1.17e-05 ours -19.23747205  mvbf(dirichlet, gaussian)
LPCHECK const 0 grad 8.25e-05 ours 108.4897988                  mvbf(logistic_normal, gaussian), repaired program
```

The logistic-normal multivariate program is brms's with two
identifiers repaired in its generated quantities (brms-C, section 8);
the model block, which `log_prob()` evaluates, is brms's own. The
helper gained `code_fix =` for that and reads `Llncor_<resp>`
(`tests/testthat/helper-brms.R`).

### refcat = NA

Under ML a common shift of every predictor leaves the likelihood
unchanged. Decided: fit it, do not refuse. Each reason measured: brms
allows it, and its own default-prior fit does not converge (R-hat
1.709, slopes 519 with sd 696, `brms-fits.txt`); frmtmb's SE check
flags every coefficient ("Standard errors are not available for 6 of
6 parameters"). CORRECTED in fix round 1: an earlier version of this
sentence said the check names the flat direction. It did not: it
listed the parameters with the generic reasons (review finding 4).
Since fix round 1 the warning names it (section 12). The category
probabilities are determined (max
|fitted(NA) - fitted(ref a)| 2.1e-06, logLik difference -3.0e-09,
`shapes-lane.txt`); a prior on each predictor gives a MAP fit with
finite SEs and a centered shift (`test-refcat.R`); and
`dirichlet(refcat = NA)` already behaves this way. A refusal would also
fire on a model whose category formulas do identify the shift (a
category without an intercept), which the SE check reads correctly.

### Shapes against brms

`shapes-lane.txt` against `brms-fits.txt` and `brms-fits2.txt`:
`fitted()` 300 x 4 x 3 `P(Y = a..c)` (categorical) and 200 x 4 x 3
`P(Y = y1..y3)` (multinomial); `fitted(scale = "linear")` 300 x 4 x 2
`eta1, eta2`, 300 x 4 x 3 under `NA`, 300 x 4 x 2 unnamed for a mixture,
300 x 4 for a two-category categorical; `predict()` multinomial
200 x 4 x 3 `y1..y3`; draws (`frm_sample()`, 150 draws):
`posterior_linpred()` 150 x 300 x 2 and 150 x 200 x 2 with no dimnames;
`posterior_epred()` multinomial 150 x 200 x 3 `y1..y3`;
`posterior_predict()` 150 x 200 x 3 `y1..y3`; `transform = TRUE`
refused with brms's sentence. All as brms. `fitted(newdata)` without
the trials column is refused by both. One pre-existing difference
stays: `posterior_epred()` of categorical draws names its observation
margin (`1, 2, ...`) where brms leaves it NULL.

`conditional_effects()` on multinomial against brms on seed 12
(`dev/catmulti-ce-compare.R`, log `ce-compare.txt`): 300 of 300 rows
matched on `x` and `cats__`; |frmtmb - brms| / brms se median 0.022,
max 0.085; se ratio median 0.990, range 0.947 to 1.039.

### logistic_normal() in mvbf()

`dev/catmulti-mvln.R lane`, log `mvln-lane.txt`, seed 13: logLik
-239.3980194005 against the sum of the univariate fits -239.3980194003,
relative difference 8.61e-13; two logistic_normal() responses 5.58e-13.
Summary rows `sigmab_Y, sigmac_Y, sigma_w, lncor_Y(b,c)`;
`lncor_Y(b,c)` 0.2797076 (univariate 0.2797077), Est.Error 0.05321804
(0.05321803). `set_prior("lkj(2)", class = "lncor", resp = "Y")` moves
it to 0.2754697; `resp = "w"` is refused ("No logistic-normal
correlation matches class=lncor, resp=w"). `frm_sample()` draws report
`lncor_Y__b__c` (posterior mean 0.2693, sd 0.0508) under the default
`lkj(1)`, disclosed as "lncor (resp = Y)". `simulate()` of any
multivariate fit stays refused ("not supported yet for multivariate
fits"), as on 0.71.0.

The frmtmb.sample half seen to fail alone (`dev/catmulti-mixarm.R`,
the lane's core with rellib-r9's frmtmb.sample, log `mixarm-base.txt`):
the draws were named `Y_lncor_raw_1` and carried no `lkj(1)` default;
`posterior_linpred()` was 5 x 200 and 4 x 200. With the lane's
frmtmb.sample: 12 and 4 passes, 0 failures (`mixarm-lane.txt`).

### The bin-2 rows (brief item 2)

`dev/brmsport-bin2-run.R` on the two generated files, rellib-r9 against
the lane library, records `bin2-<file>-<arm>.tsv`: `log_lik:504` and
`posterior_predict:349` do NOT flip. Their setup lines `log_lik:502`
and `posterior_predict:347` (`prep$family <- multinomial()`) now pass,
and the rows stop one line later, inside the bin-2 adapter:
"brms_adapt: prep$dpars has no 'mu', which frmtmb's multinomial dpar
'mu' reads". A bare `multinomial()` is the deferred placeholder whose
one dpar is `mu` until data resolve it; the adapter has to build the
family from `prep$data$Y`'s columns, or call `multinomial(K = ncat)`.
That file is lane bin2port's and was not edited. Both blocks then call
`dirichlet_multinomial()`, which frmtmb does not have, so the rows stay
blocked after that too. Both files: 52 and 54 passes, 0 failures, on
both arms; no other row moved.

### The cumulative tied-threshold NaN (coordinator item)

Cause confirmed: `ord_cumulative_logpmf()` evaluates the clamped
interior pair on every row and masks it with `(1 - i1 - iK) *`; an end
row's pair at a tie is `log(0) = -Inf`, and `0 * -Inf` is NaN on
doubles, folded away on the tape (rtmb-pitfalls item 22). Fixed with
lane bin2port's one-line change, `a + (i1 + iK)`.
`dev/catmulti-defect-cumulative.R` (copied from bin2port), logs
`defect-cumulative-base.txt` and `-lane.txt`: every NaN of the base
table becomes the exact value (thresholds 0 0 0, y = 1: -0.403186
-0.693147 -1.541008), and the sum at tau = (0, 0, 1) is -1.650688 on
the tape and on doubles. `dev/catmulti-tie-identity.R` (one process,
seed 20261014, 50 random problems per link, 40 rows, K = 5), log
`tie-identity.txt`, 0.71.0's body against the lane's at distinct
thresholds:

```
taped value identical 200 of 200; taped gradient identical 200 of 200; value on doubles identical 200 of 200
```

over logit, probit, probit_approx and cloglog. cauchit has no log-odds
form and takes the plain-difference path, which has no mask. A fit
reads only the taped value and gradient, so no fit moves.
hurdle_cumulative reads the same function, and its y = 0 rows were NaN
at a tie too (bin2port's reviewer, `dev/catmulti-bin2port-rev-cumfix.R`
copied); fixed by the same line and pinned. sratio, cratio and acat
form no difference of two distribution values, so a tie leaves them
finite (pinned, unchanged). The grouped `thres(gr = )` density
evaluates its interior pair only on interior rows, and the `cs()`
density per category; neither has the pattern. `test-ordinal-tie.R`
on rellib-r9: 31 failures in its 3 tests (all three links fail).

### Stretch items (brief item 6)

- `confint(parm = "lncor__b__c")` and `"rescor__y1__y2"` on a 2 x 2
  matrix are now aliases of the one chart coordinate, with the
  existing "internal scale" message. On rellib-r9 the advice of the
  refusal fails three ways (`dev/catmulti-profile-advice.R`, log
  `profile-advice.txt`): `hypothesis(fit, "cor_g__Intercept__x = 0",
  method = "profile")` finds no `b_cor...` without `class = NULL`, and
  with it says "is not linear in the parameters; use method = 'boot'";
  so does `sigma = 1`. The refusal now names `class = NULL` and
  `method = "boot"`; `test-portability.R` pinned the old advice and was
  updated.
- `nobs()` counting the rows `subset()` drops: NOT a defect. brms
  2.23.0's `nobs()` of `bf(y | subset(keep) ~ x)` with 1 of 300 rows
  dropped is 300, and its summary says "Number of observations: 300"
  (`dev/catmulti-brms-nobs.R`, log `brms-nobs.txt`); the review's 299
  is brms's Stan-data N, which `log_lik()` has as columns (299).
  Unchanged.
- `drop1()` on a family without `mu`: reads its scope from the first
  location predictor on the main formula, so a term leaves every
  category's logit (2 df on a three-category fit, LRT equal to the
  refit's to 1e-8 relative); refused by name when every location
  predictor has its own formula.
- `marginaleffects::predictions()` on a matrix response: NOT fixed. It
  fails inside marginaleffects's data.table step ("Column 2 ['Yc'] is
  length 600 but column 1 is length 200; malformed data.table",
  `dev/catmulti-me.R`, log `me.txt`) on both builds; on a categorical
  fit it works. A data method that hides the matrix column is its own
  change.

### Routes a multinomial mean now reaches

`dev/catmulti-mult-routes.R`, logs `mult-routes-base.txt` and
`-lane.txt`. `frm_linpred(type = "response")` gives the n x K counts;
`conditional_effects()` the per-category display (above);
`insight::get_predicted()` 600 values on 200 rows, the flattened counts,
where 0.71.0 refused (brms's own `get_predicted()` on a multinomial fit
errors: "'length = 2' in coercion to 'logical(1)'"); `emmeans(epred =
TRUE)` is refused ("not one number per grid row");
`frm_linpred(se.fit = TRUE, type = "response")` is refused ("not one
number per observation"); `dharma_residuals()` refused by name.
Unchanged: `residuals()` (refused as in brms), `predict()`,
`bayes_R2()`, `loo()`, `pp_check()`.

## 4. Tests that pin the changes, seen to fail on rellib-r9

Each run by `dev/catmulti-test1.R <arm> <pkg> <file>` (logs
`t-*-lane.txt` and `t-*-base.txt`); the base arm lists every failing
test uncapped.

```
test-refcat.R            lane pass=65 fail=0   base pass=0 fail=3 err=11 (all 11 tests)
test-refcat-brms.R       lane pass=14 fail=0   base pass=2 fail=1 err=3 (gated; the mvbf dirichlet test passes on base, as expected)
test-simplex-mvbf.R      lane pass=18 fail=0   base pass=4 fail=0 err=2 (the two logistic_normal tests)
test-ordinal-tie.R       lane pass=119 fail=0  base pass=88 fail=31 (all 3 tests)
test-corr2-alias.R       lane pass=11 fail=0   base err=3 (all 3 tests)
frmtmb.sample:
test-linpred-stack.R     lane pass=12 fail=0   base pass=2 fail=1 err=2 (dim 5 x 200 where 5 x 200 x 2)
test-simplex-mvbf-draws.R lane pass=4 fail=0   base err=1 (core refusal); lane core + base sample: fail=2 err=1
```

Behavioral failures on base, not only missing arguments: a brms
`categorical(refcat = "b")` object fitted with names `mub, muc`;
`fitted()` on `multinomial(K = 3)` refused; `fitted(scale = "linear")`
300 x 4; `posterior_linpred()` S x N; the NaNs of the tie table.

Existing tests updated for the new signatures or behavior:
`test-multinomial.R` (bare `multinomial()` and `fitted()` no longer
refused), `test-brms-families.R` (`multinomial(3)` to `multinomial(K =
3)` or `multinomial()`), `test-simulate-density.R` (`multinomial(K =
K)`), `test-numerical-robustness.R` (`categorical(K = 3)` for the
internal constructor), `test-portability.R` (the new advice),
frmtmb.sample `test-ordinal-mixture-draws.R` (a mixture's
`posterior_linpred()` is now 2 x n x 2, as brms gives it).

## 5. Test runs

Every test file of core and the seven extensions, one R process per
file, 8 at a time, `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`, lane
library ahead of rellib-r9 (`dev/catmulti-run.sh all1
dev/catmulti-jobs-all.txt`; summed from `all1.log`):

```
          frmtmb pass=20032 fail=0 err=0 skip=14 warn=0 files=240
 frmtmb.coupling pass=542 fail=0 err=0 skip=5 warn=0 files=11
      frmtmb.eam pass=1808 fail=0 err=0 skip=3 warn=0 files=30
   frmtmb.latent pass=360 fail=0 err=0 skip=2 warn=0 files=10
    frmtmb.learn pass=501 fail=0 err=0 skip=2 warn=0 files=15
      frmtmb.ode pass=549 fail=0 err=0 skip=1 warn=0 files=11
   frmtmb.sample pass=2929 fail=0 err=0 skip=1 warn=0 files=55
   frmtmb.spline pass=594 fail=0 err=0 skip=1 warn=0 files=16
RAN 388 of 388; 388 RESULT lines name the lane library
```

The skips are the environment-gated tiers: `test-scale.R` in each
package (`FRMTMB_SCALE_TESTS`), `test-drmtmb-agreement.R` (13) and
`test-fuzz.R` (1). The last two were rerun gated
(`FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true`, `extra1.log`):
131 and 2 passes, 0 failures, 0 skips. The gated brms-suite files show
no stale "now HOLDS" row: no ported row flips (section 3).

After the suite, `R/confint.R` (aliases, advice, `drop1()`),
`R/conditional-effects.R` and `R/diagnostics.R` changed, so every file
that calls `confint()`, `drop1()`, `step()`, `hypothesis()`,
`conditional_effects()` or `dharma_residuals()`, in core and every
extension, plus `test-compat.R` and `test-compat-register.R`, was
rerun (`ci2.log`, 75 files; `ce2.log`, 61 files): 0 failures but one,
`test-portability.R`, which pinned the old advice; updated, 96 passes.

OpenBLAS 0.3.26 with LAPACK routed to it (`dev/ciharden-openblas.sh
0.3.26 lapack`, `OPENBLAS_NUM_THREADS=4`, `dev/catmulti-run-ob.sh
ob326`, `ob326.log`), the files that test changed numerical code:

```
test-refcat.R 65, test-ordinal-tie.R 119, test-simplex-mvbf.R 18,
test-corr2-alias.R 11, test-multinomial.R 17, test-ordinal-fitted.R 129,
test-simulate-density.R 591, test-famgaps.R 100,
test-numerical-robustness.R 723, test-simplex-families.R 125,
frmtmb.sample test-linpred-stack.R 12, test-simplex-mvbf-draws.R 4:
every one fail=0 err=0; RAN 12 of 12
```

No test of mine decides a rounding-near branch: the tolerances are
ratios to the measured magnitude (1e-12 to 1e-4 relative) or to a
measured standard error, and the tie test compares against exact
values at points 1 to 2 units from any branch.

One session per package (`dev/catmulti-onesession.R`, `test_dir(...,
load_package = "none")`, ungated, every extension in the library
path; logs `onesession-*.txt`):

```
ONESESSION frmtmb pass=17648 fail=9 err=0 skip=190 warn=0 files=231
ONESESSION frmtmb.sample pass=2713 fail=0 err=0 skip=4 warn=0 files=50
```

The 9 core failures are 5 tests of `test-compat.R` and
`test-compat-register.R` ("every rule wins at least one pair", "rule
patterns name real kinds, groups, and features", "declared post-fit
methods are callable", "expects = accepts a name another package
already supplies", "a declared but unresolved rule side is reported,
not dropped"). The same 5 tests fail when the same worktree tests run
in one session on rellib-r9 (`onesession-frmtmb-base.txt`; there the
lane's new tests fail too, as they should). They are the known 0.71.0
one-session failures lane testiso owns; reported, not fixed.

## 6. What I did not do, and why

- The bin-2 adapter change that would let `log_lik:504` and
  `posterior_predict:349` reach frmtmb's multinomial: lane bin2port's
  file (section 3). `dirichlet_multinomial()` is not in frmtmb.
- `marginaleffects::predictions()` on a matrix response (section 3).
- A `refcat = NA` refusal: measured and decided against (section 3).
- `multinomial(K = 3)` keeps positional names `mu2, mu3` even when the
  matrix has column names; brms has no `K`. It is the spelling for a
  family built away from data (frmtmb.latent's `hmm()` components use
  it). `categorical(K = )` does the same.
- `categorical()` reads a factor's levels as given, unused levels
  included; brms's `extract_cat_names()` reads `levels(factor(y))`,
  which drops them. Pre-existing, not changed.
- A brms categorical category with its own formula cannot be checked
  against brms: brms-B (dev/famsimplex-findings.md) fires for the
  SECOND predicted category too, `bf(f ~ x, mub ~ z)` under
  `refcat = "c"`, not only the first, which corrects that record's
  "first predicted category".
- `insight::get_predicted()` on multinomial returns the flattened
  counts; brms's errors. Left as is; no frmtmb test pins it.
- `posterior_epred()` of categorical draws names the observation
  margin; brms does not. Pre-existing.

## 7. Defects found but not fixed

None in frmtmb beyond those above.

## 8. brms 2.23.0 defects met

### brms-C. logistic_normal() in mvbf() does not compile

- Seen: 2.23.0, 2026-10-14 (`dev/catmulti-brms-fits.R`, "mv-ln";
  stancode in `brms-probe.txt`).
- Repro: `brm(bf(Y ~ x, family = logistic_normal()) + bf(w ~ x,
  family = gaussian()) + set_rescor(FALSE), data = dd)` with `Y` a
  three-column composition.
- Expected: a fit. Actual: "Semantic error ... Identifier \"Llncor\"
  not in scope. Did you mean \"Lncor\"?" The generated quantities read
  `Llncor` and `nlncor` where the parameters and transformed data
  declare `Llncor_Y` and `nlncor_Y` (`stan-predictor.R`, the
  `gen_def` glue has no `{p}` suffix). With two logistic_normal()
  responses the generated quantities also declare `Lncor` and `lncor`
  twice. Renaming (`rename_family_cor_pars()`) would also give both
  responses' correlations the one name `lncor__b__c`.
- Workaround: frmtmb fits the model; the test repairs the two
  identifiers and checks brms's model block.
- Status: not reported.

## 9. Files touched

Core: `R/families.R` (categorical, multinomial, the cumulative fix,
`brmsfamily()` doc, docs section), `R/families-simplex.R` (refcat doc,
`lncor_blocks()` prefix, `lncor_values()`), `R/predict.R`
(`cat_probs()`, `fitted()` clamp and linear stacking,
`mu_stack_dpars()`, `check_mv_extra_family()`), `R/priors.R` (the
`lncor` default row and `resolve_lncor()` resp), `R/confint.R`
(aliases, advice, `drop1()`), `R/conditional-effects.R`,
`R/diagnostics.R`, `R/compat.R` (two rows), `man/frmtmb-families.Rd`,
`man/brmsfamily.Rd`, `man/drop1.frmtmb_fit.Rd`, `NEWS.md`, `README.md`,
`vignettes/brms-migration.Rmd`; tests `test-refcat.R`,
`test-refcat-brms.R`, `test-simplex-mvbf.R`, `test-ordinal-tie.R`,
`test-corr2-alias.R` (new), `helper-brms.R`, `test-multinomial.R`,
`test-brms-families.R`, `test-simulate-density.R`,
`test-numerical-robustness.R`, `test-portability.R`.

frmtmb.sample: `R/methods-draws.R`, `R/draws-brms.R`, `R/sample.R`,
`man/posterior_epred.Rd`, `NEWS.md`; tests `test-linpred-stack.R`,
`test-simplex-mvbf-draws.R` (new), `test-ordinal-mixture-draws.R`.

Scripts: `dev/catmulti-*.R` and `dev/catmulti-*.sh`.

Shared files other lanes may touch: `R/families.R`, `R/predict.R`,
`R/priors.R` (lanes prilabel and sampleprior), `R/compat.R`,
`R/confint.R`, `NEWS.md`, frmtmb.sample (`R/sample.R`: lane
sampleprior changes its default priors), `tests/testthat/helper-brms.R`.

## 10. R CMD check --as-cran

`dev/catmulti-check.sh`, built and checked in `dev/catmulti-check/`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
`R_LIBS` lane library, rellib-r9, user library (logs
`check-frmtmb.log`, `check-frmtmb.sample.log`):

- frmtmb: `Status: 1 NOTE`, the expected one ("Skipping checking math
  rendering: package 'V8' unavailable"); tests 11 min, OK.
- frmtmb.sample: `Status: OK`.

The check trees and tarballs were deleted after; the logs stay local.

## 11. Version bump and floor

Core: a minor bump (new arguments, BREAKING changes to
`multinomial()`'s signature and dpar names, new mvbf support).
frmtmb.sample: a minor bump, and it needs the new core: its
`posterior_linpred()` stacking reads the same families core's
`fitted()` stacks, and its multivariate `lncor` names and default
match core's `lncor_<resp>` and only exist where core fits
`logistic_normal()` in `mvbf()`. Its `Imports: frmtmb (>= 0.71.0)`
floor must rise to the new core version at release.

## 12. Fix round 1 (review `dev/reviews/2026-10-14-catmulti.md`)

Every new test was seen to fail on the round-1 behavior, not only on
rellib-r9: `dev/catmulti-revert-arm.R` puts each round-1 body back into
the loaded namespaces (this session only, nothing installed) and runs
the four files (log `revert-arm.txt`):

```
REVERTED test-refcat.R pass=71 fail=11 err=1
  not passing: refcat = NA counts the reference model's degrees of freedom
  not passing: conditional_effects(categorical = TRUE) reads fitted()
REVERTED test-simplex-mvbf.R pass=18 fail=2 err=1
  not passing: the correlation's name reads brms's Stan name of the response
REVERTED test-linpred-stack.R pass=13 fail=1 err=0
  not passing: transform = TRUE is refused on a two-category categorical too
REVERTED test-simplex-mvbf-draws.R pass=4 fail=1 err=0
  not passing: the draws name the correlation by brms's Stan name
```

The tests that existed before the round pass in the reverted arm, so
the new tests are what the fixes change. On rellib-r9 the new tests
fail as well (missing arguments; `t2-*-base.txt`).

1. **`lncor_<resp>` used the raw response name (MEDIUM).** Fixed: the
   prefix is `lncor_<brms_stan_name(resp)>` in `lncor_blocks()` and in
   frmtmb.sample's `draws_lncor_cols()`. `Y_a` gives `lncor_Ya__b__c`
   and `lncor_Ya(b,c)`, as brms's `Llncor_Ya`; `cbind(a, b, c)` gives
   `lncor_cbindabc__b__c`, which `hypothesis(class = NULL)` addresses
   (finite estimate and SE). `set_prior(resp = )` takes brms's Stan
   name (`Ya`) as well as frmtmb's (`Y_a`) for every class, through
   `resp_from_stan_name()` at the top of the prior resolver; the
   coordinator asked for brms's spelling, and frmtmb's own spelling
   stays valid. `default_prior()` still prints frmtmb's raw name (To
   file, below). The helper maps `Llncor_<stanname>` and every
   `<dpar>_<stanname>` suffix back to the response
   (`brms_lp_parts()`, `tests/testthat/helper-brms.R`), and a gated test
   checks the `Y_a` model against brms's repaired program.
2. **`refcat = NA` counted the flat shift as degrees of freedom
   (MEDIUM).** Fixed in `logLik()` (`loglik_df()`, `R/methods-fit.R`),
   which `AIC()`, `BIC()`, `extractAIC()`, `df.residual()`, `anova()`
   and `drop1()` read. The count is the same model's with the first
   category as reference: its frame is assembled from the fit's own
   formula and data (not fitted), and the difference in outer
   parameters is taken off; the result is cached. The families carry
   the reference twin as `refcat_first` (categorical, multinomial,
   dirichlet). With a category formula of its own the df is `NA`,
   with the reason as `attr(, "df_why")`, and `anova()` and `drop1()`
   refuse with it ("... with a formula of its own on mua the shift is
   only partly free ..."); a multivariate fit with such a response
   also gets `NA` with its reason. The review's numbers, pinned in
   `test-refcat.R` on its data (`rev_cat_data(23, 3, n = 300)` and
   `(21, 3)` with `z`): df 4 and 4 (was 6 and 4), AIC equal at
   622.6884 (was 626.69), `anova()` of the two copies Chi Df 0 (was
   2), `drop1()` Df 2 and 2 for both terms (was 3), dirichlet(NA) df
   equal to dirichlet()'s. Decided, against the review's "no prior on
   its location dpars" condition: the count applies with priors too,
   because `logLik()` is the likelihood, whose model is the reference
   model whatever penalty the fit carried.
3. **`transform = TRUE` on a two-category categorical (LOW).** Fixed:
   refused whenever `draws_mu_stack()` is non-NULL, in brms's words
   ("Valid distributional parameters are: 'mub'"); the second sentence
   now says the family has no `mu`. `dpar = "mub"` still works.
4. **The SE warning did not name the shift (LOW).** Fixed: a fit with
   a `refcat = NA` family adds "Under refcat = NA a common shift of
   every category's predictor leaves the likelihood unchanged, so the
   category probabilities are determined and the coefficients are not:
   set a reference category, or a prior on each predictor"
   (`se_lost_message()`, `R/se-check.R`), pinned by
   `allow_warnings(require = )`. The findings sentence in section 3 is
   corrected in place.
5. **`conditional_effects(categorical = TRUE)` on `categorical()`
   (LOW, pre-existing).** Done, it was small: `ce_fitted_bands()`
   adds `categorical()` to the families whose per-category band is
   read off `fitted()`, without changing its `method = "predict"`
   refusal. The test checks estimates and SEs against
   `fitted(newdata = )` (1e-10 relative). Against brms 2.23.0 on seed
   11 (`dev/catmulti-ce-compare-cat.R`, log `ce-compare-cat.txt`):
   300 of 300 rows matched, |frmtmb - brms| / brms se median 0.031,
   max 0.112; se ratio median 1.015, range 0.945 to 1.108.
6. **Clamped Wald intervals (INFO).** Documented in `?frmtmb-families`:
   `fitted()`'s intervals for `categorical()` and `multinomial()` are
   Wald intervals held inside the range, so a bound can be exactly 0.
   Filed below: an interval on the log or logit scale.
7. **Nits.** "Found by lane bin2port." removed from NEWS; the
   multinomial integer refusal says it once ("Family 'multinomial'
   requires integer responses. The counts must be non-negative whole
   numbers"); the docs say a two-category fit's linear-scale
   `fitted()` stays `n x 4`; `print()` of a family shows its reference
   ("Reference category: y2", "none (refcat = NA)"), placeholders
   included.

## 13. To file

For the consolidation to move into `dev/test-backlog.md` and
`dev/upstream-bugs.md`. Nothing is posted anywhere.

- **brms-C (upstream).** brms 2.23.0 cannot compile `logistic_normal()`
  in `mvbf()`: its generated quantities read `Llncor` and `nlncor`
  where the program declares `Llncor_<resp>` and `nlncor_<resp>`
  (`stan-predictor.R`, the `gen_def` glue has no `{p}` suffix); with
  two such responses `Lncor` and `lncor` are declared twice, and
  `rename_family_cor_pars()` would give both responses' correlations
  the one name `lncor__b__c`. Details in section 8.
- **brms-B correction (upstream).** The `use_glm_primitive_categorical()`
  failure (`invalid 'y' type in 'x && y'`) fires for any categorical
  category with a formula of its own, not only the first predicted
  one: `bf(f ~ x, mub ~ z)` under `refcat = "c"`, where `mub` is the
  second predicted category, stops `stancode()` too
  (`dev/catmulti-brms-b.R`, log `brms-b.txt`: `mua ~ z` and `mub ~ z`
  under `refcat = "c"`, and `muc ~ z` under the default reference, all
  three "invalid 'y' type in 'x && y'").
  `dev/famsimplex-findings.md` section 6 says "the FIRST predicted
  category".
- **`ranef()` of a categorical fit** is `25 x 2`; brms's is
  `25 x 4 x 2` (Estimate, Est.Error, Q2.5, Q97.5 by `mua_Intercept`,
  `muc_Intercept`). Pre-existing (review finding 8).
- **`coef()` of a categorical fit** is a list by dpar, so `coef(fit)$g`
  is `NULL`; brms's is `25 x 4 x 4`. Pre-existing.
- **`variables()` order** on a categorical fit is `b_mua_Intercept
  b_mua_x b_muc_Intercept b_muc_x`; brms's is `b_mua_Intercept
  b_muc_Intercept b_mua_x b_muc_x`. The names agree. Pre-existing.
- **`default_prior()` prints a multivariate response by its raw name**
  (`Y_a`), where brms prints its Stan name (`Ya`); `set_prior(resp = )`
  now takes both. Pre-existing for every class.
- **A log- or logit-scale interval for `fitted()`** of `categorical()`
  and `multinomial()`, so that a bound near 0 is not clamped to exactly
  0 (review finding 6: 57 of 900 lower bounds at 0 on a rare category).
- **The bin-2 adapter** for `log_lik:504` and `posterior_predict:349`
  (section 3; lane bin2port's file).
- **`marginaleffects::predictions()` on a matrix response** (section 3).

## 14. Fix round 1: runs

- Also changed by item 5: `band = "profile"` on a non-ordinal
  per-category display was refused in the ordinal family's words ("an
  ordinal category probability"); it now names the family. Two
  existing tests pinned the old wald refusal and were updated:
  `test-ce-bands.R` ("a nominal family gets bootstrap bands and honest
  refusals") and `test-brms-methods.R` ("a nominal per-category display
  needs a bootstrap band"); both fail on the fix and pass after the
  update (`all2.log`, then `ce3.log`: 170 and 1006 passes).
- Every file of core and the seven extensions, one process each,
  gated, after the fixes (`dev/catmulti-run.sh all2`, summed from
  `all2.log`; RAN 389 of 389, every RESULT line names the lane core):

```
          frmtmb pass=20084 fail=2 err=1 skip=14 warn=0 files=241
 frmtmb.coupling pass=542 fail=0 err=0 skip=5 warn=0 files=11
      frmtmb.eam pass=1808 fail=0 err=0 skip=3 warn=0 files=30
   frmtmb.latent pass=360 fail=0 err=0 skip=2 warn=0 files=10
    frmtmb.learn pass=501 fail=0 err=0 skip=2 warn=0 files=15
      frmtmb.ode pass=549 fail=0 err=0 skip=1 warn=0 files=11
   frmtmb.sample pass=2932 fail=0 err=0 skip=1 warn=0 files=55
   frmtmb.spline pass=594 fail=0 err=0 skip=1 warn=0 files=16
```

  The 2 failures and 1 error are the two pins above. After the update
  and the profile wording, every file that reaches
  `conditional_effects()` or `dharma_residuals()` in all packages,
  plus `test-compat.R` and `test-compat-register.R`, was rerun
  (`ce3.log`): RAN 61 of 61, 0 failures. Skips are the gated tiers, as
  in section 5.
- Gated `test-refcat-brms.R` (`lpreport2.txt`): 17 passes; the new
  `Y_a` check: `LPCHECK const 3.410605132e-13 grad 2.14e-05 ours
  80.36868313`.
- OpenBLAS 0.3.26 with LAPACK (`ob326b.log`): the 12 numerical files of
  section 5 plus `test-ce-bands.R` and `test-confint-anova.R`, RAN 14
  of 14, 0 failures (`test-refcat.R` 89, `test-simplex-mvbf.R` 26,
  `test-linpred-stack.R` 14, `test-simplex-mvbf-draws.R` 5).
- frmtmb.sample in one session (`onesession-frmtmb.sample-2.txt`):
  `ONESESSION frmtmb.sample pass=2716 fail=0 err=0 skip=4 warn=0
  files=50`.

Files touched in fix round 1, beyond section 9: `R/methods-fit.R`
(`loglik_df()`), `R/se-check.R`, `R/conditional-effects.R`
(`ce_fitted_bands()`, the profile wording), `R/priors.R`
(`resp_from_stan_name()`), `R/families.R`, `R/families-simplex.R`,
`R/confint.R`; frmtmb.sample `R/methods-draws.R`, `R/draws-brms.R`;
tests `test-refcat.R`, `test-simplex-mvbf.R`, `test-refcat-brms.R`,
`helper-brms.R` (`brms_lp_parts()` matches brms's Stan names),
`test-ce-bands.R`, `test-brms-methods.R`, frmtmb.sample
`test-linpred-stack.R`, `test-simplex-mvbf-draws.R`; scripts
`dev/catmulti-revert-arm.R`, `dev/catmulti-ce-compare-cat.R`,
`dev/catmulti-brms-b.R`. `R/methods-fit.R` is shared with lanes that
touch `logLik()`.
- `R CMD check --as-cran` after fix round 1 (`dev/catmulti-check.sh`,
  logs `check-frmtmb.log`, `check-frmtmb.sample.log`): frmtmb
  `Status: 1 NOTE` (the V8 math-rendering note), frmtmb.sample
  `Status: OK`.

## 15. Fix round 2 (re-check, finding 9 and the SE nit)

**9. `loglik_df()` failed on every transformed term (MEDIUM).** Round
1 rebuilt the reference-category model from `fit$frame$data_frame`,
the MODEL frame (`y`, `log(x2)`), so a term that reads a raw variable
(`poly()`, `log()`, `scale()`, `I()`, a term `update()` adds) lost its
df and AIC. Fixed: `frm()` builds the reference twin's frame while it
holds the original data, with the same `data2`, `na.action` and
`drop_unused_levels` (`subset()` is an addition term of the formula
and comes along), and stores its template lengths on the frame as
`refcat_first_outer`. `loglik_df()` reads that record and never
evaluates the formula again. `NA` with a reason stays for a category
with its own formula, a multivariate fit, and a fit built outside
`frm()` (which carries no record); each reason is worded for the user.
Ordinary fits are untouched: `refcat_first_outer()` returns `NULL` at
once, and `loglik_df()` returns `n_outer_est()` before reading
anything.

- New test (`test-refcat.R`, "refcat = NA counts the df of transformed
  terms too"), on `y ~ poly(x, 2) + log(x2)`, `y ~ scale(x) + I(x2^2)`,
  `y ~ log(x2)` and `update(f1, . ~ . + log(x2), newdata = dd)`: the
  df equals the reference fit's, AIC is finite and equal, `anova()`
  Chi Df 0. Against round 1's `loglik_df()` put back in-session
  (`dev/catmulti-revert-arm2.R`, log `revert-arm2.txt`): `REVERTED
  test-refcat.R pass=90 fail=11 err=1`, failing exactly this test and
  the SE-nit test below; 102 passes on the fix.
- Ordinary fits keep 0.71.0's numbers exactly. `dev/catmulti-df-sweep.R`
  (seed 31, n = 400, `poly(x, 2) + log(x2) + scale(x2) + I(x^2) + (1 |
  g)`) on the lane build and on rellib-r9, logs `df-sweep-lane.txt`
  and `df-sweep-base.txt`: gaussian, poisson, binomial with trials,
  Gamma (log), negbinomial, cumulative, categorical(), multinomial(K =
  3) and dirichlet() give the same df, AIC and BIC to 10 decimals
  (`diff` shows only the lane-only extra line). That extra line,
  `categorical(refcat = "b")`, has the same df (12) and an AIC 0.036
  lower (772.6221 against 772.6577). That is expected: with `(1 | g)`
  each category has its own independent group effect, so moving the
  reference changes the random-effect model, as in brms. Without
  group terms the two references give the same likelihood (section 3).
- Cost (`dev/catmulti-twin-time.R`, seed 101, 1e5 rows,
  `y ~ poly(x, 2) + log(x2)`, minimum of 3 interleaved rounds, log
  `twin-time.txt`): whole `refcat = NA` fit 9.25 s; the twin's count
  alone 0.14 s (1.5% of the fit); control, the NA model's own frame by
  `dry_run = "frame"`, 0.21 s. That is measurable but small, and only
  `refcat = NA` fits pay it; any other fit pays one `vapply()` over its
  responses.

**SE-sentence nit.** It was cheap to tell: the sentence now attaches
only when the lost parameters reach every category's predictor of the
NA family, which the common shift moves; an sd at 0 alone does not.
Pinned by a direct test of `se_lost_message()` (two of three
predictors lost: no sentence; all three: the sentence), which fails on
round 1's predicate (`revert-arm2.txt`).

**To file (added).** `drop1()` on ANY fit with a transformed term
fails, refcat or not: its refits read `fit$frame$data_frame`, the
model frame, which has no raw `x2` for `log(x2)`. On rellib-r9,
`drop1(frm(bf(y ~ x + log(x2)), data = d))` gives "The model uses `x2`,
which is not a column of `data`" (`dev/catmulti-drop1-base.R`, log
`drop1-base.txt`). Pre-existing, the same root as finding 9; not fixed
here because it needs the original data kept on the fit or a refit
route through `update(newdata = )`.

Fix round 2 runs. `frm()` changed, so every file of core and the seven
extensions ran again, one process each, gated (`dev/catmulti-run.sh
all3`, `all3.log`; RAN 389 of 389, every RESULT line on the lane core):

```
          frmtmb pass=20101 fail=0 err=0 skip=14 warn=0 files=241
 frmtmb.coupling pass=542 fail=0 err=0 skip=5 warn=0 files=11
      frmtmb.eam pass=1808 fail=0 err=0 skip=3 warn=0 files=30
   frmtmb.latent pass=360 fail=0 err=0 skip=2 warn=0 files=10
    frmtmb.learn pass=501 fail=0 err=0 skip=2 warn=0 files=15
      frmtmb.ode pass=549 fail=0 err=0 skip=1 warn=0 files=11
   frmtmb.sample pass=2932 fail=0 err=0 skip=1 warn=0 files=55
   frmtmb.spline pass=594 fail=0 err=0 skip=1 warn=0 files=16
```

OpenBLAS 0.3.26 with LAPACK (`ob326c.log`): `test-refcat.R` 102 and
`test-confint-anova.R` 35 passes, 0 failures. frmtmb.sample in one
session (`onesession-frmtmb.sample-3.txt`): pass=2716 fail=0.
`R CMD check --as-cran` for frmtmb: `Status: 1 NOTE` (V8 math
rendering). Files touched in fix round 2: `R/fit.R` (one call in
`frm()`), `R/methods-fit.R` (`loglik_df()`, `refcat_first_outer()`,
`frame_n_outer()`), `R/se-check.R`, `tests/testthat/test-refcat.R`;
scripts `dev/catmulti-revert-arm2.R`, `dev/catmulti-df-sweep.R`,
`dev/catmulti-twin-time.R`, `dev/catmulti-drop1-base.R`. `R/fit.R` is
a shared file.
