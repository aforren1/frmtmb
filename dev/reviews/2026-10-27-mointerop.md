# Review: lane mointerop (one- and two-value `mo()`, marginaleffects interop, rescor rows)

Reviewer: adversarial lane review, 2026-10-27. Tree:
`frmtmb-wt-mointerop` (branch `wt-mointerop` at 435e7ce6, uncommitted).
Library under test: `wt-mointerop-lib` (built 00:26 UTC, after the last
edit of `R/`), in front of `rellib-r14`; the "base" arm is `rellib-r14`
alone. I installed nothing. Scripts are `dev/review-mointerop-*`; their
logs are in `dev/review-mointerop-out/`. Every number below is from
those logs.

## Verdict: FIX

No blocker. The core of the lane is right: two values fit as brms's
`Jmo = 1` (equal to `lm()` on the code, also with a factor multiplier),
one value as `Jmo = 0`, fits with 3+ categories are unchanged (my rerun
of `dev/mointerop-equal.R`: 112 lines, identical on both arms),
`re_formula` and the padded-row fix work, and the gap-aware default
covariance is exact at observed categories on every structure I tried
(two gap terms, `mo(x) * z`, `(mo(x) | g)`, `cumulative()`,
`sigma ~ mo(x)`). One major: the new `vcov =` spellings that go
through a whole-outer-vector matrix (`"HC0"`, `~ g`,
`vcov(fit, full = TRUE)`) give `NaN` for every standard error, with no
warning, on exactly the fits this lane makes gap-aware or newly
accepts. Fix 1, fix or document 2 to 4, the rest at your discretion.

## Tests run (one file per process, NOT_CRAN=true)

`dev/review-mointerop-tests.sh`, logs `dev/review-mointerop-out/tests/`.
Every `lib:` line names the library of the arm.

| package | file | arm | gate | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|---|---|
| frmtmb | test-mo-coding.R | lane | | 155 | 0 | 0 | 0 | 0 |
| frmtmb | test-mo-coding.R | base | | 120 | 4 | 4 | 0 | 1 |
| frmtmb | test-marginaleffects-args.R | lane | | 29 | 0 | 0 | 0 | 0 |
| frmtmb | test-marginaleffects-args.R | base | | 10 | 16 | 1 | 0 | 0 |
| frmtmb | test-emmeans.R | lane | | 92 | 0 | 0 | 1 | 0 |
| frmtmb | test-se-check.R | lane | | 134 | 0 | 0 | 0 | 0 |
| frmtmb | test-message-uniqueness.R | lane | | 6 | 0 | 0 | 0 | 0 |
| frmtmb | test-compat.R | lane | | 767 | 0 | 0 | 0 | 0 |
| frmtmb | test-brms-suite-methods.R | lane | gated | 161 | 1 | 0 | 0 | 0 |
| frmtmb.sample | test-mo-coding.R | lane | gated | 50 | 0 | 0 | 0 | 0 |
| frmtmb.sample | test-loglik-autocor.R | lane | | 111 | 0 | 0 | 2 | 0 |

The one failure in `test-brms-suite-methods.R` is
`brmsfit-methods:924 now HOLDS but is recorded as 'divergence'`, as the
lane says. The recorded rows 921 to 949 (`dev/mointerop-out/rec-methods-
{base,lane}.tsv`) match the lane's table: 921 and 935 now run, 924
holds, 923, 927 and 937 fail on the class or the `data =` policy.

Bin-2 tier files (`dev/review-mointerop-bin2-run.R`, my copy of
`dev/brmsport-bin2-run.R` with its own staging directory and rellib-r14
behind the port library; logs `bin2-*-{lane,base}.txt`):

| file | arm | pass | fail |
|---|---|---|---|
| bin2 log_lik | lane | 97 | 4 |
| bin2 log_lik | base | 97 | 4 |
| bin2 posterior_predict | lane | 61 | 2 |
| bin2 posterior_predict | base | 61 | 2 |

Every failure is "now HOLDS but is recorded as 'unclassified'" (rows
log_lik:107, 109 and their value checks, max abs diff 8.88e-16;
posterior_predict:69, 72). The rows hold on rellib-r14 too: the port
needs no package change (see finding 8).

## Findings

### 1. major: `vcov = "HC0"`, `~ g` and a whole-outer-vector matrix give `NaN` everywhere, silently, on gap, one-value and boundary fits

The lane's default path is gap-aware (`G Vp G'` with the conditional
covariance) and zeroes a one-value coefficient. The `Vfull` path of
`interop_vcov()` uses the user's (or `vcov_cluster()`'s) matrix as both
`V` and `Vp`, and `vcov_cluster()` forms `B S'S B` with `B` holding `NaN`
in the lost rows, so every entry is `NaN`. marginaleffects then reports
`NaN` for every standard error, and nothing warns.

Evidence (`dev/review-mointerop-gap.R`, `-one.R`, `-vcov.R`, lane):

- Gap fit, seed 7, `y ~ mo(x)`, `x` in 1, 2, 5 (200 rows),
  `predictions(newdata = x = 1, 2, 5)`: default SE 0.1222210 0.1135415
  0.1212624; `vcov = 2 * get_vcov()` 0.1728466 0.1605719 0.1714909
  (sqrt 2, right); `vcov = "HC0"` NaN NaN NaN; `vcov = vcov(fit, full =
  TRUE)` NaN NaN NaN; 0 warnings. `lm(y ~ factor(x))` with `"HC0"`:
  0.1224166 0.1103619 0.1248376.
- One-value fit, seed 5, `y ~ mo(x1) + z`: `vcov = "HC0"` and `~ g`
  NaN NaN, 0 warnings (`vcov_cluster(fit, seq_len(n), "CR0")` is `NaN`
  in all 16 entries). `lm(y ~ z)` HC0: 0.08706383 0.11396959. Here
  `vcov(full = TRUE)` works (0.08701016 0.12146224) because only the
  `mox1` row is `NaN` and `zero_const()` clears it.
- RE fit at a boundary (seed 31, `y ~ x + (1 | g)`, sd at 0): `vcov =
  ~ g` NaN NaN NaN, 0 warnings.

On rellib-r14 these spellings were ignored, so this is not a regression
in numbers; it is the new interface failing on the fits that the same
lane advertises as "NaN only at an empty category". NEWS and
`?frmtmb-marginaleffects` promise `"HC0"`, `~ g` and full matrices
without this caveat.

Required fix: in `interop_user_vcov()`, for `"HC0"` and `~ g` form the
sandwich on the determined directions (bread `Vp`, the covariance
conditional on the lost directions, which `hyp_par_cov()` already
gives), then carry it through the same `G` as the default path; for a
user's full matrix with `NaN` in lost rows, either do the same
projection or refuse with a message naming the lost parameters. At the
least, warn when the returned matrix has `NaN` in a row a prediction
reads. Add a test on the gap fit (HC0 against `lm(factor(x))` HC0 at
1, 2, 5, `NaN` only at 3, 4) and on a one-value fit.

### 2. minor: one-value fits: downstream warts the lane did not test

`dev/review-mointerop-one.R`, seed 5, 200 rows, `y ~ mo(x1) + z`,
`x1 = 4`. What works: `summary()`, `fixef()`, `confint()`, `fitted()`
(SE 0.08701016 0.12146224, equal to `lm(y ~ z)` with the ML sigma),
`conditional_effects()`, `emmeans()`, marginaleffects with the default
`vcov`, a saved and reloaded fit, `update()`, `sigma ~ mo(x1)`, a
one-level ordered factor (brms's `standata()` gives `Jmo = 0` too),
and `frm_sample()`'s refusal from a fit and from a formula. What does
not:

- `predict()` (the predictive draws) drops all parameter uncertainty:
  "the covariance of the estimates could not be recovered from the
  Hessian, so the draws are taken at the estimates alone". The `NaN`
  row is a coefficient that no prediction reads; the lane already gives
  it variance 0 for marginaleffects (`interop_mo_const()`), and the same
  treatment fixes `predict()`. (Gap fits have the same fallback on both
  arms, pre-existing.)
- `logLik()` df 4 and AIC 658.5448 where `lm(y ~ z)`, and `lm()` with
  an aliased zero column, give df 3 and AIC 656.5448; `anova(f0, f1)`
  reports `Chi Df 1` for a term that cannot change the likelihood.
- `y ~ mo(x1) + z + (mo(x1) | g)` fits and prints `sd(mox1) = 1.00`
  (the start value) with `NaN` SE and two warnings. A variance of a
  column that is 0 on every row is not a parameter of the model.

Required fix: zero variance for `D = 0` coefficients in `predict()`'s
draw covariance; refuse `(mo(x) | g)` when the term has one value (or
document it); and either exclude `D = 0` coefficients from the df or
document the convention in `?frm`'s mo() section.

Is accepting one value right at all? Yes, on balance: brms builds the
model (`standata()` `Jmo = 0`, `default_prior()` lists `simo mox11`)
and fails only when Stan samples, and frmtmb already fits other
undetermined directions with `NaN` SEs and a warning. The ported row
924 depends on it. The warts above are what make it a half feature.

### 3. minor: the refusal of `"HC1"` gives a false reason; REML full matrices fail cryptically

The message (and `?frmtmb-marginaleffects`) says `"HC1"` to `"HC5"`
"need the hat matrix of a linear model". HC1 does not: it is HC0 times
`n / (n - k)`. Measured (`-vcov.R`, seed 31, 300 rows):
`sqrt(vcovHC(HC1) / vcovHC(HC0))` = 1.00335 = sqrt(300/298) for `lm`.
`"stata"`, which marginaleffects maps to HC1, gets the same message.

Also, a REML fit: `get_vcov(fit, vcov = vcov(fit, full = TRUE))` stops
with "the covariance does not line up with get_coef()", while the help
page says a full matrix is cut down. REML integrates the fixed effects
out, so its outer vector has no beta rows.

Required fix: either support HC1 (k = the number of mean-model
coefficients, as `sandwich` counts them for `lm`), or give the true
reason (the `n / (n - k)` factor needs a parameter count that a
distributional or random-effect model does not define). Name REML in
the full-matrix refusal.

### 4. minor: the help page's glm claim is wrong for non-canonical links

`?frmtmb-marginaleffects`: `~ g` "is what `sandwich::vcovCL()` gives a
glm by default". frmtmb's bread is the observed information;
`sandwich`'s glm bread is the expected information. They agree only for
canonical links. Measured (`-vcov.R` and `-misc.R`, seed 31, 300 rows,
`Gamma(link = "log")`): slope SE frmtmb HC0 0.04046535, `vcovHC(HC0)`
0.03587813 (ratio 1.128); `~ g` ratio 1.128. A hand-built observed-bread
sandwich gives 0.04046535, so the mechanism is that and not an error.
Bernoulli probit: ratio 1.001. Poisson with offset: 1.00000001. Gaussian
with weights: equal.

Required fix: say that the bread is the observed information, which
equals `sandwich`'s for gaussian and canonical-link fits and differs
for others (give the Gamma number). NEWS's own claims (gaussian HC0,
Poisson `~ g`) are accurate.

### 5. minor: an argument frmtmb ignores looks the same as `re_formula`

marginaleffects warns "not known to be supported" for every argument of
a `frmtmb_fit`, and the help page tells users the warning is expected
for `re_formula`. So a spelling frmtmb ignores is indistinguishable:
`predictions(f, newdata = nd, re.form = NA)` (lme4's spelling) returns
the conditional predictions 1.049222 1.388058 1.718826, where
`re_formula = NA` gives 0.7421286 1.1982772 1.6544258
(`dev/review-mointerop-reform.R`, seed 41, `(1 + x | g) + (1 | h)`).
brms's `sample_new_levels`, `dpar`, `resp`, `ndraws` and `draw_ids`,
which marginaleffects lists for `brmsfit`, are dropped the same way
(by reading `get_predict.frmtmb_fit()`; I ran only `re.form`).

Required fix: `get_predict.frmtmb_fit()` refuses, by name, dots it does
not honor among the names marginaleffects lists for brmsfit, lme4 and
glmmTMB (`re.form`, `allow.new.levels`, `sample_new_levels`, `dpar`,
`resp`, `ndraws`, `draw_ids`), or honors them.

### 6. minor: a gap fit with another lost parameter: `get_vcov()` is `NaN` everywhere, no warning

`dev/review-mointerop-misc.R`, seed 23, 300 rows, `y ~ mo(x) + z +
(mo(x) | g)`, `x` in 1, 2, 5: the fit stops at "singular convergence
(7)" with `cor = 1` and `theta_3` lost. Lane: every entry of
`get_vcov()` is `NaN`, 0 warnings, so `predictions()` at the observed
categories 1, 2, 5 is `NaN` with no message (the new warning fires only
when a row reads the gap). Base: the beta block is finite, the simplex
rows `NaN`, and a warning names the gap. `keep_nan` reads only the
diagonal of `pc$V`; here `Vp` itself is `NaN`. The fit is degenerate
and warns itself, so this is minor; a converged `(mo(x) | g)` gap fit
is right (seeds 51 and 52, 600 rows: finite at 1, 2, 5, `NaN` at 3,
`re_formula = NA` equal to `fitted(re_formula = NA)`).

Required fix: when `G Vp G'` has `NaN` rows the default path did not
expect, warn (or fall back to the base behavior with its warning).

### 7. nit: the gap warning counts marginaleffects' padding rows

`dev/review-mointerop-gapre.R`, 4 requested rows, `g` with 20 levels:
"1 of 24 predictions read a mo() category that no row is in". Count
rows with a positive `rowid` only.

### 8. minor: rescor rows: the port is sound, the `posterior_predict()` rows do not reach the broken path

Confirmed the lane's reported defect, identical on both arms
(`dev/review-mointerop-rescorpp.R`, seed 3, 300 rows, `rescor_matrix()`
0.820242): frmtmb.sample's `posterior_predict(ds)` is 300 x 300 (one
response, where brms returns draws x rows x responses), and the draws of
`y1` and `y2` at rows 1 to 5 correlate -0.052 -0.040 -0.005 0.082 0.026,
where core's `predict(summary = FALSE)` gives 0.795 0.809 0.861 0.831
0.818. `log_lik(ds)` is right: max relative difference to the bivariate
normal density at each draw 1.6e-15. So the defect is pre-existing,
and it is major for users of frmtmb.sample (`pp_check()` and any
predictive check of a rescor model is wrong), but it is not this lane's
regression.

The ported rows: brms's log_lik:107/109 assert a length and the value
checks compare the numbers (8.88e-16); that is a real test of
`rescor_row_loglik()`, which frmtmb.sample's `log_lik()` runs.
posterior_predict:69/72 assert only `dim(pred) == c(ns, nvars)`, as
brms's do, through core's joint draw. This is the bin-2 convention
(other adapters use core's `sim_draw()`), so the rows are not hollow by
the ledger's rules, but a pass there says nothing about frmtmb.sample.
All four rows also hold on rellib-r14 (tier runs above): they flip
because of the adapter, not a package change.

Required fix: at consolidation, the ledger note for posterior_predict:69
and 72 says the route is core's `predict()` and that frmtmb.sample's
`posterior_predict()` of a rescor fit fails the analogous shape and
joint-draw check; file the defect in `dev/test-backlog.md` as major.

### 9. nit: tier files are red until the verdicts are regenerated

Gated `test-brms-suite-methods.R` (1 fail, row 924), bin-2 log_lik (4)
and posterior_predict (2) fail in tier mode on the lane build because
the verdicts say the old outcome. Expected by process; consolidation
must regenerate `dev/brmsport-verdicts.tsv` before any tier run counts.

### 10. nit: a yardstick borrowed from another fit

`tests/testthat/test-mo-coding.R` lines 455 and 464 scale the gap fit's
error by `1e-4 * max(abs(se0 / ex - 1))`, where `se0` and `ex` belong to
the 4-category fit above. It is a measured ratio, but of a different
model. Use the gap fit's own simplex-free standard error.

## Checked and right

- Two values against brms 2.23.0 (`dev/review-mointerop-two.R`, seed
  5): `default_prior()` rows, names and order equal brms's (`b mox`,
  `simo mox1`; with `* z`, `simo mox:z1` before `simo mox1`); the prior
  strings differ by frmtmb's flat ML default, which is pre-existing.
  `conditional_effects()` grid 1, 2; `z` held at `x = 1`, and at 2 when
  the median is 2; new data 3, 0 and 1.5 refused. `emmeans` pairs
  -0.701 (SE 0.162; lm 0.163 with df). `mo(x) * f` on two values:
  logLik -322.7916 = `lm(y ~ I(x - 1) * f)`; brms 2.23.0's `standata()`
  on that formula errors ("replacement has 2 rows, data has 3",
  brms#1828). `variables()` lacks `simo_mox1[1]` (lane item 7, open).
- One value vs brms: an ordered factor with one used level of two gives
  brms `Jmo = 0` and frmtmb `D = 0`. Unused top or bottom levels of an
  ordered factor are dropped by both (`Jmo = D = 2` for levels 1:4 with
  2:4 or 1:3 used), so a gap at the top or bottom cannot occur
  (`dev/review-mointerop-brmsord.R`).
- Gap SEs (`dev/review-mointerop-gap.R`, `-gapre.R`, `-ord.R`): two
  gap terms (seed 21, `x` in 1, 2, 5 and `w` in 1, 3, 4) match
  `lm(factor(x) + factor(w) + z)` at all 9 observed cells to 6 or 7
  digits, logLik equal; `set_coef(get_coef())` round trip 0 difference;
  a perturbation of `zeta1_cum1` moves only the `x = 2` prediction;
  `mo(x) * z`, `cumulative()` and `sigma ~ mo(x)` with gaps: `NaN` only
  at gap rows (base: `NaN` everywhere).
- `vcov =`: Poisson with offset HC0 ratio 1.00000001; gaussian weights
  HC0 equal to `lm(weights =)`; distributional HC0 equals
  `vcov_cluster(seq_len(n), "CR0")`; nonlinear HC0 vs
  `sandwich(nls)` 0.999; RE fit full matrix equals the default
  (max rel 0); a prior fit is refused by `vcov_cluster()`'s reason;
  lower-case `"hc0"` accepted, as marginaleffects accepts it for `lm`.
- `re_formula` (`-reform.R`, both arms): `NA`, `~ 0`, `~ (1 | h)` and
  `~ (1 | g)` equal `fitted()` on a `(1 + x | g) + (1 | h)` fit;
  `slopes()` and `avg_comparisons()` honor it (0.456 = `fixef`);
  `allow_new_levels = TRUE` works; Poisson link and response scales
  equal `fitted()`. Base returned 63 rows of conditional predictions.
  The categorical padded-row case keeps 2 rows per category.
- NEWS: no lane names, no `dev/` paths, no em dashes; BREAKING lists the
  `mo()` change, `vcov =` and the gap `get_coef()` names.

## Pre-existing, not this lane (for the backlog)

- `fitted()` gives `Est.Error` 0, not `NaN`, at the empty category on
  the degenerate fit of finding 6 (both arms).
- `drop1()` omits `mo()` terms (both arms; also a 4-category fit).
- `predict()` of a gap fit drops parameter uncertainty (both arms).

## Short check of fix round 1 (2026-10-28)

Scripts: `dev/review-mointerop-libmatch.R`, `-rdmatch.R`, `-fix1.R`,
`-fix1b.R`, `-mv.R`, `-tests2.sh`; logs in `dev/review-mointerop-out/`
(`fix1-lane.txt`, `fix1b-base.txt`, `mv-{lane,base}.txt`, `tests2/`).
I installed nothing.

### Verdict: SHIP

Every finding of the first review is fixed or answered, and I could
not break the fixes. Three nits and one pre-existing defect are below.
None of them blocks the merge.

### The build is the tree

`wt-mointerop-lib` frmtmb was built 2026-10-28 02:38:48 UTC. A
deparse comparison of every top-level function in `R/` against the
namespace: 2115 defined, 2114 identical. The one left is
`covstruct_core`, which the namespace holds as a list, not a function.
`R/interop.R` was edited after the build (00:02 local), but only in
roxygen: the installed `frmtmb-marginaleffects` and `frm` help pages
render identically to the tree's `man/` files.

### Tests (lane library, one file per process, NOT_CRAN=true)

| file | pass | fail | error | skip | warn |
|---|---|---|---|---|---|
| test-mo-coding.R | 164 | 0 | 0 | 0 | 0 |
| test-marginaleffects-args.R | 57 | 0 | 0 | 0 | 0 |
| test-sandwich.R | 58 | 0 | 0 | 0 | 0 |
| test-interop.R | 45 | 0 | 0 | 0 | 0 |
| test-se-check.R | 134 | 0 | 0 | 0 | 0 |
| test-confint-anova.R | 35 | 0 | 0 | 0 | 0 |
| test-message-uniqueness.R | 6 | 0 | 0 | 0 | 0 |

### Verified (`dev/review-mointerop-fix1.R`)

- **Finding 1: gap fit.** Seed 7; `x` in 1, 2, 5; `g` in 20 groups of
  10 rows.

  | `vcov =` | frmtmb at x = 1, 2, 5 | `lm(factor(x))` at x = 1, 2, 5 | frmtmb at x = 3, 4 |
  |---|---|---|---|
  | `"HC0"` | 0.1224166, 0.1103619, 0.1248376 | the same | NaN |
  | `"HC1"` and `"stata"` | 0.1233452, 0.1111991, 0.1257846 | 0.1233451, 0.1111990, 0.1257846 | NaN |
  | `~ g` | 0.1082513, 0.1336043, 0.1179098 | 0.1087993, 0.1342808, 0.1185068 | NaN |

  The `~ g` ratio is sqrt(197/199), lm's `(n - 1) / (n - k)`, as the
  help page documents. The full matrix gives NaN with the new warning
  naming `zeta1_cum1`.
- **Finding 1: the conditional bread against a direct computation.**
  I took the generalized inverse of the numerical outer Hessian, which
  has two eigenvalues of 1.5e-4 against 62 to 2464. I formed
  `J ginv(H) S'S ginv(H) J'` over (beta, sigma, S_obs). Its maximum
  relative difference from `get_vcov(vcov = "HC0")` is 6.6e-8, the
  numerical-Hessian error.
- **Finding 1: one-value fit (seed 5).** HC0 is 0.08412656, 0.122077,
  equal to `lm(y ~ z)`. HC1 is 0.08455037, 0.122692, equal to lm.
  `~ g` is lm's vcovCL times sqrt(198/199). The full matrix gives
  0.08417197, 0.11864604.
- **Finding 1: boundary RE fit (seed 31).** `~ g` gives 0.0826883 and
  0.1223907, equal to `vcovCL()` of `glm(y ~ x)`, ratio 1 (it was
  NaN). HC0 and HC1 are refused because the random effect crosses the
  rows, as before.
- **Finding 1: non-gap fits are bitwise unchanged** against
  rellib-r14's `vcov_cluster()`. All 8 of 8 comparisons are
  `identical()`: gaussian HC0; gaussian, Poisson-with-offset,
  4-category `mo()` and nested `(1 | g)` with `~ g` (full matrices);
  and the lane's `get_vcov(vcov = "HC0")` and `~ g` against the base
  fixed-effect block.
- **Finding 2: one-value fits.**
  - df 3, AIC 643.2701 and BIC 653.1651, equal to `lm(y ~ z)`.
  - `anova()` gives Df 3 and 3, Chi Df 0.
  - `predict()` draws with no warning: Est.Error 1.1817 and 1.1834,
    against 1.1950 and 1.1904 for `y ~ z` (4000 draws, within Monte
    Carlo error).
  - The `(mo(x1) | g)` refusal is clear and reasonable; `(1 | g)` with
    `mo(x1)` still fits.
- **Finding 3: HC1.**
  - `k` is 3 for the gap fit, 2 for the one-value fit, 2 for `y ~ x`
    and 5 for `mo(k) + z` on 4 categories. These equal the column
    counts of `lm(factor(x))`, `lm(y ~ z)` and `lm(factor(k) + z)`.
  - The HC1/HC0 ratio on `mo(k) + z` is 1.008439, equal to lm's.
  - With `sigma ~ x`, `k` counts the mean coefficients only (2).
  - On a fit with random effects, HC1 is refused (the effect crosses
    the rows), so `n - k` never arises there. marginaleffects does the
    same or less for the other packages: glmmTMB stops with "`vcov`
    must be `TRUE`, `FALSE`, `"HC0"`", and lme4 warns and silently
    uses the default covariance for both `"HC1"` and `~ g`.
  - The REML full-matrix refusal names REML and offers `get_coef()`.
- **Findings 4 and 5: help text and `get_predict()`.**
  - The help page states the observed-information bread and gives the
    1.128 Gamma(log) ratio.
  - `predictions(dpar = "sigma")` equals `fitted(dpar = "sigma")` in
    estimate and SE (to 1e-8), and `avg_slopes(dpar = "sigma")` runs.
  - `resp = "y2"` equals `fitted(resp = "y2")`.
  - `re.form` and `ndraws` are refused with frmtmb's reason inside
    marginaleffects' wrapper.
  - `avg_comparisons(by =)`, `hypotheses()` and `avg_predictions(wts =)`
    still run.
- **Finding 6: seed 23.** The coefficient block is finite (equal to
  rellib-r14's) and a warning names the cause. Predictions at observed
  rows are still NaN, which the warning says. This is acceptable on a
  fit that stopped at singular convergence.
- **Finding 7.** The warning now reads "1 of 4 predictions" on 4 rows
  of a 20-level factor.
- **Finding 10.** The gap test now takes its yardstick from the gap fit
  itself (`frm_lp_basis(fg)` without the `zeta1` columns).

### Findings of the short check

S1. **nit:** the help page says HC1's `k` is "the number of
population-level coefficients". `interop_hc1_k()` also counts the
determined steps of each `mo()` simplex and drops one-value terms: 3
for `y ~ mo(x)` on the gap data, where there are 2 coefficients. The
code is right (it matches lm); the text should say so.

S2. **nit:** the refusal of `"HC2"` to `"HC5"` says leverage is "not
defined" for a fit with random effects or distributional parameters.
For a plain gaussian or glm fit, leverage is defined; frmtmb just does
not compute it. Say that instead.

S3. **minor, pre-existing (identical on rellib-r14), for the backlog:**
`predictions(fit)` on a multivariate fit without `resp` silently
returns the first response only. Seed 3, 200 rows,
`mvbf(y1 ~ x, y2 ~ x)`: 3 rows equal to `fitted()`'s `y1` column, with
no `group` column and no `y2` (`dev/review-mointerop-mv.R`). It should
return every response keyed by `group`, as for a categorical outcome,
or refuse without `resp`.

### The rescor record

The lane's findings file no longer records frmtmb.sample's rescor
`posterior_predict()` defect. Three things were removed:

- the subsection in section 6;
- item 1 of "Defects found, not fixed";
- the bullet in "What needs the user".

`dev/mointerop-rescor-pp.R` is deleted, and "Files touched in fix
round 1" lists it as removed. This review keeps the record: finding 8
above has the numbers, and the repro is `dev/review-mointerop-rescorpp.R`
(seed 3, 300 rows; logs `rescorpp-{lane,base}.txt`, identical on both
arms). The numbers: `posterior_predict(ds)` is 300 x 300; the
cross-response draw correlations at rows 1 to 5 are -0.052, -0.040,
-0.005, 0.082 and 0.026, against 0.795 to 0.861 from core's
`predict()`; `log_lik(ds)` is exact (1.6e-15). Consolidation should
file it in `dev/test-backlog.md` as major. I found nothing else missing
from the findings file: the remaining sections, the tables and the
other defects are intact.
