# Review: the conditional_effects lane (`wt-ce`)

Reviewer's ledger, written incrementally. Every number was reproduced
on this machine; where mine differs from the lane's, both are shown.

Environment: R 4.6.1, private library `scratchpad/rce-lib` (worktree
core 0.51.0 + `frmtmb.sample` from `extensions/`), main's core in
`scratchpad/rce-lib-main`, private Stan cache copied to
`scratchpad/rce-stan-cache`. Explicit `lib.loc=` / `--library=`
everywhere.

## Scope, and one correction to the brief

`git -C ...frmtmb-wt-ce diff --name-only 2210aa1` gives 20 tracked
files; untracked adds `dev/ce-findings.md` and
`tests/testthat/test-unpinned-seams.R`. Matches the brief.

**The brief says main is at `2210aa1`. It is at `9b9011a`.** `2210aa1`
is an ancestor, and the three commits on top touch `_pkgdown.yml`,
`docs/**` and `extensions/frmtmb.spline/tests/testthat/test-deriv.R`
only - nothing under `R/`, `tests/testthat/` or `man/`. So the merge
base is effectively unmoved for this lane and every diff below is
against `2210aa1` as instructed. Main is clean and I did not touch it.

**A second correction, to the lane's claim as the brief restates it.**
`predict_mean_se()` is NOT new: it is at `R/predict.R:1626` in
`2210aa1` and the lane's diff does not touch it. What the lane added is
the ROUTE to it - `mean_display`/`pred_dpar` at
`R/conditional-effects.R:1553-1555` and the band branch at `:1641`.
The joint delta method itself is pre-existing and was already reachable
through `predict(type = "response", se.fit = TRUE)`. This makes the fix
smaller and safer than advertised, not larger.

---

## 1. The new default estimate and band

### The identities reproduce exactly

`bf(y ~ x, zi ~ x) + zero_inflated_poisson()`, n = 300, and
`bf(y ~ x, hu ~ x) + hurdle_poisson()`, my own data and seeds:

| shape | `max|ce$estimate__ - predict(type="response")|` | ratio pt 1 | ratio pt 100 |
| --- | --- | --- | --- |
| ZIP | **0** (`identical()` TRUE) | 1 | 1 |
| hurdle | **0** | 1 | 1 |

Columns back: `x, y, cond__, effect1__, estimate__, se__, lower__,
upper__`. Every Wald lower bound strictly positive.

### The joint delta-method SE, checked independently

Hand-built, from `vcov(fit)` (4x4: `(Intercept), x, zi_(Intercept),
zi_x`) and central finite differences of
`m(theta) = (1 - plogis(b_zi0 + b_zi1 x)) * exp(b_mu0 + b_mu1 x)`,
against `predict(type = "response", se.fit = TRUE)$se.fit` over all
100 grid points:

```
max abs diff  1.379e-10
max rel diff  3.060e-10      (finite-difference precision)
mean identity max abs  0
```

And `ce$se__` is exactly `se.fit / m` (log band scale), max abs
**5.55e-17**; the reported `lower__` reproduces
`exp(log(m) - z se.fit/m)` to 1.11e-16. So `ce_band_scale()` does what
it documents.

**The "joint" is load-bearing.** Zeroing the mu/zi cross-blocks of
`vcov()` and redoing the same delta method changes the SE by up to
**68%** on this shape. A per-dpar band would have been visibly wrong.

### The 0.128 bootstrap gap is Monte Carlo noise, not a defect

I reproduced the lane's exact shape (`scratchpad/ce-after.R:64-71`:
`set.seed(11)`, `plogis(-1 + 0.5x)`, `exp(0.5 + 0.4x)`), and its Wald
widths to the digit (`0.35913956 0.33003752 2.50673325`).

Seed sweep at the lane's `boot = 200`, max relative width difference:

| seed | 4 | 7 | 13 | 17 | 23 | 31 | 42 | 99 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| max rel | **0.1278** | 0.1968 | 0.0773 | 0.0634 | 0.0899 | 0.1230 | 0.1983 | 0.1321 |

min 0.063, median 0.125, max 0.198. Seed 4 reproduces the lane's
0.1278 exactly and sits at the median, so the lane did not seed-shop.

The decisive check is the bootstrap SIZE, not the seed. Same fit,
same seed 4:

| B | max rel width diff |
| --- | --- |
| 200 | 0.1278 |
| 800 | **0.0409** |

The per-point differences shrink toward zero across the whole grid, not
just at the maximum. **0.128 is the percentile bootstrap's own sampling
error at B = 200 and it converges away.** It is the gap one expects at
that n and that B, and it is not a defect. The lane's phrasing ("within
0.128") states a one-draw measurement as if it were a bound; the
finding is sound, the wording oversells its stability.

Caveat worth recording: the gap is SHAPE-dependent, not just B- and
seed-dependent. On a steeper ZIP (mu slope 0.62, zi slope 0.82, same
n = 300) the same comparison at B = 200 gives 0.41 and 0.53 at two
seeds. See section 1b for whether that also converges.

### 1b. The tail gap is real on a steep shape and does NOT converge

Same comparison on a steeper ZIP (mu slope 0.62, zi slope 0.82, n = 300),
seed 17:

| B | 200 | 800 | 3000 |
| --- | --- | --- | --- |
| max rel width diff | 0.4065 | 0.3606 | **0.2851** |

The middle of the grid converges (0.03-0.05 at B = 3000); the two END
points do not (0.180 and 0.285). That residual is a genuine
Wald-vs-percentile difference where the estimate is poorly determined
and its bootstrap distribution is skewed: the Wald band is symmetric on
the log scale and the percentile band is not. It is in the
CONSERVATIVE direction (Wald wider), the point estimates are identical
(max abs 0), and `band = "boot"` / `band = "profile"` are documented
alternatives. Not a defect; worth one sentence in the help page's band
section, which currently implies the two agree.

---

## 2. The mixture hook

### The likelihood is untouched, measured rather than argued

The same script run against the lane's core and against main's core
(`rce-lib` vs `rce-lib-main`): a `bf(y ~ 1, theta1 ~ x) + mixture(g, g)`
fit and a three-component `theta1 ~ x, theta2 ~ x` fit.

| quantity | lane | main |
| --- | --- | --- |
| K2 `logLik` | -705.792686471644 | **identical** |
| K2 `AIC` | 1423.58537294329 | **identical** |
| K3 `logLik` | -1319.94619746563 | **identical** |
| `fitted()[1:3]` | -0.6396788770 0.8699484078 0.5954389018 | **identical** |
| `mixture_probs()` colMeans | 0.617023997657 0.382976002343 | **identical** |
| `simulate(seed = 5)[1:4]` | -1.0942054925 ... | **identical** |
| `sum(residuals())` | -0.000347025003868 | **identical** |

Only the theta REPORTING scale moves:

| quantity | lane | main |
| --- | --- | --- |
| K2 max abs resp minus plogis(eta) | **2.22e-16** | 1.68286399094 |
| K3 max abs r1 minus softmax1 | **1.11e-16** | 0.228355769690 |
| K3 max abs r2 minus softmax2 | **1.67e-16** | 0.264893995653 |
| K3 max of r1 + r2 | 0.720046264853 (< 1) | 0.459769264113 |

Claim (2) holds, and holds for K = 3 as well as the K = 2 the lane
measured.

### The hook's consumers are exactly right

`dpar_report_hook()` has exactly two call sites: `R/predict.R:1279`
(the user-facing `predict()`, and only under `type == "response"`) and
`R/conditional-effects.R:1513`. Nothing else.

The brief asks whether `fitted()` and `mixture_probs()` consult it.
They should NOT and do not, and that is correct rather than a gap:
neither reports a dpar. `fitted()` (`R/predict.R:1767`) returns the
response mean through `eval_dpars()`; `mixture_probs()`
(`R/families.R:3120`) returns the posterior class RESPONSIBILITIES
through `eval_dpars()` too. Both are density-facing, and both are
bit-identical to main in the table above.

Density-facing paths all read the natural scale:

- `predict_mean_response()` to `dpars_natural()` (`R/predict.R:811`)
- `conditional_effects(method = "predict")` to `dpars_natural()`
  (`R/conditional-effects.R:1690`)
- nonlinear `nl_dpar_refs` to `linkinv(predict(type = "link"))`
  directly (`R/predict.R:1241-1248`), which is the third path the lane
  fixed and the easiest one to have missed
- `predict_mean_se()` builds its dpars as
  `lp[["link"]]$linkinv(ed[["eta"]])` inline and never went near the
  reporting path
- `eval_dpars()`, `log_lik`, `simulate()`: untouched by the diff

### DEFECT (minor, new): the theta SE is one-predictor, so K >= 3 runs wide

`dpar_response$deriv` is `p (1 - p)`, the derivative of `p_k` with
respect to its OWN predictor, and `predict()` multiplies it by that
predictor's `se_eta` alone. For K = 2 there is no other theta and the
rule is exact. For K >= 3 the softmax also depends on the other theta
predictors (`dp1/dtheta2 = -p1 p2`) and those terms are dropped.

Measured on the K = 3 fit, `theta1`, against a finite-difference joint
delta method over both theta predictors' coefficients using `vcov(fit)`:

```
reported se vs one-predictor rule   max rel 2.64e-10   (it IS that rule)
reported se vs JOINT delta          max rel 0.1726
                                    min rel 0.0891
direction: reported is WIDER (conservative)
K = 2 reported vs exact delta       max rel 4.10e-10   (exact)
```

The reported response-scale standard error of a mixing weight runs
**9 to 17 percent wide on a three-component fit**. It is conservative,
it replaces a number that was previously meaningless, and the code
comment at `R/predict.R:1324-1327` is honest that it is "the same
one-predictor rule the link inverse gets below". But the lane measured
K = 2 only, where the shortcut is exact, so this is unmeasured and
undocumented. Punch list.

---

## 3. re_formula = NULL (sleepstudy)

`bf(Reaction ~ Days + (Days | Subject)) + gaussian()`:

- `Subject` column present and **all NA** (a factor with the fit's
  levels), so the frame says which group it is.
- point estimate `identical()` to the population curve, max abs **0**.
- band strictly wider at **every** grid point.

The band is wider by EXACTLY the random-effect contribution. Hand-built
`sqrt(x' V x + z' S z)` with `z = (1, Days)`, `V` the 2x2 fixed-effect
block of `vcov()`, and `S` the `VarCorr()` covariance
(565.512, 11.056; 11.056, 32.682):

```
hand se vs cn$se__     max abs 1.42e-14   max rel 2.66e-16
hand pop se vs cp$se__ max abs 1.95e-14
```

`se__` goes 6.63 to 24.69 at grid point 1, and 14.22 to 60.12 at point
100. Nothing but `z' S z` is added; no uncertainty in `theta` itself is
propagated, which is the standard conditional-variance convention and
matches `lp_extra_var()` (pre-existing, untouched by this lane).

`conditions =` with an existing level still conditions on that level:
`conditions = list(Subject = "310")` gives a curve equal to
`predict(re.form = NULL)` on Subject 310 to **max abs 0**, and 88.05
away from the population curve. `ce_build_nd()` tests
`nm %in% names(cset)` before the `na_vars` branch, which is why.

For contrast, main's silent default (first observed level, 308) sits
84.50 from the population curve on my fit. The lane records 86.55; the
difference is which fit, not which behavior, and the defect is the same
either way.

---

## 4. Ordinal and nominal against brms

Reproduced on the tier's own `r12a` (cumulative) and `r13` (nominal)
shapes with real `brmsfit` objects built from the private Stan cache.

`categorical = TRUE`:

| | frmtmb | brms |
| --- | --- | --- |
| key | `x:cats__` | `x:cats__` |
| columns | x, y, cond\_\_, cats\_\_, effect1\_\_, effect2\_\_, estimate\_\_, se\_\_, lower\_\_, upper\_\_ | **the same names in the same order** |
| nrow | 300 | 300 |
| effects attribute | `x:cats__` | `x:cats__` |

`categorical = FALSE` on the ordinal fit:

| | frmtmb | brms |
| --- | --- | --- |
| key | `x` | `x` |
| nrow | 100 | 100 |
| estimate 1 to 4 | 1.25552494263 1.26473640665 1.27421074528 1.28395078673 | **the same** |

Max abs 4.44e-16, max rel 2.67e-16 over the whole curve. brms emits its
documented warning ("Predictions are treated as continuous variables
... Please set 'categorical' to TRUE").

Nominal family, `categorical = FALSE`:

- brms: **errors**, "Please set 'categorical' to TRUE."
- frmtmb: errors, naming the reason ("the expected category number
  needs ORDERED categories, and a nominal family's are not ordered").

The lane's claim, that brms refuses its own default for a nominal
family and warns for an ordinal one, is exactly right.

---

## 5. ranef() keys: the behavior change, and one thing NEWS does not say

Measured against main's core on the same two models.

| | main | lane |
| --- | --- | --- |
| `names(ranef(fit))` | Days \| Subject | Subject |
| `ranef(fit)$Subject` | **NULL** | the matrix |
| `ranef(fit)[["Days \| Subject"]]` | the matrix | **NULL** |
| `attr(ranef(fit)[[1]], "term")` | absent | Days \| Subject |
| `names(coef(fit))` | Subject | Subject (unchanged) |
| `names(VarCorr(fit))` | Days \| Subject | Days \| Subject (unchanged) |
| `as.data.frame(...)$grp` | Days \| Subject | Days \| Subject (unchanged) |

Two terms on one factor, `(1 | g) + (0 + z | g)`: the lane gives two
entries both named `g`, `term` attributes `1 | g` and `0 + z | g`, and
`VarCorr()` and the `grp` column both still keyed by block. `r2$g`
returns the FIRST entry silently, which `man/ranef.Rd` documents
("index by position when a model has them").

**The NEWS bullet is accurate but incomplete.** It says the block label
"rides along in each matrix's `term` attribute, is what
`as.data.frame()` puts in its `grp` column, and is still `VarCorr()`'s
key". All three verified true. What it does not say is the thing that
breaks existing code silently: **`ranef(fit)[["Days | Subject"]]` now
returns NULL rather than a matrix.** That is the migration hazard, a
NULL flows on into arithmetic as `numeric(0)` rather than erroring, and
a BEHAVIOR CHANGE bullet should name it. Punch list.

---

## 6. Merge against the sibling lanes

Both siblings sit at `2210aa1`, the same base. Three-way `git merge-file`
(ce = ours, `2210aa1` = base, sibling = theirs):

| file | vs spline-core | vs priors |
| --- | --- | --- |
| `R/predict.R` | **0 conflicts** | 0 (priors does not touch it) |
| `R/families.R` | **0 conflicts** | 0 (priors does not touch it) |
| `R/methods-fit.R` | 0 | n/a |
| `R/conditional-effects.R` | 0 | n/a |
| `tests/testthat/test-simulate-ergonomics.R` | n/a | **0 conflicts** (both edit it) |
| `vignettes/brms-migration.Rmd` | n/a | **0 conflicts** (both edit it) |
| `NEWS.md` | **1 conflict** | **1 conflict** |

Every merged file PARSES, and so do the spline lane's own current
copies: they are not mid-edit. The only conflicts are `NEWS.md`, where
all three lanes open a `# frmtmb (development version)` section.
Trivial and expected.

**Correction to the brief:** the priors lane does not edit `R/fit.R`.
Its tracked set is `NAMESPACE`, `NEWS.md`, `R/priors.R`,
`R/simulate-new.R`, three man pages, one vignette and four test files.
It does not touch `R/predict.R` or `R/families.R` at all, so there is
no predict/families overlap with priors to report. It DOES share
`tests/testthat/test-simulate-ergonomics.R` and
`vignettes/brms-migration.Rmd` with this lane, which the brief did not
mention; both merge clean.

### Semantic overlap with spline-core: real, adjacent, compatible

`frm_lp_basis()` is at `spline/R/predict.R:2942`, inside a block
appended at the END of the file (`@@ -2733,3 +2757,412 @@`). It is
nowhere near `dpars_natural()` (`:755`) or `dpar_report_hook()`
(`:787`). No overlap there.

The real adjacency is INSIDE `predict.frmtmb_fit()`. The spline lane
replaces the inline variance assembly

```
var_eta <- pmax(rowSums((A %*% V) * A), 0)
ev <- lp_extra_var(object, ed, use_re); ... for (gv in ev$gp) ...
```

with `lp_basis_out(object, jc, ..., lp_extra_var_vec(...), ...)`, and
the ce lane inserts its hook branch in the statement IMMEDIATELY after
(`out <- if (!is.null(hook)) ...`). They merge textually clean, and the
merged result is semantically coherent:

- `lp_extra_var()` still EXISTS in the spline lane (`:1460`); the new
  `lp_extra_var_vec()` (`:2986`) wraps it, and `predict_mean_se()`
  (`:1590`) still calls `lp_extra_var()` directly.
- In the merged file the spline's `lp_basis_out()` produces `var_eta`
  and `se_eta` at lines 1340-1346, and the ce lane's hook consumes
  `se_eta` at 1352-1353. The contract between them, a length-n `se_eta`
  with NA on non-estimable rows, is preserved by both.

One coupling to flag to whoever merges SECOND, not a defect today:
this lane makes the DEFAULT `conditional_effects()` display route
through `predict_mean_se()` for zero-inflated, hurdle, `trials()` and
truncated families. Before, that function was only reachable from an
explicit `predict(type = "response", se.fit = TRUE)`. So any spline
lane change to the extra-variance path becomes visible on the default
plot for the first time. The spline lane leaves `predict_mean_se()`
alone, so nothing breaks now; it is a reason to re-run the ce band
tests after the second merge, not to hold either lane.

---

## 7. The two "unpinnable" seams: the claim is RIGHT

I built my own copy with BOTH lines deleted (`scratchpad/rce-src2` to
`rce-lib2`) and went further than the lane did. Probe family: a
`family_finalize()` that replaces the `sigma` link with a data-derived
`scaled_logit`, ALSO replaces the `mu` link (identity to log), and in a
second variant reorders `family$dpars`.

Five call paths, both variants, both libraries:

| path | with the lines | without |
| --- | --- | --- |
| `get_prior(route = "sample")` | table | **identical** |
| `get_prior(route = "fit")` | table | **identical** |
| `par_template()` no prior | template | **identical** |
| `par_template(prior = list(beta = prior_normal(0, 2)))` | template | **identical** |
| `par_template(prior = set_prior("normal(0,2)", class = "b", dpar = "sigma"))` | template | **identical** |

The last two matter: they are the ONLY paths that reach
`resolve_prior_input()`, and the lane's own test calls `par_template()`
with no prior, so the lane never exercised them either. They are still
identical.

The structural reason, measured rather than asserted. Diffing the
parsed response spec against the carried one:

```
fields differing: family, dpars
primary_dpars parsed : mu       carried: mu        <- NOT refreshed
family$dpars  parsed : mu,sigma carried: sigma,mu  <- IS refreshed
dpars[[1]]$link      : identity -> log
dpars[[2]]$link      : log -> scaled_logit(0,4.343)
```

So the carry DOES change the spec. It is unobservable because
`get_prior()` reads only `rspec$primary_dpars` (`R/priors.R:760`) and
`rspec$nlpars` (`:762`) from the response spec, and `resolve_priors()`
reads only `fit$frame[["par_template"]]`, the FRAME, never the spec's
families or links.

**Verdict: "not pinnable from outside" is correct**, and the lines are
not dead code. They are defensively right for a future stage that reads
a dpar's link from the spec. Pinning the helper's contract instead of
inventing an observable is the right call, and `test-unpinned-seams.R`
does that.

One nit: the file's comment says the measurement was made "with a
family whose `family_finalize()` replaces a link, and again with one
that swaps which dpar is primary". A finalize CANNOT swap which dpar is
primary; `primary_dpars` is never refreshed, as the same comment goes
on to say. The variant actually available is reordering
`family$dpars`, which is what I ran. Wording, not substance.

---

## 8. The remaining claims, spot-checked

All on a plain `bf(y ~ x + z) + gaussian()` unless noted.

| claim | measured |
| --- | --- |
| (4) moderator held at EXACT `mean +/- sd` | held `-1.054522591116474 -0.010683912860323 1.033154765395827`; **max abs difference from exact 0**. `signif(, 3)` would have moved it by 4.52e-3, which is the lane's own number |
| (4) label only rounded | `effect2__` levels `1.03, -0.01, -1.05`, two decimals, descending |
| (C) two-way row order | first effect varies SLOWEST (`x` repeats 3 times per value while `z` cycles), brms's order |
| (5) `int_conditions` values | `list(z = c(low=-1, mid=0, hi=1))` gives z in `-1 0 1`, `effect2__` levels `hi, mid, low` |
| (5) `int_conditions` function | `function(v) quantile(v, c(.25,.75))` reproduces the quartiles exactly |
| (5) unknown `int_conditions` name | warns "int_conditions names no variable of the model data: nope" |
| (5) unknown dots | warns "conditional_effects() is ignoring unknown argument(s): nosucharg", against the right function |
| (C) `method = "posterior_epred"` | `identical()` to `method = "epred"` |
| (C) `method = "posterior_linpred"` | refused by name, points at `dpar =` |
| (6) `cond__` always present | present with a single condition set, one level |
| (1) gaussian band bit-identical | `se__`, `lower__` and `upper__` all **max abs 0** against the plain link-scale band |
| (9) `mo()` steps by one | 4 rows at `0 1 2 3` (was 100) |
| (C) `trials()` held at 1 | message emitted once, `nt` held at 1, estimate a probability in (0.234, 0.874) |

### The exclusion table

`brms_exclusions()` at `tests/testthat/helper-brms-methods.R:551` now has
**22 rows: 0 D, 18 P, 4 C**, exactly the lane's claim. Counted by hand
from the table:

| list | rows | class |
| --- | --- | --- |
| `brms_ce_shapes` | r5, r13 | P (2) |
| `brms_ce_shapes` | r12e | C (1) |
| `brms_linpred_shapes` | r12e, r13, r17 | C (3) |
| `brms_meanlink_shapes` | r15, r16, rC16, r12a-e, r13, r17 | P (10) |
| `brms_resid_shapes` | r13, r12a-e | P (6) |

Eleven D rows are gone. Both re-added rows carry a reason that is not
the old finding: `r12e` under finding 13 (a `cs()` term stores `"csz"`,
not the variable) and `r13` under finding 9 (a nominal family has no
thresholds, so the Wald band is refused by name and its own block
covers the shape under `band = "boot"`). Both are honest
reclassifications, not laundering: I confirmed the nominal refusal
above in section 4, and `r13`'s own block does exist.

---

## 9. Defects I found that the lane did not

### D-A. `band = "boot"` silently drops the random-effect variance under `re_formula = NULL`

**The most serious finding, and it is created by this lane's item 7.**

`R/conditional-effects.R:1550` sets `na_vars` so the grid's grouping
column is `NA`; `:1562` then hands that grid to `ce_boot_draws()`,
which predicts through `ce_boot_one()` (`:540`) with
`re.form = re_formula, allow_new_levels = TRUE`. A NEW level's
conditional modes are zero in every bootstrap refit, so no refit ever
draws a group effect and the bootstrap band carries coefficient
uncertainty ONLY.

sleepstudy, `resolution = 6`, `boot = 100`, `seed = 11`:

| band | width across the 6 grid points |
| --- | --- |
| population (`re_formula = NA`), wald | 26.00 26.69 31.20 38.21 46.59 55.73 |
| `re_formula = NULL`, wald | 96.78 107.89 131.89 163.18 198.34 235.65 |
| `re_formula = NULL`, **boot** | **24.61 26.81 31.32 38.33 46.55 55.89** |

`boot / waldPop` = 0.947 1.004 1.004 1.003 0.999 1.003. **The
`re_formula = NULL` bootstrap band IS the population band**, to within
its own Monte Carlo error. It is about 4x narrower than the wald band
for the same stated quantity, and the point estimate is identical
across all three (251.40519), so nothing on the returned object says
the argument stopped mattering.

This is a regression in kind, not only in degree. On main the same call
gave 76.6-199.1 (2.9-3.6x the population band): wrong for a different
reason, but visibly not the population band. On the lane it is exactly
the population band, so `re_formula = NULL, band = "boot"` and
`re_formula = NA, band = "boot"` now return the same interval.

What makes it worse: **the package's own refusals send users down this
exact path.**

- `R/conditional-effects.R:1521` - `band = "profile"` with
  `re_formula` refuses and says "use `band = "wald"` or `"boot"` to
  condition on random effects".
- `R/conditional-effects.R:1527` - an ORDINAL per-category display with
  `re_formula` refuses `wald` outright and says "with re_formula use
  `band = "boot"`". So on an ordinal fit conditioning on a new group,
  `boot` is the ONLY band available and it is the one that drops the
  group variance. There is no way to get a correct new-group band
  there at all.

Fix: either draw a new group's effects from the fitted covariance once
per bootstrap replicate (which is also what the help page says brms
does per posterior draw), or refuse
`band = "boot"` with `re_formula = NULL` by name and say the wald band
is the one that carries the group variance. The second is a few lines
and is honest; the first is the right answer.

### D-B. `effect2__` collapses two moderator values that round to the same 2 dp

`ce_effect2()` (`R/conditional-effects.R:330`) builds the display factor
from `round(v, 2)` and takes its levels from `sort(unique(r))`. When two
grid values round together, the factor gets ONE level for two distinct
curves, and the user's own labels are dropped because the guard
`length(nm) == length(lv)` no longer holds.

```
int_conditions = list(z = c(lo = 0.001, hi = 0.002))
z values kept    : 0.001 0.002        (correct)
effect2__ levels : 0                  (ONE level for two curves)
effect2__ nlevels: 1
names lo / hi    : silently dropped
rows             : 200                (data is right; only the label is wrong)
```

`plot()` groups on `effect2__`, so the two curves are drawn as one
series. The data frame is correct, the display is not. This is newly
reachable because `int_conditions` is newly implemented: the default
`mean +/- sd` rarely collides, but a user-supplied pair of close
quantiles easily does.

Fix: make the factor's levels the ordered VALUES and use the rounded
number only as the label text, so distinct values always get distinct
levels.

### D-C. The mixing weight's response-scale SE is one-predictor (K >= 3)

Covered in section 2. `R/families.R:2865` and `:3451` (`p * (1 - p)`),
consumed at `R/predict.R:1329`. Exact at K = 2, 9-17% wide at K = 3,
conservative. Unmeasured by the lane, which tested K = 2 only.

### D-D. `ranef()`'s re-key breaks `frmtmb.sample`, which the lane did not run

See section 10. Five real assertion failures in the sibling package,
and `ranef()` on a `frmtmb_draws` object changes its keys too
(`extensions/frmtmb.sample/R/methods-draws.R:154` delegates to core's
`ranef()` per draw), with no NEWS entry on either side for the sampling
surface.

### D-E. The `frmtmb.sample` draws method is half-migrated, wider than recorded

The lane records this as "keeps the OLD column set". It is more than
that. `conditional_effects.frmtmb_draws()`
(`extensions/frmtmb.sample/R/conditional-effects-draws.R`) shares
`ce_grids_build()` at `:45`, so it SILENTLY INHERITED four grid changes:

- brms's two-way row order (`ce_build_nd()`)
- exact `mean +/- sd` moderator values instead of `signif(, 3)`
- `mo()` stepwise grids (4 points, not 100)
- `trials()` held at 1, **including the new `message()`**, which now
  fires from a sampling call that never used to emit it

while keeping the old behavior on everything the method assembles
itself:

- `:81`, `:91` `d <- g$nd[g$ev]` - old narrow columns, no held values,
  no `effect1__` / `effect2__`
- `:98` `cond__` only when `length(cond_sets) > 1L` - the old rule
- `:87` `cats__` appended last - old column order
- `:102` `ce_finalize()` called without `cats_key`, so an ordinal draws
  display is still keyed `"x"`, not `"x:cats__"`
- `:41` `ce_cats_display()` - `categorical =` is not honored and is not
  even a formal
- `:45` no `int_conditions`, no `na_vars` - so `int_conditions` is
  unimplemented there and `re_formula = NULL` still silently takes the
  first observed level, the very defect item 7 fixed in core

So the two surfaces now disagree in BOTH directions. The lane's
one-line summary understates it.

### D-A, the ordinal case: a documented argument that does literally nothing

Making the worst case concrete. `bf(y ~ x + (1 | g)) + cumulative()`,
25 groups, n = 400:

```
conditional_effects(fit, re_formula = NULL)
  ERROR: the ordinal per-category delta method is written for the
  population-level curve; with re_formula use band = "boot", or ask
  for dpar = "mu"

conditional_effects(fit, re_formula = NULL, band = "boot", boot = 80, seed = 3)
  vs the same call with re_formula = NA:
    max abs difference in WIDTH    : 0
    max abs difference in ESTIMATE : 0
    g column in the returned frame : all NA
```

So on an ordinal fit the Wald band is refused BY NAME and the user is
sent to `band = "boot"`, which returns a frame **bit-identical** to the
population one - while setting the grouping column to `NA` to say it
conditioned on a new group. The argument is accepted, documented,
recommended by the package's own error message, and is a no-op.

---

## 10. Runs

One `test_file()` per process, private library, explicit `lib.loc=`.
My harness parents the test environment on `asNamespace("frmtmb")`;
without that, internals like `ord_tau_from_raw()` are invisible and the
tier reports 9 spurious errors. Worth knowing for anyone repeating this.

### Gated tiers (`FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true`, private Stan cache)

| tier | result | lane's claim |
| --- | --- | --- |
| `test-brms-methods.R` | **47 tests, 950 assertions, 0 fail, 0 error, 0 skip**, 133 s | 47 / 950 / 0 - **exact** |
| `test-brms-likelihood.R` | **32 tests, 352 assertions, 0 fail, 0 error, 0 skip**, 91 s | 32 / 352 - **exact** |
| `test-brms-likelihood.R` at `2210aa1` (main's core, `2210aa1`'s test file) | **32 tests, 351 assertions, 0 fail** | "351 before, plus one added" - **exact** |
| row 17 (mixture) before | 2 assertions, 0 failed | at 1e-8 |
| row 17 (mixture) after | passes inside the 352 | unchanged by the diff |

Row 17's block is byte-identical in the diff (the only edit to that
file is the `ranef()` key assertion in "check C: row 6"), so the
mixture rows are the same assertions at the same tolerance on both
sides, and both pass.

### Surviving exclusion rows, by class

**22 rows: 0 D, 18 P, 4 C.** Enumerated in section 8. Matches.

### Named files

| file | tests | assertions | fail |
| --- | --- | --- | --- |
| `test-message-uniqueness.R` | 1 | 6 | 0 |
| `test-bracket-access.R` | 3 | 8 | 0 |
| `test-ce-bands.R` | 19 | 139 | 0 |
| `test-ce-facets.R` | 6 | 92 | 0 |
| `test-effects.R` | 9 | 54 | 0 |
| `test-methods-audit.R` | 10 | 57 | 0 |
| `test-simulate-density.R` | 50 | 435 | 0 |
| `test-unpinned-seams.R` | 3 | 23 | 0 |
| `test-mvn-mixture.R` | 9 | 277 | 0 |
| `test-ordinal-fitted.R` | 13 | 87 | 0 |

There is no `test-conditional-effects*.R` or `test-mixture*.R` in this
repository; the brief names files that do not exist. The ce tests are
`test-ce-bands.R`, `test-ce-facets.R` and `test-effects.R`; the mixture
coverage is `test-mvn-mixture.R` plus mixture blocks inside 16 other
files.

### Full core suite, one file per process, by name

**110 of 110 files ran. 1097 tests, 6160 assertions, 0 failures,
2 errors, 88 skips.** No file missing (`comm` against the file list is
empty).

The 2 errors are both in `test-influence-plot.R` and are MY harness,
not the lane: `local_mocked_bindings()` needs pkgload, which a bare
`test_file()` does not provide. Under `R CMD check` the same file
passes - `checking tests ... [371s] OK`.

The lane reports 6168 assertions; I get 6160. The 8-assertion gap is
exactly the two `test-influence-plot.R` tests that my harness aborted
before they finished asserting. Consistent, not a discrepancy.

### `R CMD check --as-cran`

`_R_CHECK_CRAN_INCOMING_=false`, pandoc 3.8.3 from RStudio's quarto
tools on PATH.

- Run 1 (`R CMD build --no-build-vignettes`): `checking tests ...
  [371s] OK`, examples OK, `--run-donttest` OK. Two vignette WARNINGs,
  both artifacts of my own `--no-build-vignettes` (no `inst/doc`).
- Run 2 (full build, vignettes included): see the appended result.

### Roxygen

`roxygenise(load_code = "installed")` on a clean copy: `man/` **byte
identical** to the worktree's, `NAMESPACE` identical. Idempotent, as
claimed. (`load_code = "source"` fails on `habit_prep`; use the
installed loader.)

### `frmtmb.sample` suite, one process

**136 tests, 892 assertions, 5 FAILURES, 1 error.**

The 1 error is the pkgload artifact again (`test-parallel-chains.R`).
The **5 failures are real and are caused by this lane's `ranef()`
re-key**:

| file:line | assertion | actual |
| --- | --- | --- |
| `test-sample-direct.R:67` | `expect_named(ranef(ds_form), "1 \| g")` | `"g"` |
| `test-sampling-ported.R:560` | `expect_named(re_d, "1 \| g")` | `"g"` |
| `test-sampling-ported.R:561` | `dim(re_d[["1 \| g"]])` | `NULL` |
| `test-sampling-ported.R:562` | `colnames(re_d[["1 \| g"]])` | `NULL` |
| `test-sampling-ported.R:566` | compares draws ranef to `ranef(fit)[["1 \| g"]]` | `-Inf >= -Inf` (both NULL) |

`ranef.frmtmb_draws()`
(`extensions/frmtmb.sample/R/methods-draws.R:154`) delegates to core's
`ranef()` once per draw, so the sampling package's OWN public output is
re-keyed too. Neither package's NEWS mentions the sampling surface.

**The lane did not run this suite.** Its verification table covers the
two tiers and the core suite only.

---

## 11. The blast radius of the default-display change is larger than NEWS says

The headline bullet ends: "Zero-inflated, hurdle, `trials()`, truncated
and mixture fits all move; **every other family is unchanged, bit for
bit**."

That last clause is false. I ran the default `conditional_effects()`
across seven families against main's core and diffed the printed
estimates, standard errors and lower bounds at 15 digits:

| family | lane vs main |
| --- | --- |
| gaussian, poisson, bernoulli, Gamma, negbinomial, student | **bit-identical** |
| **lognormal** | **MOVED** |

```
lognormal, main : est -0.605200987419253 ... 1.320171637097145
lognormal, lane : est  0.627163713767872 ... 4.300902819018206
```

Main was plotting `linkinv(eta_mu)` with lognormal's identity link on
`mu`, i.e. the LOG-SCALE LOCATION, as if it were the response. It goes
NEGATIVE for a strictly positive response. The lane plots
`exp(mu + sigma^2/2)`, the actual mean, which agrees with `fitted()`
and `predict(type = "response")`. **The fix is right and it is a bigger
win than the lane claims** - this was a live defect on lognormal that
the original review never listed.

The gate is `mean_is_mu()` (`R/predict.R:639`): FALSE for any family
carrying a `post$mean_fn` whose body is not literally `dpars[["mu"]]`.
Enumerating the built-in families through it, **14 move**:

```
asym_laplace, beta_binomial, binomial, cox, hurdle_gamma,
hurdle_lognormal, hurdle_poisson, lognormal, shifted_lognormal,
zero_inflated_asym_laplace, zero_inflated_beta,
zero_inflated_binomial, zero_inflated_negbinomial,
zero_inflated_poisson
```

and 23 do not (acat, bernoulli, beta, categorical, compois, cratio,
cumulative, exgaussian, exponential, Gamma, gaussian, geometric, huber,
inverse.gaussian, nbinom1, negbinomial, poisson, skew_normal, sratio,
student, tweedie, von_mises, weibull).

The NEWS enumeration covers the zero-inflated and hurdle families, and
`binomial` / `beta_binomial` under "trials()". It does NOT cover
**lognormal, shifted_lognormal, asym_laplace and cox**, and the "every
other family is unchanged, bit for bit" clause actively tells those
users their plots did not move. They did.

---

## 12. `R CMD check --as-cran`

`_R_CHECK_CRAN_INCOMING_=false`, pandoc 3.8.3 from RStudio's quarto
tools on PATH, private library.

Run 1 (`R CMD build --no-build-vignettes`) - `Status: 1 ERROR,
3 WARNINGs, 2 NOTEs`. Every one is an environment or harness gap on
THIS box, not a lane defect:

| item | cause |
| --- | --- |
| `checking files in 'vignettes' ... WARNING` | my `--no-build-vignettes`, so no `inst/doc` |
| `checking package vignettes ... WARNING` | same |
| `checking PDF version of manual ... WARNING` | **no LaTeX on this box** (`pdflatex` not on PATH, no MiKTeX) |
| `checking PDF version of manual without index ... ERROR` | same |
| `checking HTML version of manual ... NOTE` | "Skipping checking math rendering: package 'V8' unavailable" |
| `non-standard things in the check directory ... NOTE` | `frmtmb-manual.tex`, left by the failed LaTeX run |

The substantive checks all passed: **`checking tests ... [371s] OK`**,
`checking examples ... OK`, `checking examples with --run-donttest ...
OK`, Rd usage/contents/cross-references OK, dependencies OK.

Run 2 (full build, vignettes included) - `Status: 1 ERROR, 1 WARNING,
2 NOTEs`. Both vignette WARNINGs cleared:

```
* checking installed files from 'inst/doc' ... OK
* checking files in 'vignettes' ... OK
* checking package vignettes ... OK
* checking tests ... [329s] OK
* checking re-building of vignette outputs ... [276s] OK
```

The four remaining items are EXACTLY the LaTeX and V8 ones:
`checking PDF version of manual ... WARNING`,
`checking PDF version of manual without index ... ERROR`,
`checking HTML version of manual ... NOTE` (V8 unavailable), and
`non-standard things in the check directory ... NOTE`
(`frmtmb-manual.tex`, left behind by the failed LaTeX run). Nothing
that depends on the package's own code, documentation, tests or
vignettes is anything but OK.

**I cannot reproduce the lane's "Status: OK - no ERROR, no WARNING, no
NOTE" on this box, and I do not think that is a finding against the
lane**: the three residual items need LaTeX and V8, neither of which is
installed here. Everything that depends on the package's own code and
documentation is clean.

---

## 13. What the wald path gets RIGHT, for contrast

D-A is confined to `band = "boot"`. The intersection of item 1 (the new
expected-response display) and item 7 (`re_formula = NULL`) works
correctly on the wald path. `bf(y ~ x + (1|g), zi ~ x) +
zero_inflated_poisson()`, 30 groups, n = 450:

```
estimate == predict(type = "response", allow_new_levels = TRUE)  max abs 0
population se__ : 0.2141 0.1782 0.1563 0.1796 0.2815
new-group  se__ : 0.7992 0.7903 0.7857 0.7907 0.8198
wider everywhere: TRUE      estimates identical: 0
g column all NA : TRUE      all lower__ > 0: TRUE
```

`predict_mean_se()` reaches `lp_extra_var()` at `R/predict.R:1649`, so
the joint delta method over every dpar AND the new-level variance are
both in the band. That is the design working. The bootstrap path simply
never got the same treatment.

---

## 14. Item 10, second half: `se()`'s unused sigma

| | main | lane |
| --- | --- | --- |
| `predict(type = "response", dpar = "sigma")` on `y \| se(s)` | **1** | **0** |
| `predict(type = "link", dpar = "sigma")` | 0 | 0 (unchanged) |
| `sigma(fit)` | 0 | 0 (unchanged) |
| `se(s, sigma = TRUE)`, response dpar sigma | 0.00010187 | **0.00010187 (unchanged)** |

The reporting hook fires only for the mapped-out case; a model that
really does estimate the extra sigma is untouched. No collateral damage.

---

## 15. Housekeeping

Main was READ ONLY: one `git archive 2210aa1` into scratch and one
`install.packages()` from the checkout. It is clean.

**Main moved during this review, and not by me**: `9b9011a` ->
`316a28b` (`Merge branch 'wt-gddm-ref'`). Everything it added is under
`extensions/frmtmb.ddm/` and `dev/reviews/`; nothing under core `R/`,
`tests/`, `man/` or the root `NEWS.md`. So the ce lane's merge base and
every conflict result in section 6 are unaffected. Both sibling
worktrees are still at `2210aa1`.

The ce worktree still has exactly the lane's 20 tracked changes and its
2 untracked files, plus `dev/review-ce.md` which is mine.

---

## 16. Punch list

Ordered by what should block a merge.

### Must fix before merge

**P1. `band = "boot"` drops the random-effect variance under
`re_formula = NULL`.**
`R/conditional-effects.R:1550` (`na_vars`), `:1562` (`ce_boot_draws`),
`:540` (`ce_boot_one`).
On sleepstudy the band is 0.24x the wald band and 1.00x the POPULATION
band; on an ordinal fit it is bit-identical to the population frame
(max abs 0 in estimate AND width) while the grouping column is set to
`NA` to claim otherwise. The package's own refusals at `:1521` and
`:1527` route users onto this path, and at `:1527` it is the only path
allowed. Fix: draw the new group's effects from the fitted covariance
per bootstrap replicate, or refuse the combination by name and point at
`band = "wald"`. Either way `:1527`'s message has to change, because
today it recommends the broken combination.

**P2. `ranef()`'s re-key breaks the `frmtmb.sample` suite.**
`R/methods-fit.R:633-643`, consumed by
`extensions/frmtmb.sample/R/methods-draws.R:154`.
Proven by running the sibling suite against both cores:
**main's core 897 assertions / 0 failures; the lane's core 892 / 5
failures** - the 5 lost assertions are exactly the 5 failures.
`test-sample-direct.R:67`, `test-sampling-ported.R:560,561,562,566`.
Fix the five strings in this lane or land the two lanes together, and
give `frmtmb.sample` a NEWS line: its own `ranef()` output is re-keyed.

### Should fix before merge (cheap)

**P3. `effect2__` collapses distinct moderator values that round
together.** `R/conditional-effects.R:330` (`ce_effect2`). Build the
factor's levels from the ordered VALUES; use `round(v, 2)` only as the
label text. Newly reachable because `int_conditions` is newly
implemented.

**P4. The default-display NEWS bullet names the wrong set of families.**
`NEWS.md`. "Zero-inflated, hurdle, `trials()`, truncated and mixture
fits all move; every other family is unchanged, bit for bit" is false:
`mean_is_mu()` moves **14** families, and **lognormal,
shifted_lognormal, asym_laplace and cox** are not in the bullet's list.
I measured lognormal moving from -0.605..1.320 to 0.627..4.301. Name
them, and say the lognormal case was a defect of its own (main plotted
a negative curve for a strictly positive response).

**P5. The `ranef()` NEWS bullet does not say what breaks.** `NEWS.md`.
It says what is KEPT - the `term` attribute, the `grp` column,
`VarCorr()`'s key, all verified true - but not that
`ranef(fit)[["Days | Subject"]]` now returns `NULL`. That is the silent
migration hazard.

### Follow-ups, not blockers

**P6. The mixing weight's response-scale SE is one-predictor for
K >= 3.** `R/families.R:2865` and `:3451`, consumed at
`R/predict.R:1329`. 9-17% wide on a K = 3 fit, conservative, exact at
K = 2. Fix with the full softmax Jacobian or document it. The lane
measured K = 2 only.

**P7. The band help implies Wald and bootstrap agree.** They agree in
the middle of the grid; at the extremes of a steep shape the Wald band
stays ~28% wider at B = 3000 and does not converge. One sentence.

**P8. `frmtmb.sample`'s draws method is half-migrated.**
`extensions/frmtmb.sample/R/conditional-effects-draws.R:41,45,81,87,91,98,102`.
It silently inherited four grid changes through `ce_grids_build()` (row
order, exact moderator values, `mo()` steps, `trials()` held at 1 WITH
its new `message()`) while keeping the old frame shape, the old
`cond__` rule, the old ordinal key, and the old `re_formula` /
`categorical` semantics. No NEWS entry anywhere.

**P9. Three claim-wording corrections.**
- `dev/ce-findings.md`: `predict_mean_se()` is presented as part of the
  fix. It is pre-existing at `2210aa1` and untouched; only the ROUTE to
  it is new.
- `dev/ce-findings.md` / `NEWS.md`: "within 0.128 relative width of a
  200-refit bootstrap" is one draw, not a bound. Median 0.125 over 8
  seeds, range 0.063-0.198, falling to 0.041 at B = 800. "Bootstrap
  Monte Carlo error at B = 200" is the stronger and true claim.
- `tests/testthat/test-unpinned-seams.R`: a `family_finalize()` cannot
  "swap which dpar is primary"; `primary_dpars` is never refreshed. The
  available variant is reordering `family$dpars`.

---

## 17. Every edit I made

**One file, created by me, nothing else touched:**

- `C:\Users\adf44\source\r\frmtmb-wt-ce\dev\review-ce.md` (this file).

No file in the worktree was modified or deleted. No commit. The main
checkout was read only and is clean. The sibling worktrees were read
only. Everything else lives under `scratchpad/rce-*`. I killed no
processes.

---

## VERDICT: GO-WITH-FIXES

The lane's ten fixes are real, and where I could check them
independently they are better than the lane's own evidence shows. The
ZIP and hurdle identities are exact (max abs 0). The joint
delta-method SE matches a finite-difference computation from `vcov()`
to 3e-10 relative, and the cross-dpar covariance it keeps is worth 68%
of the SE - dropping it would have been visibly wrong. The gaussian
band is bit-identical (`se__`, `lower__`, `upper__` all max abs 0). The
mixture hook moves the reporting scale and NOTHING else: `logLik`,
`AIC`, `fitted()`, `mixture_probs()`, `simulate()` and `residuals()`
are bit-identical to main on both a K = 2 and a K = 3 fit, and the hook
has exactly two call sites, both user-facing. `re_formula = NULL`'s
band is wider by exactly `z' S z`, to 2.7e-16 relative, and a
user-supplied `conditions =` still pins an observed level (max abs 0
against `predict(re.form = NULL)`). The ordinal work matches brms's
keys, columns and numbers to 4e-16, and brms really does refuse
`categorical = FALSE` for a nominal family. Every count reproduces:
47/950, 352 (351 before), 110/110 files, 22 exclusion rows at
0 D / 18 P / 4 C, roxygen idempotent. The merges against both siblings
are clean and every merged file parses. The "unpinnable seams" claim is
right, and I verified it through two call paths the lane never
exercised. The lane also fixed a lognormal defect nobody had listed.

It is not a GO because of two things the lane did not measure.

**P1 is a live defect this lane created.** `re_formula = NULL` with
`band = "boot"` returns the population band - on an ordinal fit,
bit-identically so - while marking the grouping column `NA` to claim
otherwise. The package's own error messages route users there, and on
an ordinal fit it is the only band allowed, so there is no way to get a
correct new-group band at all. An interval silently 4x too narrow on a
documented argument is worse than the defect item 7 set out to fix,
because the old behavior at least produced a visibly different band.
The wald path is correct, so the fix is contained.

**P2 leaves the monorepo red.** The `ranef()` re-key breaks five
assertions in `frmtmb.sample` - proven by 897/0 against main's core
versus 892/5 against the lane's - and silently re-keys that package's
own `ranef.frmtmb_draws()` output. The lane did not run the sibling
suite; its verification table covers the two tiers and the core suite
only.

P3 through P5 are cheap and should ride along; P4 in particular,
because the headline BEHAVIOR CHANGE bullet currently tells lognormal,
shifted_lognormal, asym_laplace and cox users that nothing moved for
them, and it did.

None of this argues against the design. The default display SHOULD be
the expected response, the mixing weight SHOULD report a probability,
and `re_formula = NULL` SHOULD mean a new group. Fix P1 and P2 and this
is a GO.

### The `frmtmb.sample` draws-method gap: wait for its own lane

The COLUMN gap should wait. It needs `ce_frame()` exported through
`R/sampling-api.R`, a sibling lane's file, and the draws method's
output is self-consistent today even though it is no longer consistent
with core's. Blocking on it would couple two lanes over a cosmetic
divergence.

But keep that separate from P2, which is NOT the draws-method gap and
must not be waved through with it. P2 is five failing assertions caused
by a core change in THIS lane. A cosmetic column divergence can wait
for a lane; a red sibling test suite cannot. And P8 - the four grid
behaviors the draws method silently inherited, including a new
`message()` now firing from sampling calls - is already user-visible
and should get a `frmtmb.sample` NEWS line in whichever lane lands
first.

---

# Punch re-check, 2026-09-05

Against the lane's punch round (`dev/ce-findings.md:242`) and my verdict
at `:1003` above. Same rules: private library `scratchpad/rce-lib`
(core + `frmtmb.sample` reinstalled from the worktree), main's core in
`rce-lib-main`, private Stan cache copy, `rce-` prefix only, no commit.

Scope grew from 20 tracked files to **26**: `extensions/frmtmb.sample/`
gains `NEWS.md`, `R/methods-draws.R` and two test files;
`man/mixture.Rd` and `man/predict.frmtmb_fit.Rd` are new. Untracked is
unchanged (`dev/ce-findings.md`, `dev/review-ce.md`,
`tests/testthat/test-unpinned-seams.R`).

## P1 - FIXED, and verified against the `z'Sz` I derived

The mechanism is a placeholder level: `ce_boot_grids()`
(`R/conditional-effects.R:675`) puts an observed level in the grid so
the design maps it, and `ce_draw_new_levels()` (`:698`) overwrites that
level's coefficients per replicate.

**The decisive check.** I pulled the bootstrap replicate matrix off the
returned object and compared the per-grid-point standard deviation of
the draws to the hand-built `sqrt(x' V x + z' S z)` from section 3 of
this review - the same `V` (fixed-effect block of `vcov()`) and the
same `S` (`VarCorr()`'s `[[565.512, 11.056], [11.056, 32.682]]`).
sleepstudy, resolution 6, **seed 2024** (the lane used 11), `boot = 400`:

```
boot draw sd  : 24.3376 27.5750 33.7118 41.4805 50.1281 59.2711
hand sqrt(x'Vx + z'Sz) : 24.6880 27.5240 33.6459 41.6284 50.5985 60.1157
ratio         : 0.9858 1.0019 1.0020 0.9964 0.9907 0.9860
wald new se__ : 24.6880 27.5240 33.6459 41.6284 50.5985 60.1157
z'Sz alone    : 23.7805 26.6684 32.6906 40.4714 49.1830 58.4105
```

The bootstrap replicate spread now reproduces the wald new-level
standard error to within 1.4%, and the random-effect term dominates it
(23.78 of 24.69 at grid point 1). That is the property that was missing,
measured directly rather than inferred from a width ratio.

Band widths at my seed:

| band | grid points 1..6 |
| --- | --- |
| wald pop | 25.998 26.692 31.204 38.205 46.586 55.728 |
| wald new | 96.775 107.892 131.889 163.180 198.343 235.649 |
| boot pop | 24.423 25.563 30.044 35.929 43.478 51.997 |
| boot new | 103.018 110.528 131.630 157.219 191.687 228.063 |

`boot new / boot pop` **4.22 4.32 4.38 4.38 4.41 4.39** (was 1.00
throughout); `boot new / wald new` 1.06 1.02 1.00 0.96 0.97 0.97. The
lane reports 3.56-4.18 and 0.88-0.98 at `boot = 100`; mine at
`boot = 400` sits closer to 1, which is the expected direction.
Estimate unchanged (max abs 0 against the population curve), `Subject`
still all `NA`.

**The ordinal case**, `y ~ x + (1 | g)` + `cumulative()`, 25 groups,
n = 400, **seed 777**, `boot = 300`, 12 rows:

```
ratio new/pop : 2.93 2.98 4.62 5.94 2.95 2.44 3.06 2.47 3.96 4.91 3.52 3.12
wider at EVERY row : TRUE   (min 2.44, max 5.94)
estimate new vs pop: max abs 0
g column all NA    : TRUE
```

The lane reports 1.80-4.57; mine is 2.44-5.94 at a different seed and
B. The property that matters is that the previous result - width
**bit-identical** to the population band, max abs 0 - is gone.

### The skipped block types are exactly `lp_extra_var()`'s, structurally

Not merely the same list: the **same construction**.

| | `lp_extra_var()` | `ce_draw_new_levels()` |
| --- | --- | --- |
| skip list | `R/predict.R:1535` `c("gr_cov", "gr_prec", "car", "spde")` | `R/conditional-effects.R:665`, same four, `draw = FALSE` -> coefficients zeroed |
| covariance | `:1536` `covstruct_registry[[cs]]$vcov(th[theta_idx], bk)` | `:715` the identical call |
| Student-t | `:1542` `S * student_var_factor(dist_nu)` | `:717` the identical line |

So the two bands cannot drift apart for any registry covstruct: they
call the same `vcov()` on the same `theta` slice with the same
Student-t adjustment. Verified by reading both, and by the numeric
agreement above.

### What happens for a `gp`/`hsgp` block

`ce_group_vars()` (`R/conditional-effects.R:439`) excludes
`smooth`/`gp`/`hsgp`, so such a block's variable never enters
`na_vars`; `ce_new_level_spec()` (`:648`) skips them again. The result
is that `re_formula = NULL` is a **no-op** for a gp model, in BOTH
bands. Measured on `y ~ x + gp(xs, k = 10)` (the block comes back as
`hsgp`), `ce_group_vars()` empty:

```
wald new se__ == wald pop se__ : TRUE  (max abs 0)
boot new width == boot pop width: TRUE  (max abs 0)
```

That is the right answer and it is self-consistent: a GP's "levels" are
basis functions, so there is no unseen level to be new about. No
asymmetry between the bands.

### What happens for an `rr` block

`rr` is the one block the spec treats specially: it stores standard-
normal FACTORS (rank per level) and expands them through the loadings
in `expand_b()`, so `ce_new_level_spec()` sets `dim = rank` and
`ce_draw_new_levels()` writes `rnorm(rank)` rather than a covariance
draw (`R/conditional-effects.R:711-713`). That is correct by
construction - `expand_b()` multiplies by `L`, giving covariance
`L L'`, which is exactly `covstruct_registry[["rr"]]$vcov()`.

Verified on `y ~ x + z1 + z2 + rr(1 + z1 + z2 | g, d = 2)` (dim 3,
rank 2, 40 groups), seed 31, `boot = 300`:

```
wald new se__  : 0.3577 0.3183 0.3137 0.3452
boot new / wald new width : 1.0607 1.0652 1.0352 1.0411
g column all NA : TRUE
```

The two routes agree to ~6%, which is bootstrap noise at B = 300.

**One residual, not a defect.** For a Student-t block both paths use
`S * nu/(nu-2)` and both then treat the result as GAUSSIAN - the wald
band by construction, and now the boot band because it draws
`L %*% rnorm()`. The bootstrap could have drawn an actual t and been
more faithful than the wald band it is being made consistent with;
consistency was the stated goal, so this is a deliberate choice worth
one line in the help page rather than a defect. `R/predict.R:1538-1541`
already carries the honest comment for the wald side.

## P2 - FIXED, and the draws defect WAS created by the re-key

I reproduced the assembly without sampling, by running the old loop
(`for (tn in names(per[[1]])) { M0 <- per[[1]][[tn]]; ...; out[[tn]] <- st }`)
and the new one (by position) against each core's real `ranef()` output
for `y ~ x + z + (1 | g) + (0 + z | g)`:

| core + assembly | result |
| --- | --- |
| **main** core, OLD by-name | names `1 \| g`, `0 + z \| g`; **length 2**, values 1.106170 and 0.336090 - CORRECT |
| **lane** core, OLD by-name | names `g`, `g`; **length 1**, value 1.106170 only |
| **lane** core, NEW by-position | **length 2**, 1.106170 and 0.336090, `term` attrs `1 \| g` and `0 + z \| g` - CORRECT |

True per-block values are 1.106170 and 0.336090. So the defect is
**created by the re-key**, not pre-existing: on main the block labels
are distinct and by-name assembly is correct.

**One correction to the lane's description.** The ledger says the old
code "would have returned the FIRST block's draws twice". It does not:
`out[[tn]] <- st` writes the same name twice, so the second block is
**dropped entirely** and the list comes back with length 1, not 2.
Worse than described, and the fix handles both.

## NEW RESIDUAL (R1, blocking-class): the `boot =` reuse path reopens P1

Found while probing the fix, not claimed by the lane either way.

`conditional_effects()` returns its bootstrap on the result
(`R/conditional-effects.R:1886`) so a second call can reuse it, and the
help page recommends that. `ce_boot_draws()` guards the reuse with
`ce_boot_key()` (`:612`), which hashes the grids' `nd` columns but
**not** the new-level spec.

The two grids are indistinguishable to that hash, because two
independent choices happen to coincide:

- `ce_ref_value()` (`:63`) holds an unvaried factor at `levels(col)[1L]`
- `ce_new_level_spec()` (`:652`) uses `bk[["levels"]][1L]` as its
  placeholder

Both are `"308"` on sleepstudy, so a population call and a
`re_formula = NULL` call produce byte-identical `nd` and therefore the
same key. Measured:

```
population band width          : 22.784 24.282 32.917 47.074 59.890
new-level  band width          : 92.509 121.285 168.378 213.087 269.325

reuse a POPULATION boot object in a re_formula = NULL call
  -> ACCEPTED, width 22.784 24.282 32.917 47.074 59.890   (the POPULATION band)

reuse a NEW-LEVEL boot object in a population call
  -> ACCEPTED, width 92.509 121.285 168.378 213.087 269.325  (4x too WIDE)
```

The first line is P1 exactly as it was before this round - a
`re_formula = NULL` call returning the population interval, silently -
reachable through a documented argument the help page encourages. The
second is the converse error. Neither warns.

This is created by the P1 fix: before it, the direct path and the reuse
path agreed (both wrong), so the guard had nothing to catch. Now the
direct path is right and the reuse path is not.

Fix: put the new-level spec in the key - `length(nspec)` alone would
separate the two cases, and the block ids would be tighter. Two lines
in `ce_boot_key()` plus its call site, which already receives `nspec`.

## P3, P4, P5, P6 - all confirmed

**P3.** `ce_effect2()` now takes levels from the distinct values:

```
int_conditions = list(z = c(lo = 0.001, hi = 0.002))
  -> 2 levels, "hi" "lo", both z values kept (0.001, 0.002), 200 rows
same unnamed -> 2 levels, "0.002" "0.001"
default mean +/- sd -> levels 1.03 -0.01 -1.05, values still EXACT
                       (max abs difference from mean +/- sd = 0)
```

**P4.** NEWS names all fourteen; I enumerated them through
`mean_is_mu()` here and the list matches character for character. The
"every other family is unchanged, bit for bit" clause now sits after
the explicit list, plus `trunc()` and `mixture()`, and `lognormal` gets
a bullet of its own. Reproduced on my own data (response range 0.203 to
10.003, so different from the lane's 0.077 to 7.809):

```
old curve (mu on its identity link) : -0.6695 -0.0579 0.5536 1.1652 1.7767
new curve (the ce estimate)         :  0.6091  1.1227 2.0695 3.8147 7.0316
new == predict(type = "response")   : max abs 0
old curve negative anywhere         : TRUE (min -0.6695)
```

Same phenomenon at a different seed: the old curve is negative for a
strictly positive response.

**P5.** The bullet carries
**"What breaks: `ranef(fit)[["Days | Subject"]]` now returns `NULL`"**
in bold at `NEWS.md:159`, with the `numeric(0)` hazard and the three
ways back (factor name, `[[1]]`, the `"term"` attribute), and it points
at `frmtmb.sample`. Substance satisfied. Small correction to the relay:
it does not *lead* with that - the bullet opens with the re-key
rationale and what is kept, and the breakage lands eight lines in. It
is bold and unambiguous where it sits, so I would not hold anything for
it.

**P6.** Documented in both places and accurately:
`man/mixture.Rd` says the SE is "the delta method through its OWN
predictor", exact for two components, and "CONSERVATIVE - measured
5.5% to 26.1% wider than the joint delta method on a three-component
fit, never narrower"; `man/predict.frmtmb_fit.Rd` carries the short
form and points at `?mixture`. My own measurement last round was 8.9%
to 17.3% on a different fit; the direction (always wider) is what the
documentation asserts and it holds.

**My earlier P8 (the half-migrated draws method)** is now documented
too, in `extensions/frmtmb.sample/NEWS.md`: the four inherited grid
changes are named, including that the `trials()` message "is new here",
and so is the list of what is NOT inherited, ending "the two surfaces
disagree in both directions until the draws method is migrated". That
was the ask.

## Runs

| run | result | lane's claim |
| --- | --- | --- |
| `test-brms-methods.R`, gated, warm | **47 tests, 950 assertions, 0 fail, 0 error, 0 skip**, 136 s | 47 / 950 - exact |
| `test-ce-bands.R` | **21 tests, 159 assertions, 0 fail**, 70 s | 21 / 159 - exact |
| `test-unpinned-seams.R` | **4 tests, 28 assertions, 0 fail** | 4 / 28 - exact |
| `frmtmb.sample`, one process | **136 tests, 899 assertions, 0 fail**, 1 error | 136 / 888 / 0 |
| `test-perf.R` alone | **2 tests, 3 assertions, 0 fail**, 19.6 s | passes alone - confirmed |

Two notes on those numbers. The `frmtmb.sample` error is the same
`local_mocked_bindings()` / pkgload artifact of MY harness that
`test-parallel-chains.R` produced last round, not a lane defect. My 899
assertions against the lane's 888 is the gate: I ran with
`FRMTMB_BRMS_FIT_TESTS=true` and the Stan cache, so the lane's 2 skips
executed here and contributed the extra 11. **0 failures either way -
P2 is fixed.**

`test-perf.R` passed for me even with 36 R processes from other lanes
on the box, so I could not reproduce the contention failure the lane
reports; it passes alone, which is the claim.

The P1 pin in `test-ce-bands.R:553` is a good one: it asserts the
estimate is unchanged, the grouping column is all `NA`, `boot new >
2 * boot pop`, `boot new / wald new` in [0.6, 1.5] and `boot pop /
wald pop` in [0.8, 1.2], then repeats the shape for the ordinal case.
Loose enough to survive a seed change, tight enough that a regression
to the population band (ratio 1.0) fails the first assertion.

## Updated verdict: GO-WITH-FIXES (one blocker left, two lines)

P1's direct path, P2, P3, P4, P5 and P6 are all genuinely fixed, and P1
is fixed for the right reason rather than to a number: the bootstrap
replicate spread now reproduces `sqrt(x' V x + z' S z)` to within 1.4%,
which also validates the block indexing and the intercept-slope
correlation, not just the magnitude. The two bands draw from the same
`covstruct_registry[[cs]]$vcov()` call with the same Student-t
adjustment, so they cannot drift apart per covstruct. `gp`/`hsgp` is a
consistent no-op in both bands; `rr` draws standard-normal factors and
lands within 6% of the wald band. The `frmtmb.sample` suite is green.
The methods tier is unmoved at 47/950.

**R1 is the one thing I would hold for**, and only because it is two
lines: `ce_boot_key()` cannot distinguish a population bootstrap from a
new-group one, so `boot =` reuse hands back the population band for a
`re_formula = NULL` call - P1's exact symptom, through the path the
help page recommends. Fix the key and this is a GO.

Residuals, none blocking:

- **R2** `R/predict.R:1538-1541` and `R/conditional-effects.R:717` - a
  Student-t block is drawn GAUSSIAN with the t variance in the boot
  band, matching the wald band's own approximation. Deliberate
  (consistency was the goal) but the bootstrap could have drawn a real
  t; one line in the band help would record the choice.
- **R3** `dev/ce-findings.md` and `extensions/frmtmb.sample/NEWS.md` -
  both say the old by-name assembly "would have returned the first
  block's draws twice". It drops the second block instead: `out[[tn]]
  <- st` writes one name twice, so the list comes back length 1. Worse
  than described; the fix handles it either way.
- **R4** the P5 bullet does not lead with the breakage (`NEWS.md:151`
  opens with the rationale, the bold "What breaks" is at `:159`).
  Cosmetic.

## R1 re-check, 2026-09-05 - FIXED

`ce_boot_key()` (`R/conditional-effects.R:612`) now takes `nspec` and
folds in `ce_new_level_key()` (`:641`), which reduces the spec to
`vars, parts, idx, dim, rr, draw`. The value is kept beside the key as
`bs$ce_new` (`:842`) and compared BEFORE the grid reason (`:790`),
which it has to be: the two grids are byte identical, so the grid
comparison could never be the one that fires.

### The grids really are byte identical, so the group key is doing the work

Built both grid sets through the internals and compared:

```
population boot grid Subject : 308
new-level  boot grid Subject : 308
identical(as.list(nd_pop), as.list(nd_new))         : TRUE
keys built WITHOUT the group component, identical   : TRUE   <- the old bug
keys built WITH    the group component, identical   : FALSE  <- the fix
```

So the separation comes from the group component alone, not from any
side effect on the grid component. That also means a genuinely
different grid still reaches the grid reason, which cases 5 and 6 below
confirm.

### Keyed on what identifies the draws, not on incidentals

`ce_new_level_key()` keeps exactly:

```
$ vars : chr "Subject"     which grouping columns are blanked
$ parts: chr "308"         which placeholder level the design maps
$ idx  : int [1:2] 1 2     where the draw is written (level one's 2 coefs)
$ dim  : int 2
$ rr   : logi FALSE        which draw rule
$ draw : logi TRUE         drawn, or zeroed for gr_cov/gr_prec/car/spde
```

The block object is NOT in it - I checked: no `theta_idx`,
`group_name` or `covstruct` appears anywhere in the reduced value. So
the key turns on which blocks are redrawn, where they are written, and
whether they are drawn or zeroed, which is precisely what makes one set
of draws valid for one curve and not the other. `idx = 1 2` also
confirms the level-major layout the `z' S z` agreement implied.

Nothing incidental is in either component. `prob` is deliberately
absent - the same draws answer any coverage - and I confirmed that
reuse across a coverage change is ACCEPTED and simply narrows the band
(case 7). The grid component keys on `as.list(g$nd)`, the column values
in row order; row order is IDENTIFYING for a bootstrap band rather than
incidental, because the draws are stored per grid row and a permutation
would misalign them. Row names and the data frame class are stripped by
`as.list()`, which the comment at `:605` states and which is the right
thing to strip.

### The four cases, my seed 2024 objects, resolution 5, `boot = 150`

```
population band width : 25.2820 27.5369 34.3637 42.9691 52.0401
new-level  band width : 88.5401 92.4620 125.1045 160.3074 197.8630
ce_new entries        : population 0, new-level 1
```

| # | reuse | result |
| --- | --- | --- |
| 1 | population draws -> `re_formula = NULL` call | **REFUSED**: "predictions for a different group: this call conditions on a NEW group (re_formula = NULL), and those draws do not carry that group's effects, so their percentiles are the population band" |
| 2 | new-level draws -> population call | **REFUSED**: "those draws carry a NEW group's effects ... and this call is the population curve, so their percentiles are too wide for it" |
| 3 | population draws -> population call | **ACCEPTED**, band reproduced **max abs 0** |
| 4 | new-level draws -> `re_formula = NULL` call | **ACCEPTED**, band reproduced **max abs 0** |
| 5 | population draws -> resolution 7 | **REFUSED**, and with the GRID reason ("the effects, the resolution, the conditions or the data are not the same") |
| 6 | new-level draws -> different effects/conditions | **REFUSED**, GRID reason |
| 7 | population draws -> same call, `prob = 0.8` | **ACCEPTED**, band correctly narrower at every point |

Two distinct messages in the two directions, each naming what is wrong
with those particular draws rather than a generic mismatch. The grid
reason still reaches the calls that deserve it.

One property worth recording because it is easy to get wrong and the
lane got it right: an OLD bootstrap object saved before this change has
no `ce_new`, and `boot$ce_new %||% list()` makes it compare equal to a
population call's empty spec. So a pre-existing saved object is still
accepted for the population curve it was actually computed for, and
refused for a new-group call. Backward compatible in the safe
direction. Confirmed by stripping `ce_new` off a real object:

```
old obj -> population call  : ACCEPTED, 25.2820 27.5369 34.3637 42.9691 52.0401
                             (its own band, reproduced exactly)
old obj -> re_formula = NULL: REFUSED, "predictions for a different group"
```

### R2, R3, R4

- **R2 fixed.** `R/conditional-effects.R:1234-1238` now records that a
  Student-t block "enters either band as a GAUSSIAN with the t
  variance, `nu / (nu - 2)` times the scale matrix ... the bootstrap
  draws the same way so that the two bands stay comparable. It is the
  right variance around a heavier-tailed truth, not the right
  quantile." That is the choice, stated.
- **R3 fixed.** `dev/ce-findings.md:316` and `:533` and
  `extensions/frmtmb.sample/NEWS.md:12` all now say the second block
  was DROPPED and the list came back with one entry. **One leftover:**
  the code comment at
  `extensions/frmtmb.sample/R/methods-draws.R:167` still reads
  "silently replacing the second block's draws with the first's". Stale
  by one round; cosmetic, and the code below it is correct.
- **R4 fixed.** `NEWS.md:151` now opens
  "BEHAVIOR CHANGE, and **what breaks is
  `ranef(fit)[["Days | Subject"]]`, which now returns `NULL`**", with
  the `numeric(0)` hazard and the three ways back, before the
  rationale. It leads with the breakage as asked.

### Runs

| run | result | claim |
| --- | --- | --- |
| `test-ce-bands.R` | **22 tests, 167 assertions, 0 fail**, 42 s | 22 / 167 - exact |
| `test-boot.R` | **7 tests, 43 assertions, 0 fail**, 11 s | - |
| `test-brms-methods.R`, gated, warm | **47 tests, 950 assertions, 0 fail, 0 error, 0 skip**, 116 s | 47 / 950 / 0 - exact |

The new pin is `test-ce-bands.R:659`, "a bootstrap cannot be reused
across the new-group boundary". It asserts the new band is more than
twice the population band before testing reuse, then both refusals by
message - so it fails loudly if the key stops discriminating AND if the
underlying P1 fix regresses.

---

# FINAL VERDICT: GO

Every blocker I raised is fixed, and each was verified by reproducing
the mechanism rather than by matching a number.

- **P1** the boot band under `re_formula = NULL` now carries the group
  variance: its replicate spread reproduces the `sqrt(x' V x + z' S z)`
  I derived independently, to within 1.4% at my own seed, which
  validates the magnitude, the block indexing and the intercept-slope
  correlation together. Ordinal wider at every row. `gp`/`hsgp` a
  consistent no-op in both bands; `rr` within 6%. The two bands share
  the same `covstruct_registry[[cs]]$vcov()` call and the same
  Student-t adjustment, so they cannot drift per covstruct.
- **R1**, which that fix opened, is closed at the level of the cause:
  the grids are byte identical and the key now carries a reduced
  new-level spec that says which blocks are redrawn and where. All four
  reuse directions behave, with two distinct messages, and a different
  grid still gets the grid reason.
- **P2** `frmtmb.sample` is green (136 tests, 899 assertions, 0 fail on
  my gated run), and the `ranef.frmtmb_draws()` defect the re-key
  created is fixed by position-indexing - I confirmed on main's core
  that by-name assembly was correct there, so it was created rather
  than pre-existing.
- **P3-P6** all reproduced independently: distinct `effect2__` levels
  with names kept, all fourteen families named in NEWS with the
  lognormal defect called out, the `ranef()` bullet leading with the
  breakage, and the K >= 3 theta SE documented as conservative in two
  help pages.

Residuals, none blocking and none needing a re-check:

- `extensions/frmtmb.sample/R/methods-draws.R:167` - stale comment,
  says "replacing the second block's draws with the first's" where the
  ledger and NEWS now correctly say the second block was dropped.
- The `frmtmb.sample` draws method is still half-migrated
  (`conditional-effects-draws.R`): it inherits four grid behaviors from
  core and none of the frame shape. This is now fully documented in
  that package's NEWS, in both directions, and belongs to its own lane.
  It should NOT hold this merge - the column divergence is cosmetic and
  self-consistent, and unlike the `ranef()` breakage it leaves no test
  red.

Merge it.
