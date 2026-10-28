# Lane mointerop: one- and two-category `mo()`, marginaleffects interop, rescor rows

Base: 435e7ce6 (frmtmb 0.76.0), branch `wt-mointerop`, library
`C:/Users/adf44/source/r/wt-mointerop-lib` in front of `rellib-r14`.
Started 2026-10-27. Nothing is committed. Logs are in
`dev/mointerop-out/` and `dev/mointerop-bin2-log/` (gitignored); every
number below is pasted from them, with the script that wrote it.

Other lanes of the round: sigmainn and priordraws (dfuns and reslot were
stopped before they began).

## 1. One- and two-category `mo()` as brms 2.23.0

### brms 2.23.0, run

`dev/mointerop-brms-mo12.R` (seed 5, 200 rows, `x` in 1..2, `z`, `g` in
8 groups; log `dev/mointerop-out/brms-mo12.txt`):

- Two values: `standata()` gives `Jmo = 1`, `Xmo_1` 0 (104 rows) and 1
  (96 rows), `con_simo_1 = 1`. `stancode()` declares
  `array[Imo] int<lower=1> Jmo` and `simplex[Jmo[1]] simo_1`.
  `default_prior()` has the rows `b mox` (flat, vectorized) and
  `simo mox1 dirichlet(1)`. A two-level ordered factor gives the same
  Stan code (`identical()` TRUE) and the row `simo moxo1`.
- The fit (2 chains of 2000, seed 3): `mox` 0.693 (0.160); `variables()`
  is `b_Intercept b_z bsp_mox sigma Intercept simo_mox1[1] lprior lp__`;
  `simo_mox1[1]` has mean 1, sd 0 and no R-hat; `epred(x = 2) -
  epred(x = 1)` equals `bsp_mox` in every draw (means 0.6931478 and
  0.6931478). `conditional_effects(effects = "x")` has the grid 1, 2;
  `effects = "z"` holds `x` at 1 (`round(median())`); `z:x` moderates
  with the one value 1. New data 3, 0 and 1.5 stop with "Invalid values
  in variable 'x'".
- One value (`x1 = 4`): `standata()` gives `Jmo = 0`, `Xmo_1` 0,
  `con_simo_1` of length 0; `default_prior()` lists `b mox1` and
  `simo mox11 dirichlet(1)`; `brm()` compiles and then stops: "Jmo[1]
  is 0, but must be greater than or equal to 1.000000", 0 draws.

So brms fits two values and cannot sample one. frmtmb now fits both by
maximum likelihood; `frm_sample()` samples two values and refuses one,
with brms's reason.

### The change

- `R/frame.R`: the refusal of `D < 2` is gone; the simplex component has
  `max(D - 1, 0)` coordinates.
- `R/objective.R`: `mo_simplex()` of no coordinates is the one weight 1;
  `mo_coords()` of fewer than two weights is `numeric(0)`
  (`mo_sphere_basis(1)` would call `contr.helmert(1)`, which stops).
  `mo_frame_terms()` builds the names with `recycle0 = TRUE`, so a term
  on one value has none (without it `paste0()` gave `"simo_mox11[]"`),
  and gains the field `label`.
- `R/fit.R`: `mo_search_terms()` skips a term with `D < 2`.
- `R/priors.R`: a Dirichlet on a simplex of fewer than 2 weights adds no
  entry (its density is the constant 1); `dirichlet_logdens()` returns
  0 for it. The length check keeps brms's message ("Expected input of
  length 1").
- `R/methods-fit.R`: `summary_mo_frame()` reports a one-weight simplex as
  brms does (1, error 0, interval [1, 1]) and no row for `D = 0`.
- frmtmb.sample: `mo_sample_check_const()` refuses `D = 0`;
  `draws_mo_cols()` builds the chart column names with `recycle0 = TRUE`
  (it gave `"zeta1_"` for `D = 1`, and the `simo_` column was never
  added); `mo_simplex_nlp()` skips `D < 2`; `mo_draws_to_chart()` guards
  an empty index (`max(integer(0))`); `draws_natural_cols()` takes the
  added names with `seq_along(names) > length(internal)` instead of an
  empty negative index, which selected nothing.

Consumers checked and unchanged because they already guard `D < 2` or
an empty coordinate vector: `se_mo_gaps()` (`D < 2` next),
`se_mo_null()` (gaps only), `lp_delta_mo()` (`!length(z)` next),
`mo_indicator_basis()` (`D < 1` NULL), `interop_extra_comps()` (only
components with coordinates), `mo_codes()` (`hi = lo + D`).

### Measured (`dev/mointerop-mo12.R`, seed 5; logs `mo12-base.txt`, `mo12-lane.txt`)

Base: every model is refused ("mo() needs at least 3 ordered
categories; mo(x) has 2", "... has 1").

Lane:

| model | frmtmb logLik | reference logLik |
|---|---|---|
| `y ~ mo(x) + z` | -309.634588678 | `lm(y ~ I(x - 1) + z)` -309.634588678 |
| `y ~ mo(xo) + z` (ordered, 2 levels) | -309.634588678 | |
| `y ~ mo(x) * z` | -309.274576804 | `lm(y ~ I(x - 1) * z)` -309.274576804 |
| `y ~ mo(x) + z + (mo(x) \| g)` | -295.174298261 | |
| `y ~ mo(x1) + z` (one value) | -318.635068016 | `lm(y ~ z)` -318.635068016 |

- Two values: `summary()` prints "Monotonic Simplex Parameters: mox1[1]
  1 0 1 1"; `fitted(newdata = x 1:2, z = 0)` SE 0.1117038 0.1163428,
  exact (lm with the ML sigma) 0.1117038 0.1163428; `x = 3`, `0` refused
  ("outside the fitted range 1..2"), `1.5` refused (whole numbers);
  `conditional_effects()` grid 1, 2, `z` held at `x = 1`, moderator 1,
  as brms; `emmeans(at = list(x = 1:2))` works;
  `dirichlet(c(1, 1))` refused "Expected input of length 1 for simo
  coef "mox1"", `dirichlet(2)` gives the ML logLik (the density is
  constant).
- One value: `mox1` estimate 0 with `NaN` SE and the warning "Standard
  errors are not available for 1 of 4 parameters. mox1: the likelihood
  is flat along it ..."; `variables()` `b_Intercept b_z bsp_mox1
  sigma`; `default_prior()` lists `simo mox11` as brms does;
  `fitted(newdata = x1 = 4)` finite SE 0.08417197; `x1 = 3` refused
  ("outside the fitted range 4..4").

`frm_sample()` (`dev/mointerop-sample.R`, lane, log `sample-lane.txt`;
2 chains of 2000, seed 3, against brms's own Stan program on the same
data, compiled through rstan):

| par | frm mean | frm sd | brms mean | brms sd | z | sd ratio |
|---|---|---|---|---|---|---|
| b_Intercept | 1.8280 | 0.1103 | 1.8292 | 0.1104 | -0.28 | 0.999 |
| b_z | 0.4861 | 0.0797 | 0.4876 | 0.0824 | -0.57 | 0.967 |
| bsp_mox | 0.6968 | 0.1611 | 0.6931 | 0.1598 | 0.59 | 1.008 |
| sigma | 1.1542 | 0.0595 | 1.1512 | 0.0580 | 1.58 | 1.026 |

brms max R-hat 1.000606. `variables(ds)`: `b_Intercept b_z bsp_mox
sigma simo_mox1[1] lp__`; `simo_mox1[1]` range 1 1;
`epred(x = 2) - epred(x = 1)` mean 0.6968379 = `bsp_mox` mean
0.6968379; `log_lik()` 2000 x 200, finite. One value: refused, "the
mo() term mox1 has one category ... brms 2.23.0 cannot sample it
either: its Stan program rejects Jmo = 0". (The first run of this
script, before the refusal was wired in, sampled the one-value model
and got 37 transitions past the maximum treedepth: the improper
posterior.)

### Fits with three or more categories are unchanged

`dev/mointerop-equal.R` (seed 5, 300 rows; 8 models: `mo(k)`,
`mo(k) + z`, `mo(k) * z`, `mo(k) + (mo(k) | g)`, ordered factor, the gap
`1, 2, 5`, `sigma ~ mo(k)`, poisson), every estimate, logLik, `fixef()`
and `fitted()` with SEs at 17 significant digits: `equal-base.txt` and
`equal-lane.txt`, 112 lines each, `diff` empty (IDENTICAL).

## 2. The ported block "update has reasonable outputs" (brmsfit-methods:913 to 949)

Run gated (`FRMTMB_BRMS_FIT_TESTS=true`) in record mode on both arms
(`dev/mointerop-out/rec-methods-base.tsv`, `rec-methods-lane.tsv`), and
once as the tier runs it (`t-suite-methods-lane.txt`: pass 161, fail 1,
the stale verdict of row 924 below):

| row | rellib-r14 | lane |
|---|---|---|
| setup 921 (`update(fit1, newdata = new_data, ...)`, `Exp = 4`) | fails: "mo() needs at least 3 ordered categories; mo(Exp) has 1" | runs |
| 923 `is(up, "brmsfit")` | STALE (reads `up` of the failed 921) | fails on the class, the recorded reason (divergence / class) |
| 924 `attr(up$data, "data_name") == "new_data"` | STALE | HOLDS |
| 927 `expect_error(update(fit1, data = new_data), "use argument 'newdata'")` | the refusal of mo() | fails because frmtmb refits, the recorded reason (divergence / policy) |
| setup 935 (`formula = ~ . + I(exp(Age))`, `newdata = new_data`) | fails as 921 | runs |
| 937 `is(up, "brmsfit")` | STALE | fails on the class, the recorded reason |

The port is literal: brms's own `new_data` with `Exp = 4`; no adapted
input. Row `brmsfit-methods:924` flips from divergence / policy to pass.
The tier file says so: "brms brmsfit-methods:924 now HOLDS but is
recorded as 'divergence' ... Change its row in
dev/brmsport-verdicts.tsv and regenerate." Rows 923, 927 and 937 keep
outcome and class, and now reach the reasons their verdicts give. Not
edited here: `dev/brmsport-verdicts.tsv`, the manual verdict file, the
ledger (consolidation regenerates them).

## 3. marginaleffects `vcov =`

marginaleffects 1.0.0 calls `get_vcov(model, vcov = vcov)`. Its default
method returns NULL for `FALSE`, aligns a matrix to `get_coef()` by
names (`align_vcov_to_coef()`), calls a function on the model, and
sends strings and formulas to `insight::get_varcov()` (sandwich). For
glmmTMB it takes only TRUE, FALSE and "HC0"; for brmsfit it warns that
`vcov` is not supported. frmtmb's method had the signature
`(model, ...)` and ignored all of it.

`get_vcov.frmtmb_fit(model, vcov = NULL, ...)` now: `FALSE` NULL;
`NULL`/`TRUE` the model-based matrix; otherwise
`interop_user_vcov()`: a matrix over `get_coef()` matched by names
(unnamed: in order, with marginaleffects' own warning text), a matrix
over the whole outer vector (`vcov(full = TRUE)`,
`vcov_cluster(full = TRUE)`) cut down through `interop_vcov(Vfull =)`, a
function, `"HC0"` = `vcov_cluster(each row, "CR0")`, a formula =
`vcov_cluster(formula, "CR1")`; every other string refused by name.

`dev/mointerop-me-probe.R` (seed 5; `me-probe-base.txt`,
`me-probe-lane.txt`):

| spelling | base | lane | lm / glm control |
|---|---|---|---|
| `vcov = 2 * get_vcov()` SE ratio | 1 1 1 | 1.414214 x 3 | lm 1.414214 x 3 |
| `vcov = FALSE` | no SE | no SE | lm no SE |
| `vcov = "HC3"` | default SEs, silent | refused by name | lm: ratio 0.968 1.004 1.040 to default |
| `vcov = "HC0"` (gaussian) | default SEs | 0.11440524 0.08412656 0.12207697 | lm 0.11440524 0.08412656 0.12207697 |
| `vcov = ~ g` (gaussian) | default SEs | 0.1955346 0.1974637 0.2066839 | lm 0.1960278 0.1979617 0.2072052 (vcovCL's HC1 adds (n-1)/(n-k): ratio 1.002522 = sqrt(199/198)) |
| `get_vcov(vcov = ~ g)` SEs, poisson | | 0.05133256 0.03147581 | glm 0.05133256 0.03147580 |
| `get_vcov(vcov = "HC0")` SEs, poisson | | 0.05934874 0.06342205 | glm 0.05934875 0.06342205 |
| function, unnamed matrix | default SEs | read (unnamed: warning) | |
| `diag(2)` (wrong size) | default SEs | refused, names the 3 coefficients | lm errors too |
| `"kenward-roger"` | marginaleffects' df error | refused by name | |
| `"unconditional"` | | marginaleffects refuses the class itself | |
| `vcov = "HC0"` on `(1 \| g)` | | refused: the random effect crosses the clusters | |

## 4. Gap-aware marginaleffects standard errors

On a fit with an empty interior `mo()` category (`se_mo_gaps()`),
`get_coef()` carries the cumulative weights `S_c = w_1 + ... + w_c` at
the interior categories that rows are in (`zeta<j>_cum<c>`), not the
chart coordinates. `S_obs` moves no prediction along the undetermined
direction (`w_c - w_(c+1)` changes `S` only at the gap `c`), so its
covariance `G Vp G'` from `cov_fixed_prop` (the covariance conditional
on the lost directions, the one `fitted()` reads) is the right one.
`set_coef()` stores the full cumulative vector with `NaN` at the gap
categories (`model[["mo_cum_override"]]`), which `mo_col_values()`
reads; marginaleffects' baseline of its finite difference also goes
through `set_coef()`, so the Jacobian row of a prediction at a gap code
is `NaN` and its standard error alone is `NaN`. One warning per fit
(`mo_gap_nan_warning()`, rewritten). The estimates come from the
unperturbed fit and stay finite. A one-value term's coefficient gets
variance 0 (`interop_mo_const()`): no prediction reads it, and a NaN
there made every SE NaN.

`dev/mointerop-gap.R` (`gap-base.txt`, `gap-lane.txt`):

Case 1, seed 7, `y ~ mo(x)`, `x` in 1, 2, 5 (200 rows):

| | x = 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| exact (lm on factor(x), ML sigma) | 0.1222210 | 0.1135415 | | | 0.1212624 |
| `fitted()` | 0.1222210 | 0.1135415 | NaN | NaN | 0.1212624 |
| `predictions()` base | NaN | NaN | NaN | NaN | NaN |
| `predictions()` lane | 0.1222210 | 0.1135415 | NaN | NaN | 0.1212624 |

`avg_comparisons(x = c(1, 5))`: exact 4.1261620, SE 0.1721701; lane
4.1261620, 0.1721701; base SE NaN. `c(1, 2)`: exact 0.8899317,
0.1668222; lane the same; base NaN. `c(1, 3)` (gap): lane SE NaN.
`avg_predictions()` over the data: exact 2.62025100, SE 0.06859636;
lane the same; base NaN. `get_coef()` names: base `(Intercept) mox
sigma_(Intercept) zeta1_1 zeta1_2 zeta1_3`, lane `(Intercept) mox
sigma_(Intercept) zeta1_cum1`.

Case 2, seed 11, `y ~ mo(x) + z`, `x` in 1, 2, 3, 5, 6 (300 rows):
exact SE at 1, 2, 3, 5, 6 0.1281900 0.1359410 0.1254415 0.1292299
0.1279079; lane `predictions()` 0.1281899 0.1359409 0.1254415 NaN
0.1292300 0.1279080 (NaN at 4); base NaN at all six. `avg_slopes(z)`
SE lane 0.05512498, exact 0.05512476, base NaN. `avg_comparisons(x =
c(2, 6))` lane 2.6206850 (0.1866555), exact 2.6206830 (0.1866555).
`vcov = 2 * get_vcov()` ratio lane 1.414214 at the five observed
categories, NaN at 4.

## 5. `predictions(re_formula =)` through marginaleffects

How marginaleffects decides: `sanity_dots()` keeps a hard-coded list
of known arguments per `class(model)[1]` (`brmsfit`: `draw_ids`,
`incl_autocor`, `nlpar`, `ndraws`, `re_formula`, `allow_new_levels`,
`sample_new_levels`, `dpar`, `resp`; `glmmTMB`: `re.form`,
`allow.new.levels`, `zitype`; ...). An argument outside the list gives a
WARNING ("These arguments are not known to be supported for models of
class `frmtmb_fit`: re_formula. All arguments are still passed to the
model-specific prediction function ..."), not a refusal, and the
argument is passed on. There is no registration seam (no option, no S3
generic in the list), so frmtmb cannot add itself; that needs a request
upstream (not posted). `sanitize_model_specific()` is an internal
generic for checks only.

So the defect was frmtmb's: `get_predict.frmtmb_fit()` dropped its dots.
It now passes `re_formula` and `allow_new_levels` to `frm_linpred()`.

brms 2.23.0 through marginaleffects could not be run here: marginaleffects'
brmsfit methods need the package `collapse`, which is not in the user
library ("Package `collapse` required"; `brms-mo12.txt`). brms's own
`fitted(re_formula = NA)` on the same model gives 1.646569 2.140818
2.635067 (posterior means), frmtmb's 1.670024 2.161305 2.652586 (ML).

`dev/mointerop-me-probe.R`, `(1 | g)` fit, 3 new rows:

| | base | lane |
|---|---|---|
| `predictions(re_formula = NA)` | 1.370196 2.567631 2.212198 plus 8 padding rows (the conditional ones) | 1.670024 2.161305 2.652586 (= `fitted(re_formula = NA)`) |
| `re_formula = ~ 0` | as above | as `NA` |
| `re_formula = ~ (1 \| g)` | conditional | conditional (= `fitted()`) |

A second defect found on the way: `predictions(fit, newdata = nd)` with
a factor column that has levels no row of `nd` uses returned extra rows
(11 for 3 rows of an 8-level factor). marginaleffects pads new data with
a row per level and drops them by their negative `rowid`;
`get_predict()` numbered the rows afresh. It now keeps `newdata$rowid`.

## 6. The rescor rows: log_lik:107, log_lik:109, posterior_predict:69, posterior_predict:72

Adapter (`dev/brmsport-bin2-adapters.R`, section "Multivariate gaussian
and student-t responses with residual correlation"): brms's prep (Mu
draws x rows x responses, Sigma draws x K x K, nu per draw) becomes a
frmtmb `mvbf(y1 ~ 0 + r, y2 ~ 0 + r, y3 ~ 0 + r) + set_rescor(TRUE)`
model (`frm(dry_run = "objective")`), one coefficient per row set to
Mu, each sigma the root of Sigma's diagonal, the residual correlation
through the exported `us_theta()`, student's shared nu through its link.
Values: `rescor_row_loglik()` (what frmtmb.sample's `log_lik()` runs per
draw) and one replicate of `predict(summary = FALSE, propagate_error =
FALSE)`, core's joint draw. brms's block hands one prep (with nu) to
both helpers; `log_lik_gaussian_mv()` never reads nu, so the adapter
ignores it for gaussian. brms's log_lik prep has Y with 4 rows and Mu
with 12; brms reads row i of both, so the model's rows are Y's. The spec
(`dev/brmsport-bin2-adapters-spec.R`) lists the two blocks in
`bin2_adapted` and `gaussian_mv`, `student_mv` in the substitutions;
`dev/brmsport-bin2-gen.R`: "wrote 6 files; 242 assertion ids; 4 without
a verdict; 153 adapted; 73 value checks"; only the log_lik (+59 lines)
and posterior_predict (+48) tier files changed.

Probe (`dev/mointerop-rescor-probe.R`, seed 1, `rescor-probe-lane.txt`):
frmtmb against brms's `log_lik_gaussian_mv()` / `log_lik_student_mv()`
at every row: max abs diff 3.55e-15 (gaussian), 3.55e-15 (student).
Joint draws (4000 per helper, one Sigma): frmtmb mean 1.0153 -1.0015
0.4888, brms 1.0284 -0.9813 0.4882, target 1 -1 0.5; covariances agree
with the target to sampling error (gaussian diag 0.807 0.963 1.021
against 0.807 0.918 1.070; student 1.188 1.467 1.545 against 1.210
1.376 1.605).

Recorded on the lane build (`FRMTMB_PORT_LIB` the lane library,
`R_LIBS_SITE` rellib-r14; `dev/mointerop-bin2-log/run-*.txt`,
`rec-*.tsv`):

    RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=101 fail=0 err=0 skip=0 blocks=16
    RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=63 fail=0 err=0 skip=0 blocks=15
    assert frmtmb log_lik:107 unclassified TRUE
    value  frmtmb log_lik:107 HOLDS: 10 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 1.51e-16
    assert frmtmb log_lik:109 unclassified TRUE
    value  frmtmb log_lik:109 HOLDS: 10 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 2.11e-16
    assert frmtmb posterior_predict:69 unclassified TRUE
    assert frmtmb posterior_predict:72 unclassified TRUE

The four rows move from "cannot transfer / not-ported (cost)" to pass,
literal (not adapted-input: brms's own prep runs). The ledger, the
verdicts and `dev/brmsport-bin2-survey.tsv` are not edited
(consolidation regenerates; the survey's four "not ported (cost)"
notes are then overridden by the ported blocks).

### The route of posterior_predict:69 and 72

Per the bin-2 convention (other adapters read core's `sim_draw()`),
posterior_predict:69 and 72 pass through core's `predict()` joint
draw: brms's assertion is the shape `c(ns, nvars)`, and the adapter
answers brms's helper with one replicate of
`predict(summary = FALSE, propagate_error = FALSE)` per draw. The four
rows also hold on rellib-r14 (review, tier runs): they move because
of the adapter, not a package change.

## 7. Optional: `variables()` of an ML fit has no `simo_` entries

Not done. `variables.frmtmb_fit()` lists `hyp_env_vals()`, the
hypothesis environment, which frmtmb.sample also reads for its derived
draws columns; adding `simo_` there means adding the simplex to the
hypothesis layout and its delta-method rows. Row
`brmsfit-methods:978` stays a divergence anyway (it also asks for
`r_`, `s_`, `prior_` and `lp__`). Left open.

## Tests

New or changed:

- `tests/testthat/test-mo-coding.R`: "mo() refuses what brms refuses
  and takes one or two values" (rewritten from the old refusal),
  "mo() on two values is its coefficient times the code", "mo() on one
  value is 0 and leaves its coefficient open", "brms codes one and two
  values as frmtmb does", and the gap part of "marginaleffects'
  standard errors carry the mo() simplex" (rewritten: finite and exact
  at observed categories, NaN at the gap; was "NaN everywhere").
  Lane pass 155, fail 0; rellib-r14 pass 120, fail 4, error 4, warn 1
  (`t-test-mo-coding-{lane,base}.txt`): the four errors are the old
  refusal; the four failures are the gap SEs (NaN at x = 1, 2, 5; no
  `zeta1_cum` names; NaN contrast 1 to 5).
- `tests/testthat/test-marginaleffects-args.R` (new): lane pass 29;
  rellib-r14 pass 10, fail 16, error 1
  (`t-test-marginaleffects-args-{lane,base}.txt`): the doubled vcov's
  ratio 1, `vcov = FALSE`, unnamed and wrong-size matrices, HC0, the
  cluster sandwich, the refusals of HC3 and others, `re_formula`
  (conditional predictions returned), and 11 rows for 3.
- frmtmb.sample `tests/testthat/test-mo-coding.R`: "a two-value mo()
  is drawn with brms's constant simo_ weight", "a one-value mo() is
  refused by the sampler, as by brms's", and the gated "frm_sample() of
  a two-value mo() agrees with brms". Gated: lane pass 50; rellib-r14
  pass 29, error 3 (`ts-mocoding-gated-{lane,base}.txt`).

## Suites

Every test file of core and of the eight extensions, one file per R
process, against the lane library (frmtmb and frmtmb.sample from this
worktree, the other extensions from rellib-r14 behind it, so each loads
this core): `dev/mointerop-par.sh`, runner `dev/mointerop-run-test.R`,
6 at a time. Every `lib:` line names wt-mointerop-lib for frmtmb (456
of 456 plain, 113 of 113 gated). The plain tier ran on the build before
the last fix (a duplicated message template, found by
`test-message-uniqueness.R`, fail 2); the 45 core files that reach
interop, marginaleffects or `mo(`, and that file, were rerun on the
final build (`r2plain`). The gated tier (FRMTMB_BRMS_FIT_TESTS=true)
ran the 113 files that mention a brms fit gate, the brms suite, `mo(`,
marginaleffects or `rescor`; its one failure is row
`brmsfit-methods:924`, which now holds against its recorded verdict
(section 2). `test-drmtmb-agreement.R` skipped 13 under its own gate
and was rerun with FRMTMB_DRMTMB_FIT_TESTS=true: pass 131, fail 0. Then
`dev/release/run-sessions.sh mointerop-plain 6` with FRMTMB_RELLIB set
to the lane library (log `dev/release/mointerop-plain.log`), on the
final build. Generated by `dev/mointerop-summary.R` and grep
(`dev/mointerop-out/sum-all.txt`):

```
source: dev/mointerop-out/plain.log 
plain ran 456 of 456 
             pkg files  pass fail error skip warn
          frmtmb   282 20127    2     0  204    0
 frmtmb.coupling    11   551    0     0    5    0
      frmtmb.eam    32  1838    0     0    5    0
   frmtmb.latent    12   681    0     0    2    0
    frmtmb.learn    15   430    0     0   13    0
      frmtmb.ode    11   549    0     0    1    0
   frmtmb.sample    70  3276    0     0   15    0
  frmtmb.spatial     8   478    0     0    6    0
 frmtmb.survival    15   552    0     0    0    0
files with fail, error or warn: 1 
    pkg                      file pass fail error skip warn
 frmtmb test-message-uniqueness.R    5    2     0    0    0
files without a RESULT line: 0 
source: dev/mointerop-out/r2plain.log 
r2plain ran 45 of 45 
    pkg files pass fail error skip warn
 frmtmb    45 3781    0     0  102    0
files with fail, error or warn: 0 
files without a RESULT line: 0 
source: dev/mointerop-out/gated.log 
gated ran 113 of 113 
             pkg files pass fail error skip warn
          frmtmb    92 8709    1     0   13    0
      frmtmb.eam     1   13    0     0    0    0
   frmtmb.sample    17 1494    0     0    0    0
  frmtmb.spatial     2  193    0     0    0    0
 frmtmb.survival     1   15    0     0    0    0
files with fail, error or warn: 1 
    pkg                      file pass fail error skip warn
 frmtmb test-brms-suite-methods.R  161    1     0    0    0
files without a RESULT line: 0 
SESSION frmtmb forward files=272/282 pass=20128 fail=0 err=0 skip=204 warn=0
SESSION frmtmb reverse files=272/282 pass=20278 fail=0 err=0 skip=204 warn=0
SESSION frmtmb preload files=272/282 pass=20371 fail=0 err=0 skip=206 warn=0
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
SESSION frmtmb.sample forward files=65/70 pass=3276 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.sample reverse files=65/70 pass=3276 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.sample preload files=65/70 pass=3276 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
MOINTEROP-PLAIN ran 27 of 27
LEAK lines: 46, untagged: 0
```

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`;
Rblas 74 symbols and Rlapack 616 routed to OpenBLAS), lane build:
`test-mo-coding.R` pass 155, `test-marginaleffects-args.R` pass 29,
`test-se-check.R` pass 134, `test-interop.R` pass 45, frmtmb.sample
`test-mo-coding.R` pass 9, skip 2 (ungated); `dev/mointerop-equal.R`
base and lane under OpenBLAS: identical (`ob-equal-*.txt`). The copy
`dev/ciharden-out/` is deleted.

Test files: core 281 at 435e7ce6, 282 now (`test-marginaleffects-args.R`
added); frmtmb.sample 70 and 70. No test moved or dropped. In core's
`test-mo-coding.R` the refusal of two values is rewritten (the refusal
is gone), and the gap assertion (NaN at every category) is rewritten
to the new behavior.

## Overlap with the other lanes

Files this lane touches that the brief names as likely overlaps:
`R/frame.R` (the `mo()` loop only), `R/predict.R` (`mo_col_values()`
only), `R/interop.R`, `R/priors.R` (the simplex prior branch and
`dirichlet_logdens()`), `NEWS.md`, `man/` (`frm.Rd`,
`frmtmb-sampling-api.Rd`, new `frmtmb-marginaleffects.Rd`),
frmtmb.sample `R/sample.R` (mo helpers and one call in `frm_sample()`)
and `R/draws-brms.R` (one line), `dev/brmsport-bin2-adapters*.R` and the
two generated bin-2 tier files. NAMESPACE is unchanged.

## R CMD check --as-cran (once, at the end)

`dev/mointerop-check.sh <pkg>`, built and checked in
`dev/mointerop-check/<pkg>/`, `R_LIBS` = lane library, rellib-r14, user
library, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on
PATH. The frmtmb tarball was built at 21:39, after the last edit of
`R/` (20:19) and `NEWS.md` (21:01).

- frmtmb: `Status: 1 NOTE`, the V8 math-rendering NOTE on the HTML
  manual. Tests in one session: `[ FAIL 0 | WARN 0 | SKIP 232 | PASS
  19881 ]` (`dev/mointerop-check/frmtmb/frmtmb.Rcheck/tests/testthat.Rout`).
- frmtmb.sample: `Status: 1 NOTE`, the same V8 NOTE. Tests: `[ FAIL 0 |
  WARN 0 | SKIP 23 | PASS 3262 ]`.

## Versions and floors

- frmtmb: a minor bump (BREAKING: `mo()` on one or two values fits;
  `get_coef()` names on a gap fit; `vcov =` strings refused by name; a
  new help topic). No new export; NAMESPACE unchanged.
  `mo_frame_terms()` (exported for the extensions) gains the field
  `label`, and `mo_simplex()` / `mo_coords()` take a simplex of one
  weight.
- frmtmb.sample: a minor bump; its floor is the new frmtmb (it calls
  `mo_frame_terms()`'s `label`, and its two-value draws need
  `mo_simplex(numeric(0))`). Its NEWS says "Needs the frmtmb that fits
  `mo()` on one or two values".
- No other extension changed; they load the new core unchanged (all
  suites above).

## Ported rows this lane flips (not edited in the ledger)

- `brmsfit-methods:924`: divergence / policy to pass (literal).
- `brmsfit-methods:923`, `:927`, `:937`: outcome and class unchanged;
  they now reach the reasons their verdicts state (class name; the
  `data =` spelling) instead of "STALE" or the mo() refusal.
- Bin 2: `log_lik:107`, `log_lik:109`, `posterior_predict:69`,
  `posterior_predict:72`: cannot transfer / not ported (cost) to pass,
  literal; the two log_lik rows' value checks hold (max abs diff 8.88e-16).

Expected at consolidation: bin 1 pass 391 to 392 of 494; bins 1 and 2
619 to 624 of 823 (the four bin-2 rows leave the "not ported (cost)"
class, which stays in the rate's denominator).

## Defects found, not fixed

1. marginaleffects 1.0.0 warns on `re_formula` for any class not in its
   hard-coded list (section 5). Upstream request needed (not posted).
2. Optimizer precision, pre-existing: `frm(y ~ mo(x) + z)` on two values
   stops at "relative convergence (4)" about 2.5e-6 relative from
   `lm()`'s exact coefficients, where `frm(y ~ I(x - 1) + z)` reaches
   them to the last digit; three categories show the same gap (fitted
   0.971483844 against lm 0.971484642, `dev/mointerop-twofit.R`).
   `newton_polish()` does not move these fits. Statistically nil (1e-5
   standard errors); the tests measure against the estimates' own
   standard errors.
3. `dev/brmsport-bin2-survey.tsv` still carries the "not ported (cost)"
   notes of the four rescor rows; the ledger's run overrides them.

## What I did not do

- Item 7 (`simo_` in `variables()` of an ML fit): section 7.
- Refusal of saved 0.75.0 draws and emmeans' `round(median())`: skipped
  by the brief.
- No upstream issue was filed.

## What needs the user

- One value: frmtmb fits it by maximum likelihood (coefficient `NaN`
  SE, warning) where brms only builds it and its sampler stops. The
  brief asked for brms's behavior; brms has no fitted behavior here, so
  I matched its model (`Jmo = 0`, the term 0) for `frm()` and its
  refusal for `frm_sample()`. If the user prefers a refusal at `frm()`
  too, row `brmsfit-methods:924` goes back to a divergence.
- A marginaleffects feature request (a way for a package to declare its
  prediction arguments) would remove the `re_formula` warning.

## Files touched

Core: `R/frame.R`, `R/objective.R`, `R/fit.R`, `R/priors.R`,
`R/methods-fit.R`, `R/predict.R`, `R/interop.R`, `R/compat.R`,
`R/sampling-api.R`, `man/frm.Rd`, `man/frmtmb-sampling-api.Rd`,
`man/frmtmb-marginaleffects.Rd` (new), `_pkgdown.yml`, `NEWS.md`,
`vignettes/inputs.Rmd`, `tests/testthat/test-mo-coding.R`,
`tests/testthat/test-marginaleffects-args.R` (new).
frmtmb.sample: `R/sample.R`, `R/draws-brms.R`, `NEWS.md`,
`tests/testthat/test-mo-coding.R`.
Port: `dev/brmsport-bin2-adapters.R`, `dev/brmsport-bin2-adapters-spec.R`,
`dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R`,
`dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-posterior_predict.R`
(regenerated by `dev/brmsport-bin2-gen.R`).
Lane scripts: `dev/mointerop-*.R`, `dev/mointerop-*.sh`, this file.
`dev/mointerop-rtmb0.R` shows RTMB takes a zero-length parameter
component (the two-category simplex). `dev/mointerop-me-src*.R` print
the marginaleffects internals quoted in sections 3 and 5.

## Fix round 1 (review `dev/reviews/2026-10-27-mointerop.md`, verdict FIX)

Probe: `dev/mointerop-fix1.R` (the review's constructions: seed 7 gap
data with `g` in 20 groups of 10 rows; seed 5 one-value data; seed 31
of `dev/review-mointerop-vcov.R`; seed 41 RE data; seed 23 and seed 51
of `-misc.R` and `-gapre.R`), logs `dev/mointerop-out/fix1-lane.txt`,
`fix1-base.txt` (the base arm stops at the one-value fit, which
rellib-r14 refuses).

### 1. Sandwiches on the determined directions (major)

`cluster_sandwich_full()` (R/sandwich.R) forms `B S'S B` with the bread
`cov_fixed_prop`, the covariance conditional on the lost directions,
where the fit lost a parameter, and returns both the sandwich `Vp` and
`V` with `NaN` in the lost rows only. The cluster scores are 0 along an
exactly flat direction, so for any quantity that does not move along
one, this is the sandwich of the determined parameters.
`vcov_cluster()` reads it (before: one `NaN` row of the bread spread to
every entry, so `vcov_cluster()` was all `NaN` on any fit with a lost
parameter; base `fix1-base.txt`: every entry of the seed-7 gap fit's
`vcov_cluster(~ g, full = TRUE)` NaN). `interop_user_vcov()` passes
`Vp` to `interop_vcov(Vfull_p =)`, which carries a gap simplex through
the same `G` as the default.

| seed 7 gap fit, `predictions()` SE | x = 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| `vcov = "HC0"`, lane | 0.1224166 | 0.1103619 | NaN | NaN | 0.1248376 |
| `lm(factor(x))`, `"HC0"` | 0.1224166 | 0.1103619 | | | 0.1248376 |
| `vcov = "HC1"`, lane | 0.1233451 | 0.1111990 | NaN | NaN | 0.1257846 |
| `lm(factor(x))`, `"HC1"` | 0.1233451 | 0.1111990 | | | 0.1257846 |
| `vcov = ~ g`, lane | 0.1082512 | 0.1336043 | NaN | NaN | 0.1179098 |
| `lm(factor(x))`, `~ g` (vcovCL HC1) | 0.1087993 | 0.1342808 | | | 0.1185068 |
| base, `"HC0"` | NaN | NaN | NaN | NaN | NaN |

The `~ g` ratio 0.99496 is `sqrt((n - k) / (n - 1))` = sqrt(197/199),
lm's HC1 adjustment, documented. One-value fit (seed 5): HC0 0.07415456
0.11901178 = `lm(y ~ z)`; HC1 0.07452814 0.11961134 = lm; `~ g`
0.07551369 0.17875544 against lm's vcovCL 0.07570414 0.17920628 (the
same HC1 factor). Boundary RE fit (seed 31, sd 5.7e-5): `vcov = ~ g`
0.1708337 0.0826883 0.1202062, finite (was NaN NaN NaN).

A user's whole-outer-vector matrix with `NaN` in a lost row
(`vcov(fit, full = TRUE)` of the gap fit) cannot be projected; it now
warns, naming the row: "vcov =: the covariance has no finite value for
zeta1_cum1, which the fit does not determine, so every standard error
marginaleffects forms from it is NaN. vcov = "HC0", "HC1", ~ g and the
default (TRUE) are formed on the directions the data determine".

### 2. One-value fits

- `loglik_df()` subtracts `mo_const_df()`: df 3 and AIC 593.0153, equal
  to `lm(y ~ z)` (was df 4, AIC 2 more); `anova(y ~ z, y ~ mo(x1) + z)`
  reports Df 3 and 3, Chi Df 0.
- `predict()`: `fit_draw_space()` zeroes the one-value coefficient's
  row of the conditional covariance (`mo_const_outer_pos()`), and
  `predict_par_drawer()` draws when it is the only lost parameter
  (`mo_const_only`). Est.Error 1.044101 1.043271 at z = 0, 1 (4000
  draws, seed 1), against `y ~ z`'s 1.053868 1.049813; no fallback
  warning (was "the covariance of the estimates could not be recovered
  from the Hessian, so the draws are taken at the estimates alone").
- `(mo(x1) | g)` on one value is refused in `mo_re_link()`: "mox1 in
  mo(x1) | g: x1 has one value, so the mo() term is 0 on every row and
  its coefficient cannot vary over the levels of a group. Remove it
  from the group-level term". frmtmb fits other undetermined sds with a
  NaN SE; here the column is identically 0, so the parameter is not a
  parameter of the model at all, and a refusal is what brms's sampler
  gives too.
- Not done: `predict()` of a gap fit (pre-existing, both arms): a
  prediction at a gap code does move along the lost direction, so the
  plug-in fallback stays there. Backlog.

### 3. HC1, stata, refusals, REML

`"HC1"` and `"stata"` (marginaleffects 1.0.0 maps "stata" to HC1 and
"robust" to HC3 for lm; `dev/mointerop-stata.R`, seed 31: lm HC1 =
stata = 0.06112569 0.0653368) are HC0 times `n / (n - k)`,
`interop_hc1_k()`: the population-level coefficients plus the
determined steps of their `mo()` simplexes, less one-value terms
(`y ~ mo(x)` on the gap data counts 3, as `lm(y ~ factor(x))`). frm
`y ~ x` HC1 and stata 0.2267092 0.1153291 0.1998435 = lm. Refusals
now give each string's reason (`interop_vcov_refusal()`): HC2 to HC5
and "robust" need leverages, HAC / NeweyWest / kernHAC a time order,
satterthwaite / kenward-roger denominator df, OPG the inverse outer
product. A full matrix on a REML fit: "this fit is REML, which
integrates the fixed effects out, so that vector has no rows for them.
Give a matrix over get_coef() ...".

### 4. The help page

`?frmtmb-marginaleffects` says the bread is the observed information,
equal to sandwich's for gaussian and canonical-link fits, and gives the
Gamma(log) case: slope SE ratio 1.127856 (`fix1-lane.txt`, the review's
construction).

### 5. `get_predict()` arguments

`me_refuse_dots()` refuses `re.form`, `allow.new.levels`,
`random.only`, `include_random`, `sample_new_levels`, `ndraws`,
`draw_ids`, `incl_autocor` and `zitype` by name with a reason;
`dpar` (and brms's `nlpar`, which frm_linpred() reaches through `dpar`)
and `resp` are passed to `frm_linpred()`. marginaleffects wraps the
error but carries frmtmb's reason ("marginaleffects passed re.form =,
which frmtmb's prediction does not take: lme4's spelling; frmtmb takes
brms's re_formula ="). `predictions(fs, dpar = "sigma")` on
`bf(y ~ x, sigma ~ x)` gives 1.309742 1.278432 1.247870 =
`fitted(dpar = "sigma")`.

### 6. The degenerate `(mo(x) | g)` gap fit (seed 23)

When `G Vp G'` has `NaN` outside the rows the shown covariance lost,
`interop_vcov()` falls back to the shown covariance with the simplex
rows `NaN` (by index, not by a product, which would spread the `NaN`)
and warns once per fit (`mo_gap_lost_warning()`): the coefficient
block is finite (Intercept 1.215013e-02, z 2.783409e-03, mox
4.240265e-03, sigma 1.724623e-03 on the diagonal), `zeta1_cum1` NaN,
as rellib-r14 gave with its warning.

### 7 and 10

The gap warning counts rows with a positive `rowid`: "1 of 4
predictions" on 4 rows of a 20-level factor (was "1 of 24"). The gap
part of core `test-mo-coding.R` takes its yardstick from the gap fit's
own simplex-free standard error (`frm_lp_basis()` without the `zeta1`
columns).

### Pre-existing, for the backlog (not fixed)

- `fitted()` gives `Est.Error` 0, not `NaN`, at the empty category on
  the degenerate fit of item 6 (both arms; review).
- `drop1()` omits `mo()` terms (both arms; also a 4-category fit).
- `predict()` of a gap fit draws at the estimates alone (both arms).

### Tests seen to fail on rellib-r14

- `test-marginaleffects-args.R` (now 57 expectations): lane pass 57;
  rellib-r14 pass 15, fail 37, error 2, warn 1
  (`f1-test-marginaleffects-args-{lane,base}.txt`). New blocks: the
  gap fit's HC0/HC1/~ g and `vcov_cluster()`, the one-value fit's
  HC0/HC1, the boundary RE fit's `~ g` and `vcov_cluster()`, HC1 and
  stata against lm with the refusals' reasons and the REML refusal,
  `get_predict()`'s refusals and `dpar`, the padding count.
- `test-mo-coding.R`: lane pass 164; rellib-r14 pass 122, fail 5, error
  6, warn 2 (`f1-test-mo-coding-{lane,base}.txt`). New blocks: "a
  one-value mo() term is no degree of freedom and no variance", "a gap
  fit that lost another direction warns in get_vcov()".

### Suites after fix round 1

On the build installed at 22:38:48, after the last edit of `R/` (fit.R
22:38:17). One file per process, every test file of core and the eight
extensions (`f1plain`, 456 of 456 with `lib:` naming the lane library);
the 113 gated files (`f1gated`; its one failure is row
`brmsfit-methods:924`, still recorded as a divergence;
`test-drmtmb-agreement.R` skips under its own gate, which passed 131 in
the first round); one session per suite for core and frmtmb.sample,
three orders (`dev/release/mointerop-f1.log`); OpenBLAS 0.3.26 with
LAPACK (74 Rblas and 616 Rlapack symbols routed) for the changed
numerical files and `dev/mointerop-equal.R` (base and lane identical;
the copy is deleted). Generated (`dev/mointerop-out/sum-f1.txt`):

```
source: dev/mointerop-out/f1plain.log 
f1plain ran 456 of 456 
             pkg files  pass fail error skip warn
          frmtmb   282 20163    0     0  204    0
 frmtmb.coupling    11   551    0     0    5    0
      frmtmb.eam    32  1838    0     0    5    0
   frmtmb.latent    12   681    0     0    2    0
    frmtmb.learn    15   430    0     0   13    0
      frmtmb.ode    11   549    0     0    1    0
   frmtmb.sample    70  3276    0     0   15    0
  frmtmb.spatial     8   478    0     0    6    0
 frmtmb.survival    15   552    0     0    0    0
files with fail, error or warn: 0 
files without a RESULT line: 0 
source: dev/mointerop-out/f1gated.log 
f1gated ran 113 of 113 
             pkg files pass fail error skip warn
          frmtmb    92 8746    1     0   13    0
      frmtmb.eam     1   13    0     0    0    0
   frmtmb.sample    17 1494    0     0    0    0
  frmtmb.spatial     2  193    0     0    0    0
 frmtmb.survival     1   15    0     0    0    0
files with fail, error or warn: 1 
    pkg                      file pass fail error skip warn
 frmtmb test-brms-suite-methods.R  161    1     0    0    0
files without a RESULT line: 0 
SESSION frmtmb forward files=272/282 pass=20165 fail=0 err=0 skip=204 warn=0
SESSION frmtmb reverse files=272/282 pass=20315 fail=0 err=0 skip=204 warn=0
SESSION frmtmb preload files=272/282 pass=20408 fail=0 err=0 skip=206 warn=0
SESSION frmtmb.sample forward files=65/70 pass=3276 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.sample reverse files=65/70 pass=3276 fail=0 err=0 skip=15 warn=0
SESSION frmtmb.sample preload files=65/70 pass=3276 fail=0 err=0 skip=15 warn=0
MOINTEROP-F1 ran 6 of 6
LEAK lines: 34, untagged: 0
RESULT test-confint-anova.R pass=35 fail=0 error=0 skip=0 warn=0
RESULT test-interop.R pass=45 fail=0 error=0 skip=0 warn=0
RESULT test-marginaleffects-args.R pass=57 fail=0 error=0 skip=0 warn=0
RESULT test-mo-coding.R pass=164 fail=0 error=0 skip=0 warn=0
RESULT test-sandwich.R pass=58 fail=0 error=0 skip=0 warn=0
RESULT test-se-check.R pass=134 fail=0 error=0 skip=0 warn=0
```

R CMD check --as-cran of frmtmb (tarball built 23:35, after the last
edit): `Status: 1 NOTE` (the V8 math-rendering NOTE); tests `[ FAIL 0 |
WARN 0 | SKIP 232 | PASS 19918 ]`. frmtmb.sample did not change in this
round; its check of the first round (`Status: 1 NOTE`, V8) stands.

### Files touched in fix round 1

`R/sandwich.R` (`cluster_sandwich_full()`, `vcov_cluster()`),
`R/interop.R` (`interop_vcov(Vfull_p =)` and its fallback,
`mo_gap_lost_warning()`, `interop_user_vcov()`, `interop_vcov_checked()`,
`interop_hc1_k()`, `interop_vcov_refusal()`, `mo_const_outer_pos()`,
`get_predict()`, `me_refuse_dots()`, the help topic),
`R/methods-fit.R` (`loglik_df()`, `mo_const_df()`), `R/brms-shapes.R`
(`fit_draw_space()`), `R/predict-brms.R` (`predict_par_drawer()`),
`R/frame.R` (`mo_re_link()`), `R/fit.R` (mo() text), `man/frm.Rd`,
`man/frmtmb-marginaleffects.Rd`, `NEWS.md`,
`tests/testthat/test-mo-coding.R`,
`tests/testthat/test-marginaleffects-args.R`; dev scripts
`dev/mointerop-fix1.R`, `dev/mointerop-stata.R`,
`dev/mointerop-seed23w.R`. Removed: `dev/mointerop-rescor-pp.R`.

Bump and floor unchanged: frmtmb minor; frmtmb.sample minor with the
new frmtmb as its floor.
