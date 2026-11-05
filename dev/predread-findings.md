# Lane predread: the readers of the fitted predictor

Base d5d0056d (frmtmb 0.78.0, frmtmb.sample 0.25.0). Every number below names
its script, with its seed in the script. The logs are in `dev/predread-log/`
(gitignored); the numbers a claim rests on are pasted here.

## 1. What changed

### D4: `cs()` offsets, and every layered predictor

brms's `posterior_linpred()` with no `dpar`, run on `sratio()` with `cs(z)`
(`dev/predread-cs-layout-brms.R`, seed 13, log `cs-layout-brms.txt`):

| call | brms 2.23.0 |
|---|---|
| `posterior_linpred()` | draws x rows x thresholds, dimnames `NULL` |
| `posterior_linpred(dpar = "mu")` | draws x rows: the shared part alone |
| `posterior_linpred(transform = TRUE)` | draws x rows: the link inverse of the shared part |
| `fitted(scale = "linear")` | rows x 4 x thresholds, `eta1`, `eta2`, ... |
| `fitted(dpar = "mu", scale = "linear")` | rows x 4: the shared part |
| `emmeans(~ x)` | the mean over thresholds; `~ x / rep.meas` gives each, `"1"`, `"2"`, ... |
| mixture, `cs()` in `mu1` only | 3 + 1 layers; fitted names `eta1` to `eta4`; `dpar = "mu1"` is shared |
| mv, `y` and a `cs()` response `o` | linpred layers `y`, `o1`, `o2`; fitted `y`, `eta1`, `eta2` |

So "no `dpar` named" is what stacks, and this also holds for families with
several location predictors. `frm_linpred(type = "link")` returned the first
predictor of a `categorical()` or a `mixture()` without a word: the same defect
in another cell. One builder now serves all of them.

- `lp_layers()` decides the layers, `linpred_layers()` is `frm_linpred()`'s and
  `lp_basis_layers()` is `frm_lp_basis()`'s. `lp_basis_layer()` is one layer's
  design, and under `cs()` it holds the shared design in every threshold's rows
  plus each `cs()` column in its own threshold. `frm_linpred(se.fit = TRUE)` is
  read off that design, as before, so se.fit and the basis are one route.
- `ord_linear_per_threshold()` is removed. `fitted(scale = "linear")` of a `cs()`
  fit now reads `frm_linpred()` for the estimate and the standard error. On p013's
  model and grid the analytic standard error equals the 0.78.0 finite-difference
  one to the test's relative `1e-4` (`test-linpred-layers.R`, "se.fit has the
  matrix's shape").
- `frm_curve()`, `frm_curve_deriv()`: the grid once per layer, layer slowest, a
  new column `.eta`; the band and `"Sigma"` over all layers. The derivative in a
  `cs()` variable is `bcs_z[k]` with its own standard error (p013's fit: -0.0246
  with se 0.1353 at threshold 1, against 0 with se 0 on 0.78.0;
  `dev/predread-cscurve.R`, logs `cscurve-before.txt`, `cscurve-after.txt`).
- `frm_curve_feature()`: each layer searched on its own, through a `layer_pick`
  field that `frm_linpred()` and `frm_lp_basis()` read. `curve_extra_deriv()`
  reads the layers, and a new level's draw covaries across thresholds
  (`kronecker(J, E)`).
- `emmeans`: no `dpar` on a `cs()` fit takes the grid route and `rep.meas`, as
  brms; the design route had no column for the term.
- `get_predict()` keys the layers in `group`. `insight::get_predicted()` names the
  shared predictor, because it returns one number per row.
- frmtmb.sample's `posterior_linpred()` and `fitted(scale = "linear")`: brms's
  layouts, including `y`, `o1`, `o2` on a multivariate model (`draws_resp_stack()`
  takes `thresholds =`), compared in `dev/predread-sample-cells.R` (log
  `sample-cells.txt`): 11 layouts, each identical to brms's.

Three readers that had their own defect: `frm_linpred(dpar = "mu", type =
"response")`, `fitted(dpar = "mu")`, `conditional_effects(dpar = "mu")`,
`frm_curve(dpar = "mu", transform = TRUE)` and through them frmtmb.sample's
`posterior_epred(dpar = "mu")` and `posterior_linpred(transform = TRUE)` returned
`mu`, because an ordinal location's link slot is the identity, where brms applies
the family's link (brms at -0.9075: 0.2875; found by
`test-linpred-layers-brms.R`, relative difference 1.14; new helper
`dpar_response_link()`); `curve_cov_check()` read a row with a
standard error of exactly 0 on both routes as a disagreement of 1 (found on
`frm_curve(dpar = "mu")` at `x = 0`; fixed by comparing the difference); and
`fitted(dpar = "mu", scale = "linear")` gave threshold layers where brms gives the
shared part.

Pair p013: measured verdict "defect" on the baseline, "equivalent" now, L2 metric
Inf to 0 (section 3). **The adjudication table (`dev/equiv-data/adjudication.tsv`
line 10) still forces p013 to "defect", so `verdict` in results.tsv reads "defect"
until that row is removed at consolidation. I did not edit `dev/equiv-data/`.**

### `me()` on new data

brms 2.23.0 run at fixed parameter values (`dev/predread-me-brms-probe.R`, data
seed 41, 2000 equal draws, log `me-brms-probe.txt`). Each of these is a
measurement, not a recollection:

| case | brms |
|---|---|
| sd over draws of the epred, rows sdx 0.2, 0.4, 0.8, b = 1.013 | 0.1995, 0.4129, 0.7900 (`|b| sdx` = 0.2026, 0.4053, 0.8106) |
| `sdx = 0` in new data | "Measurement error should be positive." |
| no `sdx` column | "variables can neither be found ... 'sx'" |
| `NA` in `x` | that row `NA` with an rnorm warning |
| one latent across responses, `cor(y, y2)` by row | 0.0186, 0.0005, -0.0170: independent |
| one latent across dpars, `cor(mu, log sigma)` by row | 0.0048, 0.0001, 0.0313: independent |
| `me(gr = )`, two rows of one level | cor 1; rows of two levels -0.036 |
| `me(gr = )`, new level | refused unless `allow_new_levels = TRUE` |
| `me(gr = )`, `x` varies in a level | "Measured values and measurement error should be unique" |
| `me(gr = )`, no `g` column | "object 'g' not found" |
| correlated `me()`, `cor(Xme_x, Xme_w)` by row at corme 0.064 | 0.0133, 0.0245, 0.0013: `corme` not read |
| in sample, default `save_pars()` | redraws, with its warning; frmtmb reads the fit's latent values |

Change: `me_newdata_value()` reads `fit$me_draw`, an environment one replicate or
posterior draw holds, and then draws from `normal(x, sdx)`, once per predictor and
call. The store makes one value of one call read by every column of the predictor
(`me(x, sx) + me(x, sx):z`) and by every evaluation of the draw. `draws_fit_at()`
sets it; core's `predict()` and `simulate()` set it per replicate.

**The ML analog, chosen and measured; fix round 1 changed the summary half.** The
latent value of a new row is measurement noise in the grid, not an estimate. So
`predict()` and `simulate()` on new data draw it per replicate, also under
`propagate_error = FALSE`, because at fixed parameter values brms's spread still
holds it. `frm_linpred(se.fit = TRUE)`, `frm_curve()` and `emmeans()` report the
observed value and add its variance as a source that is not coefficient
uncertainty (`me_newdata_extra()`, `lp_extra_var_vec()`, `lp_extra_cov()`), as for
a new grouping level. `fitted()` reports the mean over the draw (section 8, M1).

Measured against brms (`dev/predread-me-compare.R`, data seed 41, sampler seed 5,
logs `me-compare-before.txt` on rellib-r16, `me-compare-after.txt`). Row sdx 0.2,
0.4, 0.8. "mean z" is `(mean_frmtmb - mean_brms) / MCSE`.

A. `predict(newdata, propagate_error = FALSE)` against `posterior_predict()` at the
frm() estimate, 4000 draws a side, sd ratio frmtmb / brms:

| model | before | after |
|---|---|---|
| `me(x, sx) + z` | 0.964 0.817 0.550 | 1.025 1.019 1.007 |
| `me(xg, sg, gr = g) + z` | 0.987 0.980 0.990 | 1.004 0.980 1.015 |
| `me(x, sx) + me(w, sw)` | 0.950 0.805 0.554 | 1.008 0.989 0.957 |
| `me(x, sx) + me(x, sx):z` | 0.936 0.815 0.548 | 0.997 0.987 1.003 |

B. The same 200 `frm_sample()` draws read by both packages (each draw one chain of
a Fixed_param brmsfit, 40 replicates a side), `posterior_epred()` sd ratio:

| model | before | after |
|---|---|---|
| `me(x, sx) + z` | 0.410 0.186 0.152 | 1.010 0.995 0.994 |
| `me(xg, sg, gr = g) + z` | 0.584 0.645 0.796 | 1.011 1.010 0.992 |
| `me(x, sx) + me(w, sw)` | 0.520 0.183 0.141 | 0.998 1.018 1.009 |
| `me(x, sx) + me(x, sx):z` | 0.402 0.157 0.156 | 1.017 1.011 0.984 |

`me(gr = )` on the same log: sd of row 1 minus row 3, the two levels, 0.5202 against
brms's 0.5201 (0.3854 before).

One more defect this found: **the emmeans design route left a `me()` term out of
the mean**, because `emm_is_design()` listed `mo()` and `mi()` and missed `me()`.
On `me-readers`'s fit, `emmeans(~ z)` read 0.9586 at z = 0, the intercept with no
`b x-bar`; `frm_linpred()` at the mean of `x` reads 0.9429, and brms's median over
its noisy draws is 0.9537 (MC error of about 0.008). brms's HPD is 0.245 to 1.756,
frmtmb's interval now 0.179 to 1.707, where the standard error of 0.051 gave 0.858
to 1.059 before (`dev/predread-me-readers.R`, log `me-readers-after.txt`).

The repro of the release review, `dev/rel078-fix2-me-newdata.R`, read 0.603 on its
own model and draws.

Divergences from brms kept on purpose, each loud:
- `NA` in `x` is refused. brms returns `NA` with an rnorm warning.
- A factor `me(gr = )` column with `NA`: brms's `gr[is.na(gr)] <- ` fails on a
  factor, warns "invalid factor level", and reads both rows as one level (cor 1
  above), against its own comment (#706) and against a character column, where
  each row gets a level of its own. frmtmb gives each row its own level. brms
  defect, candidate for `dev/upstream-bugs.md` (not edited; not reported).

Brms's independence across two predictors was matched at first; the coordinator
decided to share one latent as a stated divergence (section 8, D1).

### `resp`

brms 2.23.0 run for each case (`dev/predread-resp-brms-probe.R`, data seed 5, log
`resp-brms-probe.txt`). On a univariate model brms ignores `resp` in
`posterior_epred`, `posterior_linpred`, `posterior_predict`, `fitted`, `predict`,
`log_lik`, `bayes_R2`, `conditional_effects` and `pp_mixture`, whatever it holds:
`character(0)`, `c("y", "y")`, `"nosuch"`, `c("nosuch", "other")`, `NA_character_`,
`1L`. `resp_read()` in core and `draws_resp_read()` in frmtmb.sample do the same.
The brief named the first two; the unknown name follows, because brms answers it
the same way and one rule is simpler than a list. On a multivariate model an
unknown name stays refused.

`conditional_effects()`: names `y.y_x`, `y2.y2_x`, `y3.y3_x` without `resp` and on
`character(0)`, `y2.y2_x` on one, in the order asked on several, on both the fit and
draws (`ce_by_response()`, `ce_finalize(resp_prefix =)`). `pp_mixture()` of a
multivariate model takes one response, on both the fit and draws; none or several
get brms's "Argument 'resp' must be a single variable name". `single_response()`
reads a `resp_pick` field for it. `bayes_R2()` and `loo_R2()` read an empty `resp`
as every response (`brms_resp_index()`).

Divergence kept: a name given twice stays refused, per `draws_resp_set()` and the
brief, in `conditional_effects()`, `pp_mixture()`, `bayes_R2()` and `loo_R2()`.
brms gives two panels or two rows from `conditional_effects()` and `bayes_R2()`,
stops with an internal "no applicable method" in `posterior_epred()`, and refuses
in `pp_mixture()`.

### `update(<draws>, prior = new)`

brms's merge, through core's `update_merge_prior()` (now exported): a new row takes
its slot and the other rows stay; a prior with no rows keeps them. Test
`test-readers-resp-me.R`, "update(prior = ) merges".

### Compat rows changed (`R/compat.R`, each on its own line)

- `cs_pred()` x `frm_curve`: was the default "untested", now "works".
- `me()` x `predict`: "conditional" to "works".
- `me()` x `emmeans`: "conditional" to "works".
- `me()` x `fitted`, `me()` x `simulate`: status unchanged, note rewritten.
- `me()` x `marginaleffects`: new row, "conditional" (section 9, R2). `s()`, `t2()`,
  `gp_pred()` and `kind:covstruct` x `marginaleffects`: new rows, "conditional".
  `marginaleffects` is a new method in the vocabulary.
- The residual correlation terms x `predict`: status unchanged, note states that
  `cov = TRUE` answers on new data.

### New exports (sampling API)

`ce_by_response()`, because frmtmb.sample's method needs the joined panel object
core's method makes. `resp_refuse_single()`, brms's refusal that `pp_check_resp()`,
core's and frmtmb.sample's `pp_mixture()` share (the message-uniqueness test). And
`update_merge_prior()`, the merge core already had. frmtmb.sample also sets two
fit fields core reads, `me_draw` (as it sets `krige_draw`) and `resp_pick`.

Added by fix round 1: `ce_na_vars()`, for M4 on draws, and
`autocor_for_newdata()`, for m7 on draws. Fix round 2 took `ce_me_grouped()` back.

## 2. The audit

Which code builds the predictor, read from the source:

| builder | who reads it |
|---|---|
| `lp_eta_design()` | `frm_linpred()`, `frm_lp_basis()`, `ord_probs()`, `ord_prob_se()`, `predict_mean_se()`, `fitted()`, `predict()` off its fast path, `frm_curve*()`, emmeans's grid route, `get_predict()`, core CE, frmtmb.sample's `posterior_linpred`, `posterior_epred`, `posterior_predict` on new data or under `re_formula` |
| `eval_dpars()` | the objective's twin: in-sample `simulate()`, `predict()` fast path, `log_lik`, frmtmb.sample's `posterior_predict` fast path |
| `sim_newdata_design()`, `sim_newdata_dpars()` | `simulate(newdata = )`: own sums of the same parts; `frm_linpred()` for a body and now for `me()` |
| `emm_basis_design()` | emmeans, parametric predictors with no `re_formula`: own `model.matrix()` |
| `linpred_arma_cond()` | `frm_linpred()` on `mu` under `cov = FALSE` |
| `lp_basis_nl()`, the `nl_body` branch of `frm_linpred()` | a nonlinear body |
| `cs_offsets_add()`, `with_cs_offsets()` | `.cs` beside the shared `mu` for the densities and the simulators |

**What this table cannot see.** The review measured it against 14 mutants and it
missed 11: its "`frm_linpred(se.fit)` vs `frm_lp_basis`" cells compare two outputs
of one route, no cell compares a standard error with brms, `cs_cum`'s grid sat on
the reference level, and `me` had `sx = 0.3` and no `gr =`. Section 8, m3, has the
cells that close that. It is a record; the test files are the gate.

Each term kind by reader, at the frm() estimate against a Fixed_param brmsfit
(`dev/predread-audit.R`, data seed 101, 4000 equal draws), and on frm_sample()'s
draws (`dev/predread-audit-sample.R`, data seed 201, sampler seed 7, each draw a
chain). `A` agrees: relative difference below 1e-8, or each row's mean and
variance within the Bonferroni point at 0.001 over the rows. `x` differs. `-` no
reading: a side refuses, or brms draws where frmtmb reads a mean. Empty: the cell
does not apply. Generated by `dev/predread-audit-table.R`:

<!-- predread-audit begin -->
### core after the change

before (rellib-r16): 365 cells; 322 agree, 30 differ, 13 no reading.

after: 365 cells; 345 agree, 8 differ, 12 no reading.

| reader | plain | group | factor | offset | poly | mo | mo_inter | me | s | s_by | t2 | hsgp | ar_cov | ar_cond | nl | sigma_dpar | mixture | categorical | cumulative | cs | cs_cum | mv | mi |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| frm_linpred(link, newdata) | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |
| frm_linpred(link, in sample) | A | A | A | A | A | A | A | - | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |
| frm_curve(re_formula = NA) | A | A | A | A | A | A | A | A | A | A | A | A | A | - | A | A | A | A | A | A | A |  |  |
| frm_lp_basis(eta) | A | A | A | A | A | A | A | A | A | A | A | A | A | x | A | A | A | A | A | A | A |  |  |
| frm_linpred(se.fit) vs frm_lp_basis | A | A | A | A | A | A | A | A | A | A | A | A | A | x | - | A | A | A | A | A | A |  |  |
| fitted(scale = linear) | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |  |  |
| fitted(newdata) | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |
| fitted(in sample) | A | A | A | A | A | A | A | - | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |
| residuals(in sample) | A | A | A | A | A | A | A | - | A | A | A | A | A | A | A | A | A |  |  |  |  |  |  |
| predict(newdata, no error) | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |
| predict(in sample, no error) | A | A | A | A | A | A | A | - | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |
| simulate(newdata) | A | A | A | A | A | A | A | A | A | A | A | A | A | x | A | A | A | A | A | A | A |  |  |
| simulate(in sample, re_formula = NULL) | A | A | A | A | A | A | A | - | A | A | A | A | A | x | A | A | A | A | A | A | A |  |  |
| get_predict(link) | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |  |  |
| get_predict(response) | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A | A |  |  |
| emmeans ~ x (2 points) | A | A | A | A | A | A | A | - | A | A | A | A | A | A | x | A | - | - | A | A | A |  |  |
| conditional_effects(x) | A | A | A | A | A | A | A | - | A | A | A | A | A | A | A | A | A | - | x | x | x |  |  |

Cells that changed (before to after; metric is the relative difference, or the ratio of a z to its bound):

| shape | reader | was | now | metric was | metric now |
|---|---|---|---|---|---|
| ar_cov | predict(newdata, no error) | - | A | NA | 0.5399 |
| categorical | frm_curve(re_formula = NA) | x | A | Inf | 0 |
| categorical | frm_linpred(link, in sample) | x | A | Inf | 1.11e-16 |
| categorical | frm_linpred(link, newdata) | x | A | Inf | 0 |
| categorical | frm_lp_basis(eta) | x | A | Inf | 0 |
| categorical | get_predict(link) | x | A | Inf | 0 |
| cs | emmeans ~ x (2 points) | x | A | 0.004834 | 0 |
| cs | frm_curve(re_formula = NA) | x | A | Inf | 0 |
| cs | frm_linpred(link, in sample) | x | A | Inf | 0 |
| cs | frm_linpred(link, newdata) | x | A | Inf | 0 |
| cs | frm_lp_basis(eta) | x | A | Inf | 0 |
| cs | get_predict(link) | x | A | Inf | 0 |
| cs_cum | emmeans ~ x (2 points) | x | A | 0.08256 | 5.551e-17 |
| cs_cum | frm_curve(re_formula = NA) | x | A | Inf | 0 |
| cs_cum | frm_linpred(link, in sample) | x | A | Inf | 0 |
| cs_cum | frm_linpred(link, newdata) | x | A | Inf | 0 |
| cs_cum | frm_lp_basis(eta) | x | A | Inf | 0 |
| cs_cum | get_predict(link) | x | A | Inf | 0 |
| mixture | frm_curve(re_formula = NA) | x | A | Inf | 1.11e-16 |
| mixture | frm_linpred(link, in sample) | x | A | Inf | 2.22e-16 |
| mixture | frm_linpred(link, newdata) | x | A | Inf | 1.11e-16 |
| mixture | frm_lp_basis(eta) | x | A | Inf | 1.11e-16 |
| mixture | get_predict(link) | x | A | Inf | 1.11e-16 |

Cells that still differ:

- ar_cond, frm_lp_basis(eta): metric 0.1676 (len 5 vs 5)
- ar_cond, frm_linpred(se.fit) vs frm_lp_basis: metric 0.0317 (len 5 vs 5)
- ar_cond, simulate(newdata): metric 1.75 (z mean 6.81 var 3.05 crit 3.89)
- ar_cond, simulate(in sample, re_formula = NULL): metric 5.993 (z mean 28.22 var 3.54 crit 4.71)
- nl, emmeans ~ x (2 points): metric 0.04908 (len 2 vs 2)
- cumulative, conditional_effects(x): metric Inf (len 20 vs 5)
- cs, conditional_effects(x): metric Inf (len 20 vs 5)
- cs_cum, conditional_effects(x): metric Inf (len 20 vs 5)

### frmtmb.sample after the change

before (rellib-r16): 66 cells; 52 agree, 5 differ, 9 no reading.

after: 66 cells; 57 agree, 0 differ, 9 no reading.

| reader | group | offset | mo | me | s | hsgp | sigma_dpar | nl | categorical | cs |
|---|---|---|---|---|---|---|---|---|---|---|
| posterior_linpred(newdata) | A | A | A | A | A | A | A | A | A | A |
| posterior_epred(newdata) | A | A | A | A | A | A | A | A | A | A |
| posterior_linpred(in sample) | A | A | A |  | A | A | A | A | A | A |
| posterior_epred(in sample) | A | A | A |  | A | A | A | A | A | A |
| fitted(scale = linear, summary = FALSE) | A | A | A |  | A | A | A | A | A | A |
| posterior_predict(newdata) | A | A | A | A | A | A | A | A | A | A |
| conditional_effects(x, epred) | - | - | - |  | - | - | - | - | - | - |

Cells that changed (before to after; metric is the relative difference, or the ratio of a z to its bound):

| shape | reader | was | now | metric was | metric now |
|---|---|---|---|---|---|
| cs | fitted(scale = linear, summary = FALSE) | x | A | Inf | 0 |
| cs | posterior_linpred(in sample) | x | A | Inf | 0 |
| cs | posterior_linpred(newdata) | x | A | Inf | 0 |
| me | posterior_epred(newdata) | x | A | 0.7858 | 0.01622 |
| me | posterior_linpred(newdata) | x | A | 0.7827 | 0.02413 |

Cells that still differ:


<!-- predread-audit end -->

The cells that still differ, each read:

- **`ar_cond` x `frm_lp_basis(eta)`, `frm_linpred(se.fit)` against the basis:** the
  basis leaves the one-step part out, 0.168 relative. emmeans needs that: brms's
  `emmeans()` reads `incl_autocor = FALSE`, and that cell agrees. Now stated in
  `?frm_lp_basis`. `frm_curve()` refuses with its reason, where it refused with a
  message about a misread covariance. brms's `posterior_linpred(re_formula = NA)`
  gives the one-step mean there, and `frm_linpred()` agrees with it; so a curve of
  it is **filed: the design needs the column for the earlier residual, which
  `autocor_cond_fd_se()` differences today.**
- **`ar_cond` x `simulate()`:** by design. `simulate()` walks its own recursion, and
  `predict()`, which is brms's analog, agrees (`predict(newdata)`,
  `predict(in sample)` read A). In the code since 0.76.0.
- **`nl` x emmeans:** the divergence `?frmtmb-emmeans` states: the grid also holds
  the covariates of the parameters' formulas.
- **ordinal x conditional_effects:** the per-category default, a ledger divergence
  (`brmsfit-methods:217`).

Filed, not fixed, with numbers:

- **emmeans of a `mixture()` or `categorical()` with no `dpar` refuses**, and brms
  stops too ("Perhaps a 'data' or 'params' argument is needed") on this model, so
  this is not a difference.
- **`frm_linpred(se.fit = TRUE)` is refused on a body**, documented; `frm_lp_basis()`
  is the route.
- `me()` x `emmeans` and `conditional_effects`: brms's emmean is a median over noisy
  draws, so the cell has no exact reading (`-`). `me-readers-after.txt` has the
  numbers side by side, and the standard errors are in `test-me-newdata.R`.
- The cell `me` x `predict(newdata)` of the audit grid (`sx = 0.3`) has no power,
  because `(b sdx)^2` is small beside `sigma^2`; it read A before and after. Table A
  above is the reading.

## 3. The equivalence pairs

All 172 pairs, on the final build. Baseline
`frmtmb-wt-r16base/dev/r16base-equiv-out`, this run `dev/predread-equiv-out`, 8
processes, library `wt-predread-lib` then `rellib-r16` by absolute path
(`dev/predread-equiv-diff.R`, log `equiv-diff.txt`):

- Measured verdict changes: 1. **p013 `cs_pred()` x `frm_curve`, defect to
  equivalent, L2 metric Inf to 0.**
- `verdict` after the adjudication table: 0 changes, because the table forces p013
  (section 1). The headline therefore reads 83.6 until that row goes.
- The L2 metric of the other 171 pairs is unchanged, the simulators included: an
  unchanged seed draws the same replicates.

## 4. Tests

New: `tests/testthat/test-linpred-layers.R`, `test-linpred-layers-brms.R`,
`test-me-newdata.R`, `test-me-newdata-brms.R`, `test-resp-read.R`,
`extensions/frmtmb.sample/tests/testthat/test-readers-resp-me.R`.

Rewritten, each because brms answers otherwise (logs `resp-brms-probe.txt`,
`cs-layout-brms.txt`). No test dropped.

- `test-brms-postfit.R`: `pp_mixture(fit, resp = "z")` on a univariate fit,
  refusal to the same answer. `test-predict-lp-basis.R`: the refusal of an unknown
  `resp` moved to the multivariate case in `test-resp-read.R`. frmtmb.sample's
  `test-loo.R`: `log_lik(resp = "nope")` on a univariate model likewise.
- `test-brms-methods.R`, "a multi-column linear predictor": it recorded the
  N-vector as a divergence and now pins brms's columns, cs() included.
- `test-ce-bands.R`, `test-ce-parity.R`: panel names `y.y_x`, `x.x_w`, `y.y_w`.
- `test-ordinal-fitted.R`, `test-review-v29.R`, frmtmb.sample's
  `test-sampling-ported.R`: `dpar = "mu"` on the response scale and `transform =
  TRUE`, from the identity to the family's link inverse.
- frmtmb.sample's `test-incl-thres-draws.R` reads the layers through
  `frm_linpred()`, because `ord_linear_per_threshold()` is removed.

Rewritten by the final pass of fix round 1, each because the fix answers:
`test-compat.R` (the `marginaleffects` row, and one column per response from
`frm_linpred()` on a multivariate fit), `test-predict-newdata.R` and frmtmb.sample's
`test-loglik-autocor.R` (`cov = TRUE` on new data answers). frmtmb.sample's
`test-simulators.R` gains a guard on the branch of `log_lik()` that still names the
route as refused; it reads the same on rellib-r16. Fix round 2's rewrites are in
section 9.

Seen to fail on rellib-r16 (`dev/release/predread-seen-files/`):

<!-- predread-seen begin -->
```
test-compat.R pass=766 fail=2 err=0 skip=0 warn=0
test-dirichlet-multinomial.R pass=68 fail=0 err=0 skip=0 warn=0
test-linpred-layers.R pass=14 fail=12 err=12 skip=0 warn=0
test-me-newdata.R pass=42 fail=65 err=8 skip=0 warn=0
test-predict-newdata.R pass=12 fail=0 err=1 skip=0 warn=0
test-resp-read.R pass=4 fail=4 err=6 skip=0 warn=2
test-mixture-mvn.R pass=310 fail=1 err=0 skip=0 warn=0
test-loglik-autocor.R pass=134 fail=1 err=1 skip=0 warn=0
test-readers-resp-me.R pass=15 fail=26 err=7 skip=0 warn=0
test-linpred-layers-brms.R pass=1 fail=2 err=2 skip=0 warn=0
test-me-newdata-brms.R pass=26 fail=27 err=0 skip=0 warn=3
```
<!-- predread-seen end -->

Every suite, one file per process, gated (`FRMTMB_BRMS_FIT_TESTS=true`):

<!-- predread-suite begin -->
| package | files | pass | fail | error | skip | warn | core read from |
|---|---|---|---|---|---|---|---|
| frmtmb | 291 | 23468 | 0 | 0 | 2 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.eam | 32 | 1844 | 0 | 0 | 3 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.latent | 12 | 682 | 0 | 0 | 2 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.learn | 15 | 501 | 0 | 0 | 2 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.sample | 79 | 4224 | 0 | 0 | 1 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.spatial | 8 | 561 | 0 | 0 | 0 | 0 | C:/Users/adf44/source/r/wt-predread-lib |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 | C:/Users/adf44/source/r/wt-predread-lib |

files with no RESULT: 0 


files with a skip: 8 
- frmtmb test-curve-scale.R: skip 1
- frmtmb test-fuzz.R: skip 1
- frmtmb.coupling test-scale.R: skip 5
- frmtmb.eam test-scale.R: skip 3
- frmtmb.latent test-scale.R: skip 2
- frmtmb.learn test-scale.R: skip 2
- frmtmb.ode test-scale.R: skip 1
- frmtmb.sample test-scale.R: skip 1
<!-- predread-suite end -->

One session per suite, three orders, plain tier, on the final build
(`dev/release/predread-sess.log`): `frmtmb` forward 20799, reverse 20949, preload
21042 passes; `frmtmb.sample` 3872 in each. Fail 0, err 0 in all 27 sessions. A
session holds the files the plain tier reaches: 279 of 291 and 74 of 79 in those
two suites, and all files in the other seven. The review read the rest as the gated
`-brms` files. LEAK lines: 46, every one tagged `[with load ...]`; 0 untagged.

OpenBLAS 0.3.26 and 0.3.32, each with LAPACK, 4 threads, every gate on
(`dev/ciharden-scan.sh`, logs `dev/ciharden-log/predread-ob26.sum`,
`predread-ob32.sum`): 26 files reach the change, fail 0, err 0 in all 26 on both.
The R copies under `dev/ciharden-out/` are deleted.

`R CMD check --as-cran` in `dev/predread-check/`, then deleted: core, frmtmb.sample
and frmtmb.latent each `Status: 1 NOTE`. Each note is the V8 math-rendering one.

No test here builds its own library stack or child process.

## 5. What I did not do

- No curve under `cov = FALSE` (section 2).
- No `allow_new_levels = TRUE` on draws for a `me(gr = )` level: draws refuse that
  argument for every term, core's fit method takes it.
- No edit of `dev/equiv-data/` or `dev/brmsport-verdicts.tsv`.

## 6. The seam for lane latres, and overlaps

Latent residuals on new rows have the shape of `me()` here, and should read the
same two seams: a draw per replicate through a store on the fit (`fit[["me_draw"]]`,
set in `draws_fit_at()`, `predict_simulate()` and `sim_fit_draws()`), and a summary
at the mean whose variance goes through `lp_extra_var_vec()` and `lp_extra_cov()`
(`me_newdata_extra()` is the template, an `ed[["me_parts"]]` entry per source).

Overlaps: `R/predict.R` (`frm_linpred`, `frm_lp_basis`, `pred_design`,
`lp_eta_design`, `lp_extra_var_vec`, `lp_extra_cov`, `fitted_point`,
`single_response`; nlfix reads the nl branch moved into `lp_basis_layer()`),
`R/me.R`, `R/predict-brms.R`, `R/simulate-newdata.R`, `R/compat.R`, `NEWS.md`,
`NAMESPACE`, `man/`, frmtmb.sample `R/methods-draws.R`, `R/loo.R`,
`R/loo-predict.R`, `R/conditional-effects-draws.R`, `NEWS.md`.

## 7. Bump and floor

Minor for core and for frmtmb.sample: every change above is BREAKING in NEWS.
frmtmb.sample needs the new core, for `ce_by_response()`, `resp_refuse_single()`,
`update_merge_prior()`, `ce_na_vars()`, `autocor_for_newdata()`, and
`frm_linpred()`'s layers. frmtmb.latent changes one test, and no code or version.

## 8. Fix round 1 (review `dev/reviews/2026-11-02-predread.md`)

The review gave FIX: 0 BLOCKER, 7 MAJOR, 8 MINOR, NITs. The sections above now
describe the result. This section records what the round changed and the
evidence. Logs are in `dev/predread-log/` (gitignored). Each script names its
seed in its header.

### MAJOR

**M1. `fitted()` on a new `me()` row is the mean.** `response_mean_gh()`
(`R/me.R`) integrates on the response scale by Gauss-Hermite at 41 nodes, at the
row's total variance `sum (g sdx)^2`, with the exact remainder of `me(x):me(z)`
(m4). Its signature is that of latres's copy, so the merge keeps one. Latres's
body also guards the tails, and that body is the one to keep (S4).
`frm_linpred()` stays at `x`, where the value is exact.

`dev/review-predread-mefam.R`, seed 7971, `sdx = 1.2`, 4000 draws, at fixed
values (`mefam-fix.txt`; the review's copy is `mefam-review.tsv`):

| family | brms mean | `fitted()` | brms sd | `Est.Error` | brms 2.5% | `Q2.5` | brms 97.5% | `Q97.5` |
|---|---|---|---|---|---|---|---|---|
| gaussian | 1.0856 | 1.1017 | 0.8794 | 0.8840 | -0.6544 | -0.6308 | 2.7827 | 2.8343 |
| poisson | 4.5298 | 4.4618 | 5.0098 | 4.9549 | 0.5005 | 0.5149 | 17.7465 | 17.3268 |
| negbinomial | 4.7793 | 4.8011 | 6.0956 | 6.0421 | 0.4253 | 0.4417 | 20.9559 | 20.2285 |
| lognormal | 4.9643 | 5.0671 | 5.2056 | 5.2243 | 0.6785 | 0.6651 | 18.1292 | 18.7290 |
| gamma | 4.2090 | 4.2891 | 4.0325 | 4.3556 | 0.5837 | 0.5777 | 15.5425 | 15.6895 |
| beta | 0.7192 | 0.7202 | 0.1499 | 0.1504 | 0.3712 | 0.3701 | 0.9370 | 0.9363 |
| bernoulli | 0.7371 | 0.7388 | 0.1643 | 0.1673 | 0.3390 | 0.3357 | 0.9584 | 0.9587 |
| zero_inflated_poisson | 3.5185 | 3.4774 | 4.9741 | 4.4107 | 0.3121 | 0.3153 | 14.9249 | 14.7284 |
| hurdle_poisson | 3.5264 | 3.5350 | 4.8598 | 4.8784 | 0.8427 | 0.8405 | 15.43 | 15.5294 |
| cumulative, category 1 | 0.1966 | 0.1967 | 0.1105 | 0.1135 | 0.0503 | 0 | 0.4724 | 0.4191 |
| categorical, category 1 | 0.3112 | 0.3336 | 0.1045 | 0.1424 | 0.1157 | 0.0546 | 0.5192 | 0.6126 |

Three rows need a word:

- Cumulative: the category's interval is a Wald interval clamped to 0 and 1. It
  is in the support, and it is not a quantile. Filed below.
- Categorical: frmtmb shares one latent across the category predictors (D1), and
  brms draws one per predictor. The difference is the D1 divergence.
- Zero-inflated Poisson: `Est.Error` 4.41 against 4.97, and the quantiles agree.
  The brms sd is 4000 draws of a heavy tail; the review's table gives the same
  brms value, and I did not separate that error from the delta method's.

Against the truth, `dev/review-predread-sesim.R`, seeds 1 to 200
(`dev/predread-log/sesim-fix.rds`, read by `sesim-sum.R`): truth / estimate 0.993
to 1.002 over five rows. At `sdx = 1`: `fitted()` 2.107 against 2.11, and 6.065
against 6.02. The review read 1.28 times.

**M2. The interval stays in the support.** `fitted_response_band()`: the link band
through the inverse link. Same table. Same simulation, at `sdx = 1`: the interval
covers the true mean at the drawn latent on 95.5 and 98.0 percent of seeds,
against the review's 88 and 91 percent. Negative lower ends: 0, against 200 of
200.

**M3. The draw variance on families with their own mean.** `predict_mean_se()`
adds the term's variance, and so does the `conditional_effects()` band. See the
`Est.Error` column above, where the review read 0.1496 at every `sdx`. NEWS is
corrected: "the mean of the draw" is true of `x`, and is no longer said of the
response.

**M4. `conditional_effects()` on `me(gr =)` answers, on fits and on draws.** The
group is held at `NA` (`ce_na_vars()`), and `ce_me_grouped()` passes
`allow_new_levels`. `misc-fix.txt`: effect `xg` reads -1.4999, 0.51866, 2.5372,
4.5558 on the fit and -0.966, 0.501, 2.266, 4.738 on draws. `method = "predict"`
answers. With `conditions` on the group, both refuse, as brms does.

**M5. `simulate(newdata)` of `me()` beside `ar(cov = FALSE)` counts the term
once.** The `sim_newdata_dpars()` sum branch rebuilds `lp_eta_design()` per
replicate, so it does not read nlfix's split. `misc-fix.txt`, block G:

| build | mean by row |
|---|---|
| `rellib-r16` | 1.62, 1.61, 1.60, 1.61 |
| this lane, before | 1.63, 1.08, 1.94, 2.78 |
| now | 1.63, 1.614, 1.612, 1.613 |

The test holds two response columns identical. On the merge
(`dev/predread-merge-check.sh`, library `wt-predread-lib/merge`),
`test-me-newdata.R` has 104 passes.

`fitted()`'s `Est.Error` under `cov = FALSE` carries the term through
`linpred_arma_cond()`, against brms's sd 0.364, 0.478, 0.479, 0.476.

**M6. Two grids, same `x` and `sdx`, read one draw in `frm_curve(contrast =)`.**
Measured first, `dev/predread-me-contrast.R`, seed 7402 (`me-contrast.txt`):

| quantity | value |
|---|---|
| coefficient of `z` | 0.6856, se 0.0514 |
| `frm_curve(contrast =)` se | 0.0514 |
| frmtmb `emmeans` pairs SE | 0.9662 (estimate -0.6856) |
| brms `emmeans`, HPD | -2.668 to 1.101, width 3.77 |
| brms two calls, sd of the rows | 0.9777, 0.9565 |
| brms one call, rows 1 and 3, sd of the difference | 0.9727 |
| `me(gr =)`, `frm_curve(contrast =)` se | 0.039, the coefficient's se |
| `me(gr =)`, brms one `newdata`, sd of the difference | 0 |
| `me(gr =)`, brms two calls, sd of the difference | 0.527 |

Where the two agree and where they do not:

- `frm_curve(contrast =)` and brms disagree on a row-level `me()`. brms draws each
  row, and frmtmb reads one draw for rows with the same `x` and `sdx`. This is a
  stated divergence: a contrast of two grids that differ in `z` asks about `z`, and
  the error in `x` is not in that difference. Help of `frm_curve()`.
- On `me(gr =)` in one `newdata`, the two agree: the term adds no variance.
- frmtmb's `emmeans` agrees with brms in kind: two grid rows read two draws, SE
  0.97, against a brms HPD width of 3.77. Help of `frmtmb-emmeans`.

**M7. `marginaleffects::predictions()` refuses the default vcov, with a message.**
`get_vcov.frmtmb_fit()` refuses a `me()` fit (class `frmtmb_me_interop`) and names
`frm_linpred(se.fit = TRUE)`. `vcov = FALSE` gives the mean, through
`fitted_me_moments()`. Compat row `me()` x `marginaleffects`, status "refused";
method `marginaleffects` is new in the vocabulary.

### For the user

**D1. One latent, shared across responses and dpars on new data, on fits and on
draws.** The store key in `me_newdata_latent()` is the term and not the predictor.
This is a stated divergence from brms: the Stan program has one latent, and brms's
`prepare_predictions_sp()` draws one per predictor on new data.

`share-fix.txt`, `cor(y1, y2)` by row:

| reader | row 1 | row 2 | row 3 |
|---|---|---|---|
| `posterior_epred()` | 1 | 1 | 1 |
| `posterior_predict()` on draws | 0.7208 | 0.9148 | 0.9796 |
| `predict()` on a fit | 0.7366 | 0.9211 | 0.9781 |
| target if shared, predictive | 0.7385 | 0.9184 | 0.9782 |

Two dpars, predictive skew by row: 0.3338, 0.3637, 0.6349. The review's
`y1 - y2` widths are 15.08 (brms) and 3.18 (one latent). NEWS and `?frmtmb-me` give
the reason and the numbers. `dev/upstream-bugs.md` has a `brms-NN` entry; its
number is settled at consolidation, and nothing is posted.

### MINOR

- **m1.** `curve_spec()` and `emm_target()` read `resp_read()`.
  `frm_curve_deriv()`, `frm_curve_feature()` and `emmeans()` ignore `resp` on a
  univariate model. `test-resp-read.R`.
- **m2.** Panels are effect-first (`ce_by_response()`), and `y2.y2_f` is found.
  `misc-fix.txt`: brms, fit and draws read `y2.y2_x y3.y3_x y2.y2_z y3.y3_z
  y2.y2_f y3.y3_f`.
- **m3.** The audit has power: `dev/predread-audit2.R`, 20 cells with standard
  errors against brms on non-reference grids. Base: 20 of 20 agree. It catches 14
  of 14 mutants (table below), where the first audit missed 11 of the 14. It is a
  record. The tests are the gate.
- **m4.** `me(x):me(z)` keeps the exact product remainder, as an entry with
  `k = NA` in `me_newdata_extra()`. `test-me-newdata-brms.R`.
- **m5.** `NA` in `x` or `sdx` is refused, where brms gives `NA` for that row with
  a warning. Stated in NEWS and `?frmtmb-me`. `test-me-newdata.R`.
- **m6.** The label: `frm_curve()` prints "the band holds variance that is not
  coefficient uncertainty" through `check$extra`. The help of `fitted()`,
  `frmtmb-emmeans`, `conditional_effects()` and `frmtmb-marginaleffects` states it.
  Core NEWS states that a name given twice is refused.
- **m7.** `predict(newdata)` and `posterior_predict(newdata)` under `cov = TRUE`
  answer through `autocor_for_newdata()`, the route `simulate(newdata)` takes. Only
  a family's own simulator refuses. On the merge with latres: `test-me-newdata.R`
  passes. The final pass found that `log_lik(newdata =)` still said this route was
  refused. It now tests `draws_pp_refuses_newdata()`, the condition the route uses.
- **m8.** Three older defects:
  - Multivariate with no `resp`: FIXED. `frm_linpred()` stacks one column per
    response, as brms does. `frm_lp_basis()` and the curves refuse and name `resp`.
  - A dpar with no predictor on `scale = "linear"`: FILED, not fixed. frmtmb reads
    the link value, brms the natural one. Shape 0.776 against 2.17, `disc` 0
    against 1, 11 dpars by 5 readers, 55 cells of 683, fit and draws
    (`dev/review-predread-dpar.R`).
  - `prior_summary()` of draws marks default rows `source = "user"`: FILED, not
    fixed (`dev/review-predread-last.R`).

### NITs

- The `frm_curve(transform = TRUE)` message on a mixture no longer says `mu`.
- `test-me-newdata.R`'s header cites `test-readers-resp-me.R`.
- Comments and test headers cite scripts, not gitignored logs.
- The seen-to-fail block of section 4 is generated from the logs.
- The OpenBLAS line is corrected: both versions ran (section 4).

### The audit's mutants (`dev/predread-audit2-mut.sh`)

Cells that differ, of 20 (`mut-*-audit2.txt`):

| mutant | cells | mutant | cells |
|---|---|---|---|
| `cs_se` | 4 | `me_shared_z` | 1 |
| `cs_se_rev` | 2 | `me_gr_perrow` | 2 |
| `cs_factor_nd` | 4 | `me_store_off` | 1 |
| `me_var_double` | 3 | `emm_me_design` | 1 |
| `me_var_drop` | 3 | `link_identity` | 2 |
| `me_half_sd` | 2 | `curve_check_ratio` | 1 |
| `resp_read_off` | 1 | `emm_cs_design` | 2 |

### What the consolidation must do

This lane does not resolve these. Each needs both sides in one tree.

- **S2.** In nlfix's `frm_linpred()` body branch, change
  `out <- if (type == "response") lp[["link"]]$linkinv(eta) else eta` to use
  `dpar_response_link(rspec, lp, dpar)$linkinv(eta)`.
- **S3.** In `latres_skip_se()`, keep `dim(v)`, as `v` is a matrix under layers.
- **S4.** Keep latres's copy of `response_mean_gh()`, in `R/autocor-latent.R`, and
  delete the one in `R/me.R`. The signatures are identical. Latres's body sets aside
  a node where the family's mean fails, and gives `NaN` where a cell lost real
  weight, and the copy in `R/me.R` does not: a `me()` row under `trunc()` at a large
  `sdx` needs that guard, and no test here reaches it, so add one on the merge. Keep
  one rule:
  `fitted_me_moments()` and `me_row_loads()` take a total variance that adds the
  latent residual, and they replace the calls to `latres_fitted_mean()` in
  `fitted_point()` and `get_predict()`. `extra_var` holds both variances.
- **S5.** `R/compat.R`, `group:autocor` x `nl`: take nlfix's text.
- **Conflict hunks.** `get_predict()` in `R/interop.R` and `idx_at` in
  `methods-draws.R` keep both sides. `dev/predread-merge-check.sh` holds the
  resolution. In `get_predict()` the `dots["resp"] <- list(resp_read(...))` line
  must stay, or `test-resp-read.R` errors on `character(0)`.
- **p013.** Remove the `dev/equiv-data/adjudication.tsv` row. The measured verdict
  went from defect to equivalent (L2 Inf to 0). I do not edit that file.
- `test-latres-brms.R` and `test-nl-brms.R` have no reading on a merge, from the
  review or from this lane.
- The last merge run (`dev/release/predread-merge.log`) precedes the last edits. It
  read `test-me-newdata.R` 104, `test-me-newdata-brms.R` 53, `test-linpred-layers.R`
  83, `test-readers-resp-me.R` 100, all passes, and `test-resp-read.R` 86 passes and
  1 error. The `get_predict()` line above addresses that error, and I did not rerun
  the merge after it. The merge has also not read the `log_lik()` change or the
  rewritten pins of `test-compat.R`, `test-predict-newdata.R` and
  `test-loglik-autocor.R`.

### Filed and left

- `marginaleffects` on an exact `gp()` at an unseen position carries no kriging
  variance: the same class as M7, older in kind.
- `frm_curve()` under `cov = FALSE` refuses, and names `frm_linpred(se.fit =)`.
- The `Q` columns of an ordinal or categorical category are a clamped Wald interval.
- `allow_new_levels` on draws is refused.

### The final pass

The final pass found four failures in `dev/predread-final3.sh`, all from fixes of
this round:

- `test-compat.R`, two: the method vocabulary and multivariate `frm_linpred()`.
  Pins rewritten.
- `test-predict-newdata.R` and `test-loglik-autocor.R`: pins of the old refusal.
  Rewritten.
- The second of these showed a defect in the code: `log_lik(newdata =)` said
  `posterior_predict(newdata =)` is refused under `cov = TRUE`, where it now
  answers. Also stale: the refusal message in `draws_pp_setup()`, the help of
  `posterior_predict()`, and the compat row `ar()` x `predict`. All are corrected.

Then the whole pass ran again from the install, with no source edited after it
started. A session restart stopped one rerun in its suite step, and no result of
that rerun is used. Results: section 4 for the suites, sessions, OpenBLAS and
check; section 3 for the pairs; section 2 for the audit table.

## 9. Fix round 2 (`Short check of fix round 1` in the review)

This section records the last round. Where it differs from section 8, this section
holds: R2 replaces section 8's M7, and m5 and the cumulative item below replace
section 8's m5 and its filed `Q` line. Rulings are the coordinator's.

### R1. `fitted()`'s `Q` when more than one predictor reads the term

The one line in `me_band_end()`, as ruled. `dev/review-predread-r1-band.R`, seed
8401 (`r2-band.txt`), `Q2.5` at `sdx` 0.4, 1.2 to `Q97.5` at `sdx` 0.4, 1.2:

| model | before | now | draws |
|---|---|---|---|
| `me()` in `mu` and `zi` | 0.49, 0.33 to 2.88, 3.10 | 1.23, 0.66 to 2.11, 2.69 | 1.30, 0.69 to 2.02, 2.67 |
| mixture | 2.61, 1.51 to 7.09, 12.23 | 3.07, 1.62 to 6.03, 11.42 | 3.14, 1.62 to 5.94, 11.18 |
| poisson | 1.71, 0.55 to 5.37, 16.61 | 1.71, 0.55 to 5.37, 16.59 | 1.75, 0.55 to 5.34, 16.07 |

**The new test found one more defect, in the same function.** On another data set
of the same recipe the two ends were swapped: `Q2.5` 2.235, `Q97.5` 1.292, against
draws of 1.389 to 2.164. The sign came from `mm$grad`, which moves the first
predictor alone, and a `zi` that rises with the same term turns the mean the other
way. The sign now comes from the mean at the draw's own `+1` and `-1`. This is one
line beyond the ruling, and it changes no number in the table above.

The test is "fitted()'s interval holds when two predictors read the term",
`test-me-newdata.R`. `dev/predread-r1-mutant.R` (`r1-mutant.txt`): fix round 1's
`me_band_end()` fails 10 of 18 expectations, the one-line change with the old sign
fails 6, and the current code fails 0. I did not add an audit cell: the tests are
the gate.

### R2. `marginaleffects`: a classed warning in place of the refusal

`get_vcov()` gives a warning of class `frmtmb_interop_se_part`, and compat `me()` x
`marginaleffects` is "conditional". Coefficient, slope and contrast answers are back
(`r2-check3.txt`, seed 7402), and each matches rellib-r16:

| call | now | rellib-r16 |
|---|---|---|
| `avg_slopes(z)` | -0.4909, se 0.0485 | -0.491, se 0.0485 |
| `avg_slopes(xo)` | 1.175, se 0.0526 | 1.175, se 0.0526 |
| `avg_comparisons(f)` | 0.6427, se 0.0935 | 0.643, se 0.0934 |
| `hypotheses("b2 = 0")` | -0.4909, se 0.0485 | -0.491, se 0.0485 |
| `predictions(newdata)` | se 0.085, 0.092, with the warning | se 0.085, 0.092 |

**The same warning reaches a smooth, an exact `gp()` and an `hsgp()`.** This was
cheap: one test on `re_blocks`. Compat rows `s()`, `t2()`, `gp_pred()` x
`marginaleffects` are "conditional".

**It can truthfully reach every model with a group-level term, and it does.** Compat
row `kind:covstruct` x `marginaleffects` is "conditional". `dev/predread-mfx-se.R`,
seed 7411 (`r2-mfx.txt`), four rows, the last at an `x` far from the data:

| fit | `predictions()` se | `fitted()` `Est.Error` |
|---|---|---|
| gaussian, `(1 \| g)` | 0.2089, 0.2064, 0.2084, 0.2393 | 0.1037, 0.1205, 0.1109, 0.1582 |
| poisson, `(1 \| g)` | 0.1896, 1.4350, 0.1293, 0.5648 | 0.1240, 0.4875, 0.1490, 0.6231 |
| lognormal, `(1 \| g)` | 0.1973, 2.0073, 0.1811, 0.6728 | 0.0883, 1.0480, 0.0867, 0.4087 |
| poisson, no group | 0.0729, 0.0723, 0.1132, 1.4503 | the same to these digits |

So the error goes both ways: too wide where the level is well known, too narrow
elsewhere. The review of latres reads 0.0774 against 0.2019 on `poisson` with
`(1 | g)`, and lognormal ratios 0.48 to 0.59. These agree in kind.

Under `re_formula = NA` the two agree to these digits, and `get_vcov()` gets
`re_formula` in its dots, so that call gives no warning. The review's "`(1 | g)`
agrees to all digits" is on a group with no effect in its data
(`dev/review-predread-r1-check3.R`). I read that as a variance at its boundary, and
I did not print that variance.

**One thing the coordinator may want to weigh.** A default `avg_slopes()` on any
mixed fit now warns, though for a coefficient on a linear predictor it is right, and
`get_vcov()` cannot tell. `marginaleffects` on `glmmTMB` gives the same standard
errors with no warning (`test-v07.R` pins them to 1e-3), so this warning is stricter
than that package's convention. To drop it, delete the third clause of
`interop_warn_se_part()`, the `kind:covstruct` row, and the group-level case of its
test; the five rewrites of test files take `mfx_se_part` and need no other change.

Tests rewritten for it, each only to allow the warning: `test-marginaleffects-args.R`,
`test-famgaps.R`, `test-mo-coding.R`, `test-v07.R`. `get_vcov()` is called once per
`predictions()` or `avg_slopes()` call, measured, so the warning comes once.

### Cumulative `Q`

`prob_band_quantiles()` replaces `prob_clamp_quantiles()`: a Wald interval on the
logit, the route `conditional_effects()` has. On a new `me()` row it is centered on
the value at the observed `x`, where the logit is close to symmetric in the draw, and
its spread comes from the mean's own variance. `dev/review-predread-r1-check1.R`,
seed 7971, `x = 1`, `sdx = 1.2` (`r2-check1.txt`):

| category | was | now | `conditional_effects()` | brms, 100000 draws |
|---|---|---|---|---|
| 1 | 0 to 0.419 | 0.049 to 0.462 | 0.046 to 0.479 | 0.048 to 0.469 |
| 2 | clamped Wald | 0.096 to 0.312 | 0.084 to 0.344 | 0.067 to 0.234 |
| 3 | clamped Wald | 0.171 to 0.333 | 0.197 to 0.296 | 0.135 to 0.244 |
| 4 | 0.092 to 0.734 | 0.152 to 0.717 | 0.135 to 0.745 | 0.138 to 0.739 |

The two end categories are monotone in the draw, and the interval is close. The two
middle categories are not: a hump puts most draws near its top, and no symmetric
interval on the logit reads that. Category 3's upper end is 0.333 where brms's is
0.244, and `conditional_effects()` is off the same way. An exact reading needs the
draw's own quantile per cell, which is not this round's change. It is filed below.

Categorical category 1: 0.120 to 0.626, against frmtmb draws' 0.102 to 0.624. A
plain category probability moves with the same formula, and multinomial counts read
it on the share of the trials.

### M6

The test is now a tolerance. `me_same_value()` is `all.equal()`'s tolerance, relative
to the pair and the column's scale. `dev/review-predread-r1-check2.R`
(`r2-check2.txt`):

| grid | was | now |
|---|---|---|
| `dx = 0` | 0.0485 | 0.0485 |
| `dx = 1e-12` | 0.9979 | 0.0485 |
| `dx = 1e-6` | 0.9979 | 0.9979 |
| `sdx` 0.6 and 0.6000001 | 0.9979 | 0.9979 |

The last two are past `all.equal()`'s tolerance, 1.5e-8 relative, so they read two
draws. A jump remains at the tolerance, as any rule of one draw or two has one.

`?frm_curve` and `?frmtmb-emmeans` say what each reader answers, and which reader
gives which effect: `frm_curve(contrast =)` answers one unit with one unknown `x*`;
`emmeans` answers two new rows with their own errors, as brms's draws do; for `z`
read `frm_curve(contrast =)` or the coefficient, and for `x` read
`frm_curve_deriv()` (0.0526, where a contrast in `x` reads 0.998).

### m5

A row with no `x` or `sdx` reads `NA`, with a warning of class `frmtmb_me_missing`,
and the other rows answer. `dev/predread-me-na.R`, seed 7421 (`r2-na.txt`):

- `frm_linpred()`, `se.fit`, `fitted()`, `predict()`, `simulate()`, `frm_curve()`,
  `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`: `NA` on those
  rows, on gaussian, poisson and `me(gr = )`. A plain predictor with `NA` reads the
  same on every one of these.
- Under `me(gr = )` a row with no value does not take its level's draw from anyone,
  and its level's other rows answer.
- `fitted(fit, newdata = <the fit's own data>)` with 2 rows dropped for `NA` in `x`:
  200 rows, `NA` on rows 5 and 9.
- One warning per missing column per reader's call, however many designs the reader
  builds, in the frame of the reader. On draws, 1 warning over all draws.
- brms reads `NaN` and warns. An absent column, or a non-numeric one, still stops.
- `posterior_predict()` on draws raises one R warning per draw on such a row, "NAs
  produced", as a plain predictor's row does today. I left that as it is, in both.
- A draw stream with no `NA` is unchanged: the draw is the same `rnorm()` call.

### m8

- `frm_curve_feature()` with no `resp` on a multivariate fit asks for `resp`, through
  `curve_one_response()`, the check `frm_curve()` has.
- `insight::get_predicted()` on a multivariate fit: class `get_predicted`, columns
  `Row`, `Response`, the columns of `data`, `Predicted`, dim 6 x 5, equal to brms's
  layout and values (`r2-check5b.txt`), and its `ci_data` joins on `Row` and
  `Response`, as insight's method needs.

### Exports

`ce_me_grouped()` is internal. frmtmb.sample reads `length(ce_na_vars(x, TRUE)) > 0L`.
Five stay: `ce_by_response`, `resp_refuse_single`, `update_merge_prior`, `ce_na_vars`,
`autocor_for_newdata`. Section 1's list holds except for that name.

### NEWS and help

- The categorical mean under D1: NEWS, `?fitted`, `?frmtmb-me`. 0.3336 against 0.3106.
- The constant dpar as a divergence: NEWS and `?fitted`. 2.01 against 7.48 unmodeled,
  2.01 in both with `shape ~ 1`, and 0.776 against 2.17 on another fit.
- `se__` on lognormal and hurdle fits with no `me()`: NEWS. 0.068 against 0.519.
- `log_lik(newdata =)` still refuses where brms answers, and the predicate changes
  only the hint: frmtmb.sample's NEWS.

### Filed

- MAJOR for the next round, if the warning is not enough: `predictions()` on a fit
  with a group-level term, a smooth or a field has the wrong level's standard error,
  with numbers above. The warning says so, and it does not fix it.
- A divergence from brms: `se__` of `conditional_effects()` on a family with a link is
  on the link scale, 0.068 against 0.519 on lognormal and hurdle, 0.102 against 0.809
  on poisson. Its help says "link scale", and brms's is on the response scale.
- MINOR: a middle category's `Q` under a new `me()` row. Table above.
- `posterior_predict()` on draws, one R warning per draw on a row of `NA`, plain or
  `me()`.

### What the consolidation must do

Unchanged from section 8, with the review's updates: S2, S4a (delete `response_mean_gh()`
from `R/me.R`), S4b, S5, and S4c as ruled. The review's
`dev/review-predread-merge2-resolve.R` matches this tree: I changed no text it
searches for. Its exact matches are the `loads` lines of `fitted_me_moments()`, which
I did not edit, and the row `"Under cov = TRUE predict() draws`.

`dev/predread-merge-check.sh` on the final build, against latres and nlfix as they
stand, library `wt-predread-lib/merge`, since deleted (`merge-r2.txt`): 5 of 5 files,
512 passes, fail 0, err 0: `test-me-newdata.R` 159, `test-me-newdata-brms.R` 53,
`test-linpred-layers.R` 83, `test-resp-read.R` 111, `test-readers-resp-me.R` 106. This
keeps the older resolve of `get_predict()`, both sides, so it does not apply S4b.
`test-latres-brms.R` and `test-nl-brms.R` were not in this run; the review's merge read
them green before this round's edits.

### The final pass

The first final pass found two more effects of the cumulative change, both on
`fitted()`'s category branch:

- A `mixture_mvn()` fit reaches that branch with a matrix of means, not of
  probabilities. The logit gave 4 `NaN` warnings, seen as `warn = 4` in frmtmb.latent.
  The old clamp was wrong there too: it held a mean's `Q` between 0 and 1, so the
  interval left its own estimate. The branch now reads `fam_is_polytomous()`, and such
  a fit has a plain Wald interval. New expectation in frmtmb.latent's
  `test-mixture-mvn.R`, which fails on rellib-r16, stated in NEWS. So R CMD check ran
  on frmtmb.latent too.
- `test-dirichlet-multinomial.R` pinned `Q97.5` at exactly 50 trials, the clamp's
  value. It now holds the band inside the trials and above the estimate. This rewrite
  passes on rellib-r16 as well, as the bound is true of both.

The whole pass then ran again from the install, with no source edited after it started.
The suites, sessions, OpenBLAS and check are in section 4, the pairs in section 3, and
the audit in section 2, all generated or rewritten on this build. Every `lib:` line
reads `wt-predread-lib` by absolute path, 474 of 474 in the suites, and the seen list
reads `rellib-r16`. The logs `r2-*.txt` and `r1-mutant.txt` were taken on this build.

Coverage of the 14 mutants of section 8 was not run again. `dev/predread-audit2.R` on
this build still agrees on 20 of 20.
