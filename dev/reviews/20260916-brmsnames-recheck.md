# Recheck: lane brmsnames, plan item 2.6c, punch round 1

Adversarial recheck of `frmtmb-wt-brmsnames` (base `aa9227e`) after punch
round 1. The brief was to FALSIFY the claim that the round-0 BLOCKER and
four MAJORs are fixed.

Libraries. BASE arm: `C:/Users/adf44/source/r/rellib-r3`, read only.
LANE arm: `C:/Users/adf44/source/r/brmsnames-lib`, read only. Nothing was
installed anywhere. `dev/brmsnames-rev-verify.R` was rerun first: frmtmb
981 functions, 1 differs (the known `covstruct_registry` false positive);
frmtmb.sample 197, 0; frmtmb.learn 44, 0. No `R/` file in the worktree is
newer than the lane build. So the lane build is the worktree.

brms 2.23.0 fits used `pinlib` ahead of the user library, with
`StanHeaders == 2.32.10` asserted in each script, and compiled fresh
(`brm()` without `file =`; no `FRMTMB_STAN_CACHE` involvement). Eleven
programs were compiled, about 70 s each.

No git operation. No package file edited. Suites and `R CMD check` were
not rerun (the verifier ran them); nothing here depends on them.

Scripts added, all under `dev/`, prefix `brmsnames-rev2-`: `data.R`,
`brms.R`, `natural.R`, `prior.R`, `families.R`, `mixture.R`, `names.R`,
`b2lane.R`, `collide.R`, `collide-brms.R`, `ranef.R`, `ranef-pos.R`,
`simulate.R`, `silent.R`, `interop.R`, `hyp.R`, `er.R`, `wiener.R`. Logs
are in `dev/brmsnames-rev2-log/`. Fitted objects are in `dev/stan-cache/`
(`brmsnames-rev2-*.rds`, gitignored).

Verdict: **mergeable after punch**. The fixes for all five round-0
findings hold under every construction I made against them. Two new
defects came out of the fixes themselves. One is a silent wrong answer on
brms's own spelling, which is the class round 0 ranked BLOCKER. The other
is a regression on draws of the kind round-0 MAJOR 2 described, reached by
a route the fix does not cover. Both fixes are narrow.

---

## New findings

### BLOCKER 1. A mixture's `theta1` carries brms's name and a log ratio

Construction: `dev/brmsnames-rev2-mixture.R base|lane|brms`, data seed
64, n = 400, `bf(y ~ x)` with `mixture(gaussian, gaussian)`, 87.75% of
rows from component 1. frmtmb.sample seed 3 (R seed 5), chains 1, iter
400. The brms fit: chains 1, iter 400, seed 1.

MAJOR 3's fix marks every intercept-only distributional parameter with
no written formula as brms's natural-scale parameter
(`brms_coef_table()`, `R/brms-names.R` line 141). A mixture's `theta1`
passes that test, and its stored link is the identity, so the column is
renamed `theta1` and "inverted" by the identity. But frmtmb's `theta1`
is the log ratio of the mixing weights. In brms, `theta1` is the mixing
PROBABILITY, and `theta2` is its complement.

| quantity | base | lane | brms |
| --- | --- | --- | --- |
| name in `variables(fit)` | `theta1_Intercept` | `theta1` | `theta1`, `theta2` |
| estimate on the fit | 1.9662 (a coefficient) | 1.9662 | n/a |
| draws column mean | `theta1_Intercept` 1.9650 | `theta1` 1.9650 | `theta1` 0.8759, `theta2` 0.1241 |
| `hypothesis(fit, "theta1 = 0.5", class = NULL)` | error "object 'theta1' not found" | **1.4662** | n/a |
| `hypothesis(ds, "theta1 = 0.5", class = NULL)` | not reachable | **1.4650** | 0.3759 |
| `posterior_summary(ds, variable = "theta1")` | not reachable | **1.9650** | 0.8759 |
| `pp_mixture(ds)` share of class 1 | 0.8767 | 0.8767 | n/a |

The model is right in both arms (`pp_mixture()` agrees with the data's
0.8775). Only the name is wrong, and the name is the one brms users
write. `sigma1` and `sigma2` are right: lane 0.8283 and 1.1209, brms
0.8281 and 1.1654.

`frm_simulate(newparams = list(theta1 = 0.7))` gives a component 1 share
of 0.665, which is `plogis(0.7)`. `theta1 = qlogis(0.7)` gives 0.700. So
brms's name is read as a log ratio. That was already so at base, where
`theta1` and `theta1_Intercept` both gave 0.665. But round 1's
"brms names only" contract now promises the brms meaning.

No other family in `dev/brmsnames-rev2-families.R` (17 families, data seed
63) flags a parameter whose content differs from brms's: gaussian,
skew_normal, student, asym_laplace, exgaussian, shifted_lognormal,
lognormal, Gamma, weibull, inverse.gaussian, negbinomial,
zero_inflated_negbinomial, hurdle_lognormal, Beta, beta_binomial and
von_mises each flag the dpars brms lists under the same class names.
wiener (frmtmb.eam) flags `bs`, `ndt` and `bias`, which are brms's.

Fix: do not flag a mixture's `theta<k>` as natural. Either keep it as a
coefficient under brms's spelling for a written theta formula, or store
brms's simplex columns (`theta1` ... `thetaK`, softmax of the log ratios)
and invert that map in `draws_internal_matrix()`. The second is not
elementwise, so `draws_to_natural()` as written cannot do it. Pin it with
a share far from 0.5, as above, because at a share of 0.66 the log ratio
is 0.68 and the defect hides. My first construction did exactly that.

### MAJOR 1. Group-level names that collide still break every draws accessor

Construction: `dev/brmsnames-rev2-collide.R base|lane` and
`dev/brmsnames-rev2-collide-brms.R`, data seed 52
(`dev/brmsnames-rev2-data.R`), frmtmb.sample seed 3, chains 1, iter 100.
The brms fits: chains 1, iter 60, seed 1.

Round 1 suffixes duplicated COEFFICIENT names (`make.unique()` in
`brms_coef_table()`) and refuses duplicated effects. `brms_r_labels()`
has neither, so two `r_` labels can still be equal:

- **C4**, `y ~ x + (1 | gd)` where `gd` has the levels `lvl 1` and
  `lvl.1`. brms's `rename_re_levels()` maps both to `lvl.1`. brms suffixes
  the second, `r_gd[lvl.1,Intercept]__1`, and `as_draws_df(b)` works
  (30 x 23).
- **C1**, `y ~ x + (1 | gi:hi)` where the levels `1_2:3` and `1:2_3` both
  join to `1_2_3`. brms's `combine_groups()` pastes before it builds the
  factor, so brms fits ONE level there. The B8 fit has the two levels
  `1_2_2_3` and `1_2_3`. frmtmb has three levels and two identical labels.

| call | base | lane |
| --- | --- | --- |
| C4 `as_draws_df(ds)` | 50 x 21 | error "Duplicate variable names are not allowed" |
| C1 `as_draws_df(ds)`, `summary(ds)`, `posterior_summary(ds)`, `fixef(ds)`, `coef(ds)` | answer (for example, `summary` 4 x 7) | the same error, 5 of 5 |
| C1 `ranef(ds)` level names | `1:2_3`, `1_2:2_3`, `1_2:3` | `1_2_3`, `1_2_2_3`, `1_2_3` |
| C1 `hypothesis(ds, "r_gi:hi[1_2_3,Intercept] = 0", class = NULL)` | parse error | **-0.00898**, the first column; the second has mean -0.2507 |

This is a regression on draws: base answers all of these. The hypothesis
row is silent, because the draws method reads `x$draws[i, from_cols]` by
name, and R returns the first match. `hyp_env_vals()`'s duplicate stop
does not cover the draws columns. Interaction groups built from IDs that
contain `_` produce this class of collision.

Fix: apply brms's `make.unique(sep = "__")` to the full draws label
vector, as `repair_stanfit()` does, which settles C4 exactly as brms
does. For C1, brms fits a different model, so either refuse the merged
levels by name or suffix them and say so. Also make the draws
hypothesis refuse a name that matches two columns.

### MINOR 1. mv responses that rename to one name: brms refuses, the lane suffixes

`mvbf(bf(y_a ~ x), bf(ya ~ x))`, same scripts. brms: "Cannot use the same
response variable twice in the same model." Lane: `b_ya_Intercept`,
`b_ya_x`, `b_ya_Intercept__1`, `b_ya_x__1`, `sigma_ya`, `sigma_ya__1`. So
the claim "refuse where brms refuses" does not hold for this class. The
output is usable, which is why this is not ranked higher. C2 (levels
`a b` and `ab`) and C3 (`I(x^2)` beside `IxE2`) are refused with brms's
own message in both brms and the lane. At base, both fit.

### MINOR 2. Two names still differ from brms's fitted names

`dev/brmsnames-rev2-names.R lane` and `dev/brmsnames-rev2-b2lane.R`
against the brms fits of `dev/brmsnames-rev2-brms.R`:

- `s(x, by = f2)`: brms `bs_sx:f2u_1` and `sds_sxf2u_1`; lane `bs_sx:f2u_1`
  and **`sds_sx:f2u_1`**. brms drops the `:` from the `sds_` name only.
- `mo(xo)`: brms `bsp_moxo`; lane **`b_moxo`**. I did not check whether the
  two parameterizations hold the same quantity.

### NIT

- The round trip through a logit or `logm1` link is not bitwise exact.
  Readers on `zip`, `hurdle` and `student` differ from base by at most
  1.5e-14 relative (`dev/brmsnames-rev2-log/natural-compare.txt`). This has
  no numerical consequence, but `identical()` against a base result no
  longer holds for those families.
- `bn_exact()` and `bn_fd()` in `test-brms-names.R` scale by
  `max(abs(ref))` over the whole vector, not element by element, so a
  small element is compared at the scale of the largest one.
- Instrument note, pre-existing: `frm_sample(seed =)` does not fix the
  draws when priors are given, because the mode-anchored inits jitter with
  R's RNG. My first natural-scale harness reported every column different
  on `gauss_prior` for that reason. With `set.seed()` per model, the
  columns are identical (`dev/brmsnames-rev2-prior.R` shows the same log
  density in both arms at 7 points).
- Pre-existing and identical at base: on `laplace = TRUE` draws,
  `predictive_interval()` fails on "missing values and NaN's" and
  `pp_check()` fails on "NAs not allowed in predictions".

---

## Round-0 findings: closed or not

### BLOCKER 1 (`x:fe` read as `:`): closed

How I tried to reopen it: `dev/brmsnames-rev2-names.R lane` (data seed
52, sampler seed 3) on `y ~ x * z * f + I(x^2) + poly(z, 2, raw = TRUE) +
(1 + x | g)`, with factor levels `a b`, `c-d` and `e:f` (a level
containing `:`), against hand values from `fixef()`,
`varcorr_matrices()` and the draws matrix. The hand values use no parser.

| hypothesis | fit: hand, hypothesis | draws: hand, hypothesis |
| --- | --- | --- |
| `x:z:fe:f > 0` (three-way, level with `:`) | -0.078059, -0.078059 | -0.084547, -0.084547 |
| `x:fcMd - IxE2 = 0` (`-` and `I()` renamed) | 0.141394, 0.141394 | 0.145006, 0.145006 |
| `polyz2rawEQTRUE2 = 0` | 0.001516, 0.001516 | -0.001672, -0.001672 |
| `fe:f + x:fe:f = 0` | 0.931797, 0.931797 | 0.895331, 0.895331 |
| `b_x:z - sd_g__x > 0`, `class = NULL` (b_ with sd_) | 0.021741, 0.021741 | -0.023352, -0.023352 |

Four hypotheses in one call on the fit: all four equal. `scope = "coef",
group = "g"` with `Intercept + x:fe:f > 0`: 12 levels, max |hand -
estimate| 0. `r_g[lvl.1,x] - r_g[lvl.2,x] = 0`: equal.
`dev/brmsnames-rev2-hyp.R`: on `frm_multiple()` over two imputations,
`x:fe:f`, `x.1 - 0.5` (a dot in a name, a decimal) and `2 * x + .5 *
x:fcMd` equal the mean of the two fits to at most 5.6e-17. Two nonlinear
rows, `exp(x:fcMd)` and `fe:f^2`, differ by about 1e-4, which is Rubin
pooling of a nonlinear quantity and not the parser. `sd_g:h__Intercept`
with `class = NULL`, and `class = "sd", group = "g:h"`, both give the
hand value 0.588595. `r_g:h[lvl.1_p,Intercept] - r_g:h[lvl.1_q,Intercept]`
on draws equals the hand value, and `scope = "ranef", group = "g:h"` gives
24 rows equal to the hand value. The brms fits confirm the names these
hypotheses use (`b_x:z:fe:f`, `b_fcMd`, `b_IxE2`, `b_polyz2rawEQTRUE2`,
`sd_g__x`). The only failure found on a brms name is the mixture's
`theta1` (BLOCKER 1 above), which is a content defect, not a parser defect.

### MAJOR 1 (one renamer): closed, with MINOR 2 left

How I tried to reopen it: nine hostile models fitted in brms
(`dev/brmsnames-rev2-brms.R`, logs `brms-a.txt` to `brms-d.txt`), with
names compared to `variables(fit)` and `brms_par_labels()` in the lane.
The comparison ignores names frmtmb keeps internal by documented choice
(`theta_`, `b[i]`, thresholds, simplexes).

- B1 (three-way interaction, levels with space, hyphen and colon, `I()`,
  raw `poly`, `(1 + x | g)` with levels `lvl 1`): every lane name is in
  brms. brms has one more, `b_polyz2rawEQTRUE1`. frmtmb drops that column
  as aliased with `z`, which is my construction's collinearity and a
  pre-existing difference, not a name.
- `ns(x, df = 3)`, `bs(z, df = 4)`: `b_nsxdfEQ31` ... `b_bszdfEQ44`, the
  same as brms. `gr(g, by =)` is not supported by frmtmb, so it was
  dropped from the lane arm.
- `(1 | mm(g1, g2))`: `sd_mmg1g2__Intercept`, `r_mmg1g2[1,Intercept]`,
  the same as brms.
- Nonlinear `a1`, `b2`: `b_a1_Intercept`, `sd_g__a1_Intercept`,
  `r_g__a1[lvl.1,Intercept]`, the same as brms.
- mv response `y.a_b` with `sigma ~ z`, beside `y2`: `b_yab_x`,
  `b_sigma_yab_z`, `sd_g__yab_Intercept`,
  `cor_g__yab_Intercept__y2_Intercept`, `sigma_y2`, the same as brms.
- `t2(x, z)`: `bs_t2xz_1..3` and `sds_t2xz_1..3`, the same as brms.
  `s(x, by = f2)`: MINOR 2. `gp(z)`: brms `sdgp_gpz` and `lscale_gpz`;
  frmtmb keeps `theta_`, which is documented.
- `(1 | gi:hi)`: `sd_gi:hi__Intercept`, the same as brms. The labels
  collide: new MAJOR 1.
- `cs()` under `acat`: `tau_raw_*` and `bcs2_*`, internal, documented as
  not fixed (thresholds).
- `bf(y ~ x + (1 + sigma_Intercept | g), sigma ~ (1 | g))`: brms refuses
  ("Variable names starting with 'mu_', 'sigma_' ... are not allowed");
  the lane refuses with "Duplicated group-level effects". Both refuse.

### MAJOR 2 (collisions): NOT closed

The cases round 0 named hold: C2 and C3 are refused with brms's message,
and the round-1 pins cover K1 to K5. But `r_` labels are neither suffixed
nor refused, and on draws that reproduces the failure round 0 reported
(new MAJOR 1). Two classes are reachable, one where brms suffixes (C4)
and one where brms merges (C1). An mv response class is suffixed where
brms refuses (MINOR 1). The known exact-twin `(1|g)+(1|g)` refusal was
not ranked, as the brief says.

### MAJOR 3 (natural-scale dpars): the inversion is closed; the naming is not (BLOCKER 1)

This is where most of the effort went. How I tried to reopen it:
`dev/brmsnames-rev2-natural.R base|lane|compare`, data seed 41, R seed per
model, sampler seed 9, chains 2, iter 200. The data are chosen so that the
link and natural values differ widely: sigma 4 (log 1.4), negbinomial
shape 3, zero-inflated Poisson zi 0.3 (logit -0.85), student nu 4, Beta
phi 12, hurdle_gamma shape and hu, mv `sigma_ya` and `sigma_y2`, plus
`gauss` with a `normal(0, 1)` prior on b, with `laplace = TRUE`, and with
`reparameterize = TRUE`. The same stanfit is read in both arms.

Every stored draws column is `identical()` between arms, except the
natural columns. Those equal `exp()` or `plogis()` of base to 0 relative
difference.

Readers compared: `posterior_epred`, `posterior_linpred`, both again
with `dpar =` the natural dpar, `posterior_predict`, `log_lik`, `loo`,
`waic`, `bayes_R2`, `predictive_error`, `predictive_interval`,
`pp_check` (`dens_overlay` and `stat` data), `conditional_effects` and
`conditional_effects(dpar =)` with `robust = FALSE`, and `check_laplace`.
The result was 142 comparisons over 10 models: 125 `identical()`, and 17
differing by at most 1.5e-14 relative (the NIT). The default
`conditional_effects()` differs by up to 1.05 relative on every model,
and with `robust = FALSE` it is identical. So that difference is the
documented `robust = TRUE` default, not the natural scale.

Also identical between arms: `hypothesis(ds, "sigma > 3")` 1.423045;
`sigma - sd_g__Intercept > 0` 3.850048; the point evidence ratio with a
prior on b, `x = 0` 1.359545 and `x - 0.3 = 0` 2.611525
(`dev/brmsnames-rev2-er.R`). `VarCorr(ds)$residual__` Estimate 4.423045
is the mean of the natural `sigma` column.

`dev/brmsnames-rev2-wiener.R base|lane` (frmtmb.eam `wiener()`, whose
`ndt` link is a closure-bound scaled logit, data seed 77): `sum(log_lik)`
-12725.392684142722, `sum(posterior_epred)`, `sum(posterior_predict)` and
`loo()` agree to every printed digit between arms.

No extension outside frmtmb.sample reads a draws matrix: a grep for
`frmtmb_draws`, `$draws` and `draws_fit_at` outside frmtmb.sample finds
only documentation and refusals. Inside frmtmb.sample, every reader that
hands a draw to the model goes through `draws_fit_at()` or
`draws_internal_matrix()`. The two readers of the raw matrix are
`hypothesis()` (which should see natural values, as brms's does) and
`ranef()` (`r_` columns only).

The reader side is closed. What is open is which parameters get flagged
natural: BLOCKER 1.

### MAJOR 4 (`ranef(ds)`/`coef(ds)` on excluded blocks): closed

How I tried to reopen it: `dev/brmsnames-rev2-ranef.R` (data seed 71,
sampler seed 4, chains 1, iter 80) on 17 structures: `us`, `us` with
`reparameterize = TRUE`, `||`, `rr(d = 1)`, `ar1`, `ou`, `cs`, `homcs`,
`hetar1`, `homdiag`, `toep`, spatial `exp`, `gr(cov =)`, `gr(prec =)`,
`equalto`, two blocks on two factors, and `s(x1) + (1 | g)`. At draws 1,
20 and 40, each was compared against `ranef()` and `coef()` of the fit
set to that draw. The result: 0 NA cells in `ranef(ds)` and `coef(ds)` on
all 17, and a worst relative difference of 0 for the value sets.
`dev/brmsnames-rev2-ranef-pos.R` pairs level by coefficient positionally
at draw 20 on 10 of them, and the worst absolute difference is 0. The
smooth model's ML `ranef()` has 3 more tables than `ranef(ds)`, which
leaves the smooth out, as brms does. I did not rerun the laplace refusal;
the worker's block 10 has it.

### MINORs and the user decision

- Pins not routed through `brms_par_labels()`, and the Est.Error pin: not
  re-measured. I did not rerun the 22 mutants. Neither new defect has a
  pin (nothing in the suites builds a mixture name or an `r_` collision).
- Relative tolerances: `test-brms-names.R` and `test-brms-pins.R` carry
  no absolute `tolerance =`. See the NIT on `bn_exact()`.
- mv `residual__` and `bayes_R2`: the lane's `VarCorr(ds)` on the mv model
  has `residual__` rows for `sigma_ya` and `sigma_y2`. `bayes_R2(ds, resp
  = "ya")` is identical to base.
- `frm_simulate(newparams =)`, `dev/brmsnames-rev2-simulate.R` (data seed
  81): every old or bare name I tried is refused loudly, with its brms
  spelling where one exists. Those were `b_x` when only a covariate
  `b_x` exists (to `b_b_x`), `sigma_Intercept`, `x`, `shape_Intercept`,
  `a_Intercept`, `sd_g__(Intercept)`, and `sigma` when `sigma ~ 1` is
  written. With both `x` and `b_x` present, `b_x` sets x's coefficient
  (correlation with x 0.947), which is brms's reading. The one silent
  misread is the mixture's `theta1` (BLOCKER 1).

### Regression

- `dev/brmsnames-rev2-silent.R base|lane|compare` (round 0's script,
  data seed 31, sampler seed 5, own cache files): `fixef`, `fixef(flatten
  = TRUE)`, `ranef`, `coef`, the covariance matrices and `confint` are
  `identical()` on all 6 models. Draws columns identical: dist 39 of 39,
  mv 28 of 30, nl 14 of 15, ord 16 of 16, smooth 23 of 24, gp 204 of 205.
  Every non-identical column is an unmodeled sigma, equal to `exp()` of
  base with a relative difference of 0. The `r_` columns match `ranef()` at
  draws 1, 30 and 60: 180 checked, 0 wrong.
- `dev/brmsnames-rev2-interop.R base|lane|compare`, data seed 3:
  `insight` (`find_parameters`, `get_parameters`, `get_variance`,
  `get_statistic`, `find_random`), `emmeans`, `print` and `print(summary)`
  are `identical()`. `marginaleffects::avg_slopes` differs only in the
  attached model object, as in round 0.
- Suites: not rerun. The counts in the brief (core 10458, sample 1672,
  learn 460, five others green) are the verifier's and were not
  re-derived here.

## What to fix

1. BLOCKER 1: exclude a mixture's `theta<k>` from the natural flag, or
   store brms's simplex with a non-elementwise inverse. Pin with a mixing
   share far from 0.5, against brms's `theta1` meaning, on the fit, on
   draws and in `frm_simulate(newparams =)`.
2. MAJOR 1: `make.unique(sep = "__")` over the full draws labels, as
   brms's `repair_stanfit()` does. Refuse or suffix interaction-group
   levels that join to one string. Have the draws `hypothesis()` refuse a
   name that matches two columns. Pin C4 (`lvl 1` beside `lvl.1`) and C1.
3. MINOR 1 and 2 in the same pass, or record them under "Found and NOT
   fixed".
