# Lane famsimplex: brms's simplex-response families

Lane famsimplex, round of 2026-10-13, base d2ef6312 (frmtmb 0.70.0,
frmtmb.sample 0.18.0). Worktree `frmtmb-wt-famsimplex`, branch
`wt-famsimplex`, library `C:/Users/adf44/source/r/wt-famsimplex-lib`
ahead of `rellib-r8`. Nothing is committed.

`frm()` now fits brms 2.23.0's `dirichlet()` and `logistic_normal()`,
and brms's unexported `dirichlet2` through `brmsfamily("dirichlet2")`,
the spelling brms itself takes. Every route the brief names reaches
them. The numbers below are generated from the logs by
`dev/famsimplex-summarize.R` (blocks marked "generated"; it writes
`dev/famsimplex-summary.md` and `-summary2.md`, which were pasted here
and removed). Every script named is in `dev/`, and its log beside it
is local (gitignored).

## 1. What brms 2.23.0 does, as run

Read off brms itself, not recalled: `dev/famsimplex-brms-probe.R`
(stancode, standata, default_prior, validation) and
`dev/famsimplex-brms-fits2.R` (five fits, 4 chains x 1000 draws, and
every post-fit method on each). brms's source is the tarball tree in
`dev/brms-suite/brms`.

- `dirichlet(link = "logit", link_phi = "log", refcat = NULL)`: dpars
  `mu<cat>` (one per non-reference category) and `phi`; the only mean
  link is `logit`; `phi` takes log, identity, softplus, squareplus.
  Stan: `dirichlet_lpdf(y | softmax(mu) * phi)` with the reference
  category's predictor 0. `refcat = NA` predicts every category.
- `logistic_normal(link = "identity", link_sigma = "log", refcat =
  NULL)`: dpars `mu<cat>` and `sigma<cat>`, one per non-reference
  category, and a correlation matrix declared as `Llncor`
  (`cholesky_factor_corr[ncat-1]`), reported as `lncor__b__c` (upper
  triangle, column by column) and printed `lncor(b,c)` under Further
  Distributional Parameters. `lncor ~ x` is refused ("not a valid
  distributional or non-linear parameter"). `refcat = NA` is refused.
  The density is `multi_normal_cholesky(log(y_k / y_ref) | mu, diag(sigma)
  Lcor) - sum(log(y))`.
- `dirichlet2` exists (`.family_dirichlet2`, `brms:::dirichlet2`) but
  is not exported; `brmsfamily("dirichlet2")` builds it. Links log,
  softplus, squareplus, identity, logm1; one predictor per category
  (all K, no reference), `y ~ Dirichlet(alpha)`.
- Response: a matrix (`cbind(y1, y2, y3)` or a matrix column);
  categories are the column names (`1..K` when unnamed); dpar names
  are `make_stan_names(paste0("mu", cat))` (`cbind(y1, y2, y3)` gives
  `muy2`, `muy3`). Fewer than two columns: "At least 2 response
  categories are required." Row sums: `is_equal(rowSums(Y), 1)`,
  which is `all.equal()`, the MEAN absolute departure against 1.5e-8.
  Bounds: "Family 'dirichlet' requires response greater than 0." Rows
  with NA are dropped with a warning. Addition terms: weights, subset,
  index.
- Default priors: `student_t(3, 0, 2.5)` on each `Intercept_mu<cat>`,
  flat `b`; `gamma(0.01, 0.01)` on `phi`; `student_t(3, 0, 2.5)` on
  `sigma<cat>` (natural, lb 0) or on its intercept when modeled;
  `lkj(1)` on class `lncor`.
- Post-fit: `fitted()` is `n x 4 x K`, layers `P(Y = a)`;
  `posterior_epred()` `S x N x K`, names `a, b, c`; `predict()`
  `n x 4 x K` named `a, b, c`; `posterior_predict()` `S x N x K` with
  rows on the simplex; `log_lik()` `S x N`; `posterior_linpred()`
  without `dpar` `S x N x (K-1)`; `residuals()` and
  `predictive_error()` refused ("Predictive errors are not defined for
  ordinal or categorical models."); `conditional_effects()` needs
  `categorical = TRUE` ("Please set 'categorical' to TRUE.") and then
  draws a `cats__` column; `pp_check()` refused ("'pp_check' is not
  implemented for this family."); `bayes_R2()` refused. For
  `logistic_normal()`, `fitted()`, `posterior_epred()` and
  `conditional_effects()` are refused: "Cannot compute expected values
  of the posterior predictive distribution for family
  'logistic_normal'."
- A category's own formula: brms takes `bf(Y ~ x, mub ~ z)` for the
  FIRST predicted category of a `dirichlet()` response, and
  `bf(y ~ 1, mu1 ~ x1, mu2 ~ x1)` for a mixture
  (`dev/famsimplex-brms-first.R`).

## 2. What changed

### Core (`R/families-simplex.R`, new; wired in elsewhere)

- `dirichlet()`, `logistic_normal()` (exported, documented in
  `?frmtmb-families`, section "Compositional responses") and
  `brmsfamily("dirichlet2")`. Each constructor returns a DEFERRED
  family, as `categorical()` does: `resolve_deferred_families()` reads
  the categories off the response's columns and builds the family.
  The categorical machinery is reused rather than duplicated: the
  deferred-family seam, `primary_dpars` for the per-category
  predictors, the n x K category layer `fitted()` already builds from
  a matrix-valued mean (`brms_category_labels()`), and the matrix
  response path `simulate()` and `predict()` already had for
  `multinomial()`. They are `type = "continuous"` families with a
  matrix response and a `simplex` field (categories, reference,
  predicted categories, dpar names), because `type = "categorical"`
  means a 1..K index response throughout the package.
- Densities: brms's, written over `lgamma()` and RTMB's
  `logspace_add()`, so the softmax denominator accumulates in log
  space (as `cat_logit_lpdf()` does). `logistic_normal()` holds its
  correlations as the extra parameter `lncor_raw`, in the
  row-normalized Cholesky chart `set_rescor(TRUE)` uses
  (`us_chol_cor()`), so the LKJ prior and its Jacobian carry over;
  `RTMB::dmvnorm()` on the standardized log ratios, minus
  `sum(log(sigma))` and `sum(log(y))`.
- Response check (`simplex_valid_y()`): brms's order and sentences,
  but the row sums PER ROW at brms's own 1.5e-8 (section 6, brms-A).
- `fitted()`: the mean of `dirichlet()` and `dirichlet2` is a
  `post$mean_fn` returning the n x K category means, so `fitted()`
  gives brms's array with delta-method `Est.Error`. `logistic_normal()`
  declares a `mean_fn` that refuses in brms's words.
- `summary()`: `lncor(b,c)` rows under Further Distributional
  Parameters (`lncor_summary_rows()`), delta-method error, Wald
  interval on `atanh(r)`. `hypothesis()` and `variables()` read the
  correlations as `lncor__b__c` (`hyp_vals_only()`, `hyp_par_cov()`,
  `hyp_env_vals()` in `R/confint.R`).
- Links line: `mub = logit` on a dirichlet fit, as brms prints it
  (`family_link_str()`); `$link` is `"logit"` for `dirichlet()`,
  `"identity"` for `logistic_normal()` (`family_link_view()`, a
  `joint_link` field).
- Priors (`R/priors.R`): new class `"lncor"`, `lkj()` only, resolved
  onto `lncor_raw` like `"rescor"`; a default-prior table row.
- `residuals()` refused for every type, `pp_check()` with brms's
  sentence; `conditional_effects()` takes the per-category display by
  default and its Wald band from `fitted()`.
- `fam_is_polytomous()` includes the simplex families (brms's
  `is_polytomous()` does); `brms_no_mixture` gains them (regenerated
  `R/links-brms.R` from brms with `dev/famlink-gen-links.R`); the
  compatibility table (`R/compat.R`) has a `simplex` group and rows.
- `as_frmtmb_family()` carries a brms family object's `refcat`.
- **Parser fix (`R/parse.R`)**: a family with several location
  predictors now takes a formula for the FIRST one too. It was refused
  with "dpar(s) not available" on every such family:
  `categorical()`'s first predicted level and `mixture()`'s `mu1`
  included. This reverses one row of `dev/vigport-findings.md`'s
  divergence table ("brms's `bf(y ~ 1, mu1 ~ x1, mu2 ~ x1)` for a
  mixture: STILL FAILS") on the lane build; that row is a vignette
  port row, not a brms-suite row. Fix round 1 (section 9): the widening
  is for formulas only, and a constant on ANY location of such a
  family is now refused.

### frmtmb.sample

- `draws_lncor_cols()`: the draws column is brms's `lncor__b__c`, read
  back to the chart for every model evaluation (same mechanism as the
  ordinal blocks).
- Default priors: `lkj(1)` on class `lncor`, and brms's `sigma` default
  on each `sigma<k>` (scale 2.5, because the response is a composition
  and brms does not transform it). `phi`'s `gamma(0.01, 0.01)` is
  named by the existing disclosure ("no defaults for phi").
- `pp_check()` and `predictive_error()` refuse the simplex families in
  brms's words.

## 3. Evidence

### Log densities against brms's Stan functions

brms's own chunk files (`fun_dirichlet_logit.stan`,
`fun_logistic_normal.stan`, and Stan's `dirichlet_lpdf` for
dirichlet2) compiled through `rstan::expose_stan_functions()`; K = 2
to 5, every reference position and `refcat = NA`, 200 rows each, eta
up to +-30, phi from 0.05 to 1e4, sigma 0.05 to 5, random correlation
matrices (`dev/famsimplex-stan-lpdf.R`).

<!-- BEGIN generated: Stan lpdf grid (dev/famsimplex-stan-lpdf.R, seed 20261013) -->
```
WORST per family:
           family  max_abs  max_rel
1       dirichlet 1.92e-10 1.53e-11
2      dirichlet2 3.64e-12 6.55e-13
3 logistic_normal 4.73e-11 1.42e-14
rows compared: 7200  any non-finite: FALSE 
```
<!-- END generated -->

### brms's Stan program at frmtmb's estimates

`tests/testthat/test-simplex-brms.R`, gated, through the existing
`brms_lp_check()` harness: brms's generated program with flat priors,
`log_prob()` at frmtmb's estimates (check A), the gradient there
(check B), and for group effects the joint density (check C). The one
translation rule added is `Llncor` in `tests/testthat/helper-brms.R`.

<!-- BEGIN generated: brms_lp_check (dev/famsimplex-lpcheck.R) -->
```
row 16h: dirichlet with phi ~ z is brms's likelihood | LPCHECK const -1.705302566e-13 grad 9.71e-05 ours 476.6628169
row 16i: dirichlet with refcat and a category's own formula | LPCHECK const -3.410605132e-13 grad 3.14e-05 ours 1021.388378
row 16j: dirichlet with (1 | ID | g), on the joint density | LPCHECK const 5.684341886e-13 grad 9.41e-10 ours 656.9804067
row 16k: logistic_normal with sigma<k> ~ z and three lncor | LPCHECK const -4.547473509e-13 grad 0.000207 ours 893.4846694
row 16l: logistic_normal with refcat and (1 | ID | g) | LPCHECK const -3.97903932e-13 grad 8.22e-15 ours 478.7711228
row 16m: brmsfamily("dirichlet2") is brms's likelihood | LPCHECK const 1.136868377e-13 grad 0.00012 ours 446.8271215
rows 6; max |const| 5.68e-13; max |grad| 0.000207
```
<!-- END generated -->

The constant is an identity (0 up to rounding of a sum of 300 to 400
terms). The gradient bound of the harness is 1e-3; row 16k's 2.07e-4
is the largest.

### Fit agreement with brms fits on the same data

(posterior median - ML estimate) / posterior sd, per brms variable
(`dev/famsimplex-agree.R`, reading the brms fits of
`dev/famsimplex-brms-fits2.R`, brms seed 1; data seeds 11 to 15 of
`dev/famsimplex-sim.R`). Models: dir `bf(Y ~ x, phi ~ z)`; dirre
`bf(Y ~ x + (1 | ID | g))`, 20 groups; ln `bf(Y ~ x, sigmab ~ z)`;
lnre `bf(Y ~ x + (1 | ID | g))`; dir2 `bf(Y ~ x)`.

<!-- BEGIN generated: ML vs brms posterior (dev/famsimplex-agree.R) -->
```
SUMMARY |z| per model
  model abs(z).max abs(z).mean
1   dir     0.1283      0.0537
2  dir2     0.1535      0.0837
3 dirre     0.3568      0.1032
4    ln     0.1796      0.0610
5  lnre     0.4652      0.1322
variables compared: 38  max |z|: 0.465 
```
<!-- END generated -->

The four largest are the group-level sds (dirre 0.357 and 0.247, lnre
0.465 and 0.406): ML is below the posterior median, which is the
ordinary ML shrinkage of a variance component at 20 groups
(`sqrt(1 - 1/20) = 0.975`) plus the skew of the sd posterior, not a
disagreement in the model.

### fitted() and predict() against brms's on the same fits

Per cell of the n x K summary, 5th / 50th / 95th percentiles over
cells (`dev/famsimplex-postfit-agree.R`; frmtmb's `predict()` draws
4000 replicates after `set.seed(1)`, brms's uses its 4000 draws). The
fitted rows compare the ML mean and its delta-method `Est.Error`
with the posterior mean and sd of `posterior_epred()`; the predict
rows compare the predictive means and sds. `logistic_normal()` has no
fitted rows because both packages refuse its mean.

<!-- BEGIN generated: fitted() and predict() against brms (dev/famsimplex-postfit-agree.R) -->
```

== dir 
fitted: (brms mean - ML) / brms sd, 5/50/95%: -0.044 0.011 0.040 
fitted: ML se / brms sd, 5/50/95%: 0.976 0.986 1.001 
predict: (brms mean - ML mean) / brms sd, 5/50/95%: -0.039 -0.000 0.037 
predict: ML sd / brms sd, 5/50/95%: 0.969 0.997 1.024 
dims identical: TRUE  layer names identical: TRUE 

== dirre 
fitted: (brms mean - ML) / brms sd, 5/50/95%: -0.074 0.001 0.074 
fitted: ML se / brms sd, 5/50/95%: 0.959 0.988 1.019 
predict: (brms mean - ML mean) / brms sd, 5/50/95%: -0.040 0.001 0.038 
predict: ML sd / brms sd, 5/50/95%: 0.972 0.999 1.027 
dims identical: TRUE  layer names identical: TRUE 

== ln 
predict: (brms mean - ML mean) / brms sd, 5/50/95%: -0.037 0.000 0.038 
predict: ML sd / brms sd, 5/50/95%: 0.964 0.992 1.021 
dims identical: TRUE  layer names identical: TRUE 

== lnre 
predict: (brms mean - ML mean) / brms sd, 5/50/95%: -0.040 0.000 0.040 
predict: ML sd / brms sd, 5/50/95%: 0.966 0.996 1.024 
dims identical: TRUE  layer names identical: TRUE 

== dir2 
fitted: (brms mean - ML) / brms sd, 5/50/95%: -0.054 0.004 0.042 
fitted: ML se / brms sd, 5/50/95%: 0.973 0.989 1.026 
predict: (brms mean - ML mean) / brms sd, 5/50/95%: -0.037 0.000 0.034 
predict: ML sd / brms sd, 5/50/95%: 0.970 0.997 1.024 
dims identical: TRUE  layer names identical: TRUE 
```
<!-- END generated -->

`as_tmbstan()` and `frm_sample()` on both families, `loo()`,
`conditional_effects()` on draws, `hypothesis()`, `mvbf()` with a
`dirichlet()` response and `emmeans(dpar = "mub")` were run once in
`dev/famsimplex-try2.R` and `dev/famsimplex-misc.R` (logs local); the
tests cover the sampler routes in `test-simplex-draws.R`.

### Recovery on simulated data, and the SE check

100 seeds each (1001 to 1100), `dev/famsimplex-recovery.R`. Dirichlet
data by gamma normalization, logistic-normal data by an MVN draw on
the additive log-ratio scale (`dev/famsimplex-sim.R`). "fits with a
warning" counts every warning `frm()` raised, so it is the SE check's
false-alarm count (its default is `check_se = "warning"`). The SE check
does fire where it should: `dirichlet(refcat = NA)` leaves a common
shift unidentified and the fit names all six coefficients
(`dev/famsimplex-edge.log`; pinned in `test-simplex-families.R`).

<!-- BEGIN generated: recovery dir (dev/famsimplex-recovery.R dir 1001 100) -->
```
config dir : seeds 1001 to 1100 ; fitted 100 of 100 ; errors 0 
nonzero convergence codes: 0 
fits with a warning: 0 
                truth mean_est      bias    mc_se  sd_est mean_se cover95
b_mub_Intercept   0.4   0.4011  0.001067 0.004691 0.04691 0.04452    0.96
b_mub_x           0.6   0.6078  0.007753 0.005140 0.05140 0.04779    0.92
b_muc_Intercept  -0.3  -0.2911  0.008940 0.005280 0.05280 0.05084    0.92
b_muc_x          -0.5  -0.4949  0.005052 0.004458 0.04458 0.05169    0.99
b_phi_Intercept   2.0   2.0107  0.010673 0.005662 0.05662 0.05464    0.93
b_phi_z           0.5   0.4993 -0.000710 0.005755 0.05755 0.05241    0.92
median seconds per fit: 0.05 
```
<!-- END generated -->

<!-- BEGIN generated: recovery dirre (dev/famsimplex-recovery.R dirre 1001 100) -->
```
config dirre : seeds 1001 to 1100 ; fitted 100 of 100 ; errors 0 
nonzero convergence codes: 0 
fits with a warning: 0 
                                     truth mean_est       bias    mc_se  sd_est
b_mub_Intercept                      0.400   0.4040  0.0040109 0.011798 0.11798
b_mub_x                              0.600   0.6009  0.0008871 0.004218 0.04218
b_muc_Intercept                     -0.300  -0.2986  0.0014341 0.009100 0.09100
b_muc_x                             -0.500  -0.5008 -0.0007918 0.005012 0.05012
phi                                  7.389   7.4031  0.0140826 0.035012 0.35012
sd_g__mub_Intercept                  0.500   0.4738 -0.0262296 0.009089 0.09089
sd_g__muc_Intercept                  0.300   0.2830 -0.0170124 0.007188 0.07188
cor_g__mub_Intercept__muc_Intercept  0.400   0.3654 -0.0346211 0.026905 0.26905
                                    mean_se cover95
b_mub_Intercept                     0.11394    0.94
b_mub_x                             0.04454    0.97
b_muc_Intercept                     0.07925    0.90
b_muc_x                             0.04776    0.93
phi                                 0.36633    0.96
sd_g__mub_Intercept                      NA      NA
sd_g__muc_Intercept                      NA      NA
cor_g__mub_Intercept__muc_Intercept      NA      NA
median seconds per fit: 1.035 
```
<!-- END generated -->

<!-- BEGIN generated: recovery ln (dev/famsimplex-recovery.R ln 1001 100) -->
```
config ln : seeds 1001 to 1100 ; fitted 100 of 100 ; errors 0 
nonzero convergence codes: 0 
fits with a warning: 0 
                     truth mean_est       bias    mc_se  sd_est mean_se cover95
b_mub_Intercept     0.4000   0.4017  0.0016946 0.003118 0.03118 0.02972    0.93
b_mub_x             0.6000   0.6015  0.0015360 0.002528 0.02528 0.02982    0.98
b_muc_Intercept    -0.3000  -0.2936  0.0063724 0.004841 0.04841 0.04703    0.91
b_muc_x            -0.5000  -0.4938  0.0062330 0.004128 0.04128 0.04726    0.97
b_sigmab_Intercept -0.5000  -0.5022 -0.0022041 0.003747 0.03747 0.04089    0.97
b_sigmab_z          0.4000   0.4024  0.0024288 0.003623 0.03623 0.04073    0.97
sigmac              0.8187   0.8194  0.0007168 0.003172 0.03172 0.03346    0.96
lncor__b__c         0.3000   0.3029  0.0028872 0.005232 0.05232 0.05237    0.93
median seconds per fit: 0.03 
```
<!-- END generated -->

<!-- BEGIN generated: recovery lnre (dev/famsimplex-recovery.R lnre 1001 100) -->
```
config lnre : seeds 1001 to 1100 ; fitted 100 of 100 ; errors 0 
nonzero convergence codes: 0 
fits with a warning: 0 
                                      truth mean_est       bias    mc_se
b_mub_Intercept                      0.4000   0.4013  1.283e-03 0.012027
b_mub_x                              0.6000   0.6013  1.336e-03 0.003142
b_muc_Intercept                     -0.3000  -0.2939  6.107e-03 0.008886
b_muc_x                             -0.5000  -0.5012 -1.209e-03 0.004832
sigmab                               0.6065   0.6066  5.054e-05 0.002389
sigmac                               0.8187   0.8170 -1.742e-03 0.002791
lncor__b__c                          0.3000   0.3024  2.431e-03 0.004080
sd_g__mub_Intercept                  0.5000   0.4790 -2.101e-02 0.007984
sd_g__muc_Intercept                  0.3000   0.2979 -2.061e-03 0.006968
cor_g__mub_Intercept__muc_Intercept  0.4000   0.3622 -3.779e-02 0.025772
                                     sd_est mean_se cover95
b_mub_Intercept                     0.12027 0.11145    0.95
b_mub_x                             0.03142 0.03116    0.95
b_muc_Intercept                     0.08886 0.07864    0.91
b_muc_x                             0.04832 0.04178    0.89
sigmab                              0.02389 0.02200    0.93
sigmac                              0.02791 0.02964    0.95
lncor__b__c                         0.04080 0.04652    0.97
sd_g__mub_Intercept                 0.07984      NA      NA
sd_g__muc_Intercept                 0.06968      NA      NA
cor_g__mub_Intercept__muc_Intercept 0.25772      NA      NA
median seconds per fit: 0.13 
```
<!-- END generated -->

<!-- BEGIN generated: recovery dir2 (dev/famsimplex-recovery.R dir2 1001 100) -->
```
config dir2 : seeds 1001 to 1100 ; fitted 100 of 100 ; errors 0 
nonzero convergence codes: 0 
fits with a warning: 0 
                truth mean_est       bias    mc_se  sd_est mean_se cover95
b_mua_Intercept   1.0   1.0097  0.0097305 0.006172 0.06172 0.05964    0.92
b_mua_x           0.3   0.3003  0.0002732 0.006456 0.06456 0.05986    0.96
b_mub_Intercept   1.5   1.5132  0.0132245 0.006025 0.06025 0.06012    0.94
b_mub_x          -0.4  -0.4046 -0.0046380 0.006636 0.06636 0.06067    0.93
b_muc_Intercept   0.8   0.8102  0.0102059 0.005752 0.05752 0.05944    0.95
b_muc_x           0.2   0.1962 -0.0037955 0.006997 0.06997 0.05992    0.90
median seconds per fit: 0.05 
```
<!-- END generated -->

Reading. 500 fits, 0 errors, 0 nonzero convergence codes, 0 warnings.
Every coefficient and natural-scale dpar bias is within 2.2 Monte
Carlo standard errors (largest: dirichlet2 `b_mub_Intercept`,
0.0132 / 0.00603 = 2.19). The three dirichlet2 intercepts are all
positive, +0.010 to +0.013 (1.6 to 2.2 MC SE each): the log
concentrations share one upward finite-sample bias, as the Dirichlet
precision MLE has at n = 300; the slopes, which a common scale does
not move, are within 0.7 MC SE. The group-level sds sit below truth
(0.474 against 0.5 in dirre: -2.9 MC SE from 0.5, -1.5 from the ML
shrinkage target 0.5 x sqrt(1 - 1/20) = 0.487), as ML variance
components do. Wald
coverage at 100 replicates ranges from 0.89 to 0.99 against 0.95;
the binomial standard error at 100 is 0.022, so 0.89 is 2.7 SE low
(lnre `b_muc_x`, whose mean SE 0.0418 is 0.86 of its sd 0.0483: the
Wald error of a slope in a random-effects fit at 20 groups, not a
simplex effect). The sd and correlation rows carry no SE here because
the summary's random-effect table is not read by this script.

### Simulators against their own densities

`tests/testthat/test-simulate-density.R` now has the three families in
its agreement tier: on K = 2 a chi-square over 8 cells and the first
two moments of the second proportion, on K = 3 the mean and variance
of each proportion from a 400 x 400 triangle grid of `fam$lpdf`, and
the logistic-normal log ratios' correlation against the chart's 0.4,
4000 draws per cell at the file's alpha 1e-6. 56 statistical
assertions were added (the file header's count is updated to 307).

### Behavior the tests pin, seen to fail on rellib-r8

Most pins fail on rellib-r8 in the weak form (no `dirichlet()`). The
behavioral form, through routes that need no new symbol:

<!-- BEGIN generated: pins, before (dev/famsimplex-pins.R) -->
```
lib: C:/Users/adf44/source/r/rellib-r8/frmtmb 
family = brms::dirichlet() : ERROR: dirichlet is not a supported family. Supported families are: gaussian, poisson, binomial, Gamma, lognormal, studen
family = brms::logistic_normal() : ERROR: logistic_normal is not a supported family. Supported families are: gaussian, poisson, binomial, Gamma, logno
family = "dirichlet" : ERROR: dirichlet is not a supported family. Supported families are: gaussian, poisson, binomial, Gamma, lognormal, student, neg
brmsfamily("dirichlet2") : ERROR: dirichlet2 is not a supported family. Supported families are: gaussian, poisson, binomial, Gamma, lognormal, student
```
<!-- END generated -->

<!-- BEGIN generated: pins, after (dev/famsimplex-pins.R) -->
```
lib: C:/Users/adf44/source/r/wt-famsimplex-lib/frmtmb 
family = brms::dirichlet() : fits; logLik 456.8262864 
family = brms::logistic_normal() : fits; logLik 276.0969825 
family = "dirichlet" : fits; logLik 456.8262864 
brmsfamily("dirichlet2") : fits; logLik 452.5082448 
```
<!-- END generated -->

The parser fix, behaviorally, on an existing family:
`test-simplex-families.R`'s "categorical() takes a formula for its
first predicted category" fails on rellib-r8 with "dpar(s) not
available for family 'categorical': mub (available: muc)"
(`dev/famsimplex-tests/before2/`), and the vigport mixture spelling:

<!-- BEGIN generated: mixture mu1 ~ x (dev/famsimplex-mixmu1.R) -->
```
lib: C:/Users/adf44/source/r/rellib-r8/frmtmb 
[1] "ERROR: dpar(s) not available for family 'mixture(gaussian, gaussian)': mu1 (available: sigma1, mu2, sigma2, theta1)"
lib: C:/Users/adf44/source/r/wt-famsimplex-lib/frmtmb 
logLik_mu1_mu2    logLik_main 
     -788.8647      -788.8647 
```
<!-- END generated -->

(The two models are the same model, so the equal logLik is an
identity check of the new route, not a measurement.)

Test runs of the new files, before and after:

<!-- BEGIN generated: test run before1 (dev/famsimplex-run.sh) -->
```
RESULT test-simplex-families.R pass=0 fail=0 err=17 skip=0 warn=0 | core C:/Users/adf44/source/r/rellib-r8
RESULT test-simplex-brms.R pass=0 fail=0 err=6 skip=0 warn=0 | core C:/Users/adf44/source/r/rellib-r8
RESULT test-simulate-density.R pass=492 fail=1 err=1 skip=0 warn=0 | core C:/Users/adf44/source/r/rellib-r8
RESULT test-simplex-draws.R pass=0 fail=0 err=4 skip=0 warn=0 | core C:/Users/adf44/source/r/rellib-r8
RAN 4 of 4
```
<!-- END generated -->

<!-- BEGIN generated: test run before2 (dev/famsimplex-run.sh) -->
```
RESULT test-simplex-families.R pass=0 fail=0 err=19 skip=0 warn=0 | core C:/Users/adf44/source/r/rellib-r8
RAN 1 of 1
```
<!-- END generated -->

<!-- BEGIN generated: test run new3 (dev/famsimplex-run.sh) -->
```
RESULT test-simplex-families.R pass=97 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RESULT test-simplex-brms.R pass=15 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RESULT test-simplex-draws.R pass=25 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RAN 3 of 3
```
<!-- END generated -->

<!-- BEGIN generated: test run new4 (dev/famsimplex-run.sh) -->
```
RESULT test-simplex-families.R pass=100 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RAN 1 of 1
```
<!-- END generated -->

## 4. Ported brms-suite bin-1 rows

None flip. No bin-1 block of brms's suite fits these families: their
rows are in `tests.stancode.R`, `tests.log_lik.R`,
`tests.posterior_epred.R` and `tests.posterior_predict.R`, which test
brms's Stan code and its internal `prep` functions, outside bin 1;
`grep -i "dirichlet|logistic_normal|simplex"` over
`dev/brmsport-verdicts.tsv` and `dev/brmsport-ledger.tsv` matches
nothing. The gated run of the generated `test-brms-suite-*.R` files
is in section 7.

## 5. What I did not do, and why

- `logistic_normal()` in a multivariate model (`mvbf()`): core refuses
  family extras outside the ordinal families in a multivariate frame
  (`check_mv_extra_family()`), so `lncor` there is refused by name.
  `dirichlet()` in `mvbf()` is not tested.
- `posterior_linpred()` without `dpar` stays one predictor's
  `S x N` matrix; brms returns `S x N x (K - 1)`. This is the existing
  documented convention for every multi-predictor family in
  frmtmb.sample (`?posterior_epred`, "Categorical outcomes"), not a
  simplex decision, and changing it is a change for categorical and
  multinomial too.
- `categorical()` has no `refcat`; brms's `categorical(refcat = )`
  passed as a brms object loses it silently, because
  `as_frmtmb_family()` carries `refcat` only to a constructor that
  takes it. Not this lane's family.
- `multinomial()`'s `fitted()` still refuses ("declares no mean"), where
  brms returns an `n x 4 x K` array; seen while reading the
  categorical machinery (`dev/famsimplex-explore-multi.log`), not
  fixed.
- emmeans on these families was run once (`emmeans(fit, ~ 1, dpar =
  "mub")`, `dev/famsimplex-misc.log`) and has no test.
- `confint()` lists the logistic-normal correlations as their chart
  coordinates, `lncor_raw_1`, as it lists every family extra; the
  natural-scale value and interval are in `summary()` and
  `hypothesis()`.
- The simplex families' `phi` keeps no sampler default. CORRECTED in
  fix round 1: an earlier version of this line said brms's
  `gamma(0.01, 0.01)` is a density `set_prior()` does not carry. That
  was false: `set_prior("gamma(0.01, 0.01)", class = "phi")` works in
  `frm()` (MAP phi 6.570) and in `frm_sample()`
  (`dev/famsimplex-rev-sample.log`). The same false sentence was in
  frmtmb.sample's own `?frm_sample` and source comment, and is
  corrected there too. Should frmtmb.sample default to brms's
  `gamma(0.01, 0.01)` on an intercept-only `phi`, and to
  `student_t(3, 0, 2.5)` on its intercept when `phi` has a formula, as
  brms does? Yes, as brms parity. But the gap is the same for every
  family with a `phi`, `shape` or `nu` (`Beta()`, `xbeta()`,
  `beta_binomial()`, `Gamma()`, `negbinomial()`, `weibull()`,
  `student()`'s `gamma(2, 0.1)`, ...), and the disclosure tests of
  those families pin the present default. Doing it for `dirichlet()`
  alone would make the simplex families the one exception, so it is
  left for one change across families, and not done here.

## 6. brms 2.23.0 defects met (upstream-bugs.md style)

### brms-A. The simplex response check is a mean over rows

- Seen: 2.23.0, 2026-10-13.
- Repro: a 300-row `dirichlet()` data set with row 1 off the simplex
  by 1e-7 (`dev/famsimplex-brms-tol.R`, `dev/famsimplex-brms-tol2.R`).
  `brms:::data_response()` tests `is_equal(rowSums(Y), 1)`, which is
  `all.equal()`: the mean absolute departure, here 3.3e-10, against
  1.5e-8.
- Expected: refuse the row, "Response values in simplex models must
  sum to 1." Actual: `standata()` accepts it. `dirichlet()` then fails
  inside Stan with an unrelated message and `brm()` returns a fit with
  no draws; `logistic_normal()` samples, scoring a row that is not a
  composition:

<!-- BEGIN generated: brms dirichlet, one row off the simplex (dev/famsimplex-brms-tol2.R) -->
```
== row 1 off by 0
draws: 200
== row 1 off by 1e-09
draws: 200
== row 1 off by 1e-07
[1] "Error : Exception: no more scalars to read (in 'string', line 42, column 2 to column 23)"
[1] "error occurred during calling the sampler; sampling not done"
draws: 0
```
<!-- END generated -->

- Workaround: frmtmb tests each row at brms's 1.5e-8. It accepts no
  data set brms refuses (a per-row bound implies the mean bound), and
  refuses such a row by number (`dev/famsimplex-edge.log`). The
  divergence is deliberate under the tiebreaker: a per-row constraint
  tested on average is clearly wrong.
- Accepted divergence, recorded in fix round 1: in one case the
  per-row test refuses data brms FITS. With one row of 300 off by 1e-7,
  brms's `logistic_normal()` samples (the log-ratio density is defined
  off the simplex), and frmtmb refuses the row. For `dirichlet()` brms
  accepts the data and then fails in Stan, as above.
- Status: not reported.

### brms-B. A categorical category's own formula stops stancode()

- Seen: 2.23.0, 2026-10-13 (`dev/famsimplex-brms-first.R`).
- Repro: `stancode(bf(f ~ x, mub ~ z), data = d, family =
  categorical())`, `f` a three-level factor.
- Expected: Stan code with `X_mub` on `z`, as brms writes for the same
  formula on `dirichlet()`. Actual: "Error in use_glm_primitive(tmp) &&
  all.equal(tmp$dpars$mu$fe, bterms$dpars[[1]]$fe) : invalid 'y' type
  in 'x && y'" from `use_glm_primitive_categorical()`: `all.equal()`
  returns a string when the two designs differ, and `&&` refuses it.
- Workaround: frmtmb fits it (section 2, parser fix).
- Status: not reported.

## 7. Test runs

Every test file of core and of the seven extensions, one R process per
file, 8 at a time, `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`, through
`dev/release/run-tests.R` with the lane library ahead of rellib-r8
(`dev/famsimplex-run.sh all1`, job list
`dev/famsimplex-tests-all.jobs`):

<!-- BEGIN generated: whole suite, gated (dev/famsimplex-run.sh all1) -->
```
              pk  pass fail err skip warn files
          frmtmb 19094    2   0   14    0   231
 frmtmb.coupling   542    0   0    5    0    11
      frmtmb.eam  1747    0   0    3    0    29
   frmtmb.latent   360    0   0    2    0    10
    frmtmb.learn   501    0   0    2    0    15
      frmtmb.ode   549    0   0    1    0    11
   frmtmb.sample  2828    0   0    1    0    52
   frmtmb.spline   590    0   0    1    0    15
RAN 374 of 374
RESULT lines naming the lane library: 374
RESULT test-brms-families.R pass=1169 fail=2 err=0 skip=0 warn=0
```
<!-- END generated -->

- The 2 failures were a test of mine to update, not a defect:
  `test-brms-families.R` "link and link_<dpar> agree with links on
  every registry family" expected `$link` to be the `mu` link for
  every family not typed categorical; a simplex family's `$link` is
  its joint link, as brms's is. The test now reads `joint_link`;
  rerun below, 0 failures.
- The 29 skips are environment-gated tiers: 15 `FRMTMB_SCALE_TESTS`
  (core and extension `test-scale.R`), 13 `FRMTMB_DRMTMB_FIT_TESTS` and
  1 `FRMTMB_FUZZ`. The last two reach the parser change, so they were
  rerun gated (below). 0 escaped warnings.
- The gated `test-brms-suite-*.R` files (10 in core, 5 in
  frmtmb.sample) have 0 failures: no stale "now HOLDS" row, which
  agrees with section 4.

<!-- BEGIN generated: test run extra1 (dev/famsimplex-run.sh / -ob.sh) -->
```
RESULT test-fuzz.R pass=2 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RESULT test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RESULT test-brms-families.R pass=1169 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RAN 3 of 3
```
<!-- END generated -->

OpenBLAS 0.3.26 with LAPACK routed to it (`dev/ciharden-openblas.sh
0.3.26 lapack`, `OPENBLAS_NUM_THREADS=4`), the four files that test
the new numerical code (`dev/famsimplex-run-ob.sh ob326`):

<!-- BEGIN generated: test run ob326 (dev/famsimplex-run.sh / -ob.sh) -->
```
RESULT test-simplex-families.R pass=100 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RESULT test-simplex-brms.R pass=15 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RESULT test-simulate-density.R pass=561 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RESULT test-simplex-draws.R pass=25 fail=0 err=0 skip=0 warn=0 | core C:/Users/adf44/source/r/wt-famsimplex-lib
RAN 4 of 4
```
<!-- END generated -->

No test of mine decides a rounding-near branch: the density tests are
ratios of 1e3 ulps to the reference's magnitude, the fit tests ratios
to the measured logLik or standard error, and the simulation tests the
file's alpha of 1e-6.

R CMD check --as-cran (`dev/famsimplex-check.sh`, in
`dev/famsimplex-check/`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, R_LIBS
lane library, rellib-r8, user library):

- frmtmb.sample: `Status: OK` (tests 218 s).
- frmtmb: `Status: 1 NOTE`, the expected one ("Skipping checking math
  rendering: package 'V8' unavailable"); tests 12 min, OK. The first
  attempt failed at `R CMD build`, which copies the worktree before
  applying `.Rbuildignore` and met a file the concurrent
  frmtmb.sample check held open; it was rerun alone.
- After the check, one string in `R/families.R` (the mixture refusal)
  was rewrapped, a split with the same text. That rewrap left the next
  line at 88 columns, which the review found; fix round 1 rewraps it to
  80 and reruns the check (section 9).

## 8. Pre-existing, not this lane

Found by the review (`dev/reviews/2026-10-13-famsimplex.md`, R2, R6,
R7) and filed here; none is fixed in this lane.

- `confint(parm = "lncor__b__c")` and `confint(parm =
  "rescor__y1__y2")` are refused with "does not stand for a single
  internal one here". That is false when the matrix is 2 x 2, where the
  one correlation IS one chart coordinate. The refusal then advises
  `hypothesis(method = "profile")`, which fails with "not linear; use
  method = 'boot'". The same holds for rescor on rellib-r8.
- `conditional_effects(dpar = )` reports `se__` on the LINK scale while
  `estimate__` is on the response scale: `dpar = "phi"` gives se 0.054
  where brms's sd is 0.383. A gaussian `sigma ~ z` fit does the same on
  rellib-r8.
- `nobs()` counts the rows `subset()` drops: 300 where brms's N is 299.
  gaussian does the same on rellib-r8.
- frmtmb.sample draws of a group block are named `theta_1 ...`, not
  brms's `sd_`/`cor_`; identical for a categorical fit on rellib-r8
  (`dev/famsimplex-rev-sample-base.log`).
- `drop1()` on a matrix response ("No linear predictor for dpar 'mu'",
  as for `categorical()` on rellib-r8) and
  `marginaleffects::predictions()` on one (dies inside data.table on the
  matrix column; a multinomial fit fails there on rellib-r8 too,
  `dev/famsimplex-rev-basecat.log`).
- brms's `validate_prior()` takes the Stan names `"Llncor"`,
  `"Lrescor"` and `"L"` as classes; frmtmb takes only `"lncor"`,
  `"rescor"` and `"cor"`, which is the existing convention.

## 9. Fix round 1

Review `dev/reviews/2026-10-13-famsimplex.md`, verdict MERGEABLE
AFTER FIXES. What changed, item by item:

1. **A constant on a location is refused** (`R/parse.R`). The parser
   widening stays for formulas. A constant (`mu1 = 0`, `mub = 0`,
   `muc = 0`) on any location of a family with several locations is
   now refused by name. Before, the location was built from a formula
   anyway and the constant had no effect. brms 2.23.0 ignores it the
   same way, and the tiebreaker applies: a value the user wrote,
   silently replaced by an estimate, is clearly wrong. No family can
   honor the constant on this path, because a location here has no
   fixed-value form, so the refusal is the fix. Equating two locations
   (`mu2 = "mu1"`) keeps brms's own refusal. A constant on a dpar that
   is not a location (`phi = 5`, `sigmac = 1`) is held, and the test
   checks `phi = 5` against the hand density. This is marked BREAKING
   in NEWS for `categorical()`, `multinomial()` and `mixture()`, where
   a non-first location's constant used to fit.
2. **`dharma_residuals()` refuses the simplex families** in words
   (`R/diagnostics.R`), where it died on the matrix response.
   `residuals()` on a simplex fit now reaches one refusal
   (`simplex_residuals_refusal()`) for every type and for new data,
   before the polytomous sentence that pointed to `type = "osa"` and
   `dharma_residuals()`. It points to `predict()`.
3. **`conditional_effects(method = "predict")`** keeps its refusal, as
   brms refuses it, with the display's reason: one simulated
   composition is a row over every category. On `dirichlet()` it
   advises `method = "epred"` or `dpar = "mub"`; on `logistic_normal()`
   it advises `dpar = "mub"` only.
4. **Density tolerances** in `test-simplex-families.R` and the
   `log_lik()` check of `test-simplex-draws.R` are now the error over
   the reference's sum of absolute terms, in double epsilons, bound 64.
   The references carry that sum as attribute `scale`.
5. Small items:
   - `R/families.R`'s mixture refusal is wrapped to 80 columns, and
     `R/compat.R`'s `covs <- c(` is restored.
   - Section 5 is corrected on `gamma(0.01, 0.01)`, with the answer
     on a `phi` default: yes for brms parity, but for every family at
     once, so not here. The same false sentence is corrected in
     frmtmb.sample's `?frm_sample` and source comment.
   - The accepted divergence on `logistic_normal()` with one row off
     by 1e-7 is recorded in section 6, brms-A.
   - `dirichlet(refcat = 2)` takes the number as `"2"`, as brms's
     `as_one_character()` does.
   - `?confint` names `thetar_*` and `lncor_raw_*` as chart
     coordinates and says where the natural values are.
6. The pre-existing items are in section 8.
7. NEWS: the constant refusal is a Breaking-changes entry, and the
   sentences are split to 25 words or fewer, in core and frmtmb.sample.

The reviewer's script on the fixed build:

<!-- BEGIN generated: constants on locations, fixed build (dev/famsimplex-rev-pfix.R lane) -->
```
--- mixture y ~ x 
logLik -501.1224413371
--- mixture y ~ x, mu2 = 0 
ERROR: A location parameter of family 'mixture(gaussian, gaussian)' cannot be held at a value: mu2 = 0. Each l
--- mixture y ~ x, mu1 = 0 
ERROR: A location parameter of family 'mixture(gaussian, gaussian)' cannot be held at a value: mu1 = 0. Each l
--- mixture y ~ x, mu2 = "mu1" 
ERROR: Equating parameters of class 'mu' is not allowed. 
--- categorical f ~ x 
logLik -311.5006016534
--- categorical f ~ x, muc = 0 
ERROR: A location parameter of family 'categorical' cannot be held at a value: muc = 0. Each location of this 
--- categorical f ~ x, mub = 0 
ERROR: A location parameter of family 'categorical' cannot be held at a value: mub = 0. Each location of this 
--- dirichlet Y ~ x 
logLik 456.8262864196
--- dirichlet Y ~ x, muc = 0 
ERROR: A location parameter of family 'dirichlet' cannot be held at a value: muc = 0. Each location of this fa
--- dirichlet Y ~ x, mub = 0 
ERROR: A location parameter of family 'dirichlet' cannot be held at a value: mub = 0. Each location of this fa
--- dirichlet Y ~ x, phi = 5 
logLik 450.8987887798
--- logistic_normal Y ~ x, sigmac = 1 
logLik 519.0634227082
--- logistic_normal Y ~ x, muc = 0 
ERROR: A location parameter of family 'logistic_normal' cannot be held at a value: muc = 0. Each location of t
```
<!-- END generated -->

Seen to fail before the fix: `dev/famsimplex-prefix.R off` switches
the four guards off in the source, the lane library is rebuilt and
`test-simplex-families.R` run, and `dev/famsimplex-prefix.R on`
restores the files (checked byte for byte with `cmp`):

<!-- BEGIN generated: pre-fix variant (dev/famsimplex-prefix.R off, dev/famsimplex-run.sh prefix2) -->
```
RESULT test-simplex-families.R pass=108 fail=9 err=2 skip=0 warn=0
RAN 1 of 1
FAIL [a constant on a location of a family with several is refused] Expected `frm(bf(Y ~ x, muc = 0), family = dirichlet(), data = dd)` to throw a err
FAIL [a constant on a location of a family with several is refused] Expected `frm(bf(Y ~ x, mub = 0), family = dirichlet(), data = dd)` to throw a err
FAIL [a constant on a location of a family with several is refused] Expected `frm(bf(Y ~ x, muc = 0), family = logistic_normal(), data = dl)` to throw
FAIL [categorical() and mixture() refuse a constant on a location] Expected mub to throw a error.
FAIL [categorical() and mixture() refuse a constant on a location] Expected muc to throw a error.
FAIL [categorical() and mixture() refuse a constant on a location] Expected `frm(bf(y ~ x, mu1 = 0), family = mix, data = d)` to throw a error.
FAIL [categorical() and mixture() refuse a constant on a location] Expected `frm(bf(y ~ x, mu2 = 0), family = mix, data = d)` to throw a error.
ERROR [the routes a composition has no answer on say so] <frmtmb_error/error/condition> | Error: residuals(type = "response") is not defined for the '
ERROR [dharma_residuals() refuses a composition in words] Error in `sims[keep, , drop = FALSE]`: (subscript) logical subscript too long
FAIL [conditional_effects(method = "predict") says why it is refused] Expected `e` to match regexp "no closed form". | Actual text: | x | method = "pr
FAIL [conditional_effects(method = "predict") says why it is refused] Expected `grepl("epred", e, fixed = TRUE)` to be FALSE. | Differences: | `actual
```
<!-- END generated -->

The scaled error, on the tests' own draws and on 400 wider cases
("old" is the former measure, over `max(1, |ll|)`):

<!-- BEGIN generated: scaled density error (dev/famsimplex-ulps.R) -->
```

== the tests' own draws
dirichlet  worst scaled 3.68  worst old 224 
logistic_normal worst scaled 1.83  worst old 8 

== 400 wider random cases
dirichlet  worst scaled 5.38  worst old 14080 
logistic_normal worst scaled 19.02  worst old 96 
test bound: 64
```
<!-- END generated -->

Reruns. `fix1` ran every file that names `categorical()`,
`mixture()`, `multinomial()`, `mvbf()` or a simplex family, and
`test-parse.R`, in core and the extensions, gated. It ran on the
build before the refusal's last rewording and the dharma test's
split. `fix2` and the OpenBLAS run are on the final build:

<!-- BEGIN generated: parser-reaching files (dev/famsimplex-run.sh fix1) -->
```
            pk pass fail err skip warn files
        frmtmb 8293    0   0   13    0    76
    frmtmb.eam  316    0   0    0    0     4
 frmtmb.latent  224    0   0    0    0     2
 frmtmb.sample 1358    0   0    0    0    16
RAN 98 of 98
RESULT test-drmtmb-agreement.R pass=0 fail=0 err=0 skip=13 warn=0
```
<!-- END generated -->

<!-- BEGIN generated: final build, simplex and gated parser tiers (dev/famsimplex-run.sh fix2) -->
```
RESULT test-simplex-families.R pass=122 fail=0 err=0 skip=0 warn=0
RESULT test-parse.R pass=43 fail=0 err=0 skip=0 warn=0
RESULT test-fuzz.R pass=2 fail=0 err=0 skip=0 warn=0
RESULT test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0
RESULT test-simplex-brms.R pass=15 fail=0 err=0 skip=0 warn=0
RESULT test-simplex-draws.R pass=25 fail=0 err=0 skip=0 warn=0
RAN 6 of 6
```
<!-- END generated -->

<!-- BEGIN generated: OpenBLAS 0.3.26 LAPACK (dev/famsimplex-run-ob.sh ob326b) -->
```
RESULT test-simplex-families.R pass=122 fail=0 err=0 skip=0 warn=0
RESULT test-simplex-brms.R pass=15 fail=0 err=0 skip=0 warn=0
RESULT test-simulate-density.R pass=561 fail=0 err=0 skip=0 warn=0
RESULT test-simplex-draws.R pass=25 fail=0 err=0 skip=0 warn=0
RAN 4 of 4
```
<!-- END generated -->

<!-- BEGIN generated: R CMD check frmtmb (dev/famsimplex-check.sh) -->
```
  Running 'testthat.R' [20m]
* checking HTML version of manual ... [32s] NOTE
Status: 1 NOTE
```
<!-- END generated -->

<!-- BEGIN generated: R CMD check frmtmb.sample (dev/famsimplex-check.sh) -->
```
  Running 'testthat.R' [180s]
Status: OK
```
<!-- END generated -->


### Fix round 2

The refusal reaches the extension families with several locations.
The review showed it for `lba()`, `rdm()` and `hmm()`. `lca()` was
checked here (`dev/famsimplex-lca.R`): on rellib-r8, `theta2 = 0.5`
gives the free fit's logLik and `opt$par` to the last printed digit,
and `theta1 = 0` is refused as "not available"; the lane refuses
both. An offset-only formula holds the location on both builds.

<!-- BEGIN generated: lca constants (dev/famsimplex-lca.R base, lane) -->
```
arm: base  lib: C:/Users/adf44/source/r/rellib-r8/frmtmb 
primary dpars: theta1 theta2 
lca(3) free : logLik -754.7820174645; par -0.330835 0.537468 -2.373153 -0.061715 2.946487 -17.539111 0.080219 1.716812 -
lca(3) theta1 = 0 : ERROR: dpar(s) not available for family 'lca(K = 3)': theta1 (available: theta2) 
lca(3) theta2 = 0.5 : logLik -754.7820174645; par -0.330835 0.537468 -2.373153 -0.061715 2.946487 -17.539111 0.080219 1.
lca(3) theta2 ~ 0 + offset(o) : logLik -754.8737731387; par -0.561678 -2.499968 -0.174659 2.317463 -21.327451 -0.062834 
arm: lane  lib: C:/Users/adf44/source/r/wt-famsimplex-lib/frmtmb 
primary dpars: theta1 theta2 
lca(3) free : logLik -754.7820174645; par -0.330835 0.537468 -2.373153 -0.061715 2.946487 -17.539111 0.080219 1.716812 -
lca(3) theta1 = 0 : ERROR: A location parameter of family 'lca(K = 3)' cannot be held at a value: theta1 = 0. Each locat
lca(3) theta2 = 0.5 : ERROR: A location parameter of family 'lca(K = 3)' cannot be held at a value: theta2 = 0.5. Each l
lca(3) theta2 ~ 0 + offset(o) : logLik -754.8737731387; par -0.561678 -2.499968 -0.174659 2.317463 -21.327451 -0.062834 
```
<!-- END generated -->

Core NEWS names the four extension families, and frmtmb.eam's and
frmtmb.latent's NEWS have a Breaking-changes entry under the
development heading; their code is unchanged. The refusal now names
the offset-only route (`muc ~ 0 + offset(o)`), and the test checks the
clause and that such a fit reads the offset exactly. No Rd changed,
so R CMD check was not rerun.

<!-- BEGIN generated: round-2 runs (dev/famsimplex-run.sh fix3) -->
```
RESULT test-simplex-families.R pass=125 fail=0 err=0 skip=0 warn=0
RESULT test-parse.R pass=43 fail=0 err=0 skip=0 warn=0
RESULT test-lba.R pass=109 fail=0 err=0 skip=0 warn=0
RESULT test-rdm-gng.R pass=221 fail=0 err=0 skip=0 warn=0
RESULT test-hmm.R pass=119 fail=0 err=0 skip=0 warn=0
RESULT test-hmm-starts.R pass=84 fail=0 err=0 skip=0 warn=0
RESULT test-lca.R pass=105 fail=0 err=0 skip=0 warn=0
RAN 7 of 7
```
<!-- END generated -->

## 10. Files touched

Core: `R/families-simplex.R` (new), `R/families.R`, `R/parse.R`,
`R/predict.R`, `R/conditional-effects.R`, `R/confint.R`,
`R/methods-fit.R`, `R/priors.R`, `R/compat.R`, `R/links.R`,
`R/diagnostics.R` (fix round 1), `man/dharma_residuals.Rd`,
`man/confint.frmtmb_fit.Rd`,
`R/links-brms.R` (regenerated), `NAMESPACE`, `man/frmtmb-families.Rd`,
`man/brmsfamily.Rd`, `man/frmtmb-links.Rd`, `man/set_prior.Rd`,
`NEWS.md`, `README.md`, `vignettes/brms-migration.Rmd`,
`dev/famlink-gen-links.R`;
tests `tests/testthat/test-simplex-families.R` (new),
`tests/testthat/test-simplex-brms.R` (new),
`tests/testthat/test-simulate-density.R`,
`tests/testthat/test-brms-families.R`, `tests/testthat/helper-brms.R`.

Lane scripts (new, `dev/famsimplex-*`): `brms-probe.R`,
`brms-fits.R`, `brms-fits2.R`, `brms-tol.R`, `brms-tol2.R`,
`brms-first.R`, `sim.R`, `stan-lpdf.R`, `agree.R`,
`postfit-agree.R`, `recovery.R`, `lpcheck.R`, `edge.R`, `pins.R`,
`mixmu1.R`, `misc.R`, `try.R`, `try2.R`, `explore-multi.R`,
`explore-prior.R`, `tpl.R`, `dbg.R`, `summarize.R`, `build.sh`,
`run.sh`, `run-ob.sh`, `check.sh`, and the job lists `tests-*.jobs`;
in fix round 1 `prefix.R`, `ulps.R`, `summarize-fix1.R` and
`fix1-splice.R`; in fix round 2 `lca.R`. The
reviewer's `dev/famsimplex-rev-*` files are the review's, untouched.

frmtmb.sample: `R/draws-brms.R`, `R/methods-draws.R`, `R/sample.R`,
`man/frm_sample.Rd`, `NEWS.md`,
`tests/testthat/test-simplex-draws.R` (new).

Likely overlaps with other lanes: `R/families.R` (registry entries,
`fam_is_polytomous()`, `family_joint_link()`, `family_link_view()`,
`family_link_str()`, the mixture refusal text, `as_frmtmb_family()`),
`R/links-brms.R` (lane famuni regenerates it too: rerun
`dev/famlink-gen-links.R` after merging both `src` lists),
`R/parse.R`, `R/predict.R`, `R/priors.R`, `R/compat.R`, `NEWS.md`,
frmtmb.sample.

Version: a minor bump of frmtmb (new exported families and a prior
class). frmtmb.sample's code runs on frmtmb 0.70.0 (it reads only a
family field and `us_chol_cor()`), but its new test file needs the new
core, so its floor is the new frmtmb version.
