## Filed at the 0.71.0 release (2026-10-14)

Seven lanes off 0.70.0 (`dev/round-20261014.md`). What each lane filed
and did not fix, from its findings, with the record to read. The items
this round closed are marked "CLOSED at 0.71.0" in place above. Lane
famuni's own section follows this one as the lane wrote it.

### Found at the consolidation

See `dev/round-20261014.md`, "Found at the consolidation".

- **The count header of `test-simulate-density.R` was short at 0.70.0**
  (low). It said 251 statistical assertions; the file made 289 on
  rellib-r8 (`dev/rel071-simcount.R`, `dev/rel071-log/simcount-head-r8.txt`:
  the `P(Y = 1)` atoms and ten cells were not counted). The header now
  states the count measured on the merged build, 363. CLOSED at 0.71.0.
- **illse's test of an sd held near zero lost its premise under
  optstart's least-squares start** (`test-se-conv.R`; 5 of its
  expectations failed on the merged tree). The first replacement
  (every parameter started at 0) did not pin illse's R2-1 fix: the
  intercept was short too, so the test passed with the fix removed
  (the release review, 6a). The test now starts at `lm()`'s maximum
  given sd = 0 with `iter.max = 0`, so only the sd is short; it fails
  on 0.70.0, fails with `se_sd_gain_up()` stubbed to 0, and passes on
  the merged build (`dev/rel071-seconv3.R`). CLOSED at 0.71.0.
- **`_pkgdown.yml` lacked `frmtmb-extreme-families`** (lane famuni's
  new topic; `check_pkgdown()` stopped). CLOSED at 0.71.0.
- **24 bin-2 rows of the ported suite carried a false verdict** once
  lanes famuni's and famsimplex's families were in one build: 5 now
  pass, `log_lik:349` is a divergence (domain), and 18 move from
  "absent" to "not ported (cost)". Each has a true verdict and reason
  (`dev/rel071-bin2-moves.R`, `dev/rel071-log/bin2-diff.txt`). The
  ledger had stopped on 14 rows (lane postfix 9, lane famuni 5).
  CLOSED at 0.71.0.
- **The families instrument had no data generator for the five new
  families** (`dev/brms-port/families-coverage.R`): added, each on its
  own seed; 39 of 40 fitted by core, 40 of 40 with frmtmb.eam.
  CLOSED at 0.71.0.

### Filed by the release review (`dev/reviews/2026-10-14-release.md`)

#### Open - low

- **A correlation's boundary verdict depends on where the stop lands
  on a flat ridge** (section 2, `dev/relrev071-rssd.R`). On the
  false-alarm battery's `rs_sd0.02` seed 4 the merged build stops at
  t = 19.66 (rho 0.99871, code 1, the correlation's SE 1.86e4, "at
  the maximum"), and 0.70.0 at t = 23.75 (rho 0.99911, code 0, a
  boundary message); both are at the same maximum to 6e-8, and lme4
  gives rho 0.99926. ci701's `se_edge_cor()` acts only past t = 70.7.
  It goes with ci701's deferred general case, a pull in correlation
  space (lane ci701's section, "The warning side of a correlation
  past the clamp"). Pre-existing; the merge makes it visible.
- The false-alarm battery's new warnings on the merged build are 7 of
  4000 fits (reference BLAS) and 6 (OpenBLAS 0.3.26) over 100 seeds
  per model, all `rs_sd`, each within twice the base SE error; fewer
  warnings in 19 and 22 (section 2). The record's "1 new convergence
  warning" was the count over 20 seeds.

### Filed by lane dblbar (`dev/dblbar-findings.md`, "Not done")

#### Open - medium

- **sd priors per coefficient.** brms's `default_prior()` lists an sd
  row per coefficient and honors `set_prior(class = "sd", coef = )`;
  frmtmb lists class and group rows only and refuses `coef` on class
  sd. A `||` block could honor it. A priors lane.
- **Nested and multi-membership prior labels.** frmtmb labels
  `(1 + x | g1/g2)` `g2:g1` and `g1`, and `mm(m1, m2)` as written; brms
  writes `g1:g2` and `mmm1m2`, so `set_prior(group = "g1:g2")` is
  refused. Pre-existing for `|` terms (`dev/dblbar-preexist.R`).
- **brms shares one Student-t latent scale across all the Student-t
  terms of a grouping factor**; frmtmb gives each block its own.
  `?frmtmb-student-re` says to write one term to share the scale.
- **`frm_sample(seed = )` does not fix every chain** (frmtmb.sample):
  two calls with `seed = 3` and 2 chains differ unless `set.seed()`
  comes first, on rellib-r8 too (`dev/dblbar-rev2-seed2.R`).

#### Open - low

- `variables()` lists `sd_g1:g2__*` before `sd_g1__*` (brms the
  reverse); `gr(g, id = )` and group-level `cs()` refused;
  `re_formula = ~ diag(x | g)` names no term. Pre-existing.
- Messages and tables that read a block's `term_label` still say
  `1 + f | g [diag]` for a `||` term (the class "cor" prior refusal,
  `confint_varcorr()$block`); `term_label` is a matching key, so the
  change is not cosmetic.
- A diagonal within-level covariance for `gr(cov = )` / `gr(prec = )`
  would let a factor sit in a `||` term over a known matrix; now
  refused.
- frmtmb.sample's draws carry `theta_*`, not brms's `sd_` names, in
  `variables(draws)` and `as_draws_df()`.

### Filed by lane famsimplex (`dev/famsimplex-findings.md`, sections 5 and 8)

#### Open - low

- `logistic_normal()` in `mvbf()` is refused (family extras in a
  multivariate frame); `dirichlet()` in `mvbf()` is untested.
- `posterior_linpred()` without `dpar` is one predictor's `S x N`
  matrix for every multi-predictor family; brms gives `S x N x (K - 1)`.
- `categorical()` has no `refcat`, so a brms `categorical(refcat = )`
  object loses it.
- `multinomial()`'s `fitted()` refuses; brms returns `n x 4 x K`.
- The simplex families' `phi` (and every family's `phi`, `shape`,
  `nu`) keep no frmtmb.sample default where brms has one; one change
  across families.
- Pre-existing, from the review: `confint(parm = "lncor__b__c")` and
  `"rescor__y1__y2"` refused on a 2 x 2 matrix with advice that fails;
  `conditional_effects(dpar = )` reports `se__` on the link scale
  beside a response-scale estimate (gaussian `sigma ~ z` too); `nobs()`
  counts rows `subset()` drops; `drop1()` and
  `marginaleffects::predictions()` on a matrix response fail.

### Filed by lane postfix (`dev/postfix-findings.md`, "Not done", "Defects")

#### Open - medium

- **The review's Emax model M6 converges silently 114.85 below its
  maximum** (`dev/postfix-rev-nlscan.R`, pre-existing).
- `conditional_effects(<frm_multiple()>, int_conditions = )` stops with
  "matched by multiple actual arguments" (pre-existing).

#### Open - low

- `predict()` and `conditional_effects(method = "predict")` of a
  `frm_multiple()` result stay refused (a mixture of simulated
  replicates is not Rubin's rules); `update()` of one is refused
  (`frm_multiple()` keeps neither its call nor the data list).
- `fitted(ndraws = )` stays a Wald summary; `summary = FALSE` has no
  answer.
- brms's "first imputed data set" warning on in-sample methods of
  pooled draws is not copied (documented).

### Filed by lane optstart (`dev/optstart-findings.md`, "Filed")

#### Open - medium

- **At a condition number of 2e4 to 1e5 with `(1 | g)`**, where
  `newton_polish()` does not run, 7 of 40 fits stop short of lme4 at
  code 0 with no warning (up to 0.143 on the lane, 0.038 on rellib-r8;
  `dev/optstart-rev-collin.R`). Pre-existing; the start moves which
  seeds.
- **exgaussian's density** (`RTMBdist::dexgauss`) ends at "NA/NaN
  gradient" on a group sd 30 times sigma, 5 of 10 fits on both builds
  and from the true values too.

#### Open - low

- huber() reached +1.39e5 on 400 rows on rellib-r8 (a degenerate stop
  of its density); fuzz specs with `weights(w)` and `ar1()` on beta and
  Gamma reach objectives near -1e15 on rellib-r8.
- 0.70.0's initializers read an `mi()` response's imputed zeros outside
  the least-squares rule (`bf(y | mi() ~ 1)`).
- The wide design ends at "false convergence (8)" at lme4's maximum on
  1 of 3 seeds from any start; an anova start halves its time.
- `fit_assembled()` solves the unit's least squares a second time.
- The gr(g, by = f) fixture of `test-se-check.R` and `test-se-tier.R`
  pins SE verdicts at a point 1.44e-4 below the maximum.
- `fd_hessian()` keeps a step that leaves the matrix positive definite
  however large its rise.
- Not done: a within-group or anova start for sigma and the group sd;
  families outside the rule (tweedie, compois, nbinom1, beta_binomial,
  the zero-inflated and hurdle families, von_mises, inverse.gaussian,
  Gamma's inverse link, cox) keep their starts; nonlinear formulas keep
  theirs.

### Filed by lane remeasure (`dev/remeasure-findings.md`, "Defects found")

- R1 (two `mi()` fits lose their SEs at 0.70.0): CLOSED at 0.71.0 by
  lane optstart (`fd_hessian()`'s retake; `dev/optstart-findings.md`,
  "Item R1"). R2 to R4 (draws' `loo(newdata = )`, `logLik()`,
  `bayes_R2(newdata = )`): CLOSED at 0.71.0 by lane postfix (item h).
- **R5** (low): frmtmb.eam's default `ndt` link is not brms's; lane
  famuni added `wiener(link_ndt = )` with brms's links and kept the
  bounded default, and `family = brms::wiener()` now reaches frmtmb.eam
  with brms's link. The default remains a documented difference.
- **R6** (decided): `com_poisson`'s `mu` is the mean in frmtmb; brms
  2.23.0 fits the rate and post-processes another pmf (brms-24 in
  `dev/upstream-bugs.md`). The difference is documented in
  `?brmsfamily` and `?frmtmb-families` at 0.71.0.
- **R7** (low): the ported-suite harness lets staleness lapse through a
  setup line that reads a stale name (`helper-brms-suite.R`); 0 held
  assertions read one today (`dev/remeasure-stale-reach.R`).
- **R8** (low): a family's `$lcdf` slot returns the CDF, by documented
  convention; the name misleads.
- The bin-2 port leaves 49 rows "not ported (cost)" and 41 gaps in its
  denominator at 0.71.0 (`dev/brmsport-bin2-summary.md`). 18 of the 49
  are rows of families this round added that the thin adapter or the
  G1 shims do not map yet (dirichlet, dirichlet2 and logistic_normal's
  per-category `mu` and `lncor`; brms's d/p/r functions of frechet and
  gen_extreme_value): the cheapest bin-2 rows to port next.
NA
  denominator (`dev/brmsport-bin2-summary.md`).
