# Known bugs in packages that frmtmb uses or compares against

This file lists bugs and documented defects in other packages that the
frmtmb records name. Use it to file reports upstream. It does not list
frmtmb's own defects or the deliberate divergences from brms. Each
entry gives the version, a minimal repro, the expected result (and
what happens), the frmtmb workaround, the source record and the report
status. "Confirmed 2026-09-29" means that the repro was run again on
that date (R 4.6.1, Windows, `rellib-r3` and the user library). Paths
that start with `frmtmb-wt-fams2/` are in the uncommitted fams2
worktree, `C:/Users/adf44/source/r/frmtmb-wt-fams2`.

Installed on 2026-09-29: brms 2.23.0, RTMB 2.0, TMB 1.9.25, tmbstan
1.2.1, rstan 2.32.7, StanHeaders 2.39.1, RTMBode 1.0, deSolve 1.42,
mgcv 1.9.4, reformulas 0.4.4, lme4 2.0.6, knitr 1.52, typetracer 0.2.5,
autotest 0.2.0, emmeans 2.0.4, drmTMB 0.7.0, hmmTMB 1.1.2, ordinal
2026.7.26, MASS 7.3.66. gratia is not installed.

## brms

### brms-1. `hurdle_cumulative` logit density tests the wrong top category

- Seen: 2.23.0. Confirmed 2026-09-29 (source read).
- Repro: `stancode(bf(y ~ x, disc ~ 0 + x), data = d, family = hurdle_cumulative("logit"))`
  on a 0..4 response. `brms:::stan_hurdle_ordinal_lpmf`, `inv_logit`
  branch, emits `y == nthres + 2`; the probit branch emits `nthres + 1`.
- Expected: test `y == nthres + 1`. Actual: the top category reads
  `thres[nthres + 1]` and Stan stops, "index 4 out of range". brms uses
  this function whenever `disc` is modeled or `cs()` is present.
- Workaround: none needed. Row 16g uses probit for this reason
  (`frmtmb-wt-fams2/tests/testthat/test-brms-likelihood.R:463`).
- Record: `frmtmb-wt-fams2/dev/reviews/2026-09-29-fams2.md` claim 4;
  `frmtmb-wt-fams2/dev/fams2-findings.md` "brms defect"; repro
  `frmtmb-wt-fams2/dev/fams2-brms-hc-bug.R`. Met again 2026-10-05 by
  lane ordmix on the `cs()` path: `yh ~ cs(x)` with
  `hurdle_cumulative()` (logit) stops at "index 4 out of range"
  (`dev/ordmix-findings.md`, "brms 2.23.0 defects met").
- Status: not reported.

### brms-2. The softit link does not compile

- Seen: 2.23.0. Confirmed 2026-09-29.
- Repro: `stancode(y ~ x, data = d, family = bernoulli("softit"))`
  emits `return log(expm1(-p / (p - 1)));` and
  `return log1p_exp(y) / (1 + log1p_exp(y));` on vectors.
- Expected: `./` for element-wise division. Actual: Stan rejects the
  program; no softit model compiles.
- Workaround: none needed. Row 22c asserts the defect in text and fails
  when brms fixes it (`tests/testthat/test-brms-likelihood.R:1377`).
- Record: `dev/brms-likelihood-tests.md`, row 22c. Met again on
  2026-09-30 by lane ordinal: `acat("softit")` does not compile for the
  same reason, so acat softit is checked against brms's R side
  (`brms:::dacat()`) instead (`dev/ordinal-findings.md`, "brms 2.23.0
  defects met", item 3).
- Status: not reported.

### brms-3. `posterior_predict()` for `hurdle_negbinomial` does not draw the zero-truncated NB

- Seen: 2.23.0. Confirmed 2026-09-29 (source read).
- Repro: `brms:::posterior_predict_hurdle_negbinomial` draws
  `rnbinom(mu = mu - t, size = shape) + 1`, `t` a truncated exponential;
  `dev/fams-brms-pp-check.R`.
- Expected: zero-truncated NB draws. Actual: at `mu = 1.5`,
  `shape = 0.8`, P(Y = 1 | Y > 0) is 0.567 (true 0.393, z = 505). The
  shift trick is exact for the Poisson only.
- Workaround: frmtmb's `sim` draws the exact distribution
  (`R/families.R:2378`).
- Record: `dev/fams-findings.md` "Defects found and not fixed" item 1;
  `dev/round-handoff.md` "Upstream reports".
- Status: not reported.

### brms-4. The expected ordinal category of `hurdle_cumulative` is off by one

- Seen: 2.23.0. Confirmed 2026-09-29 (source read).
- Repro: `conditional_effects(fit, categorical = FALSE)` on a
  `hurdle_cumulative` fit. `brms:::ordinal_probs_continuous()`
  multiplies column `k` by `k`, and column 1 is category 0.
- Expected: weight each column by its code (0..K). Actual: E[Y] + 1.
- Workaround: frmtmb weights by the codes
  (`frmtmb-wt-fams2/R/conditional-effects.R:1936`).
- Record: `frmtmb-wt-fams2/dev/fams2-findings.md`, "DIVERGENCE FROM
  BRMS" and "Defects found and not fixed" item 3. The lane records it
  as a divergence; the user decides whether to report it.
- Status: not reported.

### brms-5. A nonlinear body that names `mu` is silently discarded

- Seen: 2.23.0. Confirmed 2026-09-29 (`stancode()`).
- Repro: `stancode(bf(y ~ mu - 2 * x, mu ~ 1, nl = TRUE), data = d)`
  gives `vector[N] mu = rep_vector(0.0, N); mu += Intercept;`.
- Expected: an error, as for `sigma` ("No non-linear parameters
  specified"). Actual: an intercept-only model, no warning.
- Workaround: frmtmb refuses the collision by name at parse time
  (`R/parse.R:1895`).
- Record: `dev/api-findings.md` (the `mu` collision); `NEWS.md` near
  line 2555.
- Status: not reported.

### brms-6. `pp_check(re.form = NA)` warns that it ignores `re.form`, then honors it

- Seen: 2.23.0. Not rerun (needs a fitted brmsfit).
- Repro: `pp_check(fit, re.form = NA)` against `pp_check(fit)` and
  `pp_check(fit, re_formula = NA)`, same seed.
- Expected: the warning and the result agree. Actual: the warning says
  "unrecognized and ignored", but the value reaches
  `posterior_predict()` and changes the plot data.
- Workaround: frmtmb refuses `re.form` on `pp_check()`
  (`R/conditional-effects.R:2539`).
- Record: `dev/reviews/20260916-argspell.md`, near line 692.
- Status: not reported. Minor.

### brms-7. ARMA predictions and the ARMA likelihood use different scales under a non-identity link

- Seen: 2.23.0. Not rerun.
- Repro: `stancode(bf(y ~ x + arma(p = 1, q = 1)), data = d, family = gaussian("log"))`
  applies `mu = exp(mu)` before the ARMA loop (response scale);
  `brms:::.predictor_arma()` adds the term to `eta` (link scale).
- Expected: one scale in both. Actual: brms's predictions disagree
  with its own likelihood.
- Workaround: frmtmb follows the likelihood everywhere.
- Record: `dev/arcov-findings.md`, near line 42.
- Status: not reported.

### brms-8. `probit_approx` is two different functions inside brms

- Seen: 2.23.0. Not rerun.
- Repro: the Stan program uses `Phi_approx()`
  (`inv_logit(0.07056 x^3 + 1.5976 x)`); `brms:::inv_link(x, "probit_approx")`
  uses `pnorm()`.
- Expected: one function per link name.
- Workaround: frmtmb uses the Stan form; the test asserts it
  (`tests/testthat/test-brms-likelihood.R:1326`).
- Record: `dev/links-findings.md`, near line 79. Met again on
  2026-09-30 by lane ordinal on `acat("probit_approx")`:
  `brms:::inv_link()` reads it as `pnorm()`, the Stan program as
  `Phi_approx()`. frmtmb follows the likelihood, and
  `test-ordinal-disc-thres.R` checks its acat against brms's formula
  with `Phi_approx` (`dev/ordinal-findings.md`, "brms 2.23.0 defects
  met", item 5).
- Status: not reported.

### brms-9. The comment on `cholesky_cor_ar1()` misdescribes what it returns

- Seen: 2.23.0. Not rerun.
- Repro: the generated function ends
  `return cholesky_decompose(mat ./ (1 - ar^2));`; its comment says
  "the cholesky factor of an AR1 correlation matrix".
- Expected: a comment that says it is the stationary covariance with
  unit innovation variance, so `sigma` is the innovation SD.
- Workaround: frmtmb rescales by `sqrt(1 - ar^2)`
  (`tests/testthat/test-brms-likelihood.R:1037`).
- Record: `dev/brms-likelihood-tests.md`, near line 938.
- Status: not reported. Documentation only.

### brms-10. `me(x, sx):f` with a factor multiplier emits both columns against one `Csp_1`

- Seen: 2.23.0. Not rerun; the repro script was in `/tmp` and is gone.
- Repro: `stancode()` of a model with `me(x, sx):f`, `f` a factor.
- Expected: one coefficient per column. Check the generated code
  before you file.
- Workaround: frmtmb refuses factor multipliers.
- Record: `dev/me-findings.md`, near line 37.
- Status: not reported.

### brms-11. The `hurdle_negbinomial` normalizer loses precision at a small mean

- Seen: 2.23.0. Not rerun.
- Repro: `log1m((shape / (mu + shape))^shape)` at `mu = exp(-25)`.
- Expected: `pnbinom(0, lower.tail = FALSE, log.p = TRUE)` accuracy.
  Actual: relative error 1.53e-7.
- Workaround: `RTMB::logspace_sub()` (`R/families.R:2351`).
- Record: `dev/fams-findings.md`, near line 27.
- Status: not reported. Minor.

### brms-12. The `brms_multilevel` vignette reads a dead UCLA data URL

- Seen: 2.23.0 vignette sources.
- Expected: a live URL or packaged data.
- Workaround: the port harness shims a mirror.
- Record: `dev/brms-vignette-port.md`, near line 40.
- Status: not reported.

### brms-13. `conditional_effects()` on `(1 | gr(g, by = f))` stops when `g` is unset

- Seen: 2.23.0. Found 2026-09-29 (lane postfit2, punch round 2); not
  rerun independently.
- Repro: fit `y ~ x + (1 | gr(g, by = f))`, then
  `conditional_effects(fit, "x", conditions = data.frame(f = "a"),
  re_formula = NULL)`. brms stops with "missing sd_g__Intercept:f".
  With `g = "99"` (an unseen level) it answers.
- Expected: a new level of `g` drawn with the SD of the row's `f`
  level, as for `g = "99"`.
- Workaround: frmtmb answers the unset case as it answers `"99"`.
- Record: `frmtmb-wt-postfit2/dev/postfit2-findings.md` section 7d;
  `dev/postfit2-p2-brms.R`.
- Status: not reported.

### brms-14. Two different unseen `mm()` members share one new-level draw

- Seen: 2.23.0 source. Found 2026-09-29 (postfit2 final check).
- Where: `data_gr_local()` numbers each member's unseen values
  separately (`match(new_gdata, new_levels) + length(levels)`, per
  member), so `"99"` in `g1` and `"98"` in `g2` both become
  `nlevels + 1` and read one draw. For a single grouping factor, two
  different unseen values get separate draws.
- Expected: distinct unseen values get distinct draws, as for `(1 | g)`.
- Workaround: frmtmb draws them independently; `?conditional_effects`
  names the difference.
- Record: `frmtmb-wt-postfit2/dev/reviews/2026-09-29-postfit2.md`,
  final check, item 2.
- Status: not reported. Likely low impact.

### brms-15. `cratio("cloglog")` with a modeled `disc` has a NaN gradient on some data

- Seen: 2.23.0. Found 2026-09-30 (lane ordinal).
- Repro: `dev/ordinal-brms-defects.R`, item 1. Fit a `cratio("cloglog")`
  model with `disc ~ x` and evaluate brms's `log_prob` and its gradient
  at frmtmb's optimum. The value equals frmtmb's logLik
  (-441.1943809182 in both), and 4 of 6 gradient entries are NaN.
- Expected: a finite gradient, as frmtmb's (1.9e-4 there). The NaN
  depends on the data: it comes when a row has `disc * (mu - thres_k)`
  above about 6.6 below its category. brms's line is `q[k] =
  log1m_exp(-exp(disc * (mu - thres[k])))`, which loses the value in
  that range. At other points (disc slope 0, or the x slope halved) the
  gradient is finite.
- Workaround: none needed. The gated Stan identity row for cratio
  cloglog with disc uses seed 20260931, whose largest such value is
  4.70 (`dev/ordinal-p1-cloglog.R`, `dev/ordinal-p1-log-cloglog.txt`;
  seed 20260930 has 8.75).
- Record: `dev/ordinal-findings.md`, "brms 2.23.0 defects met", item 1,
  and punch round 1, m3; output `dev/ordinal-log-brms-defects.txt`.
- Status: not reported.

### brms-16. `cumulative(threshold = "sum_to_zero")` with the logit link does not compile

- Seen: 2.23.0. Found 2026-09-30 (lane ordinal).
- Repro: `stancode()` of a `cumulative("logit", threshold =
  "sum_to_zero")` model with `disc` held at 1 emits
  `ordered_logistic_glm_lpmf(Y | X, b, 0)`, with the literal 0 as the
  thresholds.
- Expected: the sum-to-zero threshold vector. Actual: stanc refuses the
  program. The probit link and a modeled disc take other code paths and
  compile.
- Workaround: none needed. The likelihood rows use probit and a
  modeled disc for this structure.
- Record: `dev/ordinal-findings.md`, "brms 2.23.0 defects met", item 2;
  `dev/ordinal-brms-defects.R`.
- Status: not reported.

### brms-17. Equidistant thresholds in a multivariate or mixture model declare `delta` twice

- Seen: 2.23.0. Found 2026-09-30 (lane ordinal); widened by its review.
- Repro: `stancode()` of a multivariate model with two ordinal
  responses and `threshold = "equidistant"` (`dev/ordinal-mv.R`).
  `stan_thres()` suffixes the prior with the group only, so `delta` is
  declared twice, while the transformed parameters use `delta_y` and
  `delta_y2`. A multivariate model with ONE equidistant response also
  fails ("delta_y not in scope"), and so does an equidistant ordinal
  mixture (`dev/ordinal-rev-log-defects.txt`).
- Expected: one `delta_<suffix>` per response or component. Actual:
  stanc refuses the program.
- Workaround: none needed. frmtmb names the parameters as brms's
  transformed-parameter code does, `delta_y` and `delta_y2`; frmtmb
  fits ordinal mixtures since 0.68.0 and names one delta per
  component, `delta_mu1`, as brms's transformed parameters read them;
  `dev/ordmix-lpcheck.R` compiles brms's program with the declarations
  renamed that way (its `patch_delta`).
- Record: `dev/ordinal-findings.md`, "brms 2.23.0 defects met", item 4,
  and punch round 1, m4; `dev/ordinal-log-brms-code.txt`. Confirmed
  for mixtures 2026-10-05 by lane ordmix (`dev/ordmix-findings.md`,
  "brms 2.23.0 defects met": one bare `delta` per component, read as
  `delta_mu<k>`, and `delta_mu1_1` read beside a declared `delta_1`
  under `thres(gr = )`).
- Status: not reported.

### brms-18. A new level of an `mm(by = )` term cannot be drawn

- Seen: 2.23.0. Found 2026-09-30 (lane ceplot, punch round 1).
- Repro: fit `y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2)))`, then ask
  `conditional_effects()` for a row with a new member. brms stops with
  one of three errors, by condition
  (`dev/ceplot-rev-log/mmby-brms.txt`, seeds 1 to 3 each). Both members
  unset: "The following variables are missing in the draws object:
  {'sd_mmg1g2__Intercept:cbind(f1, f2)'}". One member seen and one
  unset: "all(bylevels %in% reframe$bylevels[[1]]) is not TRUE". Both
  unset at different by-levels: "Some levels of 'g1', 'g2' correspond
  to multiple levels of 'cbind(f1, f2)'". With both members observed,
  brms answers and equals frmtmb to 0.
- Expected: a new member drawn with the SD of its own by-level block,
  as for `(1 | gr(g, by = f))`.
- Workaround: frmtmb draws each new member in its own by-level block;
  the check is against the Wald band (`dev/ceplot-mmby.R`).
- Record: `dev/ceplot-findings.md` sections 1.7 and 10 (m9).
- Status: not reported.

### brms-19. An all-`TRUE` bernoulli response is coded as all failures

- Seen: 2.23.0. Found 2026-09-30 (lane formrobust; its review agreed).
- Repro: `standata(y ~ 1, data.frame(y = rep(TRUE, 5)), bernoulli())`.
  `data_response.brmsframe()` takes the levels from
  `levels(as.factor(Y))`, which has one level, `"TRUE"`, and a
  one-valued response becomes `c(0, value)`, so `TRUE` is coded 0
  (`dev/formrobust-brms-binary.R`,
  `dev/formrobust-log/brms-binary.txt`).
- Expected: `TRUE` coded 1, as in a response that holds both values.
  Actual: every row is a failure.
- Workaround: frmtmb reads a logical as its number, so all-`TRUE`
  codes 1. Recorded as a deliberate divergence in NEWS and
  `dev/formrobust-findings.md` section 5.
- Record: `dev/formrobust-findings.md` section 5, "Deliberate
  divergences" (a); `dev/reviews/2026-09-30-formrobust.md`, minor 2.
- Status: not reported. Low impact: a response that holds one value
  is rare.

### brms-20. An ordinal mixture's gradient is NaN where its density is finite

- Seen: 2.23.0. Found 2026-10-05 (lane ordmix).
- Repro: `dev/ordmix-lpcheck.R acat_probit`:
  `mixture(acat("probit"), cumulative("cloglog"))` on `y ~ x`, brms's
  compiled program with flat priors, `rstan::grad_log_prob()` at
  frmtmb's optimum. Each component alone (`acat_probit_1`,
  `cum_cloglog_1`) has a finite gradient (0.00137 and 0.000163).
- Expected: a finite gradient; brms's log density there equals
  frmtmb's to 2.5 ulp, and frmtmb's gradient is 2.8e-4. Actual: NaN.
  It comes from the mixture's `log_sum_exp` over a component term that
  underflows on some row, the pattern of lane ordinal's defect 1.
- Workaround: none needed; frmtmb's densities are written in log
  space.
- Record: `dev/ordmix-findings.md`, "brms 2.23.0 defects met".
- Status: not reported.

### brms-21. `posterior_linpred(incl_thres = TRUE)` on `hurdle_cumulative()` returns no linear predictor

- Seen: 2.23.0. Found 2026-10-05 (lane fixes and its review).
- Repro: `dev/fixes-rev-thres2.R`. brms runs
  `posterior_epred_hurdle_cumulative()` at the identity link: column
  `0` is `hu`, and columns 1 to K are `(1 - hu) * (tau_k - eta)`
  (`all.equal` TRUE on 220 rows).
- Expected: the threshold predictors `disc * (tau_k - eta)`, as for
  `cumulative()`, or a refusal. Actual: a probability beside linear
  predictors scaled by `1 - hu`, which is the predictor of nothing.
- Workaround: frmtmb refuses the call and says why, a divergence the
  user decided on 2026-10-06.
- Record: `dev/fixes-findings.md`, punch round 1, m9;
  `dev/reviews/2026-10-05-fixes.md`, m9.
- Status: not reported.

### brms-22. Crossing `cs()` thresholds give a negative probability and a draw from another distribution

- Seen: 2.23.0. Found 2026-10-05 (the ordmix review).
- Repro: `dev/ordmix-rev-crossing.R`. `brms:::dcumulative()` at the
  thresholds 0, 1, 0.5 returns 0.5, 0.2311, -0.1086, 0.3775, and the
  draw rule of `posterior_predict()` gives the shares 0.505, 0.222, 0,
  0.274.
- Expected: no probability below zero; at least, draws from the
  distribution `posterior_epred()` reports. Actual: a negative
  `posterior_epred()` entry, and draws from a different distribution.
  The Stan density is undefined there (`log_inv_logit_diff()`), so a
  sampler rejects the point.
- Workaround: frmtmb gives `NaN` for such a row's `fitted()`
  probabilities and `NA` for its draw, on `hurdle_cumulative()` and,
  since 0.68.0, on `cumulative()`.
- Record: `dev/reviews/2026-10-05-ordmix.md`, "Post-fit";
  `dev/round-20261005.md`, small change a.
- Status: not reported.

### brms-23. `brm(threshold = )` is accepted and dropped

- Seen: 2.23.0. Found 2026-10-05 (lane vigport, repro R7).
- Repro: `dev/vigport-repros.R` R7: the JSS paper's
  `brm(..., family = sratio(), threshold = "equidistant")`. brm()
  says "passing unknown arguments: threshold" and fits flexible
  thresholds; `stancode()` drops the argument with no message
  ("brms: delta in program, threshold = argument FALSE").
- Expected: a refusal, or the threshold structure. Actual: the
  argument is ignored, silently in `stancode()`.
- Workaround: in frmtmb `threshold` is an argument of the family
  constructor, `sratio(threshold = "equidistant")`, as in current brms.
- Record: `dev/vigport-findings.md`, "Defects found, not fixed" (R7)
  and the ranked list, item 6.
- Status: not reported.

### brms-24. `com_poisson`: the R-side pmf is not the pmf the Stan program fits

- Seen: 2.23.0. Confirmed 2026-10-13 (source read and a fit; lane
  remeasure, R6).
- Repro: brms's Stan function `com_poisson_log_lpmf(y | log_mu, nu)`
  (`chunks/fun_com_poisson.stan`) is `y * log_mu - nu * lgamma(y + 1)
  - log_Z_com_poisson(log_mu, nu)`, with series terms
  `(k - 1) * log_mu - nu * lgamma(k)`: the pmf is proportional to
  `mu^y / (y!)^nu`, so the fitted `mu` is Shmueli's rate lambda.
  brms's R side, `brms:::dcom_poisson(x, mu, shape)` (used by
  `log_lik()`, `posterior_predict()` and `loo()`), is
  `shape * (x * log(mu) - lgamma(x + 1)) - log_Z`, with
  `brms:::log_Z_com_poisson()` summing `shape * (k * log_mu - lfac)`:
  proportional to `(mu^x / x!)^shape`, so there `mu = lambda^(1/shape)`.
  At `mu = 3, shape = 2` the two pmfs have different means: 1.453549
  on the Stan side, 2.737078 on the R side (exact sums over 0:200,
  `dev/remeasure-compois-means.R`).
- Expected: one pmf. Actual: post-processing of a `com_poisson` fit
  evaluates and simulates a different pmf from the one sampled, for
  every `shape != 1`; `shape = 1` agrees (Poisson). The R-side
  normalizer is also approximate in its asymptotic branch (mass
  1.00014 at `mu = 3, shape = 2`).
- Workaround: none needed in frmtmb, which calls neither. Its own `mu`
  is the mean (`RTMB::dcompois2()`), a third parameterization, and
  `?brmsfamily` and `?frmtmb-families` say so.
- Record: `dev/remeasure-findings.md` R6; `dev/reviews/2026-10-13-remeasure.md`
  section C, `dev/remeasure-rev-compois.R` (brms fit: 2 chains of 1000,
  seed 20261013).
- Status: not reported.

### brms-25. The simplex response check is a mean over rows

- Seen: 2.23.0, 2026-10-13 (lane famsimplex, brms-A).
- Repro: a 300-row `dirichlet()` data set with row 1 off the simplex
  by 1e-7 (`dev/famsimplex-brms-tol.R`, `dev/famsimplex-brms-tol2.R`).
  `brms:::data_response()` tests `is_equal(rowSums(Y), 1)`, which is
  `all.equal()`: the mean absolute departure, here 3.3e-10, against
  1.5e-8.
- Expected: refuse the row, "Response values in simplex models must
  sum to 1." Actual: `standata()` accepts it. `dirichlet()` then fails
  inside Stan with an unrelated message ("Exception: no more scalars
  to read") and `brm()` returns a fit with no draws;
  `logistic_normal()` samples, scoring a row that is not a
  composition.
- Workaround: frmtmb tests each row at brms's 1.5e-8 and refuses such a
  row by number. In one case this refuses data brms fits: with one row
  of 300 off by 1e-7, brms's `logistic_normal()` samples.
- Record: `dev/famsimplex-findings.md`, section 6, brms-A.
- Status: not reported.

### brms-26. A categorical category's own formula stops `stancode()`

- Seen: 2.23.0, 2026-10-13 (lane famsimplex, brms-B;
  `dev/famsimplex-brms-first.R`).
- Repro: `stancode(bf(f ~ x, mub ~ z), data = d, family =
  categorical())`, `f` a three-level factor.
- Expected: Stan code with `X_mub` on `z`, as brms writes for the same
  formula on `dirichlet()`. Actual: "Error in use_glm_primitive(tmp) &&
  all.equal(tmp$dpars$mu$fe, bterms$dpars[[1]]$fe) : invalid 'y' type
  in 'x && y'" from `use_glm_primitive_categorical()`: `all.equal()`
  returns a string when the two designs differ, and `&&` refuses it.
- Workaround: frmtmb fits it.
- Correction (lane catmulti, 2026-10-14; `dev/catmulti-brms-b.R`, log
  `brms-b.txt`): it fires for any category with a formula of its own,
  not only the first predicted one. `mua ~ z` and `mub ~ z` under
  `refcat = "c"`, and `muc ~ z` under the default reference, all stop
  with the same error. `dev/famsimplex-findings.md`, section 6, says
  "the FIRST predicted category".
- Record: `dev/famsimplex-findings.md`, section 6, brms-B;
  `dev/catmulti-findings.md`, section 13.
- Status: not reported.

### brms-27. `posterior_epred()` for `wiener()` is right only at `bs = 1`

- Seen: 2.23.0, 2026-10-13 (lane famuni; `dev/famuni-wiener-findings.md`,
  section 5).
- Repro: `brms:::posterior_epred_wiener()` is
  `ndt - bias / mu + bs / mu * (exp(-2 mu bias) - 1) /
  (exp(-2 mu bs) - 1)`. On a 48-cell grid against the exact
  unconditional mean decision time (frmtmb.eam's conditional means
  weighted by the boundary probabilities) its relative error is at
  most 4.07e-16 at `bs = 1` and up to 0.637 at `bs != 1` (`bs = 2.5`,
  `mu = 1.5`, `bias = 0.7`). With the absolute starting point
  `z = bias * bs` in place of `bias` the formula agrees to 1.04e-15. A
  20000-draw RWiener Monte Carlo agrees with the exact value on four
  cells (0.44867 against 0.44902, 0.79461 against 0.79217). The
  formula is also `NaN` at `mu = 0`.
- Expected: the mean of the Wiener first-passage time. Actual: a mean
  for `bs = 1` with the starting point read as relative, so every fit
  with a boundary separation other than 1 gets a wrong expected
  response time.
- Workaround: frmtmb.eam's `fitted()` is conditional on the observed
  boundary and does not use this formula; `?wiener` documents the
  difference.
- Record: `dev/famuni-wiener-findings.md`, sections 5 and 9.
- Status: not reported.

### brms-28. `logistic_normal()` in `mvbf()` does not compile

- Seen: 2.23.0, 2026-10-14 (lane catmulti, brms-C;
  `dev/catmulti-brms-fits.R`, case "mv-ln"; the Stan code in its log
  `brms-probe.txt`).
- Repro: `brm(bf(Y ~ x, family = logistic_normal()) + bf(w ~ x, family
  = gaussian()) + set_rescor(FALSE), data = dd)`, `Y` a three-column
  composition.
- Expected: a fit. Actual: "Semantic error ... Identifier "Llncor" not
  in scope. Did you mean "Lncor"?". The generated quantities read
  `Llncor` and `nlncor`, where the parameters and the transformed data
  declare `Llncor_Y` and `nlncor_Y` (`stan-predictor.R`: the glue of
  `gen_def` has no `{p}` suffix). With two `logistic_normal()`
  responses the generated quantities also declare `Lncor` and `lncor`
  twice, and `rename_family_cor_pars()` would give both responses'
  correlations the one name `lncor__b__c`.
- Workaround: frmtmb fits the model. Its test repairs the two
  identifiers and checks brms's model block (lines 135 and 136 of 145
  change).
- Record: `dev/catmulti-findings.md`, sections 8 and 13.
- Status: not reported.

### brms-29. `update()` without new data computes a data-dependent basis again on fewer rows

- Seen: 2.23.0, 2026-10-19 (lane drop1raw, review finding 2;
  `dev/drop1raw-rev-brms2.R`, log `brms2.txt`, run on an empty fit with
  a stub stanfit, so `update()` runs without a compile).
- Repro: a fit of `z ~ poly(x2, 2) + x3` on 40 rows (seed 5), 2 of
  them `NA` in `x3`. The fit's `poly()` column is computed on all 40
  rows before `na.action`; `object$data` keeps the 38 rows after it, and
  `update(fit)` computes `poly()` again on those 38. The script prints
  "fit column == update column: FALSE".
- Expected: an update with no change reproduces the fit's design.
  Actual: a different basis, so the refit is a different model (in
  frmtmb's measurement of the same defect, `ns(x2, 3) + x3` moved the
  log-likelihood by 1.5e-2).
- Workaround: frmtmb keeps the fit's bases (`predvar_map`) on every
  refit from the fit's own data (`update()` without `newdata`,
  `drop1()`, `influence()`), and so does not follow brms here.
- Record: `dev/reviews/2026-10-19-drop1raw.md`, finding 2;
  `dev/drop1raw-findings.md`, fix round 1.
- Status: not reported.

### brms-30. `loo_epred(ndraws =)` weights one subsample's prediction by another's likelihood

- Seen: 2.23.0, 2026-10-20 (lane loodraws; `dev/loodraws-brms-ndraws.R`,
  log `dev/loodraws-log/brms-ndraws.txt`).
- Repro: a gaussian `y ~ x` fit, 60 rows (seed 1), 4 chains of 1000.
  `set.seed(42); loo_epred(fit, ndraws = 100)`. A trace on
  `brms:::validate_draw_ids()` shows two calls, one from `log_lik()`
  inside `loo()` and one from `posterior_epred()`, each with its own
  `sample()`: the two subsamples share 2 of 100 draws. The result equals
  `loo::E_loo()` of `posterior_epred()` at the second subsample under
  the PSIS weights of `log_lik()` at the first (max abs diff 0).
- Expected: one subsample for both halves. Actual: the weights belong to
  other draws than the prediction they weight; the same-draw answer
  differs by up to 8.7 percent of its largest value. `loo_linpred()`,
  `loo_predict()` and `loo_R2()` pass `...` to both halves the same way.
- Workaround: frmtmb.sample's `loo_*()` functions choose the draws once
  (`ndraws` is an evenly spaced subsample) and give the same rows to
  `log_lik()` and the prediction.
- Record: `dev/loodraws-findings.md`.
- Status: not reported.

### brms-31. `zero_inflated_asym_laplace_lccdf()` counts the point mass at `y = 0`

- Seen: 2.23.0, 2026-10-20 (lane cdfs, source read:
  `inst/chunks/fun_zero_inflated_asym_laplace.stan`).
- Repro: the Stan function takes the `y <= 0` branch at `y = 0`,
  `log_sum_exp(log(zi), log(1 - zi) + asym_laplace_lccdf(0))`, which is
  `P(Y >= 0)`; `zero_inflated_asym_laplace_lcdf(0)` is `P(Y <= 0)`, so
  the two sum to `1 + zi` at 0. brms's own R function
  `pzero_inflated_asym_laplace(0, lower.tail = FALSE)` gives
  `P(Y > 0)`, without the mass.
- Expected: `P(Y > y)` at every `y`, as every other `_lccdf` is. A row
  right-censored at exactly 0 is scored `zi` too high.
- Workaround: frmtmb's `zi_asym_laplace_cdf_pair()` is `P(Y > y)`.
- Record: `dev/cdfs-findings.md`.
- Status: not reported.

### brms-32. A truncated count's `posterior_epred()` drops the lower bound

- Seen: 2.23.0, 2026-10-20 (lane cdfs, `dev/cdfs-brms-tmean.R`, log
  `dev/cdfs-out/brms-tmean.txt`).
- Repro: `brms:::posterior_epred_trunc_discrete("pois",
  list(lambda = matrix(3, 2, 2)), lb = matrix(2, 2, 2), ub = matrix(6,
  2, 2))` is 3.9587628866, the mean over `(2, 6]`; the likelihood brms
  writes for `trunc(lb = 2, ub = 6)` is `poisson_lcdf(ub) -
  poisson_lcdf(lb - 1)`, the window `[2, 6]`, whose mean is
  3.3868613139.
- Expected: the mean of the distribution the model fits. The same sum
  serves binomial, negbinomial, negbinomial2 and geometric.
- Workaround: frmtmb's truncated means use the inclusive window
  (`nbinom_trunc_mean()`, `binomial_trunc_mean()`, poisson's
  `trunc_mean_fn`).
- Record: `dev/cdfs-findings.md`.
- Status: not reported.

### brms-33. `inv_gaussian_lccdf()` adds the second term it should subtract

- Seen: 2.23.0, 2026-10-20 (lane cdfs, source read of
  `inst/chunks/fun_inv_gaussian.stan`; `dev/cdfs-brms-ig.R`, log
  `dev/cdfs-out/brms-ig.txt`).
- Repro: the function is `log1m(Phi(a) - exp(2 shape / mu) Phi(-b))`,
  which is `log(Phi(-a) + exp(...) Phi(-b))`; the survivor is
  `Phi(-a) - exp(...) Phi(-b)`. At `y = 2`, `mu = 2`, `shape = 4` it
  returns -0.4657 where `statmod::pinvgauss(lower.tail = FALSE,
  log.p = TRUE)` and `log(1 - exp(inv_gaussian_lcdf()))` give -0.9880.
- Expected: `log(1 - F)`. Every right-censored row and every
  `trunc(lb = )` without an upper bound of an inverse.gaussian model is
  scored with the wrong survivor; frmtmb's log-density identity test
  differed from brms's program by 8.99 nats on 150 rows until its
  right-censored rows were dropped.
- Workaround: none needed; frmtmb's survivor is right
  (`inverse_gaussian_cdf_log()`), so a right-censored inverse.gaussian
  fit differs from brms's.
- Record: `dev/cdfs-findings.md`, section 6;
  `tests/testthat/test-brms-lik-cdfs.R`.
- Status: not reported.

### brms-34. A censored row under `trunc()` is the event over the whole line, divided by the window

- Seen: 2.23.0, 2026-10-20 (the review of lane cdfs,
  `dev/cdfs-rev-ctcode.R`, `dev/cdfs-rev-out/r1/ctcode.txt` and
  `conv.txt`; `stancode()` of `y | cens(cc) + trunc(lb = 2)` in lane
  cdfs, fix round 2).
- Repro: `brms::stancode(bf(y | cens(cc) + trunc(lb = 2, ub = 9) ~ x),
  family = poisson())` scores a right-censored row as
  `poisson_lccdf(Y | mu) - log_diff_exp(poisson_lcdf(ub | mu),
  poisson_lcdf(lb - 1 | mu))` and a left-censored one as
  `poisson_lcdf(Y | mu) - (the same)`; the gaussian is written the
  same way, and with `trunc(lb = )` alone a left-censored row is
  `poisson_lcdf(Y | mu) - poisson_lccdf(lb - 1 | mu)`.
- Expected: the probability of the event inside the window,
  `P(y < Y <= ub)` and `P(lb <= Y <= y)`, over the window's mass. As
  written the numerator includes mass the truncation excludes, so a
  censored row's "probability" can exceed one. With censored rows in
  `trunc(lb = 2, ub = 9)`, brms's log density differs from frmtmb's by
  6.16 on a poisson and 14.85 on a negative binomial. Right censoring
  with no upper bound, left censoring with no lower bound, and interval
  rows are the same in both.
- Workaround: none needed; frmtmb's composition is the correct one
  (`row_lpdf()`, `log_window()`), so such a fit differs from brms's.
- Record: `dev/cdfs-findings.md`, Fix round 2;
  `tests/testthat/test-cens-trunc.R`.
- Status: not reported.

### brms-35. `rdirichlet()` underflows at a small concentration

- Seen: 2.23.0, 2026-10-20 (`dev/nb2dm-brms-rdirichlet.R`, log
  `dev/nb2dm-out/brms-rdirichlet.txt`).
- Repro: `rdirichlet(20000, alpha = phi * p)` with
  `p = softmax(c(0, 0.3, -0.2))`, seed 5: 18 rows NaN at phi = 0.01 and
  9570 at phi = 1e-3. `posterior_predict_dirichlet_multinomial()` on a
  prep with phi = 1e-3 stops with "NA in probability vector";
  `posterior_predict_dirichlet()` returns 973 NaN rows of 2000.
- Expected: compositions. Actual: every gamma variate of a row
  underflows to 0 and the row is 0 / 0. A small phi is an ordinary
  estimate for strongly overdispersed counts.
- Workaround: frmtmb draws the gamma variates on the log scale,
  `log Gamma(a + 1) + log(U) / a`, then a softmax
  (`simplex_rdirichlet()`), for dirichlet(), dirichlet2 and
  dirichlet_multinomial().
- Record: `dev/reviews/2026-10-20-nb2dm.md`, finding 1;
  `dev/nb2dm-findings.md`, fix round 1.
- Status: not reported.

### brms-36. An interval-censored `cox()` row scores the interval as probability 0

- Seen: 2.23.0, 2026-10-24 (`dev/bhaz-brms-icens.R`, log
  `dev/bhaz-log/brms-icens.txt`; no compile).
- Repro: `stancode(bf(t | cens(c2, u) + bhaz(df = 6) ~ x), family =
  cox())` on 30 rows (seed 1), half with `c2 = 2`. The interval term is
  `log_diff_exp(cox_log_lcdf(rcens[n] | mu[n], bhaz[n], cbhaz[n]),
  cox_log_lcdf(Y[n] | mu[n], bhaz[n], cbhaz[n]))`. `cox_log_lcdf()`
  reads only `cbhaz`, and `cbhaz[n]` is the cumulative baseline hazard
  at `Y`: `standata()` passes `Zcbhaz` for `Y` alone, with no basis at
  `rcens`. On row 2 (Y 1.2816, rcens 1.9225, equal weights) both ends
  give -1.019042, and the interval's log probability is -Inf.
- Expected: the cumulative baseline hazard at `rcens` for the upper
  end. Actual: every interval-censored row has log density -Inf, so
  brms cannot fit an interval-censored Cox model.
- Workaround: none needed. frmtmb.survival's `cox()` builds the
  I-spline basis at each censoring and truncation bound
  (`cox_cb_entry()`), and its interval rows are scored at both ends.
- Record: `dev/reviews/2026-10-24-bhaz.md`, "Outside this lane";
  `dev/bhaz-findings.md`, fix round 1.
- Status: not reported.

### brms-37. `sar()` passes complex eigenvalues of a row-standardized matrix to Stan

- Seen: 2.23.0, 2026-10-24 (`dev/sarfcor-upstream-eigen.R`, log
  `dev/sarfcor-out/upstream-eigen.txt`).
- Repro: an 8 x 8 rook lattice, `W <- A / rowSums(A)`,
  `brm(y ~ x + sar(W), data2 = list(W = W))`. `eigen(W)` returns its
  values as complex, max |Im| 7.97e-17, and `standata()` keeps
  `eigenMsar` complex.
- Expected: the real parts. `W` is similar to the symmetric
  `D^-1/2 A D^-1/2`, so its eigenvalues are real. Actual: Stan stops in
  data initialization, "variable does not exist ... variable
  name=eigenMsar; base type=double", and `brm()` returns a fit with no
  draws. A row-standardized matrix is spdep's default and the usual
  `sar()` input.
- Workaround: frmtmb.spatial takes the real part when the imaginary
  part is below 1e-8 of the largest modulus, and refuses a matrix whose
  eigenvalues are complex beyond that. The brms rows of frmtmb.spatial's
  tests use the symmetric normalization.
- Record: `dev/sarfcor-findings.md`; `dev/reviews/2026-10-24-sarfcor.md`,
  finding 7.
- Status: not reported.

### brms-38. `sar()` with `weights()`, `se()`, `cens()` or `trunc()` writes Stan code that does not parse

- Seen: 2.23.0, 2026-10-24 (`dev/sarfcor-upstream-stanc.R`, log
  `dev/sarfcor-out/upstream-stanc.txt`).
- Repro: `stancode(y | weights(w) ~ x + sar(W), ...)` and the same with
  `se(se, sigma = TRUE)`, `cens(cen)` or `trunc(lb = lb)`, on a
  row-standardized 4 x 4 lattice, then
  `rstan::stanc(model_code = code)`.
- Expected: a refusal from brms, as for the other residual correlation
  terms, or a program that compiles. Actual: stanc gives "Semantic
  error" at the `normal_lagsar_lpdf` call for each of the four; the
  plain `y ~ x + sar(W)` type-checks. The user sees a Stan compiler
  error, not a message about the model.
- Workaround: frmtmb refuses each of the four with the term and the
  reason (the density is joint over the rows).
- Record: `dev/sarfcor-findings.md`; `dev/reviews/2026-10-24-sarfcor.md`,
  finding 7.
- Status: not reported.

### brms-39. `posterior_epred()` and `fitted()` fail on every `sar(type = "error")` fit

- Seen: 2.23.0, 2026-10-24 (`dev/sarfcor-upstream-epred.R`, log
  `dev/sarfcor-out/upstream-epred.txt`).
- Repro: a 5 x 5 lattice with the symmetric normalization,
  `fit <- brm(y ~ x + sar(Ws, type = "error"), ...)`, then
  `posterior_epred(fit)` or `fitted(fit)`.
- Expected: `mu`, since the error form does not move the mean. Actual:
  "'a' is 0-diml". `prepare_predictions_ac()` reads `lagsar` with
  `prepare_draws()`, which gives a matrix with no columns (200 x 0)
  rather than `NULL` when the fit has no `lagsar`, so
  `posterior_epred_gaussian()` takes the lag branch and `solve()` stops.
  `posterior_predict()` and `log_lik()` work.
- Workaround: none needed; frmtmb's `fitted()` gives `mu` for the error
  form. frmtmb.spatial's `posterior_epred()` comparison runs on the lag
  form only.
- Record: `dev/sarfcor-findings.md`, fix round 1 (F3).
- Status: not reported.

## RTMB

### RTMB-1. `pbeta()`, `pbinom()` and `pnbinom()` have non-finite third derivatives

- Seen: 2.0. Confirmed 2026-09-29.
- Repro:

  ```r
  F <- RTMB::MakeTape(function(p) log(RTMB::pbinom(4, 12, plogis(p[1]))), 0)
  F$jacfun()$jacfun()$jacfun()(0)   # -Inf
  ```

  Sweep (`frmtmb-wt-fams2/dev/fams2-rev-pthird.R`, seed 7): 47/2000
  pbeta, 96/472 pbinom, 196/2000 pnbinom points non-finite.
  Also: `pbeta(x, a, b)` at `x == a/(a+b)` exactly gives a NaN
  GRADIENT (and Hessian, third derivative), e.g. x = 0.3, a = 300,
  b = 700; finite 1e-12 away.
- Expected: finite derivatives (a Laplace fit needs up to the third).
- Workaround: xbeta uses `log_ibeta_half()`
  (`frmtmb-wt-fams2/R/families.R`), which calls `pbeta` only near the
  mean at large shapes, blended with the shifted form
  I_x(a+1, b) + x^a(1-x)^b / (a B(a, b)) whose tie sits at
  (a+1)/(a+b+1) (`log_pbeta_ad()`, fams2 punch round 2); no family
  calls the three otherwise.
- Record: `frmtmb-wt-fams2/dev/reviews/2026-09-29-fams2.md` claim 2 and
  the re-check's n1.
- Status: not reported.

### RTMB-2. `dbeta(log = TRUE)` has a NaN gradient at large shapes

- Seen: 2.0. Confirmed 2026-09-29.
- Repro:

  ```r
  D <- RTMB::MakeTape(function(p) RTMB::dbeta(p[1], exp(p[2]), exp(p[3]), log = TRUE), c(0.3, 0, 0))
  D$jacobian(c(0.3, log(300), log(700)))   # finite
  D$jacobian(c(0.3, log(3e3), log(7e3)))   # NaN NaN NaN
  ```

- Expected: a finite gradient (the value is right). Under `MakeTape`
  a constant `x` fails too; under `MakeADFun` only a taped `x` fails.
- Workaround: xbeta writes the log density out with `lbeta_ad()`
  (`frmtmb-wt-fams2/R/families.R:2976`).
- Record: `frmtmb-wt-fams2/dev/reviews/2026-09-29-fams2.md` B1;
  `frmtmb-wt-fams2/dev/fams2-rev-nan.R` section 3.
- Status: not reported.

### RTMB-3. Derivatives of `pnorm(x, log.p = TRUE)` lose accuracy in the lower tail

- Seen: 2.0 (TMB 1.9.25). Confirmed 2026-09-29.
- Repro: `MakeTape(function(x) pnorm(x, log.p = TRUE), 0)`; see
  `dev/phase3b-rtmb-report-pnorm.md`, "Reproduction".
- Expected: second derivative -1 at `x = -1e5`. Actual: 3169.6; at
  `-1e7` the first derivative is 9990923 (true 1e7).
- Workaround: `ddm_rt_u_floor <- 1e-10`
  (`extensions/frmtmb.eam/R/wiener-rtcdf.R:68`).
- Also (lane optima, 2026-10-07, the review's
  `dev/optima-rev-links.R`): the FIRST derivative of `pnorm(x, log.p =
  TRUE)` against the Mills-ratio series degrades smoothly from about
  `|x| = 1e3`: relative error 4.8e-11 at -1e3, 1.3e-9 at -1e4, 2.0e-5
  at -1e6, 1.5e-3 at -9.4e6, 0.34 at -1e8, and the derivative is
  exactly 0 from -1e9 on; it is NaN on 137 of 2001 points within 0.1
  percent of -4.46e9 (`dev/optima-pnorm-scan.R`). The value is exact
  to 1e154. Workaround: the probit's `logit_eta` reads the logs at
  `|eta| <= 1000` only and continues with the tail's leading term
  (`R/links.R`).
- Record: `dev/phase3b-rtmb-report-pnorm.md` (filing draft);
  `dev/optima-findings.md`.
- Status: not reported (draft ready).

### RTMB-4. `dnbinom_robust()` gradient in the size is wrong at very large sizes

- Seen: 2.0. Confirmed 2026-09-29.
- Repro:
  `G <- MakeTape(function(p) dnbinom_robust(1, 0, -p[1], log = TRUE), 30)`;
  `G$jacobian(30)` gives 1.2e-2 (25: -5.6e-5).
- Expected: about 1e-13 (the Poisson limit).
- Workaround: none; the robustness test does not sweep `shape`
  (`tests/testthat/test-numerical-robustness.R:136`).
- Record: `dev/fams-findings.md` "Defects found and not fixed" item 2.
- Status: not reported.

### RTMB-5. `OBS()` in a loop silently swaps data between responses

- Seen: 2.0. Confirmed 2026-09-29.
- Repro:

  ```r
  y <- list(c(1, 2, 3), c(10, 20, 30))
  f <- function(p) { nll <- 0
    for (r in 1:2) nll <- nll - sum(dnorm(OBS(y[[r]]), p$mu[r], 1, log = TRUE)); nll }
  obj <- MakeADFun(f, list(mu = c(0, 0)), silent = TRUE)
  nlminb(obj$par, obj$fn, obj$gr)$par   # 2 2; without OBS() 2 20
  ```

- Expected: one observation per call, or an error on a duplicate key.
  Actual: `OBS()` keys by the deparsed expression (`"y[[r]]"`).
- Workaround: `OBS()` only for a univariate, non-matrix response,
  through a fixed symbol (`R/objective.R:486`).
- Record: `dev/rtmb-pitfalls.md` item 8.
- Status: not reported.

### RTMB-6. The AD branch of `dt()` loses its digits at large `nu`

- Seen: 2.0. Confirmed 2026-09-29.
- Repro: `T <- MakeTape(function(p) RTMB::dt(0.5, exp(p[1]), log = TRUE), 0)`;
  `T$jacobian(log(1e6))` gives +3.6e-7 (true about -3.6e-7);
  `T(log(1e10))` is -1.043945 against `stats::dt()` -1.043939.
- Expected: the accuracy of the double branch (`stats::dt()`). The AD
  branch forms `lgamma((nu + 1) / 2) - lgamma(nu / 2)` as written.
- Workaround: `dt_stable()` with `lgamma_shift_diff()`
  (`R/families.R:1786`).
- Record: `dev/remlopt-findings.md`, near line 8; `NEWS.md` near 2048.
- Status: not reported.

### RTMB-7. `ppois(lower.tail = FALSE, log.p = TRUE)` is not taped

- Seen: 2.0. Confirmed 2026-09-29.
- Repro: `MakeTape(function(p) RTMB::ppois(3, exp(p[1]), lower.tail = FALSE, log.p = TRUE), 0)`:
  "Non-numeric argument to mathematical function".
- Expected: a taped upper tail, as for the lower one.
- Workaround: through 0.72.0 `poisson()` declared no `lccdf`; since
  lane cdfs (round of 2026-10-20) its survivor is the lower incomplete
  gamma, `log_igamma_two()` in `R/families-cdf.R`. `pgamma()`'s upper
  tail and `log.p` do not tape either, the same error (lane cdfs,
  `dev/cdfs-probe1.R`, log `dev/cdfs-out/probe1.txt`).
- Record: `dev/cdfs-findings.md`, section 3.
- Status: not reported.

### RTMB-8. `pbeta()`'s derivative in a shape is 0 or NaN at a whole shape

- Seen: 2.0, TMB 1.9.25, 2026-10-20 (lane cdfs, `dev/cdfs-probe2.R`).
- Repro:

  ```r
  F <- RTMB::MakeTape(function(p) RTMB::pbeta(0.3, p[1], 5), 2)
  F$jacobian(2)          # 0
  F$jacobian(2.0001)     # -0.2524
  ```

  The same at `b = 30` in the second shape, and through
  `pnbinom(prob =)` in `size`, which is `pbeta()` underneath. The
  review of 2026-10-20 (`dev/cdfs-rev-upstream.R`,
  `dev/cdfs-rev-out/upstream.txt`) mapped the cases at `x = 0.3` and
  `x = 0.4`, AD against finite differences in `(a, b)`:
  - a whole first shape: the derivative in it is 0 and the other is
    right, `(2, 5.5)` gives `0, 0.0907` against `-0.2425, 0.0907`;
  - a whole second shape: the reverse, `(2.5, 5)` gives `-0.2285, 0`
    against `-0.2285, 0.1026`;
  - both whole and small: 0 and right, `(2, 5)` gives `0, 0.0994`
    against `-0.2524, 0.0994`;
  - both whole and large: `NaN, NaN`, `(200, 300)` at `x = 0.4` against
    `-0.02187, 0.01458`.
- Expected: the derivative, -0.25242 at `a = 2` by finite differences
  of `stats::pbeta()`. The value is right.
- Workaround: frmtmb's CDFs read the incomplete beta through
  `log_ibeta_half()` (`R/families.R`), which calls `pbeta()` only at
  shapes clamped to 150 or more and blended; a whole `nu / 2` is
  frmtmb's own start for the student-t (`nu = 4`).
- Record: `dev/cdfs-findings.md`, section 3.
- Status: not reported.

### RTMB-9. `pgamma()`'s third derivatives are not finite at half of all points

- Seen: 2.0, 2026-10-20 (lane cdfs, `dev/cdfs-dbg9.R`).
- Repro: `T3 <- MakeTape(function(p) log(RTMB::pgamma(exp(p[1]),
  shape = exp(p[2]))), c(0, 0))$jacfun()$jacfun()`; non-finite at 960
  of 2000 uniform points of `log x` in (-3, 4) and `log shape` in
  (-2, 6), mostly in the lower tail. The value and first two
  derivatives are right; `ppois()`'s third derivatives are finite at
  500 of 500 points (`dev/cdfs-dbg13.R`). The review of 2026-10-20
  found the SECOND derivatives non-finite too, at 52 of 500 points,
  all in the lower half (`P < 1/2`), and the third at 238 of 500
  (`dev/cdfs-rev-rtmb9.R`, `dev/cdfs-rev-out/rtmb9.txt`), so a Hessian
  and its standard errors are reached, not only a Laplace fit.
- Expected: finite third derivatives (a Laplace fit needs them).
- Workaround: `log_igamma_two()` (`R/families-cdf.R`) does not call
  `pgamma()`: a series, a continued fraction and Gauss-Legendre
  quadrature, blended.
- Record: `dev/cdfs-findings.md`, section 3.
- Status: not reported.

## RTMBdist

### RTMBdist-1. The inverse-Gaussian upper tail is computed on the probability scale

- Seen: 1.1.0. Not rerun.
- Repro: the log upper tail reaches `-Inf` at log S = -34, the same as
  `log(1 - F)`.
- Expected: a log-scale upper tail.
- Workaround: through 0.72.0 `inverse.gaussian()` declared no `lccdf`;
  since lane cdfs both tails are formed in logs by
  `inverse_gaussian_cdf_log()` (`R/families-cdf.R`).
- Record: `dev/cdfs-findings.md`; `NEWS.md` near line 3383.
- Status: not reported.

### RTMBdist-2. `pinvgauss()` is NaN where `exp(2 shape / mean)` overflows

- Seen: 1.1.0, 2026-10-20 (lane cdfs, `dev/cdfs-ig.R`).
- Repro: `RTMBdist::pinvgauss(c(0.5, 1, 2), mean = 1, shape = 1000)` is
  `NaN NaN NaN`; `statmod::pinvgauss()` gives 6.34e-111, 0.506, 1.
- Expected: the CDF. The factor is `Inf` and the normal tail it
  multiplies is 0, from `shape / mean = 355`.
- Workaround: `inverse_gaussian_cdf_log()` adds the two terms in logs.
- Record: `dev/cdfs-findings.md`.
- Status: not reported.

### RTMBdist-3. `pskewnorm()` and `pvm()` do not tape

- Seen: 1.1.0, 2026-10-20 (lane cdfs, `dev/cdfs-probe1.R`).
- Repro: `MakeTape(function(p) RTMBdist::pskewnorm(p[1], 0, exp(p[2]),
  p[3]), c(0.4, 0.1, 2))`: "Can only determine finite status of
  constant expressions"; `RTMBdist::pvm(1, p[1], exp(p[2]))` on the
  tape: "Comparison is generally unsafe for AD types".
- Expected: taped CDFs, as their densities are.
- Workaround: `skew_normal_cdf_log()` (Owen's T as an angle integral)
  and `von_mises_cdf_log()` (the tail as a quadrature of the density)
  in `R/families-cdf.R`.
- Record: `dev/cdfs-findings.md`, section 3.
- Status: not reported.

### RTMBdist-4. `dexgauss()` is a gaussian approximation from `lambda sigma` of about 20

- Seen: 1.1.0, 2026-10-20 (the review of lane cdfs,
  `dev/cdfs-rev-exgdens2.R`; frmtmb's exposure,
  `dev/cdfs-fix1-exgfit.R`).
- Repro: `RTMBdist::dexgauss(0, mu = 0, sigma = 1, lambda = 20, log =
  TRUE)` is -0.9189385, the standard normal's log density at 0; the
  closed form `log(lambda) + lambda^2 sigma^2 / 2 - lambda y +
  log Phi(y / sigma - lambda sigma)` is -0.9214231. Against
  `brms::dexgaussian()` the largest error of the log density on a grid
  of `y` is 2.8e-14 at `sigma / beta = 15` and 0.162 at 20.
- Expected: the density. The switch also puts a kink in the
  likelihood at `lambda sigma = 20`.
- Effect in frmtmb through 0.72.0, whose `exgaussian()` density was
  `dexgauss()`: on 400 gaussian rows the objective is 0.63 off at
  `sigma / beta = 50`, and three of four near-gaussian fits stopped at
  `sigma / beta = 20` with nlminb's false convergence, where the
  closed form's maximum is at 420 to 700 and 0.0008 to 0.0024 higher
  (`dev/cdfs-out/fix1-exgfit-r10.txt`).
- Workaround: `exgauss_lx()` (`R/families-cdf.R`), the closed form as
  `log phi(z) + log M(u - z)` with Mills's ratio from a continued
  fraction, for the density and the CDF.
- Record: `dev/cdfs-findings.md`, Fix round 1.
- Status: not reported.

## TMB

### TMB-1. The taped second derivative of `tanh(x)` is NaN for |x| >= 711

- Seen: 1.9.25 (through RTMB 2.0). Confirmed 2026-09-29.
- Repro: `MakeTape(function(x) tanh(x), 0)$jacfun()$jacobian(711)`
  gives NaN (710 gives 0). Source: `TMB/include/TMBad/global.hpp`,
  `struct TanhOp`, `reverse()`, line 3534: `1 / (cosh(x) * cosh(x))`.
- Expected: 0. Fix: `Type(1.) - args.y(0) * args.y(0)`.
- Workaround: `ddm_tanh_s()` clamps at 40
  (`extensions/frmtmb.eam/R/wiener-rtcdf.R:393`).
- Record: `dev/phase3b-rtmb-report-tanh.md` (filing draft).
- Status: not reported (draft ready).

### TMB-2. The CRAN macOS binaries of TMB and RTMB fail to load on the GitHub macOS runner

- Seen: CRAN binaries after RTMB 2.0 reached CRAN (2026-09-17). Not
  checked since.
- Repro: `macos-latest`, `r-lib/actions/setup-r@v2`, binary install,
  `library(RTMB)`: "symbol not found in flat namespace
  '___kmpc_for_static_fini'".
- Expected: the binary loads. `setup-r` removes `-fopenmp`, and the
  binaries were built with it. Find out whether CRAN or `setup-r` is
  at fault before you file.
- Workaround: the macOS job builds TMB and RTMB from source
  (`.github/workflows/R-CMD-check.yaml:36`).
- Record: that workflow comment; `dev/round-handoff.md` "Upstream
  reports" (no draft text exists).
- Status: not reported.

### TMB-3. `marginal_gk` is one-dimensional and freezes its rescaling at tape time

- Seen: 1.9.25 (TMBad). Not rerun.
- Repro: a Gauss-Kronrod marginal (`quadrature = TRUE` in frmtmb) over
  a random-effect block of dimension 2 returns `NaN`; the `dim` in the
  spec is never read. Taped at a cold start, the frozen rescaling gave
  NA modes and "NA/NaN gradient" crashes.
- Expected: `dim` honored or refused; rescaling that follows the
  parameters, or a documented warning.
- Workaround: tape at the Laplace optimum; the importance-sampling
  correction for dimension 2 and up (`R/importance.R:5`).
- Record: `NEWS.md` near lines 3699 and 5824; `dev/test-backlog.md`
  near line 505.
- Status: not reported.

### TMB-4. The ADFun finalizer fails while `isNullPointer` is a lazy-load promise

- Seen: TMB 1.9.25 with RTMB 2.0, 2026-10-19 (lane illse2 and its
  review; `dev/illse2-rev-tmbfin.R`).
- Repro: the finalizer that `TMB:::registerFinalizer()` puts on an
  ADFun is `function(ptr) if (!isNullPointer(ptr)) .Call(...)`.
  `isNullPointer` is lazy-loaded. When a garbage collection runs the
  first finalizer of a session, that finalizer forces the promise, and
  every other ADFun finalizer pending in the same collection meets the
  promise under evaluation. Each prints "Error in (function (ptr) :
  promise already under evaluation" to the console and does not free
  its tape. It happens at most once per session. In frmtmb it showed on
  a 159-parameter fit that makes many tapes (`dev/illse2-p4wide.R 150
  2000 1:2`: three lines during the second fit).
- Expected: no error, every tape freed. The finalizer could call
  `.Call("isNullPointer", ...)` directly, or TMB could force the
  binding when it loads.
- Workaround: frmtmb's `.onLoad()` forces `TMB:::isNullPointer` before
  any tape exists (`R/zzz.R`); with that, the same fit prints nothing
  (`dev/round-20261020.md`, consolidation item 2).
- Record: `dev/illse2-findings.md`, section 4 and "Defects found, not
  fixed"; `dev/reviews/2026-10-14-illse2.md`, item 6 of the re-check.
- Status: not reported.

## tmbstan

### tmbstan-1. tmbstan 1.2.0 built against StanHeaders 2.39 samples a standard normal

- Seen: 1.2.0 with StanHeaders 2.39.1. Fixed in 1.2.1 (`089bc18`,
  `tools/autogen.R`).
- Repro: `TMB:::runExample("simple")`, then `tmbstan(obj)`; the
  generated `model.hpp` has `std_normal_lpdf<propto__>(y)`.
- Expected: draws from the model. Actual: N(0, 1), silently.
- Workaround: `tmbstan_build_broken()`, `check_tmbstan_build()`
  (`extensions/frmtmb.sample/R/sample.R:2483`); CI needs tmbstan >= 1.2.1.
- Record: `dev/prior-dropping-investigation.md`;
  `dev/tmbstan121-findings.md`.
- Status: reported by the user, kaskr/tmbstan#33, closed. Nothing to file.

## rstan and StanHeaders

### rstan-1. rstan 2.32.7 cannot compile a Stan program against StanHeaders 2.39.1

- Seen: rstan 2.32.7, StanHeaders 2.39.1 (the pair CRAN and RSPM
  install). Check on the next rstan release.
- Repro: `rstan::stan_model()` with no `CXX17FLAGS` line in Makevars:
  "'void_t' is not a member of 'std'".
- Expected: C++17. Actual: rstan emits `Rcpp::plugins(cpp14)`, whose
  `-std=c++1y` lands last. rstan declares `StanHeaders (>= 2.32.0)`
  with no upper bound.
- Workaround: `CXX17FLAGS = -O2 -Wall -std=gnu++17` in Makevars
  (`.github/workflows/brms-likelihood.yaml:125`; `dev/lane-rules.md`).
- Record: `dev/brms-likelihood-tests.md` (rstan compile section);
  `dev/machine-library.md` "The StanHeaders trap".
- Status: not reported.

### rstan-2. `sampling()` returns an empty stanfit on an unknown `control` name

- Seen: 2.32.7 (through `tmbstan()`). Not rerun.
- Repro: `sampling(..., control = list(nonsense_option = 1))` prints
  "error in specifying arguments; sampling not done" and returns.
- Expected: an error. A message request.
- Workaround: `check_stan_control()`
  (`extensions/frmtmb.sample/R/sample.R:74`).
- Record: `dev/stanctl-findings.md` section 1.5.
- Status: not reported.

### rstan-3. The gradient of `gamma_lcdf()` and `gamma_lccdf()` in the shape is wrong

- Seen: rstan 2.32.7 with StanHeaders 2.39.1 (Stan Math's), 2026-10-20
  (lane cdfs, `dev/cdfs-dbg16.R`).
- Repro: brms's `hurdle_gamma()` program with `cens(cc, y2)` on 200
  rows (the data of `tests/testthat/test-brms-lik-cdfs.R`, seed 104 and
  12), at frmtmb's estimates: `rstan::grad_log_prob()` in the
  unconstrained log shape is -0.1823, central differences of
  `rstan::log_prob()` give 0.000203, and frmtmb's AD and finite
  differences agree with the latter. With one censoring code at a time
  the right-censored rows give -0.00705 against -0.00027 and the
  interval rows 0.0215 against -0.000012; the left-censored rows are
  close (0.000229 against 0.000019). The values agree with frmtmb's to
  4e-13.
- Expected: the derivative of the value Stan computes.
- Not shown: that it biases a fit. The review of 2026-10-20 did not
  reproduce the gradient itself, and its censored Gamma brms fit (2
  chains of 1000 draws) agreed with frmtmb to 0.08 posterior SD with
  Rhat 1.002 (`dev/cdfs-rev-brmsfit.R`).
- Workaround: frmtmb's own CDF; frmtmb's identity test takes that
  gradient by finite differences of brms's log density.
- Record: `dev/cdfs-findings.md`, section 6.
- Status: not reported.

### Stan-1. `neg_binomial_2_log_glm_lpmf()` loses the density at a large precision

- Seen: StanHeaders 2.39.1 (rstan 2.32.7). Confirmed 2026-10-20.
- Repro: one observation, `X` a 1 x 1 zero matrix,
  `neg_binomial_2_log_glm_lpmf({1} | X, 3, [0]', phi)` through
  `rstan::expose_stan_functions()`: -17.0859375 at phi = 1e12 and 0 at
  phi >= 1e16, where `neg_binomial_2_log_lpmf({1} | 3, phi)` and
  `dnbinom()` give -17.0855369232. Already off by 1.8e-7 at phi = 1e8,
  and by 1.2e-9 relative at phi = 1e6 (y 0, eta -1; review of lane
  nb2dm, `dev/nb2dm-rev-stan1.R`).
  On 400 rows at b_Intercept -45.75, b_x 15.82, phi = 1.1e23 it returns
  +2.06e13 where the true value is -37774.33.
- Expected: the value `neg_binomial_2_log_lpmf()` gives. Actual: a
  cancellation in the GLM form; a positive log density of 2e13 is an
  attractor that one of four brms chains entered and could not leave
  (Rhat 1.56, 139 divergences).
- Reaches: brms's `negbinomial()` and `brmsfamily("negbinomial2")`
  whenever the shape (or sigma) has no formula, which is when brms
  writes the GLM primitive. frmtmb does not use it. Workaround:
  `init = 0`, under which the same fit mixed (4 of 4 chains, Rhat
  1.002).
- Record: `dev/nb2dm-findings.md` "brms 2.23.0 / Stan defect met";
  repros `dev/nb2dm-stan-glm-bug.R` (log
  `dev/nb2dm-out/stan-glm-bug.txt`) and `dev/nb2dm-brms-stuck.R`.
- Status: not reported.

## RTMBode

All six reproduce without frmtmb. Analysis, repro scripts, a
three-patch series and issue text are in
`extensions/frmtmb.ode/dev/upstream/` (`rtmbode-issues.md`). Seen on
RTMBode 1.0 (`5242257`) with RTMB 1.9 and deSolve 1.42. Check on RTMB
2.0 before you file.

### RTMBode-1. A failed deSolve call escapes the adjoint node as an R error

- Repro: `repro/03-infparm.R` (`obj$fn()` at an extreme parameter).
- Expected: a `NaN` solution. Actual: "illegal input detected..." or
  "Wrong output length"; under tmbstan the chain aborts at warmup
  iteration 1 and rstan's AD arena stays broken for the session.
- Workaround: `frm_sample()` and `as_tmbstan()` refuse `frm_ode()`
  (`extensions/frmtmb.ode/R/zzz.R:45`).
- Record: `rtmbode-issues.md` sections 1, 5.1; `dev/stanctl-findings.md` 6.1.
- Status: not reported (text and patch 01 ready).

### RTMBode-2. `events` fail on the AD path because the state vector has no names

- Repro: `repro/06-events.R`.
- Expected: events work. Actual: "too many state variables in 'event';
  should be < 0" (`MakeTape(...)$par()` drops the names).
- Workaround: `frm_ode()` splits the solve at event times
  (`extensions/frmtmb.ode/R/ode.R:203`).
- Record: `rtmbode-issues.md` 2.1, 5.2; `dev/ode-feasibility.md` 9.2.
- Status: not reported (patch 02 ready).

### RTMBode-3. `replace` and `multiply` events give wrong gradients

- Repro: `repro/06-events.R`; reachable once RTMBode-2 is fixed.
- Expected: correct gradients or a refusal. Actual: 42 % and 59 %
  gradient error; values are right. `replace` is deSolve's default.
- Workaround: the segmented solve (`extensions/frmtmb.ode/R/ode.R:203`).
- Record: `rtmbode-issues.md` 2.2; `dev/ode-feasibility.md` 9.3.
- Status: `NEWS.md` near line 5340 says "reported upstream", but
  `dev/stanctl-findings.md` section 5 says nothing was filed. Check.

### RTMBode-4. Event times off the output grid break the adjoint node

- Repro: `repro/11-event-times.R` (doses at 5, 11, 17; times every 2).
- Expected: the solution on the requested times. Actual: 16 rows for
  13 times; "EvalOp: Function must return 'real' or 'integer'".
- Workaround: the segmented solve.
- Record: `rtmbode-issues.md` 2.4.
- Status: not reported (patch 01 ready).

### RTMBode-5. The augmented system outgrows lsoda and the gradient is NaN with no message

- Repro: `repro/07-state-ceiling.R`; `dev/ode/probeE8-laplace-limit.R`.
- Expected: a message that names the limit. Actual: at 10 states the
  order-3 system has 84210 equations, lsoda's work array overflows
  (deSolve-1), and `gr` is NaN.
- Workaround: none in code; solve one small system per group.
- Record: `rtmbode-issues.md` 3, 5.3.
- Status: not reported (patch 03 ready).

### RTMBode-6. `forcings` are unusable and `approxfun()` in the dynamics is frozen at t = 0

- Repro: `dev/ode/probeH3-forcings-infusion.R`.
- Expected: forcings supported, or a refusal. Actual: "initforc should
  be loaded..."; `approxfun(...)(t)` gives 464.63 where the truth is
  278.66, because `func2tape()` tapes once at t = 0.
- Workaround: an infusion rides as extra parameters on a segment
  (`extensions/frmtmb.ode/R/ode.R:1845`).
- Record: `rtmbode-issues.md` section 4; `dev/ode-feasibility.md` 9.5.
- Status: not reported. No issue text or patch.

## deSolve

### deSolve-1. lsoda's work-array size overflows R integers with a nonsense message

- Seen: 1.42.
- Repro: `deSolve::ode(rep(1, 46340), c(0, 1), function(t, y, p) list(-y), NULL)`.
- Expected: a refusal that names the limit and the banded or sparse
  options. Actual: "cannot allocate memory block of size 134217728 Tb".
  A message request; the limit itself is structural.
- Workaround: none needed.
- Record: `extensions/frmtmb.ode/dev/upstream/rtmbode-issues.md` 5.4.
- Status: not reported (text ready).

## mgcv

### mgcv-1. `PredictMat()` on an `re` smooth returns the wrong width at an unseen level

- Seen: 1.9.4. Confirmed 2026-09-29.
- Repro:

  ```r
  g <- factor(rep(letters[1:10], each = 10))
  sm <- mgcv::smoothCon(mgcv::s(g, bs = "re"), data = data.frame(g))[[1]]
  nd <- data.frame(g = factor(c("a", "zz"), levels = c(levels(g), "zz")))
  dim(mgcv::PredictMat(sm, nd))   # 2 11; the fit has 10 columns
  ```

- Expected: an error. Actual: a wrong-width matrix, silently.
- Workaround: frmtmb refuses an unseen level before `PredictMat()`.
- Record: `dev/resmooth-findings.md` section 2 and section 8.
- Status: not reported.

### mgcv-2. `PredictMat()` on a `t2()` with an `re` margin gives an unclear error at an unseen level

- Seen: 1.9.4. Confirmed 2026-09-29.
- Repro: as mgcv-1 with `t2(x, g, bs = c("cr", "re"))`.
- Expected: an error that names the level. Actual: "non-conformable
  arguments".
- Workaround: frmtmb refuses first.
- Record: `dev/resmooth-findings.md` section 8.
- Status: not reported. Low priority.

## reformulas and lme4

### reformulas-1. `splitForm()` silently drops fixed terms that mention a special name

- Seen: 0.4.4. Confirmed 2026-09-29.
- Repro: `reformulas::splitForm(y ~ x + I(exp(x2)) + (1 | g), specials = "exp")$fixedFormula`
  gives `y ~ x`. `log(exp(x2))` is dropped too.
- Expected: only a top-level call to a special is a special term.
- Workaround: `sub_specials()` aliases special-named calls
  (`R/parse.R:1338`).
- Record: `dev/rtmb-pitfalls.md` item 9.
- Status: not reported.

### reformulas-2. `x * (1 | g)` is silently fitted as `x + (1 | g)`

- Seen: 0.4.4. Confirmed 2026-09-29.
- Repro: `reformulas::splitForm(y ~ x * (1 | g))` gives `y ~ x` and
  `1 | g`.
- Expected: an error or a warning.
- Workaround: frmtmb refuses the spelling (`R/parse.R:1175`).
- Record: `tests/testthat/test-open-issues.R:22`.
- Status: reported, lme4#196 (harvested from lme4's open issues).

## knitr

### knitr-1. `tangle_block()` cuts a nested `read_chunk()` call at the first `)`

- Seen: under pkgdown 2.2; the regex is in knitr 1.52 (source read
  2026-09-29).
- Repro: a chunk with `knitr::read_chunk(some_fn("file.R"))`, then
  `knitr::purl()`. `tangle_block` uses `"read_chunk\\(([^)]+)\\)"`.
- Expected: the call is parsed. Actual: "unexpected end of input";
  pkgdown's article build dies while `rmarkdown::render()` passes.
- Workaround: `purl = FALSE` on the chunk
  (`vignettes/bayesian-cognitive-modeling.Rmd:10`).
- Record: `dev/rtmb-pitfalls.md` item 19.
- Status: not reported.

## typetracer

Seen on 0.2.5. All four are patched with `assignInNamespace()` in
`dev/autotest-run.R`. The source of 1, 2 and 4 was read again on
2026-09-29 and is unchanged. Record: `dev/autotest-triage.md`
"Upstream workarounds". Status: none reported.

- **typetracer-1.** `insert_counters_in_tests()` pastes Windows paths
  into string literals: "'\U' used without hex digits". Expected:
  forward slashes. Workaround: `dev/autotest-run.R:53`.
- **typetracer-2.** `reload_pkg()` runs `grepl(tempdir(), lib_path)`:
  "Invalid back reference" on Windows. Expected: `fixed = TRUE`.
  Workaround: `dev/autotest-run.R:96`.
- **typetracer-3.** `reload_pkg()` deletes the traced build, so every
  later autotest mutation runs against the installed version, with no
  sign in the output. Expected: keep the build or report the version.
  Workaround: `dev/autotest-run.R:96`.
- **typetracer-4.** `trace_package_tests()` lets one failed test stop
  the whole run (`testthat::test_package()` throws). Expected: record
  and continue. Workaround: `dev/autotest-run.R:123`.

## autotest

### autotest-1. `pass_one_rect_as_other()` stops the run on a matrix parameter

- Seen: 0.1.1.010; the call is still in 0.2.0 (source read 2026-09-29).
- Repro: `do.call(data.frame, x$params[[x$i]], quote = TRUE)` on a
  parameter classified as rectangular that is a matrix.
- Expected: convert or skip. Actual: "second argument must be a list",
  and the run stops.
- Workaround: `dev/autotest-run.R:145`.
- Record: `dev/autotest-triage.md` item 5.
- Status: not reported.

## emmeans

### emmeans-1. `ref_grid()` hides the error a `recover_data()` method raises

- Seen: 2.0.4.
- Repro: `ref_grid()` runs `try(recover_data(...))` and, on an error,
  stops with "Perhaps a 'data' or 'params' argument is needed"; the
  real text goes to stderr only.
- Expected: the method's own message. A message request: emmeans's
  convention is that a method returns a string to refuse.
- Workaround: `recover_data.frmtmb_fit()` catches its errors and
  returns the text (`R/interop.R:709`).
- Record: `dev/emm-findings.md` "How the refusal message was lost";
  `dev/brmsport-findings.md` P2.
- Status: not reported.

## gratia

### gratia-1. `derivatives()` uses one fixed `eps = 1e-7` for the second derivative

- Seen: 0.11.2. Not rerun (gratia is not installed).
- Repro: `gratia::derivatives(m, order = 2, type = "central")` on an
  mgcv fit; `dev/spline-findings.md` section 4 sweeps `eps`.
- Expected: a step near 1e-4 for order 2 (error 2.2e-7). Actual: error
  8.2e-2 at 1e-7, from cancellation that grows as `u / eps^2`.
- Workaround: frmtmb.spline sets `eps` per order and scales it by the
  covariate range (`extensions/frmtmb.spline/R/curve-deriv.R:8`).
- Record: `dev/spline-findings.md` section 4.
- Status: not reported.

## drmTMB

### drmTMB-1. The REML criterion changes when `sigma` gets a random effect

- Seen: 0.7.0.
- Repro: the draft issue in `dev/drmtmb-findings.md` section 2.
- Expected: one REML criterion for nested models. Actual:
  `sigma ~ z + (1 | g)` also integrates `beta_sigma`; the larger model
  reports a REML logLik 4.67 to 5.51 below the smaller one (10 of 10
  seeds).
- Workaround: none needed.
- Record: `dev/drmtmb-findings.md` section 2.
- Status: not reported (draft ready). Check on 0.7.1 or later.

### drmTMB-2. `ar(p = 1)` in a formula fails with an internal error

- Seen: 0.7.0.
- Repro: `drmTMB(bf(y ~ x + ar(p = 1)), data = d)`: "argument "x" is
  missing, with no default" from `stats::ar()`.
- Expected: a named refusal. Workaround: none needed.
- Record: `dev/drmtmb-findings.md` section 2, last paragraph.
- Status: not reported. Add it to the drmTMB-1 report.

## hmmTMB

### hmmTMB-1. A data column named `state` is read as known states, silently

- Seen: 1.1.2.
- Repro: `dev/hmm/probeD3-hmmtmb-knownstate.R`.
- Expected: a message, or a clear statement in the documentation.
  Actual: `llk()` -1247.248 with the column, -1216.403 without.
- Workaround: frmtmb's comparisons drop the column.
- Record: `dev/hmm-feasibility.md` "hmmTMB sharp edge";
  `dev/latent-findings.md` item 3.
- Status: not reported (a documentation issue).

### hmmTMB-2. `MarkovChain$new(ref =)` cannot be used with a formula matrix

- Seen: 1.1.2.
- Repro: `MarkovChain$new(ref = rep(1, K), formula = <matrix>)`.
- Expected: `ref` works, as documented. Actual: `check_args()` refuses
  any formula matrix whose diagonal is not `"."`.
- Workaround: frmtmb compares transition matrices instead.
- Record: `dev/latent-findings.md` near lines 1163 and 1359.
- Status: not reported.

## ordinal and MASS

### ordinal-1. `clm(link = "cauchit")` reports a logLik that differs from its own likelihood

- Seen: 2026.7.26.
- Repro: `frmtmb-wt-fams2/dev/fams2-validate.R` section 6b.
- Expected: equal values. Actual: reported -653.93698; the likelihood
  at clm's estimates is -653.93406.
- Workaround: none needed.
- Record: `frmtmb-wt-fams2/dev/fams2-findings.md` 6b, "Defects found" item 4.
- Status: not reported. Not investigated; confirm the cause first.

### MASS-1. `polr(method = "cauchit")` stops 2.9 nats short and reports convergence

- Seen: 7.3.66.
- Repro: as ordinal-1.
- Expected: the optimum or a nonzero code. Actual: code 0, 2.9 below
  the optimum that `clm()` and frmtmb reach.
- Workaround: none needed.
- Record: `frmtmb-wt-fams2/dev/fams2-findings.md` sections 6, 6b.
- Status: not reported. Check the starting values first.

## rtdists

### rtdists-1. `dlba_norm()` and `plba_norm()` subtract two normal CDFs

- Seen: 0.11-6. Not rerun (a test pins it).
- Repro: `Phi(g) - Phi(h)` in the LBA density and CDF; compare with a
  200-bit Rmpfr reference.
- Expected: a log-space difference. Actual: exactly 0 on 68 of 80
  tail rows, and 2.7e-3 off on 8 rows outside the tail.
- Workaround: `exp(la) * -expm1(lb - la)`
  (`extensions/frmtmb.eam/R/lba.R:411`); test
  `extensions/frmtmb.eam/tests/testthat/test-lba.R:55`.
- Record: `extensions/frmtmb.eam/NEWS.md` near line 1091.
- Status: not reported.

## EMC2

### EMC2-1. `dWald()` and `pWald()` subtract two normal CDFs

- Seen: 3.5.0. Not rerun.
- Repro: `EMC2:::dWald()` and `1 - EMC2:::pWald()` on the grid in
  `dev/rdm-gng-findings.md`.
- Expected: no cancellation. Actual: `dWald()` 7.95e-4 relative off;
  `1 - pWald()` exactly 0 on 35 of 315 rows.
- Workaround: frmtmb.eam's own Wald functions.
- Record: `dev/rdm-gng-findings.md` "Against EMC2";
  `extensions/frmtmb.eam/NEWS.md` near line 833.
- Status: not reported.

## mice

### mice-1. `D1()` with the Reiter df returns NaN at a small number of imputations

- Seen: mice as installed at the time. Not rerun.
- Repro: `mice::D1()` with a finite complete-data df and t < 4, where
  the Reiter term `1 / (t - 4)` turns negative.
- Expected: a valid df (for example the large-sample form).
- Workaround: fall back to the large-sample df (`R/multiple.R:482`).
- Record: `NEWS.md` near line 5719.
- Status: not reported. Confirm the repro in mice first.

## hBayesDM and CmdStan

### hBayesDM-1. hBayesDM 1.2.1 does not install against StanHeaders 2.39.1

- Seen: 1.2.1. Not rerun.
- Actual: its Stan programs use pre-2.33 array syntax, which the
  current Stan rejects. hBayesDM 2.0.0 needs cmdstanr instead.
- Workaround: frmtmb.learn checks against its own Stan programs.
- Record: `dev/learn2-findings.md`, near line 741.
- Status: not reported. Probably fixed by 2.0.0; check.

### CmdStan-1. CmdStan 2.39.0's vendored TBB does not compile under GCC 14.3.0

- Seen: CmdStan 2.39.0, Rtools 4.5, Windows. Not rerun.
- Repro: `cmdstanr::install_cmdstan()`.
- Workaround: none; the hBayesDM cross-check was skipped.
- Record: `dev/learn-findings.md`, near line 712.
- Status: not reported. Check for a known issue first.

## Already reported upstream

The tests harvest these issues from the brms, lme4 and glmmTMB
trackers (`tests/testthat/test-open-issues.R` header). Each is either a
defect frmtmb shared and fixed, or a trap frmtmb pins. They are not
re-verified here, and some may be feature requests. Nothing to file.

| issue | topic | frmtmb test |
|---|---|---|
| brms#391, #734 | nlpar that collides with a column; failed nl fit | `test-diagnostics-ux.R:316`, `:334` |
| brms#494 | single-row newdata | `test-edgecases.R:87` |
| brms#1070 | interval censoring | `test-v07.R:26` |
| brms#1716 | `cens()` with dpar formulas | `test-backlog.R:71` |
| brms#1747 | row permutation and relevel for cov structures | `test-edgecases.R:121` |
| brms#1828 | `mo()`/`mi()` multipliers must be numeric | `test-open-issues.R:58` |
| brms#1903, #1923 | discrete truncation with F(lb - 1) | `test-open-issues.R:173` |
| lme4#144 | rank-deficient designs | `test-edgecases.R:96` |
| lme4#164, #819 | nonlinear SEs; nlmer newdata | `test-open-issues.R:187` |
| lme4#180, #682 | non-integer responses | `test-open-issues.R:161` |
| lme4#196 | `x * (1 | g)` (reformulas-2) | `test-open-issues.R:22` |
| lme4#234 | `(1 | a * b)` expansion | `test-open-issues.R:121` |
| lme4#303 | prediction at an aliased cell | `test-aliased-grouping.R:16` |
| lme4#156, #464 | grouping factors written as calls | `test-aliased-grouping.R:98` |
| lme4#616 | grouping levels, numeric vs character newdata | `test-open-issues.R:133` |
| lme4#622 | `anova()` with different n | `test-open-issues.R:86` |
| lme4#737 | ordinal `simulate()` type | `test-diagnostics-ux.R:247` |
| lme4#818 | `||` over a factor | `test-open-issues.R:245` |
| lme4#880 | zero prior weights | `test-open-issues.R:100` |
| glmmTMB#402 | data-dependent bases at predict time | `test-edgecases.R:76` |
| glmmTMB#439 | non-default contrasts in prediction | `test-backlog.R:3` |
| glmmTMB#625 | offsets in dpar formulas silently dropped | `test-diagnostics-ux.R:391` |
| glmmTMB#773, #937 | matrix-attribute responses and offsets | `test-edgecases.R:112` |
| glmmTMB#776 | REML `anova()` | `test-diagnostics-ux.R:353` |
| glmmTMB#923 | `re_formula = NA` without grouping columns | `test-backlog.R:25` |
| glmmTMB#983, #1143 | REML predictions | `test-open-issues.R:146` |
| glmmTMB#1082 | tensor-product smooths | `test-open-issues.R:287` |
| glmmTMB#1278 | `ar1()` with gapped levels | `test-open-issues.R:203` |
| glmmTMB#1317, #1319, #1325 | `cbind()` responses; no free parameters | `test-diagnostics-ux.R:16`, `:27`, `:61` |
| brms#674, #779; glmmTMB#634, #798, #847, #873, #1120 | regression tags in `dev/test-backlog.md` (near lines 235, 345, 359, 576) | none |
| kaskr/tmbstan#27 | parallel chains on Windows | `NEWS.md` near line 5643; `dev/rtmb-pitfalls.md` item 17 says they now work |

## Considered and excluded

- RTMB `matrix()`, `c()`, `rep(each =)`, `rowSums()` strip the advector
  class (`dev/rtmb-pitfalls.md` 1 to 3): missing overloads, not defects.
- RTMB `qgamma()` takes a rate positionally (pitfall 10): same formals
  as `stats::qgamma()` (checked 2026-09-29).
- RTMB `abs()` derivative 1 at 0; no comparison on the tape; no
  `pmax()`; no OpenMP `parallel_accumulator`; `obj$simulate()` without
  sparse `Z %*% simref`: conventions or limitations.
- RTMB `lbeta()` cancellation (`frmtmb-wt-fams2/dev/fams2-findings.md`
  B1): not confirmed; at (0.01, 5e6) value and gradient match `stats`
  and digamma on 2026-09-29.
- Matern `besselK` NaN gradient, glmmTMB `mat` the same (pitfall 12):
  an unbounded transform reaches `0 * Inf`; robustness, not a defect.
- tmbstan under `devtools::load_all()` on PSOCK workers (pitfall 17):
  expected.
- RTMBode `?ode` example gradient off by 10 %: unexplained, probably
  ill conditioning (`rtmbode-issues.md` section 4).
- Fixed-step ODE integrators give a wrong likelihood: step size.
- RTMBdist `rvm()` scalar-only and `dvm()` OSA refusal: limitations.
- WienR `1 - pWDM()` 4.4 % off at 1.2e-13 (`dev/rdm-gng-findings.md`):
  the complement is taken by the caller, so the relative error is
  cancellation, not a WienR defect.
- brms right and interval censoring of counts as `P(Y > y)`, against
  its inclusive truncation (`NEWS.md` near line 2347): recorded as a
  deliberate divergence, though the records call brms inconsistent.
- brms `re_formula` silently drops an unmatched term
  (`dev/reunc-findings.md`, `dev/simnewdata-findings.md`): a
  divergence decided by the user.
- brms `+ family` on a multivariate formula replaces every family
  (`dev/mv-findings.md`): a recorded divergence.
- brms `dpar = "theta2"` returns the reference component's fixed zero
  (`dev/brms-methods-tests.md`): a convention question.
- brms `ar()` with `cov = FALSE`, response validation in
  `default_prior()`, the `cs()` unused-level rule, `re.form` over
  `re_formula` in `posterior_epred()`: divergences or shared behavior.
- loo does not export its default `loo_compare` method (`R/loo.R:99`):
  normal S3 registration.
- posterior `rhat.default` and gratia `posterior_samples()` reached
  after a later attach (`extensions/frmtmb.sample/NEWS.md`): S3 masking
  order, fixed in frmtmb.
- `insight::get_residuals()` wrong on frmtmb fits
  (`dev/reviews/20260918-shapes.md`): frmtmb lacked a method.
- brms `fcor()` checks that `M` is positive definite with
  `min(eigen(M)$values <= 0)`, which passes a matrix with a negative
  eigenvalue (`dev/review-sarfcor-brms-nocompile.R`). Not a wrong
  model: the Stan program takes `cholesky_decompose(Mfcor)` in
  transformed data, and the fit fails there. brms fails late, where
  frmtmb refuses at the start.
- brms `sar()` bounds `rho` by the eigenvalues of `M` themselves. For a
  matrix that is not row-standardized (a 0/1 adjacency), that interval
  reaches past `1 / max eigenvalue`, where `log1m(rho * eigenMsar)` is
  not a number and Stan rejects the proposal, so brms's support is
  frmtmb's interval, at the cost of rejected proposals and initial
  values. For a row-standardized `M` the interval is narrower than the
  valid one: a deliberate divergence, recorded in frmtmb.spatial's
  NEWS and `?sar`.
- System GSL for rtdists on CI: a system requirement.
- pkgdown 2.2 tangles every article: the fault is knitr-1.
- `pkgcheck` container results: a symptom of tmbstan-1.
- drmTMB 0.7.0 `anova()` refuses every call; website functions absent
  from CRAN: release scope.
- tmbstan autogen `file.copy()` unchecked
  (`dev/reviews/2026-09-09-tmbstan.md`): a hypothetical path, never
  observed.
- RTMB `ifelse()` on an advector test: the error on 2026-09-29 is the
  documented "Comparison is generally unsafe for AD types"; the
  "argument e2 is missing" error in `dev/reviews/2026-09-08-remlopt.md`
  was not reproduced.
- `TMB::sdreport()` with no outer parameter returns the joint precision
  without dimnames (`dev/predfix-findings.md`): a quirk.
- TMB's unexported `TransformADFunObject` failed on CRAN TMB 1.9.25
  (`dev/benchmarks.md`): unexported API.
- An ODE derivative written `0 * y[k]` fails with "Tape has
  Non-consecutive outputs" (`dev/nss-findings.md`): recorded as a
  frmtmb validation defect; the upstream fault is not established.
- nloptr rejects callbacks that declare `...` (`NEWS.md` near 5989):
  probably design.
- glmmTMB `homtoep` NA logLik, `toep` false convergence, `matern`
  non-convergence, NA deviance residuals for beta and tweedie: limits
  or convergence cases on specific designs, not established defects.
- brms `predict()` with a missing grouping column and no
  `allow_new_levels` stops with base R's "object 'visit' not found"
  (`dev/adefects-findings.md` D6): a message request only.
- PyDDM snaps the bound to the dx grid; drmTMB caps correlations at
  0.999999; gratia standardizes by the full-predictor SE: design.
- Recorded as divergences, although the records call brms wrong or
  inconsistent: ICAR over a disconnected graph
  (`dev/brms-likelihood-tests.md`); `re.form` preferred over
  `re_formula` when both are given; unknown arguments swallowed in
  `posterior_*()`; the Savage-Dickey prior density by KDE
  (`dev/sample-ce-findings.md`); `resp` required on dpar priors in a
  multivariate model; the latent AR switch for non-gaussian families.
- frmtmb's own defects (for example `skew_normal()` at alpha = 0,
  `tmbstan_build_broken()` failing open, `frm_sample(laplace = TRUE)`
  NaN draws, the `curry(tv)(zv)` parse defect).
