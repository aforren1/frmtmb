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
  `frmtmb-wt-fams2/dev/fams2-brms-hc-bug.R`.
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
  refuses ordinal mixtures for other reasons.
- Record: `dev/ordinal-findings.md`, "brms 2.23.0 defects met", item 4,
  and punch round 1, m4; `dev/ordinal-log-brms-code.txt`.
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
- Record: `dev/phase3b-rtmb-report-pnorm.md` (filing draft).
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
- Workaround: `poisson()` declares no `lccdf`
  (`tests/testthat/test-cens-lccdf.R:83`).
- Record: that test comment.
- Status: not reported.

## RTMBdist

### RTMBdist-1. The inverse-Gaussian upper tail is computed on the probability scale

- Seen: 1.1.0. Not rerun.
- Repro: the log upper tail reaches `-Inf` at log S = -34, the same as
  `log(1 - F)`.
- Expected: a log-scale upper tail.
- Workaround: `inverse.gaussian()` declares no `lccdf`
  (`tests/testthat/test-cens-lccdf.R:86`).
- Record: that test comment; `NEWS.md` near line 3383.
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
