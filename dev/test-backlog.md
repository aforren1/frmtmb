# Test backlog mined from lme4 / glmmTMB / brms

Edge cases harvested from the reference packages' test suites and
issue/PR history (2026-08-31). Status: DONE items have tests in
tests/testthat/ (mostly test-edgecases.R); the rest are open work.

## Triage, 2026-09-07 (frmtmb 0.53.0)

Every entry under an `Open` heading was re-measured against the tree at
0.53.0, by running it rather than by reading for it. Of the 17, six
closed as done, one closed as moot, one halved, and nine remain; two of
the nine are restated below, because what they predicted is not what
now happens. The Done, Fixed and Verified-immune sections were not
re-measured and stay as the archive they are.

An `Open` heading is therefore live: an entry under one is work still
to do, and every closure carries the measurement that closed it.

## Done

- Frozen data-dependent bases (poly/ns/scale) via predvars stored from
  the fit-time model frame, applied to every linear predictor and RE
  component at prediction; single-row newdata. [glmmTMB#402, #512, #853;
  lme4 predict_basis; brms#494]
- Rank-deficient X: drop aliased columns per linear predictor with a
  message naming the component; frozen column set reused for newdata.
  [lme4#144; glmmTMB test-checkRank.R]
- scale()-style n x 1 matrix responses dropped to vectors. [glmmTMB#937]
- Row-permutation invariance for us/ar1; relevel invariance for us.
  [brms#1747 - silent bias with sorted rows; glmmTMB test-varstruc.R]
- Numeric/character/factor grouping equivalence. [lme4 test-factors.R]
- Duplicate multivariate responses rejected. [brms tests.standata.R]
- trials() validation: y > trials, non-integer y, constant literal.
  [brms data-response.R taxonomy]
- Gaussian rescaling equivariance of coefficients and logLik.
- NA rows dropped across response/predictors/grouping jointly.
- Pooled model comparison over imputations: `anova()` on
  `frmtmb_multiple` with the D1, D2 and D3 rules. Checked against
  `mice::D1`/`D2` on a poisson GLM, where our ML covariance and
  `df.residual` equal `glm()`'s, and against the Meng-Rubin formula
  written out by hand for D3; plus the degenerate identity (identical
  imputations collapse D3 and D2(likelihood) to the single-fit LRT)
  and a GLMM smoke case with no reference implementation, where the
  three rules have to agree with the per-imputation LRTs in rejection
  direction. tests/testthat/test-pooled-anova.R

## Addressed in v0.6 (moved up from Open)

- Non-default contrasts in prediction (global option and per-factor).
- re.form = NA without grouping columns in newdata.
- Formula-environment robustness: combined model frame stored on the
  fit; model.frame/predict/emmeans survive the calling env vanishing.
- na.exclude: fitted/residuals/predict padded via napredict; na.pass
  errors informatively; Inf responses rejected.
- Slash grouping (1|a/b) equivalence with (1|a) + (1|a:b).
- cens() x dpar-formula cross-product (vs hand-rolled reference).
- Frequency-weight equivalence for aggregated poisson data.
- dpar/nlpar names with dots or underscores rejected.
- predict() warns on unknown arguments.

## Open - high priority

- DONE, lane wt-reunc (`dev/reunc-findings.md`): `predict()` draws the
  group effects at a known level jointly from their conditional law,
  `fitted()`'s ordinal and categorical route carries them, and a
  partial `re_formula` keeps only its own terms. The entry as decided:
  DECIDED 2026-09-22 (user): match brms. `predict()` AND `fitted()` carry
  the group effects' uncertainty at a level the fit saw; `re_formula = NA`
  adds none, and a partial `re_formula` adds only its own terms, so this
  goes in one lane with the silent partial-`re_formula` defect. The lane
  must state which coverage the interval claims (conditional on each
  group's true effect, or averaged over groups) and measure that one.
  The record as filed: should `predict()` carry `Var(b | y)` at a grouping
  level the fit saw? Today it draws conditional on the modes, so its
  interval at a known level of a mixed fit is narrow: 0.9445 out of
  sample on `dev/shapes-coverage.R`'s mixed design, 0.9336 on the
  reviewer's, and adding exactly the analytic conditional variance
  moves it to 0.9609 and 0.9477. brms's draws carry the posterior of
  `r_g` at a known level, so leaving the term out DIVERGES from brms;
  the case for leaving it out is only that every other method here
  (`fitted()`, `frm_linpred()`, `residuals()`) is conditional on the
  modes. Carrying it means drawing each replicate's modes from their
  joint conditional law (the `sdreport()` joint precision has it),
  not adding a per-row variance, so that `summary = FALSE` stays
  jointly right across rows of one group. Filed by lane wt-shapes,
  punch round 2 (`dev/shapes-findings.md` section 3).

- `predict()` on a draw that is FINITE but absurd. Masking a
  non-finite cell (punch round 2) does nothing for a draw that
  overflows to a huge finite number: a poisson row whose linear
  predictor sits just under the overflow point summarizes to an
  Estimate of 1.4e146 in the recheck. The parametric-bootstrap law is
  honest about a barely identified parameter and the summary is
  useless. Candidates: report the median and quantiles only when the
  mean is dominated by a few draws, or warn when the draw law puts
  mass far outside the data's range. Not fixed; filed by lane
  wt-shapes.

- DONE, lane wt-reunc: such a term is refused now, classed and named.
  brms does NOT refuse it: `check_re_formula()` drops an unmatched term
  silently, so `~(1 | nosuch)` there equals `re_formula = NA`
  (`dev/reunc-log/brms-lane.txt`). frmtmb departs from brms here on
  purpose; `dev/reunc-findings.md` names it for the user. The entry:
  `re_formula` naming a grouping factor the model does not have is
  silently treated as `NULL`. Measured by the review of items 2.6d and
  2.6f (`dev/reviews/20260918-shapes.md`, m5): `predict(fit,
  re_formula = ~(1 | nosuch))` returns the CONDITIONAL prediction, bit
  for bit what `re_formula = NULL` returns, with no warning, on a
  grouping factor the model never had. The same holds for `fitted()`,
  `frm_linpred()` and `conditional_effects()`. `check_re_form()`
  accepts any formula and `re_form_keeps()` only asks whether it is
  `NA` or `~0`, so a typo in a grouping-factor name reads as "keep
  everything". This is the same family as the partial-`re_formula`
  defect (a formula naming SOME of the model's blocks is also read as
  all of them) that another lane filed, and it wants one fix for both:
  resolve the formula's terms against the fit's blocks, keep exactly
  those, and error on a name the fit does not have. Pre-existing; the
  shapes lane did not touch `check_re_form`.

- Constant-weight (non-)invariance for gaussian documented and tested.
  The behavior is already right and only unpinned: at 0.53.0 a constant
  `weights(w)` of 2 leaves the coefficients invariant to 5.9e-06 and
  exactly doubles the logLik (-347.5174 against -173.7587). What is
  missing is the regression test and the sentence saying so.
  [lme4 priorWeights.R]

- `simulate(re_formula = )` still reads any formula as "condition on
  every group effect". Its switch is lme4's, `NA` redraws the effects
  and anything else keeps the modes, so `~0` and a partial formula such
  as `~ (1 | g)` both simulate conditionally with nothing said. Lane
  wt-reunc fixed `re_formula` in `predict()`, `fitted()`,
  `frm_linpred()` and `frm_lp_basis()` and left `simulate()` alone,
  because brms has no `simulate()` and what a partial formula should
  MEAN there (redraw the dropped terms, as `NA` does for all of them?)
  is a decision. `dharma_residuals()` and `pp_check()` on a fit reach
  it. Filed by lane wt-reunc. RESOLVED by lane wt-simnewdata: the user
  decided on 2026-09-23 that simulate() reads re_formula as predict()
  does, and a term that is not kept is redrawn
  (`dev/simnewdata-findings.md`).

- DECIDED, waiting for a lane (user, 2026-09-24): `re_formula = NA`
  must keep EVERY smooth, as brms keeps every smooth under any
  `re_formula`. Today `predict(re_formula = NA)`, and so
  `simulate(re_formula = NA)`, drops three kinds of smooth that are
  indexed by a grouping factor: `s(g, bs = "re")`, a factor smooth
  `s(x, g, bs = "fs")`, and a `t2()` with an `re` margin. Measured by
  the review of lane wt-simnewdata (`dev/simnewdata-review/rv-smooth.R`,
  log `dev/simnewdata-review/log/smooth.txt`): of the constructions
  that fit, `predict(NA)` keeps every smooth block of `s(x)`,
  `s(x, by = f)`, `t2(x, z)`, `gp()`, `hsgp`, `sigma ~ s(x)` and
  `s(x) + (1 | g)`; it keeps block 1 of `s(g, bs = "re")` and drops
  block 2, which `simulate(NA)` then redraws (row sd over sigma 1.260);
  it keeps none of `s(x, g, bs = "fs")` (redraws 1, 2, 3; 1.741) and
  none of `t2(x, g, re)` (redraws 1, 2; 1.677). The rule to change is
  `smooth_group_block_ids()` and its use in `lp_eta_design()`;
  `sim_group_block_ids()` follows it. The refusal of a partial formula
  beside a factor smooth (`re_fit_components()`'s "unnamed" content)
  has to be revisited with it, because its stated reason is that `NA`
  drops that content. Recorded, not changed, by lane wt-simnewdata.

- The uncertainty in the VARIANCE PARAMETERS is not in any interval
  here. Measured by lane wt-reunc with an exact positive control
  (`dev/reunc-findings.md` section 3): with sigma and tau known, the
  prediction interval covers 0.9523 and 0.9513 on two designs; with the
  fit's own ML estimates in the same formula it covers 0.9403 and
  0.9337, and `predict()` and `fitted()` sit on those numbers. REML
  covers better in every arm (0.9437, 0.9380), which is the same cause
  seen from the other side. `fitted()`'s interval is the one where it
  shows most: 0.9170 and 0.8977 against the oracle's 0.9480. brms
  carries this term because its posterior includes tau. Candidates: a
  t-like widening with the profile curvature of tau, or a bootstrap
  interval. Filed by lane wt-reunc.

- A quadrature fit's SCALAR standard errors leave out the group-effect
  term and say nothing. On `y ~ x + (1 | g)` with `bernoulli()` and
  `quadrature = TRUE`, `fitted()` reports Est.Error 0.09407 and 0.08913
  at a known level against 0.09340 and 0.09095 at `re_formula = NA`, so
  the term is missing; a Laplace fit of the same data reports 0.15115
  against 0.09158. The objective marginalizes the random effects, so
  the joint covariance of such a fit carries `beta` and `theta` and no
  `b`, and `lp_delta_A()`'s guard does not fire because the design adds
  no b columns to pair. Pre-existing in 0.61.0, and more visible since
  lane wt-reunc made the finite-difference route warn in the same
  situation. Candidates: warn on the scalar route too, or refuse
  `fitted()` at a known level on a quadrature fit. Filed by lane
  wt-reunc, punch round 1.

### Closed at the 2026-09-07 triage (was: Open - high priority)

- DONE. Interval censoring (cens code 2 + y2), NA y2 off the interval
  rows. `cens_code_map` carries `interval = 2` (R/frame.R:531), an NA
  y2 is allowed off interval rows (:1086) and refused on them (:1383),
  and a gaussian fit over all four codes converges. [brms#1070]
- DONE. Discrete truncation normalizes with F(lb - 1)
  (R/objective.R:71-77): a poisson `trunc(lb = 2)` fit's logLik
  matches a hand-rolled
  `ppois(lb - 1)` reference to 1.1e-13. Duplicate of the
  "Verified immune" entry below, which already said so. [brms#1903]
- DONE. `predict(newdata =, allow_new_levels = TRUE, se.fit = TRUE)` at
  a new level of a `(1 | ID | g)` block spanning mu and sigma returns a
  finite estimate and standard error, and a `student()` fit with
  `sigma` and `nu` swapped in `bf()` gives a bit-identical logLik
  (difference 0.0e+00). [brms#779, #674]
- DONE. tests/testthat/test-aliased-grouping.R:144 pins the label order
  against `levels(droplevels(a:b))`, and `(1 | a * b)` is covered under
  "Verified immune" below. [lme4#635/#636/#945]

## Open - medium

- DONE 2026-09-22 (lane `wt-correct`, `dev/correct-findings.md`
  section 7). `?frm` says it under `REML`, and `test-smooths.R`
  measures it: over seeds 44 to 48, `log(1 / sd^2) - log(sp)` is
  -0.011 to -0.005 for the `mu` smooth under REML against 0.129 to
  0.262 under ML, while the `sigma` smooth's stays where ML left it
  (0.142 to 0.170 against 0.138 to 0.165, ratio 1.026 to 1.037). The
  test pins those RATIOS. A first design put the `sigma` smooth on its
  boundary in both packages, where the difference measures nothing;
  the recorded one has both smooths interior. The entry as filed:
  REML against mgcv on a location-scale smooth (filed 2026-09-22).
  frmtmb's `REML = TRUE` integrates the `mu` coefficients and keeps the
  distributional coefficients outer, which is the double-GLM REML
  (Smyth and Verbyla; `nlme` varFunc; pinned against
  `gls(method = "REML")` at 1e-5). mgcv's LAML integrates every
  coefficient of every linear predictor, so for `sigma ~ s(z)` the two
  `REML` criteria differ by design, and nothing says so: the only
  frmtmb-versus-`gaulss` comparison (`test-smooths.R`) runs under ML.
  Two parts. (1) One sentence in `?frm` under `REML` naming the
  difference, so a ported mgcv model compared under REML is not read as
  a defect. (2) A test fitting `bf(y ~ s(x), sigma ~ s(z))` with
  `REML = TRUE` against `gam(list(y ~ s(x), ~ s(z)), family =
  gaulss(b = 0), method = "REML")` that MEASURES the gap in the fitted
  curves and smoothing parameters and pins its size, rather than
  asserting agreement. Integrating the distributional coefficients is
  not proposed: it would turn a gaussian model's sigma back into the
  ML estimate, and the inner problem can be unbounded as sigma goes to 0.

- Draws hint that leads to a refusal (found 2026-09-22 by the round-3
  recheck of 2.6d/2.6f). On a draws object, a `newdata` holding an
  unseen level WITHOUT `allow_new_levels` gets core's error, whose hint
  says "Use allow_new_levels = TRUE". Following the hint reaches
  frmtmb.sample's refusal, since new levels are not implemented on
  draws. The same happens with `sample_new_levels` alone. Either
  implement new levels on draws or have the draws methods replace the
  hint with the refusal.

- Singular-fit detection: the isSingular verdict is done (see the
  diagnostics/UX cluster below); what is left is profile CIs on
  boundary parameters. RESTATED 2026-09-07, because it does not hang,
  and it fails in TWO ways rather than one. Sweeping 17 gaussian fits
  whose grouping factor carries no signal, all returning in well under
  a second, `confint(parm = "theta_1", method = "profile")`:
  - dies inside `approx()` with "need at least two non-NA values to
    interpolate", a raw internal message (10 of 17 here, 6 of 17 on a
    reviewer's sweep; the split is seed-dependent, both modes always
    appear);
  - or returns a CI with `lwr = NA` and a finite `upr`, warning only
    "NA/NaN function evaluation", which names no component (7 of 17
    here, 11 of 17 on the reviewer's). This is the worse mode: a
    plausible-looking half-answer rather than a refusal, and the first
    pass of this triage recorded only the error.

  One fix covers both: catch it and return NA bounds with a warning
  naming the component, which is what the Wald path already does.
  [lme4 test-isSingular.R, #660]
- The offset ARGUMENT form: `frm(..., offset = )` is still an
  "unused argument" error. HALVED 2026-09-07: the offset-only model is
  done, `y ~ 0 + offset(o)` fits, and offset inside a dpar formula was
  already done. [glmmTMB test-offset.R, #625, #286]
- gam-style exclude= for zeroing individual smooths in prediction.
  Confirmed 2026-09-07: `predict(exclude = "s(z)")` still warns
  "ignoring unknown arguments to predict(): exclude". `re.form = NA`
  keeps smooths (implemented; keep the regression test).
  [glmmTMB test-smooths.R; mgcv semantics]
- confint/vcov excluding mapped parameters everywhere. Constant dpars
  are covered and nothing else maps yet, so this is a rule to apply
  when a new map appears rather than work waiting. [glmmTMB#1120]

### Closed at the 2026-09-07 triage (was: Open - medium)

- DONE. predict type grid for zi families: on a
  `zero_inflated_poisson()` fit, `type = "response"` equals
  `(1 - zi) * conditional` exactly, maximum difference 0.0e+00. The
  third part of the same entry, truncated conditional means, is done
  too and recorded under the open-issue sweep's "Fixed": `fitted()`,
  `predict(type = "response")` and `residuals()` report
  E[Y | lb <= Y <= ub].
- MOOT. "zprob on a non-zi model returns 0 not garbage": it now refuses
  by name instead, `Unknown dpar: 'zi' for response 'y'. Available:
  mu`, which is the better answer and retires the item as written.
  [glmmTMB#798, #873, #634]
- DONE, and a duplicate. dpar/nlpar names with dots or underscores are
  refused at R/bf.R:157 with the message this entry asks for; the same
  item already sits under "Addressed in v0.6" above.
  [brms tests.brmsformula.R]

## Diagnostics and UX cluster (2026-09-01)

The porting papercuts: what a user moving from glm, lme4, glmmTMB or
brms meets before anything statistical goes wrong. Regression tests in
tests/testthat/test-diagnostics-ux.R.

### Fixed

- `cbind(successes, failures)` binomial responses are ACCEPTED, which
  is what every reference package does. The spelling is rewritten at
  parse time to the internal `successes | trials(successes + failures)`
  form, so frame assembly, `valid_y`, `fitted()`, `predict()` and
  `simulate()` need no matrix branch of their own and the two spellings
  give bit-identical fits. Verified against `glm(cbind(s, f) ~ x)`:
  coefficients, logLik and standard errors all agree to 1e-5 or better.
  `cbind(s, f) | trials(n)` (two contradictory sources of trials), a
  three-column `cbind()`, and a fractional failure column are refused
  with messages that name the actual problem; the old message spoke
  about `trials` without saying what the user had written wrong.
  Matrix-response families (multinomial) keep their own spelling.
  [glmmTMB#1319, #1325]
- A model with zero free outer parameters (`y | trials(n) ~ 0`) fits
  degenerately instead of dying inside nlminb with "'d' must be a
  nonempty numeric (double) vector". `optimize_obj()` evaluates the
  template once and reports convergence 0 with the message "no free
  parameters (degenerate model)"; `logLik` is the template's, with
  `df = 0`. `par_est_se()` and `diagnose()` were fixed alongside: an
  empty sdreport summary has no rows for `as.list(what = "Std. Error")`
  to reshape, and `check_convergence()` took `max()` of an empty
  gradient. [glmmTMB#1325, #1317]
- One-level grouping factors and gaussian OLRE now have lme4's
  three-way control vocabulary: `frmtmb_control(check_nlev_1 =)` and
  `check_olre =`, each `"warning"` (default) / `"ignore"` / `"stop"`.
  Both previously fit silently. The nlev check fires only on SCALAR
  blocks - a structured block over several terms per level (ar1, us,
  the spatial covstructs) is one realization of a field, where a single
  grouping level is the normal spelling. The OLRE check fires only when
  sigma is free to absorb the variance: `se()` and a constant sigma
  both identify the split, which is exactly the random-effects
  meta-analysis. [lme4 lmerControl checks]
- `diagnose()` gained three checks and lost a crash. It errored on
  ANY fit without random effects, because `theta` is NULL there and
  `abs(NULL)` is an error rather than an empty result. The new checks:
  a complete-separation heuristic for binomial-type fits (|coef| > 10
  on the link scale WITH a standard error to match - a genuinely large
  effect on a well-populated cell keeps a small se); predictor-scale
  warnings (|log10 sd(x)| > 3, pointing at `autoscale`); and an
  isSingular-style verdict listing every variance component on the
  boundary of its parameter space (sd at zero, |cor| at one), read off
  the estimates so it stands independently of `pdHess`.
  [glmmTMB diagnose(), lme4 isSingular]
- `simulate()` returns the response's own type. The four ordinal
  families had NO simulator at all ("Family 'cumulative' has no
  simulator yet"); they now have one, and it hands back an ordered
  factor carrying the original levels rather than 1..K codes. The
  category distributions are built in plain doubles, one branch per
  family, and reproduce `exp(lpdf)` at the observed category to 1e-15
  for cumulative, sratio, cratio and acat (the taped lpdf only scores
  the observed category, so the whole distribution is new code and is
  pinned against it). `multinomial()` gained a simulator too, returning
  an n x K count matrix with the response's column names in a data
  frame of matrix columns (the lme4 convention). `simulate()` now
  respects `na.exclude` padding like `fitted()` and `residuals()` do;
  the internal consumers that work in fitted-row space
  (`frm_bootstrap()`, `dharma_residuals()`, `pp_check()`) unpad first.
  `dharma_residuals()` refuses ordinal fits outright: DHARMa's rank
  transform needs a numeric scale, and an ordinal response has only an
  order. [glmmTMB test-simulate.R; lme4#737]
- An nlpar whose name is also a data column is refused. The nonlinear
  body resolved the name to the PARAMETER and silently ignored the
  column, so the fit ran and reported numbers for a model the user did
  not write. Silent precedence was the other option and is worse: no
  one would remember it. [brms#391]
- A nonlinear fit that dies from the default zero starting values now
  names `start=` and the template layout instead of surfacing RTMB's
  own error; a non-convergence warning on an nl fit says the same. The
  raw "NA/NaN function evaluation" from nlminb is muffled on that path,
  since it is the optimizer noticing the same undefined objective the
  error already names. [brms#734 doctrine]
- `anova()` no longer refuses every REML fit. A REML likelihood is a
  likelihood for the error contrasts of its fixed-effect design, so two
  of them are on a common scale exactly when the designs span the same
  column space - which is the usual REML comparison, of
  variance-component structures with the fixed effects held fixed. The
  check projects each design onto the other's column space, so a
  reordered but equivalent design passes (which also settles the "REML
  logLik invariance to fixed-effect term order" item). Different column
  spaces are refused with the reason, and mixing REML with ML fits is
  refused separately. [glmmTMB#776]

### Already worked (regression tests added)

- `offset()` inside a dpar formula (`sigma ~ x + offset(o)`) reaches
  the likelihood exactly - the fitted logLik matches a hand-rolled
  `dnorm(y, mu, exp(Xb + offset))` to 1e-13 - and follows through to
  `predict(dpar = "sigma")` both in sample and on newdata. Pinned down
  so it cannot regress into the silent drop glmmTMB#625 describes.
  [glmmTMB test-offset.R]

## Open-issue sweep (2026-09-01)

Mined from the *currently open* trackers of brms (145), lme4 (191) and
glmmTMB (223); 34 shortlisted and run against frmtmb. Fixed items have
regression tests in tests/testthat/test-open-issues.R, except lme4#303
and lme4#464/#156, which live in
tests/testthat/test-aliased-grouping.R.

### Fixed

- Random-effect terms crossed with `*` or `:` (`y ~ x * (1 | g)`) now
  error instead of being silently refit as `+`. [lme4#196]
- `mo()`/`mi()` interaction multipliers: the numeric type gate never
  fired for character vectors (`is.numeric(as.numeric("a"))` is TRUE),
  so the column went all-NA and the fit died at "NA/NaN gradient
  evaluation". [brms#1828]
- `anova()` rejected fits with different `nobs`; previously it compared
  likelihoods across data sets and returned a negative Chisq.
  [lme4#622]
- `||` over a factor produced a fully correlated `us` block (|cor| up to
  0.99), the opposite of what the syntax promises. `parse_linpred()` now
  expands each `||` term itself and tags every piece `diag()`, so a
  factor's levels get independent variances; the fit is identical to an
  explicit `diag(f | g)`. Numeric double bars keep lme4's block split
  and their old estimates, because `diag` and `us` coincide at dimension
  one. [lme4#818]
- Prediction at an aliased cell of a rank-deficient design returned the
  random-effect contribution alone. The fit now freezes the null space
  of the fixed-effect design, and `predict(newdata =)` returns NA (with
  NA standard errors) plus one warning naming the dropped columns for
  the rows that load on it. Rows that merely restate a kept column
  (`x2 = 2 * x`) stay exact, and the in-sample paths are untouched.
  lme4 still returns the partial sum silently. [lme4#303]
- Grouping factors written as calls - `(1 | factor(x))`,
  `(1 | interaction(a, b))` - now fit. reformulas re-evaluates the bar
  RHS inside the model frame, where the call's own arguments are not
  columns, and died with an error raised several frames down
  ("unique() applies only to vectors"). Frame assembly now points the
  bar at the frame column the expression already produced; the original
  expression is kept for labels and for newdata prediction, and `:`
  / `/` groupings keep their reformulas expansion. [lme4#464, #156]
- The quadrature defect cluster (dev/fuzz-findings.md N1-N4). TMBad's
  `marginal_gk` rescales each integrand ONCE, at whatever parameter
  values the template holds when `MakeADFun` tapes it, and freezes that
  `(mu, sigma)` pair. One cold calibration explained all of it: every
  conditional mode but the first came back `NA` (the marginalized
  objective carries none, and `parList()` slid an outer value into the
  first slot), and poisson, Gamma and Beta over nested scalar blocks -
  Beta over a single one - died at a bare `NA/NaN gradient evaluation`.
  `frm()` now fits the plain Laplace objective first, tapes the
  marginalized one at that optimum, and reads the modes back from the
  inner Newton solve there. They match `glmer(nAGQ = 25)`'s `ranef()`
  to 3e-05. `quadrature` crossed with `trunc()` is refused instead:
  the normalizer `log(F(ub) - F(lb))` underflows at the Gauss-Kronrod
  nodes, so the objective is `-Inf` even at the Laplace optimum, and
  the fit used to report `logLik = +Inf` as converged. What survives is
  a runtime limitation on hard likelihoods (singular variance
  components, nested blocks on thin data), reported as an error naming
  `quadrature` rather than an RTMB string. Regression tests in
  tests/testthat/test-quadrature-defects.R.
- `mixture()` under `REML = TRUE` or `frmtmb_control(profile = TRUE)`
  is refused. Both integrate the fixed effects out with a Laplace
  approximation about a single inner mode, and a mixture likelihood is
  invariant to permuting its components, so it is multimodal in exactly
  those coefficients. The fits used to stop at `NA/NaN gradient
  evaluation` or report a gradient near 1e9 with no guard.
  `quadrature = TRUE` stays allowed - it marginalizes the random
  effects, not the coefficients - and test-v19.R pins down that it is
  exact when the per-group integrand is univariate.
- Truncation reached only the likelihood, never the post-fit surface.
  `fitted()`, `predict(type = "response")` and `residuals()` now report
  E[Y | lb <= Y <= ub] (closed forms for every family with an `lcdf`:
  gaussian, lognormal, poisson, exponential, weibull, inverse.gaussian;
  the poisson form reuses the objective's F(lb-1) convention), and
  `simulate()`/`posterior_predict()` reject out-of-bounds draws instead
  of sampling the untruncated distribution (23% of draws used to land
  outside [lb, ub] on a Poisson trunc(2, 6) fit, which silently
  invalidated DHARMa and pp_check). Dpar-scale predictions stay
  untruncated on purpose. `residuals(type = "osa")` was wrong too: the
  taped density integrates to 1 only over [lb, ub], so the conditional
  CDF is now built on that domain (and never with `fullGaussian`, since
  a truncated gaussian is not gaussian). [brms#1923, #1903]

### Mitigated

- `ar1()`/`hetar1()` over an ordering factor with gaps still treat level
  POSITION as time - dropping times 7-9 from a 1..10 series makes
  cor(t6, t10) come back as rho, not rho^4, and biases rho itself. That
  reading is glmmTMB's, and the glmmTMB agreement tests pin it down, so
  the likelihood is unchanged and frame assembly warns instead: when the
  ordering levels are whole numbers but not consecutive, the warning
  names the gap and points at `ou()` over `num_factor()`, which is the
  correct spelling for irregular spacing. Non-integer labels stay
  silent, since position is then the only available meaning.
  [glmmTMB#1278]

### Open - medium

- `quadrature = TRUE` still breaks down on hard likelihoods: a variance
  component near zero (the fuzzer found `sd = 8.8e-05`) defeats the
  finite-difference curvature estimate `marginal_gk` calibrates with
  (`dx = 1`), and nested scalar blocks on thin data make the outer
  integrand the output of a frozen inner rescaling. `quad_fit()` tries
  three calibration points and then errors, naming `quadrature`. A real
  fix wants an integrator that recalibrates per evaluation; TMBad's
  `adaptive = TRUE` is meant to be that and is measurably worse, so it
  would have to be built rather than switched on.
  Confirmed 2026-09-07: unchanged at 0.53.0, and the fix is still an
  integrator to build rather than a flag to set.
- No `link_zi` / `link_hu`: the zero-inflation and hurdle parts are
  logit-only, as in glmmTMB. Confirmed 2026-09-07:
  `zero_inflated_poisson(link_zi = "probit")` is an unused-argument
  error. [glmmTMB#847]

(`cbind(successes, failures)` and the zero-free-parameter model moved
to the diagnostics/UX cluster below; both are fixed.)

### Verified immune (regression tests added)

- Zero prior weights are exactly equivalent to subsetting. [lme4#880]
- `(1 | a * b)` expands to a + b + a:b; lme4 still cannot. [lme4#234]
- Nonlinear fixed-effect SEs match nlme; nlmer is ~100x too small on
  the same fit, and `predict(newdata =)` works. [lme4#819, #164]
- REML predictions agree with `fixef()`. [glmmTMB#1143, #983]
- Numeric vs character grouping levels in newdata. [lme4#616]
- Non-integer binomial/poisson responses rejected. [lme4#682, #180]
- `t2()` matches mgcv; `te()`/`ti()` refused clearly. [glmmTMB#1082]
- Discrete truncation normalizes with F(lb-1). [brms#1903, #1923]

## OSA and inference surface (2026-09-01)

From the compatibility-registry probes and the grammar fuzzer.
Regression tests in tests/testthat/test-osa-inference.R.

### Fixed

- `cens()` x `residuals(type = "osa")`. A censored row contributes a
  probability MASS, and in the tape that contribution does not depend
  on the observation at all, so `fullGaussian` inverted an exactly
  singular Hessian block and `oneStepGeneric` integrated a flat slice
  to infinity and returned NaN. What is well defined is the CDF of the
  UNCENSORED rows conditional on the censoring events, so those rows go
  in `subset` and the censored rows in `conditional`; an uncensored row
  is a draw that landed inside the censoring window, so its integration
  domain is that window, exactly as a `trunc()` fit's is. Censored rows
  return NA. Verified against the analytic conditional PIT
  `qnorm(F(y) / F(c))` to 7e-09 and by KS uniformity. Row-varying
  censoring points (the distribution of an uncensored response is then
  not identified without a model for the censoring process) and
  interval censoring are refused with a message that also says why
  `dharma_residuals()` is not the fallback: `simulate()` draws the
  LATENT uncensored response, so its draws are not comparable with the
  observed censored values.
- `cens()` x `simulate()` semantics. Drawing the latent response is
  CORRECT and is what brms does: no `posterior_predict_*` method in
  brms 2.23.0 reads the censoring column (only `log_lik_censor` does),
  and brms's `pp_check()` therefore drops the censored rows outright
  ("Censored responses are not included"). The model describes the
  latent distribution; censoring belongs to the observation process.
  Documented as the default, with `simulate(censored = TRUE)` as the
  opt-in that applies the mechanism: every draw is recorded at the
  edge of the observation window, so the draws become comparable with
  the observed data. That needs one censoring point per side (type-I
  censoring), because an uncensored row's censoring point is unknown
  when the times vary by row; row-varying times and interval censoring
  are refused.
- Ordinal x `residuals(type = "osa")`. `oneStepPredict` re-tapes with
  the response promoted to a parameter, and the ordinal lpdfs index and
  compare with it (`y == K`, `tau[pmin(y, K1)]`), which no advector
  supports. The four ordinal families now carry an OSA branch that
  selects the category with a Lagrange basis over 1..K - exact in
  floating point at integer y - and applies the `@keep` data-term
  indicator that RTMB's own densities get from `dGenericOSA`. The
  residuals match the analytic randomized quantile residual to 4e-14
  and are uniform under KS for cumulative, sratio, cratio and acat,
  with or without random effects.
- Raw LAPACK error from `solve(jointPrecision)` under REML or
  `profile = TRUE`. See fuzz finding N4.

### Open - medium

- OSA under `cens()` covers the uncensored rows only. A residual for a
  censored row would have to be a randomized quantile inside the
  censoring interval, which `oneStepPredict` cannot produce; doing it
  would mean computing the conditional CDF at the censoring point
  outside TMB. Confirmed 2026-09-07: a standing design limit, not a
  defect that a fix would close.
- `dharma_residuals()` on a censored fit compares latent draws with
  observed censored values and is not valid. It neither warns nor
  refuses. Confirmed 2026-09-07: `dharma_residuals()` on a gaussian
  `cens()` fit returns normally, silently. `simulate(censored = TRUE)` makes the draws comparable, but
  the point mass it puts at each censoring point is not a distribution
  DHARMa's rank transform can use, so it is not a drop-in fix.

## Code-review fixes, v0.22-v0.25 (2026-09-01)

Confirmed defects from a review of the v0.22-v0.25 work. Regression
tests in tests/testthat/test-review-v25.R.

### Fixed

- `cens()` x `trunc()` was a silently wrong likelihood. Censoring
  contributed `log(1 - F(y))` / `log(F(y))` against the whole line and
  truncation then subtracted `log(F(ub) - F(lb))`, so only the
  normalizer was windowed: the model censored the UNTRUNCATED
  variable. Under truncation a right-censored row observes
  `y < Y <= ub` and a left-censored one `lb <= Y <= y`, so the
  numerators are `(F(ub) - F(y)) / Z` and `(F(y) - F(lb)) / Z`, with
  `Z` the window mass and the discrete inclusive lower bound still
  `F(lb - 1)`. An interval-censored row was already a windowed
  difference of CDFs and only needed the division by `Z`, which it had.
  Verified against hand-rolled likelihoods: gaussian on `[-1, 3]`
  right-censored at 1.5 matches to 1e-8, and the residual sd drops
  from 1.021 (old composition) to 0.896 against a data-generating 0.9;
  the discrete composed form matches a `ppois` reference to 1e-10.
  The corrected likelihood is also the one `simulate(censored = TRUE)`
  draws from (truncate by rejection, then clamp to the censoring
  window), which the old form was not. The OSA path is unchanged and
  still matches the analytic PIT `(F(y) - F(lb)) / (F(c) - F(lb))` on
  the uncensored rows to 5e-14, so `osa_cens_domain`'s window
  intersection and the likelihood still agree.
- `residuals(type = "deviance")` ignored `se()`. The gaussian unit
  deviance `(y - mu)^2` assumes one dispersion, which `se()` removes:
  two rows with `se` 0.1 and 0.5 got the same scale. The known
  variance now enters as a glm prior weight `sigma^2 / s_i^2`, which
  is the familiar `1 / se_i^2` when `se()` maps `sigma` out and 1
  (hence exact `glm()` agreement) when there is no `se()`.
- `quadrature = TRUE` x `predict(se.fit = TRUE)` died non-conformable.
  A marginalized objective has no `b` rows in its covariance, so the
  Z block was cbind'd against an empty column-position vector. The RE
  columns are dropped and the standard error is reported conditional
  on the modes, with a warning saying the random-effect uncertainty is
  not in it.
- `frm_sample()` mode inits could sit outside `lower`/`upper`. Stan
  turns a bound into a constrained transform, so an init at or past it
  has no preimage and rstan reports an unrecoverable initialization
  failure naming neither parameter nor bound. Every chain's init is
  clamped strictly inside the box (jittered chains included), and a
  bound that excludes the ML mode itself warns.
- `aterms_for_newdata()` dropped `vint()`/`vreal()`. A custom family
  reading `aterms$vreal1` in its `mean_fn` then returned a LENGTH-0
  prediction with no message. Those payloads are now required on
  newdata (the error names the missing column) and the general
  fallback warns rather than omitting silently.
- `quad_fit()` kept the first non-stationary candidate instead of the
  one with the lowest objective.
- `get_joint_cov()` inverted the joint precision by hand, so a
  singular REML/profile fit threw a raw LAPACK message from inside
  `predict(se.fit = TRUE)`. It routes through
  `solve_joint_precision()` and degrades to NaN plus one warning, as
  `vcov()` does.
- `decode_cens()` rejected the character codes `"0"`, `"1"` and
  `"-1"` that its own error message advertises; numeric-looking
  strings are coerced before prefix matching.

### Open - medium

- `conditional_effects(method = "predict")` still drops an unevaluable
  `vint()`/`vreal()` payload silently (`ce_aterms()`), where
  `aterms_for_newdata()` now errors or warns. The grid case is
  deliberately laxer, but a custom family that needs the payload gets
  the same length-0 mean there. Confirmed 2026-09-07 at
  R/conditional-effects.R:14: `ce_aterms()`'s `strict` set is trials,
  se, trunc_lb and trunc_ub, and vint/vreal are not in it.
- `cens()` stays refused for discrete families, so the discrete
  censored-truncated branch is reachable only through
  `build_objective()` and is tested that way. Confirmed 2026-09-07:
  a poisson `cens()` fit refuses with "cens() is not supported for
  discrete families yet (truncation is)".

## Recorded by the loose-ends lane, 2026-09-07

Found while clearing other debt; none of it is done. The first entry is a SILENT WRONG ANSWER, which this package treats as the worst category: the feature search stays quiet in exactly the case that matters and speaks in one that does not.

- `frm_curve_feature()`'s stencil span check is incoherent, and both
  main and every branch share it. When the five-point stencil leaves a
  `ps()` term's span, the guard at `curve-feature.R:235` is true and
  the code then re-asks about the WHOLE grid: it refuses when some
  unrelated grid row leaves some other term's span (evidence the search
  never touched), and stays silent in the case that actually matters,
  where the stencil left the span but the grid is clean. `parts$span`
  holds the right answer at `:224` and is discarded. `frm_curve()` and
  `frm_curve_deriv()` both warn about their own excursions; the feature
  path, which `sp_span_stop()`'s own roxygen at `curve-cov.R:50-56`
  argues is the one that most needs to speak, says nothing. Fixing it
  means changing what a refusal keys on, so it is a behavior decision
  rather than debt. Covered both ways by test-span.R's two-ps() test.
- `gp()` has no Rd topic. It is a formula special, parsed at
  R/parse.R:960 and registered at R/covstruct.R:1587, documented only
  in the vignettes and not exported, so an unqualified roxygen
  `[gp()]` resolved against installed brms until this lane unlinked
  it. A topic for it is a documentation feature rather than debt:
  `ps()` has one because `ps()` is called directly and `gp()` never is.
- `local_mocked_bindings()` without `.package` at
  extensions/frmtmb.ode/tests/testthat/test-ode.R:128 and :132. The
  same defect was cleared in core and in frmtmb.sample; ode was left
  because verifying a change there means checking a package this lane
  does not otherwise touch.
- Several test files run only under `pkgload::load_all()`, because they
  reach package internals by bare name, and `test_check()` hides it.
  Measured under a bare `test_file()` against the installed package:
  tests/testthat/test-compat-register.R fails all 27 tests at
  `frmtmb_compat_contrib`; tests/testthat/test-structure.R takes 13
  errors at `fam_structure`, `structure_unit`, `structure_gate`,
  `structure_allows`, `structure_group_codes` and
  `check_structure_block`; and frmtmb.sample's test-loo.R takes 6 at
  `mixture`, `mvbf`, `cumulative`, `rescor_matrix` and
  `expose_functions`, none of which the file attaches. Qualifying the
  calls would make each file runnable on its own.

## Recorded by the skew-normal start lane (wt-skewinit), 2026-09-22

Found while fixing the `alpha = 0` stall. Neither item is done, and
neither is the defect that lane fixed.

- **A `sigma` start taken from `sd(y)` ignores the `mu` predictor, in
  every family that does it.** `sigma` in these families is the
  CONDITIONAL standard deviation, so the marginal `sd(y)` overstates it
  by whatever the mean structure explains. On the wt-skewinit design
  the start was `log(sd(y)) = 0.662949` against a fitted
  `sigma_(Intercept)` of `-0.014760`, a factor of two. `skew_normal()`
  now starts from the residual sd; `gaussian()`, `student()`,
  `lognormal()`, `exgaussian()`, `asym_laplace()` and every other
  family written as `sigma = function(y, aterms) stats::sd(y)` in
  `R/families.R` still do not. For those families the optimizer
  recovers without help, so this is cost and conditioning rather than a
  wrong answer, and changing it moves the optimizer path of nearly
  every fit in the package. `dev/skewinit-starts.R` measures what the
  change is worth on skew_normal, where it does change the answer:
  arm `res_sdy` (residual alpha, marginal sigma) stalls on 1 of 40
  seeds of the dead-draw stream, arm `res_sdr` (residual sigma too) on
  0 of 40.

- **`skew_normal()` has a second local optimum away from `alpha = 0`,
  which the stationary-point escape does not reach.** Construction, in
  `dev/skewinit-falsealarm.R`, arm `mild_a05` seed 22: true `alpha`
  0.5, n = 200, the data from that script's `make_mild(22, 0.5, 200)`.
  The fit converges at `alpha = -0.551280` with logLik -351.81810,
  while `sn::selm` reports -351.79426 at `alpha = 0.672091`, a
  shortfall of 0.023841. The escape does not fire because `|alpha|` is
  0.55, far outside the 0.05 stationary-point window, and it should
  not: this is ordinary multimodality, not the singular point, and a
  window wide enough to catch it would refit on genuine optima. It is
  1 of 240 fits over six arms. Filed rather than fixed because the
  remedy is a multi-start policy, which is a decision about every
  family rather than about this one.

- **`diagnose()` says nothing about a dpar whose information is
  singular.** A `skew_normal()` fit can converge cleanly with `alpha`
  0.2477 and a standard error of 5.09, and every check passes:
  convergence 0, `pdHess = TRUE`, smallest covariance eigenvalue
  0.0045, max abs gradient 4.3e-05. Construction in
  `dev/skewinit-claims.R`, seed 3 of that script's `make_sym()` at
  n = 50. The only signal the user gets is the standard error. When
  `alpha` runs away instead, `diagnose()` does speak (convergence 1,
  `pdHess = FALSE`, eigenvalue -0.378), so the gap is the FINITE
  weakly-identified case. Whether `diagnose()` should carry a
  "this parameter is barely identified" line is a design question about
  every family, not about skew_normal, which is why the wt-skewinit
  lane did not add one.

## Added by wt-skewinit after punch round 1, 2026-09-22

- **`quadrature = TRUE` calibrates the Gauss-Kronrod tape at the
  Laplace optimum and never recalibrates it.** When that Laplace
  optimum is a stall, the stationary-point escape re-optimizes a tape
  frozen at the wrong point: on acceptance seed 1 with `(1 | g)` it
  recovers 24.262 of 24.272 units and still lands 1.003e-02 below the
  unstalled fit, against 9.5e-09 for the same model under Laplace.
  Reproduction: `dev/skewinit-punch1c.R`. The fix is to escape BEFORE
  `quad_fit()` calibrates, or to recalibrate after it, both of which
  are changes to `quad_fit()` rather than to the escape.

- **nlminb reports "false convergence (8)" on fits that satisfy
  frmtmb's own gradient criterion.** Construction in
  `dev/skewinit-m3.R`: `set.seed(8)`, n = 200, `x ~ N(0, 1e4^2)`,
  `y = 0.001 x + rskewnorm2(n, 0, 1.5, 4)`. The fit reaches
  logLik -342.388188574, which equals `sn::selm` to 2.3e-10, at a
  finite alpha of 5.79, with max abs gradient 4.81e-04 against a
  `grad_tol` of 1e-3, and warns. Measured causes ruled OUT: autoscale
  (`par_units` is NULL, it never engaged) and the escape (it never
  fired); the cause is which start the optimizer arrives from, and
  five extra restarts do not clear it. frmtmb reports the optimizer's
  code verbatim, so the question is whether a nonzero PORT code should
  still warn when the package's own convergence criterion is met. That
  is a policy decision about every family and every fit mode, which is
  why wt-skewinit did not change it.

- **`hmm()` does not carry a component's `post$stationary`.**
  `extensions/frmtmb.latent/R/hmm.R` copies `comp$init_dpars[[dp]]`
  into the composed family but nothing copies the stationary
  declaration, so `hmm(K, skew_normal())` keeps the alpha = 0 stall
  unguarded. `mixture()` was fixed in core in the same round and is
  the model to copy: build the declaration alongside `init_fns`, keyed
  by the same `paste0(dp, k)` name. Left to the extension's own lane.

## Added by wt-skewinit after punch round 2, 2026-09-22

- **A least-squares start for an intercept-less design helps several
  families and is NOT safe to ship as a general rule.** wt-skewinit
  briefly placed one for every family and had to withdraw it. The
  upside is real: on ordinary cell-means factor models with no scaling
  pathology the reviewer measured gaussian **+263.96**, Gamma
  **+624.94** and exgaussian **+16.85** log-likelihood units, and
  seven jobs in six families improved.

  The counterexample is why it is filed rather than shipped.
  `ypois ~ 0 + xt` with `xt` on a 1e-6 scale, n = 250, seed 23:
  the placement starts the coefficient at 1.6e6, nlminb reports
  "X-convergence (3)" before the fit moves because the relative step
  is negligible at that magnitude, and the answer is **1967.03 units**
  below `stats::glm()`, with convergence 0 and no warning. A six-start
  sweep from 0 to 1e7 reaches glm's value, so it is the start and not
  multimodality. A scale sweep puts the crossover between 1e-4 (matches
  glm) and 1e-6. The value is right on the PREDICTOR scale and ruinous
  on the PARAMETER scale.

  A safe version needs one of: a guard on the resulting coefficient
  magnitude, or keeping the new start only when its objective beats the
  old one, which costs an extra evaluation but cannot lose. Either is a
  change to `make_start()` for every family and wants its own lane.
  Batteries to reuse: `dev/skewinit-noint.R` and its comparator.

- **frmtmb sits below `glm()` on a no-intercept poisson with a
  1e-6-scaled covariate, on the UNCHANGED build.** Same construction as
  above: `dev/skewinit-noint.R` job `pois_s6` gives -512.543671042038
  against glm's -449.965387300331, a gap of **62.578284**, and
  `rellib-r3` and the lane build agree to the last bit, so this is not
  wt-skewinit's. At scales 1, 1e-2 and 1e-4 frmtmb matches glm exactly.

  **`frmtmb_control(autoscale = TRUE)` fixes it completely.** An
  earlier version of this entry said autoscale does not rescue it,
  which was WRONG and would have sent the next reader away from the one
  remedy that already works; the claim is withdrawn here rather than
  edited out. Measured in `dev/skewinit-punch3.R`, log
  `skewinit-log/punch3-lane.txt`, autoscale off against on, both
  compared to `glm()`:

  | covariate scale | autoscale off | autoscale ON |
  |---|---|---|
  | 1 | 0.000000 | 0.000000 |
  | 1e-2 | 0.000000 | 0.000000 |
  | 1e-4 | 0.000000 | 0.000000 |
  | 1e-6 | **-62.578284** | 0.000000 |
  | 1e-8 | **-62.578284** | 0.000000 |

  So the open question is narrower than it looked: not "why is the fit
  wrong" but "why does the default not turn autoscale on here", since
  the default leaves the 62.58 gap in place. The fit still reports
  convergence 0 and says nothing, which is the part worth fixing.

- **The same scale mechanism survives INSIDE a dpar that declares a
  stationary point, which is the one place wt-skewinit still places a
  start on an intercept-less design.** Reproduction in
  `dev/skewinit-punch3.R` part 2: `bf(y ~ xr, sigma ~ 1,
  alpha ~ 0 + xs)` with `xs` the covariate rescaled, n = 250, seed 1.
  At scale 1e-6 the placement starts the alpha coefficient at
  -6.769e+05 and the fit stops at **-346.8767477** with
  "X-convergence (3)", convergence 0 and no warning, against
  **-340.9903976** from a sweep of predictor-scale starts: **5.886
  units short, silently**. At 1e-9 the same, starting at -6.769e+08.
  At scales 1 and 1e-3 the default equals the sweep exactly.

  Two facts that keep this in proportion, and both belong here. Base
  is **-356.8824988** on the same design at every scale, so the lane
  is **9.65 to 18.91 units BETTER** than base, not worse. And over 48
  paired fits (12 seeds x scales 1, 1e-3, 1e-6, 1e-9,
  `dev/skewinit-alphascale-cmp.R`) the lane is better on **48 of 48**,
  worst gain +3.70194, best +26.444, mean better at every scale, with
  0 nonzero convergence codes on either build. The reviewer's own
  construction gave 33 better and 15 worse with a worst deficit of
  -0.012426; the designs differ, so both are recorded rather than
  reconciled. The defect is the shortfall against a predictor-scale
  sweep, NOT a regression against base.

  Why no test catches it: the lane's complement test uses a
  well-conditioned cell-means factor, where every column is 0 or 1 and
  the placement cannot produce a large coefficient. A test for this
  needs an intercept-less alpha design on a CONTINUOUS covariate whose
  scale is swept, asserting the default fit against the best of a
  predictor-scale start sweep on the same data. Filed beside the
  poisson entry because they share one mechanism: a start that is
  right on the predictor scale and too large on the parameter scale.

## Reference

Full agent report with per-item repro sketches and issue links:
sourced from glmmTMB/lme4/brms GitHub test suites and NEWS, 2026-08-31.
Key files to mirror: glmmTMB test-predict.R, test-formulas.R,
test-varstruc.R, test-offset.R, test-weight.R, test-NAhandling.R (lme4),
tests.brmsformula.R (brms).
