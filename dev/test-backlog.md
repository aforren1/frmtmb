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

- DECIDED NOT TO DETECT, lane wt-predfix (`dev/predfix-findings.md`
  item 2): no statistic calibrated here separates a flat direction
  from a genuine heavy tail or from extrapolation. The share of a row's
  absolute deviation carried by its single largest draw is 0.553 to
  1.000 on flat directions and 0.517 to 0.600 on well-identified
  lognormal fits with sigma 5.5 and 6.5. The link-scale se.fit is 10.0
  to 11.3 on the near-collinear 1e-2 rows and 16.4 to 20.1 on a
  legitimate poisson extrapolation to x = 300. Those collinear rows
  ARE design-space extrapolation (leverage 3.8e4 to 4.5e4); the only
  flat case that is not, an all-zero poisson cell, already raises the
  non-finite warning. The entry as filed:
  `predict()` on a draw that is FINITE but absurd. Masking a
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

- DONE, lane wt-resmooth (`dev/resmooth-findings.md`). `re_formula = NA`
  keeps every smooth. brms's `posterior_epred(re_formula = NA)` is
  BITWISE its `re_formula = NULL` on all three constructions and differs
  on `s(x) + (1 | g)` (`dev/resmooth-brms.txt`), and frmtmb now agrees on
  both: the row sd over sigma is 1.000 for each of the three, against
  1.260, 1.741 and 1.919, and 1.258 for `s(x) + (1 | g)` before and
  after (`dev/resmooth-before.txt`, `dev/resmooth-after.txt`). brms
  REFUSES a new level of the grouping factor of a smooth whatever
  `allow_new_levels` says, so frmtmb refuses too, and a missing grouping
  column with it: `re_formula = NA` is no longer a way out of either.
  Two open items go with it. The partial-formula refusal beside a factor
  smooth is now conservative rather than forced, and section 4 of the
  findings PROPOSES dropping it, for the user to decide. And the lane
  found and fixed a separate defect on the way: the finite-difference
  `Est.Error` of a category probability differenced no smooth
  coefficients at any `re_formula`, which put it 74 percent from a Monte
  Carlo reference on an ordinal `s(x)` fit (findings section 5).
  The entry as filed:
  `re_formula = NA`
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

- DONE, lane wt-predfix: the scalar route now warns, classed
  `frmtmb_modes_conditional_se`, as the finite-difference route and
  `predict()` do. The entry as filed:
  A quadrature fit's SCALAR standard errors leave out the group-effect
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

- The finite-difference `Est.Error` route costs one pair of model
  evaluations per differenced coefficient, and a smooth or `gp()` block
  can have as many coefficients as there are observations. Filed by lane
  wt-resmooth after its nits round fixed the part that was fixable.
  Counted, not timed (`dev/resmooth-batchcost.R`, log
  `dev/resmooth-batchcost-after.txt`; wall clock on three arms in
  `dev/resmooth-cost3-*.txt`):
  * `fitted(newdata = 3 rows, re_formula = NA)` on a `cumulative()`
    `gp(x)` fit with 160 coefficients makes 329 model evaluations and
    takes 1.56 s against 0.0497 s on 0.64.0, about 31x. The three rows
    are OFF the fitted `gp()` positions, so every kriging row loads all
    160 columns and no batch of two columns can be attributed by row.
  * `fitted()` on `s(x, k = 8) + (1 | g)` at 40 levels is 1.70x at
    `re_formula = NULL` and 1.48x at `NA` for the same reason at a
    smaller size: 25 and 23 evaluations against 11 and 9. Every row
    loads every basis column of a smooth, so the block cannot batch.
  * A CONTROL with no smooth is 11 evaluations on both arms and moves
    0.0731 s to 0.0537 s, which bounds the wall-clock noise at about
    1.36x and is why the counts are the instrument here.
  The route that would fix it: a category probability depends on `b`
  only through `eta`, and `eta` is LINEAR in `b`, so
  `dp/db = (dp/deta) Z` needs one evaluation pair per ROW rather than
  per coefficient, and `Z` is already built by `lp_delta_A()`. That
  needs a seam `fit_fd_se()` does not have, because it differences the
  composite `f` and cannot perturb `eta` alone. A cheaper partial
  measure this round did not take: an `fs` block's columns partition by
  LEVEL, so a row loads only its own level's basis columns and its 50
  columns could be 5 batches rather than 50. `re_b_batches()`'s
  whole-block test cannot see that; a greedy grouping on the column
  conflict graph would, at the cost of an `O(ncol^2)` sparse product
  that has to stay off the path a fit with thousands of levels takes.


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

- DONE, lane wt-predfix: the draws methods refuse an unseen level with
  or without the flag, and say the flag is refused too; core's
  new-level errors are classed `frmtmb_new_levels` so the draws side
  can catch them. The entry as filed:
  Draws hint that leads to a refusal (found 2026-09-22 by the round-3
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

- DECIDED 2026-09-24 (user): KEEP REPORTING nonzero optimizer
  convergence codes, "I'd rather people think about their fits". The
  measured case for each option is in `dev/predfix-findings.md` item 7;
  nothing changed. The entry as filed:
  **nlminb reports "false convergence (8)" on fits that satisfy
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

- DONE, lane wt-predfix: `frmtmb_control(autoscale = )` defaults to
  `NULL`, which engages the pre-fit when a qualifying column's sd is
  below 1e-3 and its coefficient is an outer parameter. The job
  below reaches glm() exactly. The entry as filed:
  **frmtmb sits below `glm()` on a no-intercept poisson with a
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

- DONE by the same default, lane wt-predfix: at 1e-6 and 1e-9 the
  default now engages and beats the predictor-scale sweep on 24 of 24
  fits, by 0.115 to 7.14 units (`dev/predfix-skewalpha.R`). The entry
  as filed:
  **The same scale mechanism survives INSIDE a dpar that declares a
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

## Recorded by the phase-3a lane (wt-phase3a), 2026-09-24

- **frmtmb.spline: a penalty on negative hazards, as rstpm2 has.** Since
  the user decision of 2026-09-24, `rp_floored()` refuses a
  non-positive `d(eta)/d(log t)` on an EVENT row and only warns on a
  CENSORED row, which is where a cure-fraction design with
  `gamma1 ~ arm` puts it: 27 of 80 such fits flag rows past an arm's
  last event, with the fitted survival rising by 5.2e-04 to 0.13
  (`dev/phase3a-findings.md`). rstpm2 prevents the rise during the fit
  instead, with an adaptive quadratic penalty on negative hazards
  (`kappa.init = 1`, raised up to `maxkappa = 1000`). A comparable
  option on `royston_parmar()` would remove the warning at its source,
  at the price of a penalized likelihood that `logLik()` would have to
  report as such. Not started; the size of the change to the density
  and to `logLik()` is the first thing to measure.
- **frmtmb.spline: `frm_curve(dpar = "mu")` refuses on a
  `gamma1 ~ x` fit.** On the cure design above (seed 20260932,
  `df = 4`), `frm_curve(fit, dpar = "mu")` stops with "the assembled
  covariance of this grid disagrees with frm_linpred(se.fit = TRUE) by 1
  relative", on the released 0.7.0 build as well. On that one fit
  `dpar = "gamma1"` answers, and `dpar = "mu"` on a proportional-hazards
  fit of the same data answers; but on the reviewer's cure design
  `dpar = "gamma1"` refuses too, on 3 of 4 fits, identically on the
  released build. So the refusal reaches both dpars of a fit with a
  covariate on `gamma1`. It refuses rather than answering wrong, so it
  is filed rather than fixed; the cause is not yet known.

## Recorded by lane wt-predfix, 2026-09-24

- DONE: `fitted()` on a multivariate fit refused while `predict()`
  answered (`dev/shapes-findings.md`, "fitted(mv)"). It returns brms's
  `n x 4 x nresp` now.
- DONE: `vcov()` on a `REML = TRUE` or `profile = TRUE` fit with no
  random effect and no free dispersion died in `vcov_estimated()`
  (`dev/correct-findings.md`, "Found and NOT fixed" 1). TMB leaves the
  joint precision unnamed when there is no outer parameter.
- DONE in punch round 1: a random slope on a badly scaled column,
  found by the reviewer (`y ~ x + (1 + x | g)`, gaussian, scale-1
  optimum -311.327238; at sd 1.86e-3 seed 1 was 3.18 units short with
  code 0, and even with autoscale engaged at 1e-6 it was 4.87 and 1.79
  short on two seeds; base identical). Autoscale rescaled X but not Z,
  so the slope's log sd started where the scale-1 fit starts it. The
  pre-fit now rescales the Z column of `us()`, `diag()`, `us_t()` and
  `diag_t()` slopes too, and the default engages for a slope column
  below sd 0.05 (`dev/predfix-p1-slopecal*.R`, 736 fits).
- OPEN, from the same calibration: a random slope whose covariance
  structure is not `us()`, `diag()` or a Student-t block (`cs()`,
  `ar1()` and the rest), or whose column is not also a fixed effect
  (`y ~ 1 + (0 + x | g)`), is not rescaled, so it keeps the stall. And
  just above the 0.05 threshold the plain fit fell short by at most
  2.2e-4 units (spread 0.06, 64 fits). Filed by lane wt-predfix.
- OPEN, pre-existing, found by the wt-predfix reviewer: `frm_sample()`
  on a one-parameter model (`y ~ 0 + x`) fails inside rstan with "no
  more scalars to read", on the 0.62.0 base as well. Not reproduced by
  the lane, which only files it.
- DONE in punch round 1 (assigned by the organizer, found by the
  wt-simnewdata reviewer): `predict()` and draws `posterior_predict()`
  drew every `cs()` model as if the term were absent, in sample and at
  newdata, silently (max |z| 1216.7). `dev/predfix-findings.md`, the
  cs() section.
- DONE by lane wt-csfactor (`dev/csfactor-findings.md`), A SILENT WRONG
  ANSWER, pre-existing (0.62.0 through 0.64.0): `cs()` on a factor is
  fitted and predicted on the factor's INTEGER CODES.
  `ord_cs_values()` calls `as.numeric()` on the factor,
  so the model is linear in the codes (logLik -446.9233 against
  -446.6361 with the two dummy columns written out), and a `newdata`
  factor is re-coded from its own levels: a single row
  `factor("c")` gets level a's probabilities (0.1876, 0.6725, 0.1399)
  instead of (0.7430, 0.0883, 0.1687), from `fitted()` and `predict()`
  alike. Found by the wt-predfix reviewer, punch round 2
  (`dev/predfix-review2/r2-cs2.R`, seed 405). brms 2.23.0 builds
  treatment-contrast dummy columns, `fcb` and `fcc`, with one `b` row
  each (`dev/predfix-brms-csfactor.R`, `dev/predfix-log/brms-csfactor.txt`),
  so the fix is to build `cs()` from `model.matrix()` with the fit's
  contrasts and levels, in the frame and at `newdata`. That is what was
  done, one `bcs<j>` per design COLUMN; `tests/testthat/test-cs-factor.R`
  pins it and was seen to fail 27 of 46 with 3 errors on the 0.64.0
  reference build.
- The frmtmb.sample floor must move to the core that carries
  `frmtmb_new_levels` and exports `cs_offsets_add()`: against 0.62.0
  core the draws check falls back to its flagged path, core's old hint
  reaches the caller again, and `posterior_predict()` cannot find
  `cs_offsets_add()`.

## Filed by wt-mvprior after punch round 1, 2026-09-24

Found by the punch-round reviewer, reproduced on the lane build and on
brms 2.23.0's `stancode()` with the same data
(`dev/mvprior-filed.R`, `dev/mvprior-log/filed.txt`). All five are
pre-existing and none is a wrong answer: each is a loud refusal where
brms accepts the model. Not fixed by that lane.

### Closed at 0.64.0

All five shipped in the brms-parity round. The entries are kept with
what closed each one; NEWS.md's 0.64.0 section has the detail.

- **Three `bf()` summed with `+` fail at construction.**
  `bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ w) + set_rescor(FALSE)` warns
  `Incompatible methods ("+.frmtmb_mvformula", "+.frmtmb_formula")`
  and stops with R's "non-numeric argument to binary operator". Two
  `bf()` work, and `mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ w))` works.
  brms accepts the three-term sum.
  SHIPPED: any number of `bf()` can be summed. The formula and the
  multivariate formula had two `+` methods and now have one, and `lf()`
  and `nlf()` take brms's `resp =` (lane mv, `dev/mv-findings.md`).
- **`cumulative()` in a multivariate model** is refused: "Families with
  extra parameters ('cumulative') are not supported in multivariate
  fits yet". brms accepts `bf(y1 ~ x) + bf(o ~ x, family =
  cumulative())`.
  SHIPPED: `cumulative()`, `sratio()`, `cratio()` and `acat()`
  responses keep their own thresholds and can share `|ID|` group
  effects. `cox()` and `mixture_mvn()` stay refused, now with the
  reason (lanes mv and thres).
- **`me()`** is not a function: `bf(y1 ~ me(xe, sde))` stops with
  `could not find function "me"`, R's message rather than a designed
  refusal. brms accepts it.
  SHIPPED: `me(x, sdx)` and `me(x, sdx, gr = g)` with brms's meaning
  and brms's names, agreeing with brms 2.23.0's Stan program to 1e-13
  on five shapes. See `?frmtmb-me` (lane me, `dev/me-findings.md`).
- **`0 + Intercept`** reads `Intercept` as a data column: "The model uses
  `Intercept`, which is not a column of `data`". brms accepts
  `y1 ~ 0 + Intercept + x` as its uncentered intercept.
  SHIPPED: brms's reserved intercept, as an uncentered class `"b"`
  coefficient, in the location, distributional and nonlinear-parameter
  formulas. `bf(center = FALSE)` is the same mechanism, and ordinal
  families refuse it as brms does (lane icpt0,
  `dev/icpt0-findings.md`).
- **`student()` with `set_rescor(TRUE)`** is refused: "rescor = TRUE
  requires all responses to be gaussian (got: student)". brms accepts
  it (a multivariate student-t).
  SHIPPED: brms's multivariate Student-t, one shared `nu`, a sigma per
  response, checked against `mvtnorm::dmvt()` and brms's own
  `log_lik()` (lane mv, `dev/mv-findings.md`).

## Filed by wt-mvprior after punch round 2, 2026-09-24

Reported by the punch-round-2 reviewer with probes in
`dev/mvprior-review2/` (`r2-probe-*.txt`, `r2-brms-defaults.txt`,
`r2-sample-*.txt`), measured on base 0.62.0 as well. Filed as reported;
the lane did not re-measure them. None is a wrong answer from the fit.

### Open - medium

- **frm_sample() gives a mixture's `sigma1`, `sigma2` and `theta` no
  default**, where brms uses `student_t(3, 0, mad)` on each sigma and
  `dirichlet(1)` on the mixing proportions, and the announcement does
  not list them among the slots it leaves flat.
- **Three spellings frmtmb accepts and brms refuses:** `cor` with
  `dpar = "mub"`; a bound (`lb`/`ub`) with `coef`; and class `b` with
  `coef = "Intercept"` on `y ~ 1`.
- **One spelling brms accepts and frmtmb refuses:** `Intercept` with
  `coef = "2"` on `cumulative()`, brms's per-threshold intercept.
- **`cratio()` with `disc ~ z`** is unsupported.
- **default_prior() rows brms has and frmtmb lacks:** `sd` rows with
  `coef = "Intercept"`, the per-threshold ordinal `Intercept` rows
  (`coef = "1"`, `"2"`, ...), and a `theta1` row on a mixture whose
  sigma is modeled.

## Filed by wt-phase3b after punch round 2, 2026-09-24

- **Core seam: a log-difference slot for interval and truncation
  masses.** frmtmb forms an interval-censored row as
  `log(F(y2) - F(y))` and a truncation window as `log(F(ub) - F(lb))`,
  both on the probability scale from one `lcdf` slot. That loses every
  digit when the two values agree to their error, and holds the mass at
  1e-300 or above. frmtmb.eam works around it for `wiener()` by
  recognizing the upper edge inside its own `lcdf` and returning the
  mass there (`ddm_rt_linterval_b()`, 4.9e-11 on the log against 751
  Rmpfr intervals), which is why left or interval censoring cannot be
  combined with `trunc()` on that family. A family slot
  `linterval(lo, hi, dpars, aterms)` returning the log mass directly
  would remove both the workaround and the refusal. Test: an interval
  holding 1e-40 of the mass before it, against Rmpfr, through the core
  path.
- **Core seam: refuse a missing `dec()` on a censored row by name.**
  For `wiener()` a left- or interval-censored row reads `dec()` (the
  boundary is known), a right-censored row does not. frmtmb's
  `na.action` drops a row with an NA `dec()` before the family sees it,
  with "1 row removed because of missing values"; the family cannot
  refuse it or keep a right-censored row whose `dec()` it would not
  read. Needs a hook that lets a family say which addition terms a row
  needs, by censoring code. Test: a right-censored row with NA `dec()`
  kept, a left-censored one refused by name.
- **Under-coverage of log sd(log bs | s) on the base build.**
  `wiener()` at 30 subjects by 400 trials, no censoring, no
  contaminant, fresh seeds 2001 to 2070 on frmtmb.eam 0.10.0: the Wald
  interval for the log subject sd of log boundary separation covers on
  60 of 70, 85.7 percent, Wilson [75.7, 92.1], mean SE 0.1307 against
  an sd of 0.1426 (SE/sd 0.917). Every other quantity covers on 65 or
  66 of 70. From the review of punch round 1
  (`dev/phase3b-review2/cov-summary-fresh.log`). Not caused by items
  3.4 or 3.5; not investigated. Test: coverage of that quantity on
  independent seeds, with enough replicates for the Wilson interval to
  exclude or include 95.

## Filed at the 0.64.0 release (2026-09-25)

The brms-parity round's "Left open" list, from
`dev/parity-round-20260925.md`. Each bullet names the lane that found
the item and the findings file with its measurement. Six lanes started
on 2026-09-28; the items they closed at 0.65.0 say so.

### Open - high priority

- **The convergence check fires on correct fits.** "Large maximum
  absolute gradient" appears on ordinal and multivariate fits whose
  likelihood identities hold to 1e-12, and on fits with an active
  bound. The 1e-3 threshold is absolute, so it does not scale with the
  data or the parameter count. Found by lanes thres, mv and arcov
  (`dev/thres-findings.md`, `dev/mv-findings.md`,
  `dev/arcov-findings.md`). Done at 0.65.0, lane wt-gradcheck.
- **Refits can lose a threshold.** `frm_bootstrap()` and other refits
  recount the ordinal thresholds on each simulated data set, so a
  category that is empty in one draw silently changes the model that
  draw fits. Lane thres (`dev/thres-findings.md`). Done at 0.65.0,
  lane wt-thresrefit: the bootstrap never recounted; `influence()` did.
- **`y ~ x + cs(x)` is not identified and is not refused.** The
  category-specific effect and the population effect of the same
  column carry the same information. Lane sratio
  (`dev/sratio-findings.md`). Done at 0.65.0, lane wt-csfactor.
- **`draw_prior_entry()` on an unordered threshold vector** may copy
  one draw into every threshold. Not known to be reachable, so prove
  reachability or unreachability by construction before fixing it.
  Lane sratio (`dev/sratio-findings.md`). Done at 0.65.0: unreachable,
  and now refused, lane wt-thresrefit.

### Open - medium

- **Adding a family to a multivariate formula.** In
  `bf(o ~ x) + cumulative() + bf(y ~ x) + gaussian()`, frmtmb fills
  only the responses that have no family, so `o` stays ordinal, while
  brms gives the last family to every response. frmtmb behaved this way
  before the round and the behavior is now documented. CLOSED
  2026-09-29: the user keeps frmtmb's rule, because brms's silently
  overwrites an ordinal response's family. Lane mv
  (`dev/mv-findings.md`).
- **Priors on `me()` hyperparameters.** brms's classes `meanme`, `sdme`
  and `corme` are refused by name, not implemented. Lane me
  (`dev/me-findings.md`).
- **REML with `me()` in `mu` is approximate**, and is registered as
  conditional. The same argument applies to the existing `mi()`
  predictor, which is registered as working, so one of the two
  registrations is wrong. Lane me (`dev/me-findings.md`).
- **frmtmb.sample `log_lik()` and `loo()` for `cov = FALSE` ARMA** are
  not built. Lane arcov (`dev/arcov-findings.md`). Done at 0.65.0,
  lane wt-arcovsample.
- **brms's latent-residual AR for non-gaussian families** is not built;
  the 0.64.0 form is brms's residual-regression one. Lane arcov
  (`dev/arcov-findings.md`).
- **`gr(g, by = f, cov = A)`** is refused by name, because brms
  correlates the by-levels through `A`, which is not a by-split. Lane
  grby (`dev/grby-findings.md`).

### Upstream

- **brms `posterior_predict_hurdle_negbinomial()` does not draw the
  zero-truncated negative binomial.** Lane fams
  (`dev/fams-findings.md`). Drafted nowhere yet; see
  `dev/round-handoff.md` for the other unfiled upstream reports.

### Closed at 0.64.0 by lane wt-records, 2026-09-28

The round also listed two older record inconsistencies, both closed:

- frmtmb.eam's tests called `frmtmb.sample::` while its Suggests was
  empty of it. `frmtmb.sample` is now suggested, and
  `.github/workflows/check-frmtmb-eam.yaml` installs it from the
  checkout and lists `extensions/frmtmb.sample/**` in `paths:`, which
  `tests/testthat/test-ci-siblings.R` asserts.
- `codemeta.json` said 0.50.0 and lacked `ordinal`. Regenerated with
  `codemetar::write_codemeta()` at 0.64.0.


## Filed by wt-thresrefit, 2026-09-29

### Open - medium

- **A character-coded ordinal response with non-numeric labels dies in
  an internal error.** `frm(bf(y ~ x), family = cumulative())` on a
  character `y` whose values are not numeric text stops with "missing
  value where TRUE/FALSE needed", which names neither the response nor
  the requirement. Reported by the punch-round-1 reviewer as nit 3, on
  the ground that this lane's new `?influence.frmtmb_fit` and `?frm` text
  tells readers a character coding gives the same influence table as
  integers, which points them at the refusal. Re-measured rather than
  filed as reported (`dev/thresrefit-p2-charresp.R`, seed 2501, n = 60,
  four categories): with labels `"none" < "mild" < "moderate" <
  "severe"` the message above, on the lane build and on base 0.64.0
  alike, so it is pre-existing and was not fixed here; with labels
  `"1"..."4"` the fit runs; and `factor(..., levels = , ordered = TRUE)`
  over the SAME four non-numeric labels fits with `n_tau = 3` on both
  arms. So the capability is there and only the spelling fails, loudly
  but anonymously. Test: the refusal names the response and says an
  ordinal response is positive integers or an ordered factor, and the
  ordered factor over the same labels still fits.
## Filed by lane wt-csfactor after punch round 1, 2026-09-29

Both found by the wt-csfactor reviewer
(`dev/reviews/2026-09-29-csfactor.md`), both pre-existing on the 0.64.0
reference build, neither this lane's.

### Open - medium

- **`fitted(mv, newdata = )` errors on ANY multivariate ordinal fit**
  with "values must be length 1, but FUN(X[[1]]) result is length 0",
  R's own `vapply()` message rather than a designed refusal. Measured
  with `cs()` and without it, same message either way, so it is not
  about `cs()`; in-sample `fitted(mv)` returns the `n x 4 x K` array
  both ways (300 x 4 x 6), and `predict(mv, newdata = )` refuses with a
  designed message. Construction:
  `frm(bf(yo ~ x + cs(fc)) + bf(yo2 ~ z + cs(fc)), sratio())` on seed 77,
  n = 300, then `fitted(mv, newdata = data.frame(x = 0, z = 0, fc =
  factor("a", levels = c("a","b","c"))))`
  (`dev/csfactor-rev-docs.R`, `dev/csfactor-rev-log/rev-docs-lane.log`).
- **`y ~ mo(m) + m` is not identified and is fitted in silence**, with no
  `cs()` term anywhere. A monotone function of `m`'s categories lies in
  the span of `m`'s own dummies, so the two coefficient blocks share a
  flat direction. Measured on both builds, seed 1907, n = 300, `m` an
  ordered factor with 4 levels (`dev/csfactor-p2.R`,
  `dev/csfactor-log/p2.txt`): `mo(m) + m` gives df 8 and logLik
  -328.534572784 where `m` ALONE gives df 5 and -328.534572790, so three
  extra parameters buy 6e-09 of log likelihood; the covariance matrix has
  a smallest eigenvalue of -6441 and 6 of 6 standard errors are
  non-finite, with only R's own `In sqrt(diag(V)) : NaNs produced` to show
  for it. Written with `m` as an UNORDERED factor the eigenvalue is -6724
  and the standard errors happen to stay finite, which is worse: nothing
  at all is visible. Lane wt-csfactor's `check_cs_identified()`
  deliberately does NOT refuse this, being about `mo()` rather than about
  `cs()`, and its comment points here. The fix belongs with a rank check
  over the FILLED design rather than the assembly-time one, since the
  `mo()` column is a zero placeholder at assembly.

### Open - minor

- **`variables()` on a `frm_sample()` draws object of an ORDINAL fit
  lists the raw internal names.** It gives `tau_raw_1`, `tau_raw_2` where
  `variables(fit)` gives `b_Intercept[1]`, `b_Intercept[2]`, and
  `bcs2_1`, `bcs2_2`, `bcs3_1`, `bcs3_2` where `variables(fit)` gives
  `bcs_fcb[1]` and the rest. It does the `tau_raw` half for an ordinal
  fit with NO `cs()` term, so the defect is about ordinal draws and not
  about `cs()`: `variables(ds)` for `sratio, yo ~ x` gives
  `b_x tau_raw_1 tau_raw_2 lp__`, while `fixef(ds)` for the same fit
  gives `Intercept[1] Intercept[2] x`. A gaussian fit's draws are
  correct (`b_Intercept b_x sigma lp__`). So `fixef()` on the draws
  object goes through `draws_fixef_ordinal()` and `variables()` does
  not. `?variables` now records the exception rather than claiming the
  draw columns follow the same convention
  (`dev/csfactor-rev-docs2.R`, `dev/csfactor-log/sample.txt`). Done at
  0.66.0, lane wt-sampfix: the draws carry brms's names; the covariance
  parameters (`theta_1`) are the one gap left, filed below.
- **The covariance machinery is not bound-aware.** At a constrained
  optimum the UNCONSTRAINED Hessian can be indefinite, because the fit
  is only a minimum along the feasible directions. `sdreport()` reads
  the full Hessian, so such a fit reports `pdHess = FALSE` and NaN
  standard errors while being exactly right. Constructed
  (`dev/gradcheck-01-construct.R`, seed 101,
  `dev/gradcheck-rev-10-flip.R` for the cap sweep): gaussian
  `y ~ x` with a true slope of 2 under
  `set_prior("", class = "b", ub = 0.1)` stops with `x` on the bound;
  the full Hessian's eigenvalues are 693.40, 73.10 and -29.63, and the
  same Hessian restricted to the two parameters no bound holds is
  positive definite, which is why the convergence check can measure a
  headroom there. Whether it bites depends on how far the bound is from
  the unconstrained optimum, not on bound-awareness: over caps 0.1 to
  1.99 on that design the bound holds `x` at every cap while `pdHess`
  is FALSE at 0.1 and 0.5 and TRUE at 1 and above. Fix: restrict the
  reported covariance to the free subspace, and say in the report that
  a bound-held parameter has no standard error rather than returning
  NaN for every parameter. Test: on the `ub = 0.1` fit, the free-set
  Hessian is positive definite while `pdHess` is FALSE, and the free
  parameters get finite standard errors.
## Filed by wt-arcovsample after punch round 1, 2026-09-29

Found by the punch-round reviewer
(`dev/reviews/2026-09-29-arcovsample.md`, section 9) and by this lane
while probing `log_lik()` for brms's `cov = FALSE` ARMA. Every item was
reproduced on the LANE build and on the base build `rellib-r3`, so none
is caused by that change. Item 1 is a wrong ANSWER and item 3 fails on
every model, which is what makes them worth tests; item 2 is a bare
internal error. Not fixed: all four are outside that lane's gap.

### Open - high

- **`posterior_predict()` and `posterior_epred()` on
  `frm_sample(laplace = TRUE)` draws return NaN silently.** The draws
  hold only the outer parameters, because the Laplace route integrates
  the random effects out, and the predictive methods evaluate the model
  as if `b` were there. `dev/arcovsample-rev-11-laplace.R` (seed 1212,
  N = 30 in 6 groups of 5), re-run by the lane into
  `dev/arcovsample-log/punch1-laplace-{lane,ref}.txt`, byte-identical on
  the two builds: `y ~ x + ar(t, g) + (1 | g)` gives **4500 non-finite
  cells of 4500**, `y ~ x + (1 | g)` gives **3000 of 4500**, and the
  only signal is a repeated base warning "NAs produced". `log_lik()` on
  the same object already refuses with a designed message
  ("needs draws of the random effects ... Resample without
  laplace = TRUE", `draws_require_b()`); the two predictive methods
  should reach the same refusal. Test: assert the refusal on a laplace
  draws object for `posterior_predict()` and `posterior_epred()`, and
  assert no NaN is returned. Done at 0.66.0, lane wt-sampfix: laplace
  draws are read in their own layout, and a call that reads the
  integrated values is refused by name (`test-laplace-draws.R`).

### Open - medium

- **`frm_sample(laplace = TRUE)` on a model with NO random effect dies
  with a bare internal error.** Same script, same on both builds:
  `Error in -obj$env$random : invalid argument to unary operator`. With
  nothing to marginalize, `obj$env$random` is `NULL` and the negation
  is R's error rather than a designed one. Either accept the
  combination as a no-op (there is nothing to integrate out, so the
  Laplace route and the full route are the same model) or refuse it by
  name. Test: construct `frm_sample(bf(y ~ x), family = gaussian(),
  laplace = TRUE)` and assert whichever is chosen, not the unary-minus
  error. Done at 0.66.0, lane wt-sampfix: accepted, with a message, and
  the draws equal those of the call without `laplace`.
- **`pp_check()`'s four `loo_*` types fail on every model.**
  `pp_check.frmtmb_draws` resolves the bayesplot function and calls it
  with the response and the predictions, and nothing in it computes
  `lw` or a `psis_object`, which `bayesplot::ppc_loo_*` require and
  which `brms:::pp_check.brmsfit` builds from `log_lik()`.
  `dev/arcovsample-rev-07-ppcheck.R` (seed 31) over THREE models, a
  plain `y ~ x` gaussian fit, one with `(1 | g)` and one with
  `ar(t, g)`: **4 of 6 types tried, 3 of 3 models, 2 of 2 builds**,
  character for character identical.

  ```
  [dens_overlay]    OK ggplot2::ggplot
  [stat]            OK ggplot2::ggplot
  [loo_pit_overlay] rlang_error: One of 'lw' and 'psis_object' must be
                    specified.
  [loo_pit]         getvarError: argument "lw" is missing, with no
                    default
  [loo_intervals]   getvarError: argument "psis_object" is missing,
                    with no default
  [loo_ribbon]      getvarError: argument "psis_object" is missing,
                    with no default
  ```

  The lane measured the same thing independently
  (`dev/arcovsample-ppcheck.R`, one model, both builds). Test: one
  block per `loo_*` type asserting a plot object, which will fail until
  `pp_check()` passes the PSIS weights. Done at 0.66.0, lane
  wt-sampfix: built as brms builds them, bitwise equal to brms 2.23.0
  on brms's own draws (`test-ppcheck-loo.R`).

### Open - low

- **`nchains.frmtmb_draws()` and `draws_derived_matrix()` share the
  `NULL`-`stanfit` bug that `draws_chain_id()` had.**
  `x$stanfit@sim$chains %||% 1L` cannot guard a `NULL` `stanfit`,
  because `@` on `NULL` is an error and never reaches `%||%`.
  `draws_chain_id()` was fixed in the arcovsample round, and the
  reviewer saw the unfixed form fail behaviorally
  (`log_lik -> draws_chain_id -> %||%` on the base build,
  `dev/arcovsample-rev-log/12-nan-ref.txt`). The remaining two show up
  as `print(VarCorr(ds))` dying in
  `nchains.frmtmb_draws -> %||%` on a draws object built with
  `stanfit = NULL`, which is how a test supplies a chosen parameter
  vector without a sampler. A sweep of every `@` read behind a `%||%`
  is its own change. Test: `nchains()`, `ndraws()` and `VarCorr()` on a
  `stanfit = NULL` draws object. Done at 0.66.0, lane wt-sampfix: every
  read goes through `draws_nchains()` (`test-draws-no-stanfit.R`).

## Filed by lane splinecurve at the 0.65.0 consolidation, 2026-09-29

### Open - medium

- **An exact `gp()` read off its observed positions leaves out the
  kriging variance in three frmtmb.spline routes.** Under the default
  `allow_new_levels = FALSE`, `frm_curve(simultaneous = TRUE)`,
  `frm_curve_deriv()` and `frm_curve_feature()` build their draws or
  their delta method from `A V A'`, while the pointwise `frm_curve()`
  band adds `frm_lp_basis()`'s per-row `extra_var`. Pre-existing, and
  small where measured: `extra_var` 7.3e-07 to 8.5e-07 on `y ~ gp(x)`
  (`dev/splinecurve-findings.md`, "For the consolidating session to
  file"). Under `allow_new_levels = TRUE` the 0.9.0 refusal already
  covers it. The fix that makes every route right is a core seam that
  returns the between-row covariance block, not only its diagonal.
  Test: at grid points off the observed `x`, the simultaneous critical
  value from the full covariance against the one `frm_curve()` reports.

## Filed at the 0.66.0 release (2026-09-29)

What the six lanes of the parity round found and did not fix. Each
item names its lane and the findings section with the measurement;
`dev/round-20260929b.md` is the round's record. Items a lane chose not
to build, with the reason, stay in that lane's "decided not to do"
section and are listed here only where a user could hit them.

### Open - high

- **Laplace-draws probe: a NaN can pass.** On
  `bf(yn ~ log(c1) * x + exp(a)^k, c1 ~ 1, a ~ 1 + (1 | g), k ~ 1,
  nl = TRUE)` with `c1 = -1` and `k = 0` on draw 1 only, draw 1 is NaN
  for its own reason at both fills, so the watch takes the all-`NA`
  pattern as its reference and `posterior_epred()` returns 480 of 480
  cells non-finite where full draws give 80. Needs two exact values at
  one draw; the result is NaN, never a finite wrong number. Fix:
  compare each draw with itself at fill 0, at one more evaluation per
  draw. Lane sampfix, R2 (`dev/sampfix-rev-r2-01-watch.R`).
- **`predict(newdata = )` with NA responses under `ar()`/`arma()` with
  `cov = FALSE`** is refused by name. brms fills a missing response with
  its predicted draws and runs the recursion over them;
  `simulate(newdata = )` does that here. Ledger row
  `brmsfit-methods:747`. Lane defects, section 6.
- **bernoulli on a two-valued response that is not 0/1** (`-1`, `-2`)
  is refused; brms codes it by level order. The coding has to be stored
  with the fit and applied on every refit and every response read from
  newdata, or a subset holding one of the two values is coded wrong in
  silence. Ledger row `standata:75`. Lane defects, section 6.

### Open - medium

- **New levels in `conditional_effects(re_formula = NULL)` follow
  brms's `sample_new_levels = "gaussian"`**, not brms's default
  `"uncertainty"`, which lends an observed level's draw and gives
  narrower bands (brms 0.60, 0.51, 0.75 against frmtmb 0.90, 0.95,
  1.05 at 40 draws). On the fit and on draws. Lanes sampfix (M7) and
  postfit2 (section 8).
- **`conditional_effects()` refuses two displays brms answers:**
  crossed `(1 | g) + (1 | h) + (1 | g:h)` with `g` and `h` at observed
  levels never seen together, on `band = "boot"` and draws (the Wald
  band answers), and an `mm()` term with a `by` variable at a new
  level. Both refused by name. Lane postfit2, section 8.
- **`conditional_effects()` of a model with `mi(x, idx = )`** stops
  with "Could not match all indices" on its grid, whose `idx` and
  `index` values are reference values. Lane aterms2, section 6.
- **`conditional_effects()` and `emmeans()` on
  `y ~ x + offset(log(time))`** fail with "non-numeric argument to
  mathematical function" and "undefined columns selected". Pre-existing
  on 0.65.0. Lane aterms2, section 7.
- **An addition term given an expression**, `weights(wt * 2)` or
  `rate(time * 2)`, fails with "invalid model formula in ExtractVars";
  brms accepts it. Pre-existing for `weights()`. Lane aterms2,
  section 7.
- **`rate()` with `cens()` or `trunc()` on `poisson()`** runs brms's
  CDF at `mu * d`, but was never compared with brms. The compat row
  says "untested". Lane aterms2, section 6.
- **frmtmb.sample puts no default prior on a mixture's `sigma1` and
  `sigma2`**: `default_priors_for()` matches `dpar == "sigma"` only, so
  a mixture samples them flat where brms uses `student_t(3, 0, 2.6)`.
  Pre-existing. Lane formula2, section 4.
- **`thres(gr = )`, `cens()`/`trunc()` and mixtures** are refused on
  `hurdle_cumulative()`, `zero_inflated_beta_binomial()` and in a
  mixture respectively; brms fits each. `disc` on the four older
  ordinal families, and `threshold = "equidistant"` or `"sum_to_zero"`
  on any of them, are not built (ledger row `priors:14`). Lane fams2,
  "Decided not to do".
- **The covariance parameters in draws keep internal names**
  (`theta_1` where brms has `sd_g__Intercept`). `VarCorr()` and
  `hypothesis()` compute the brms quantities; `variables(ds)` and
  `variables(fit)` still differ there. Needs a per-structure map.
  Lane sampfix, "Found and not fixed" 2.
- **`plot()` of a conditional-effects object and of a hypothesis**
  refuses brms's `plot`, `rug`, `stype` and `ignore_prior`. frmtmb
  draws with base graphics where brms returns ggplot objects; whether
  to return plot objects is a decision for the user. Ledger rows
  `brmsfit-methods:154`, `:162`, `:164`, `:169`, `:391`, `:396`.
- **`parnames()` on a fit** (ledger row `brmsfit-methods:995`) moves
  the generic between the two packages' owner tables; and
  `nsamples(incl_warmup = TRUE)` (`:595`) needs stored warmup. Lane
  defects, section 6.

### Open - low

- **A harmless false alarm on laplace draws**: `(0 + x | g)` at
  `newdata` with `x = 0` is refused, because the probe's `NA` times 0
  is `NA`. Lane sampfix, M2 (`dev/sampfix-rev-01-probe.R`).
- **REML with `laplace = TRUE`** samples an objective rebuilt without
  the REML integral over `beta`, so the draws are the ML Laplace
  posterior. Lane sampfix, item 3.
- **`?frm_sample`'s `prior` says list names are "parameter names as in
  the draws"**; named-list priors resolve against the internal outer
  names (`tau_raw`, `beta`). Lane sampfix, item 4.
- **`update(fit, bf(y ~ x2))` with a complete `bf()` replaces the
  model**, where brms pools the old parameter formulas into it; and
  brms replaces an equation with a later formula on the same parameter
  where frmtmb refuses it. Both loud. Lane formula2, section 4.
- **`summary()`'s Formula line** does not print `sigma1 = sigma2`, or
  any parameter formula. Lane formula2, section 4.
- **A `y ~ . - w` fit refit on data without `w`** is refused ("The
  model uses `w`"), as `y ~ x1 + w - w` is. Lane formula2, 9.2.
- **A refusal from frmtmb.sample's draws accessors names the internal
  function**, "draws_accessor_args() has no argument `chains`", where
  the caller typed `as.array()`. Lane defects, section 7.
- **`summary()`'s `Data:` line is empty** for a fit whose data came by
  value through `do.call()`; brms records a deparse cut at 50
  characters. Lane defects, section 7.
- **`xbeta()` third derivatives inside the tie blend of
  `log_pbeta_ad()`** are off by up to 0.12 relative at shapes (1e6,
  3e6), against 3e-9 for plain `log(RTMB::pbeta())`. Values and first
  derivatives are unaffected; few rows land in the band. And `xbeta()`
  at `phi` 2e4 can stop 2.4e-5 short of the reference optimum without
  a warning, on a flat ridge. Lane fams2, m1 and n5.
- **`posterior_average()`'s refusals** quote the variable without
  brms's quotes, and `pp_average()` and `model_weights()` are not
  exported. Lane postfit2, section 8.
- **A multivariate model's `logLik()` `nobs` attribute** counts rows
  that no response's `subset()` uses (80 where 72 rows carry data).
  `BIC()` reads it. Lane aterms2, re-check.

### Upstream, listed in `dev/upstream-bugs.md` for the user to file

- brms 2.23.0 `hurdle_cumulative_logit_lpmf()` tests `y == nthres + 2`
  where the top category is `nthres + 1` (lane fams2).
- RTMB 2.0: `pbeta()`, `pbinom()` and `pnbinom()` have non-finite third
  derivatives at ordinary points, and `dbeta()` a `NaN` gradient past a
  shape sum of about 1e3 (lane fams2).

## Filed at the 0.67.0 release (2026-09-30)

What the three lanes of the round of 2026-09-30 (`ordinal`, `ceplot`,
`formrobust`) found and did not fix, and what the consolidation found.
Each item names its lane and the findings section with the
measurement; `dev/round-20260930.md` is the round's record. Deliberate
divergences from brms are not listed here; each lane's findings list
them, and the ledger records the ones the ported suite reaches.

### Closed by this round (were open at 0.66.0)

- `predict(newdata = )` with NA responses under `cov = FALSE` ARMA
  (formrobust, `brmsfit-methods:747`); bernoulli on two values that are
  not 0 and 1 (formrobust, `standata:75`).
- `conditional_effects()` on crossed `(1 | g) + (1 | h) + (1 | g:h)` at
  an unseen combination and on a new `mm(by = )` member (ceplot); with
  `mi(x, idx = )` it now stops with the reason, as brms stops.
- `conditional_effects()` and `emmeans()` on an `offset()` model, and
  an addition term given an expression (formrobust).
- `disc`, `threshold = "equidistant"` and `"sum_to_zero"` on the four
  older ordinal families (ordinal, `priors:14`).
- `plot()` of conditional effects and of a hypothesis takes brms's
  arguments and returns plot objects; the six ledger rows are a
  divergence, since the objects are not ggplot objects (ceplot).
- `parnames()` on a fit and `nsamples(incl_warmup = TRUE)` (ceplot).
- `update(fit, bf(y ~ x2))` pools the stored parameter formulas, as
  brms does (formrobust).

### Open - medium

- **`emmeans()` on a fit with a transformed predictor stops** with
  "undefined columns selected": `poly(z, 2)`, `log(abs(z) + 1)` and
  `scale(z)`, with or without an offset. The reference grid is built
  from the model frame, which holds the transformed columns; lane
  ceplot's raw-variable frame (`ce_base_frame()`) solved the same
  problem for `conditional_effects()`. Pre-existing on 0.66.0 and on
  both lane builds. Found at consolidation (`dev/rel067-emmoffset.R`,
  `dev/rel067-log/emmoffset.txt`).
- **`update(fit2, formula. = bf(count ~ a + b, nl = TRUE))` stops in
  the refit** from the default starting values ("NA/NaN gradient
  evaluation"), on the lane build and the release build alike. The
  update itself is brms's now. Ledger row `brmsfit-methods:955`, a
  defect; the lane had recorded it as a class-name divergence. Found at
  consolidation (`dev/rel067-update955.R`,
  `dev/rel067-log/update955.txt`).
- **Ordinal mixtures stay refused.** brms fits them with `disc1`,
  `disc2` and thresholds fixed across components (its equidistant
  mixture does not compile, brms-17). Lane ordinal, "Not done".
- **`thres(gr = )` and `cs()` on `hurdle_cumulative()`** stay refused;
  brms fits both. Lane ordinal, "Not done".
- **`posterior_linpred(incl_thres = TRUE)`** in frmtmb.sample stays
  refused; brms returns `disc * (thres - mu)` per threshold for every
  ordinal family. Lane ordinal, "Not done".
- **`conditional_effects()` refuses `trunc()` and `se()` terms whose
  variables are not pinned in `conditions`** (`ce_aterms()`); brms
  holds such a variable at its mean, and a `min(y) - 1` bound at
  `mean(y) - 1`. Pre-existing. Lane formrobust, section 9.

### Open - low

- **A fixed `disc` in `variables()`**: brms lists `disc` (a transformed
  parameter, 1) on every ordinal fit; frmtmb lists no fixed dpar for any
  family. Lane ordinal, "Not done". Closed at 0.68.0 (lane fixes).
- **Class `"Intercept"` with `coef`** does not address one ordinal
  threshold; brms lists per-threshold rows under flexible and
  sum-to-zero thresholds. Pre-existing. Lane ordinal, "Not done".
  Closed at 0.68.0 (lane fixes).
- **`confint()` names the ordinal internal parameters `tau_raw_k`**;
  under equidistant thresholds `tau_raw_2` is `log(delta)`.
  Pre-existing. Lane ordinal, "Not done". Closed at 0.68.0 (lane
  fixes).
- **The compatibility table (`R/compat.R`)** has no rows for the
  threshold structures or `disc`. frmtmb.eam's tests read that table.
  Lane ordinal, "Not done". Closed at 0.68.0 (lane fixes).
- **frmtmb.sample's default prior on the `disc` intercept under
  `link_disc = "identity"`**: brms uses `lognormal(0, 1)`, which
  `set_prior()` does not carry, and the default-prior message names the
  gap. Lane ordinal (frmtmb.sample NEWS).
- **`fitted(ndraws = )` on a fit** is refused, so ledger row
  `brmsfit-methods:314` does not hold although `fitted()` takes
  `sample_new_levels = "old_levels"` now. Lane ceplot, section 7.
- **`frm_sample()` takes no `drop_unused_levels`**; fit with
  `frm(drop_unused_levels = FALSE)` and sample the fit. Lane
  formrobust, section 9.
- **`refit(newresp = )` of a recoded bernoulli fit** takes the 0/1
  codes that `simulate()` returns and refuses the original values by
  name. Lane formrobust, section 9.

### Upstream, listed in `dev/upstream-bugs.md` for the user to file

- brms 2.23.0: brms-15 (`cratio("cloglog")` with `disc`, NaN gradient
  on some data), brms-16 (`sum_to_zero` with the logit link does not
  compile), brms-17 (equidistant `delta` declared twice in a
  multivariate or mixture model), all lane ordinal; brms-18 (a new
  `mm(by = )` member cannot be drawn), lane ceplot; brms-19 (an
  all-`TRUE` bernoulli response coded as all failures), lane
  formrobust. Lane ordinal met brms-2 (softit) and brms-8
  (`probit_approx`) again.

## Filed at the 0.68.0 release (2026-10-06)

What the four lanes of the round of 2026-10-05 (`vigport`, `fixes`,
`gpby`, `ordmix`) found and did not fix, and what the consolidation
found. Each item names its lane and the record with the measurement;
`dev/round-20261005.md` is the round's record. Lane nanse (silent NaN
standard errors, the starts of `frm_allfit()`) was still in work at the
release and owns vigport's defects 8 and 9 and fixes' consolidation
items 1 and 3; they are listed here so that nothing is lost if that
lane's merge is delayed. Lane ordmix's own section below ("Filed by
lane ordmix after punch round 1") holds the probit underflow and the
multi-start option. Deliberate divergences are in the lanes' findings
and the ledger, and upstream defects in `dev/upstream-bugs.md`.

### Closed by this round (were open at 0.67.0)

- Class `"Intercept"` with `coef = "1"`, ... addresses one ordinal
  threshold, and `default_prior()` lists the per-threshold rows, also
  under `dpar = "mu<k>"` in an ordinal mixture (fixes, ordmix and the
  merge).
- `confint()` and `vcov(full = TRUE)` name the ordinal internal
  parameters by what they are (fixes; ordinal mixtures at the merge).
- The compatibility table has rows for `disc` and the threshold
  structures (fixes; their mixture cells at the merge).
- A fixed `disc` in `variables()` (fixes).
- `emmeans()` on a transformed predictor (fixes).
- `conditional_effects(method = "predict")`'s `estimate__` is the
  median of its draws, as brms's is (the user's decision of
  2026-10-06, at the merge; lane fixes' m8).
- `cs()` on `cumulative()` fits, as in brms (the user's decision of
  2026-10-06, at the merge).
- vigport's defect 1 (the `s()` null-space scale; fixes) and defect 6
  (the migration vignette on `hurdle_cumulative()`; at the merge).

### Open - medium

- **vigport defect 2: `stancode()`, `standata()` and `pp_mixture()` on
  a `frmtmb_fit` have no refusal.** Without brms loaded the error is
  "no applicable method"; with the brms namespace loaded they dispatch
  to brms's default and say "Data must be specified using the 'data'
  argument". Repro `dev/vigport-repros.R` R1. Lane vigport.
- **vigport defect 4: `update()` keeps a prior the new formula cannot
  use, and refuses.** brms's `update()` of `brms_overview`'s `fit2`
  recompiles and drops the stale `lkj(2)` prior without a message
  (the reviewer's `dev/vigport-rev-brms-update.R`). Repro
  `dev/vigport-repros.R` R3. Lane vigport.
- **vigport defect 7: `brms_multilevel`'s `fit_loss2` does not
  converge reliably.** nlminb reports false convergence; eight refits
  over four optimizers and two starts span -358.52 to -359.26 (0.74),
  the best `nloptr_lbfgs` at -358.5206 with code 0, 0.025 above
  nlminb's code-1 value (`dev/vigport-pr1-conv.R`). brms samples it
  with Rhat 1.01. Whether a better optimum exists is not settled.
  Lane vigport.
- **vigport defects 8 and 9, owned by lane nanse**: `ls ~ mo(income) *
  age` has all-NaN standard errors on 72 of 200 of `brms_monotonic`'s
  data sets, 55 of them silent (`dev/vigport-pr1-mo.R`); and
  `frm_allfit()` refits a nonlinear model with `start = NULL`
  (`R/allfit.R:93`; `dev/vigport-rev2-allfit.R`).
- **fixes, consolidation item 1: a linear ridge next to a coefficient
  near an undefined region gives NaN standard errors and no warning.**
  `a + b + log(c0)` with `a, b ~ 1 + x`, `c0` estimated at 5.9e-5: the
  differenced Hessian of the flat check is not finite, so the check is
  silent, and the gradient is finite, so no non-finite-gradient
  warning either (the fixes final review, item 7). Owned by lane
  nanse with the item below.
- **`fixef()` and `summary()` report NaN standard errors without a
  warning** when the outer Hessian is singular, as at a smoothing SD on
  its boundary; only `vcov()` warns. 3 of 140 gamSim smooth fits on
  0.67.0, 5 of 140 on lane fixes' build (`dev/fixes-sx-conv3.R`,
  `-conv5.R`). Lane fixes, punch round 1.
  Owned by lane nanse.
- **`fitted()` on an ordinal fit whose `disc` predictor has no fixed
  column stops** with "requires numeric/complex matrix/vector
  arguments" in `lp_eta_design()` at `X %*% est[[lp$par]][lp$idx]`:
  `bf(y ~ z, disc ~ 0 + gp(x, k = 6))` and `disc ~ 0 + (1 | g)` on
  `cumulative()`, in sample and on newdata. `est$betad` is `NULL` when
  no `disc` column is fixed. `disc ~ 0 + s(x)` (which keeps a fixed
  linear column), `sigma ~ 0 + gp(x)` and `sigma ~ 0 + (1 | g)` work,
  so the gap is the `betad` slot, not gp(). Pre-existing on 0.67.0 and
  on the gpby lane build. Found by the gpby review
  (`dev/reviews/2026-10-05-gpby.md`); scope in `dev/gpby-p1-disc2.R`,
  `dev/gpby-p1-disc2.txt`.
- **`frm_sample(fit)` on an exact `y ~ gp(x)` fit does not move**:
  stepsize NaN, acceptance 0, all 300 post-warmup transitions
  divergent, every draw at one point, on rellib-r5 and on the gpby lane
  build alike (60 points, data seed 5, `chains = 1, iter = 600,
  seed = 4`). It is not a flat prior: since frmtmb.sample 0.43.0 the
  fit route carries brms's defaults, and `prior_summary()` is the same
  four rows on both routes. The formula route on the same data and
  seed samples (acceptance 0.97 lane, 0.92 base). The fit route's
  start or its first gradient is the suspect. Repro
  `dev/gpby-p1-m1b.R lane|base`, logs `dev/gpby-p1-m1b-*.txt`,
  `dev/gpby-p1-m1.txt`.
- **`frm_sample()` repeats a warning the frame build gives**: on
  `y ~ z + cs(x)` with `cumulative()` it gives brms's experimental
  warning 6 times where `frm()` and brms's `brm()` give it once
  (`dev/rel068-cs-probe.R`, `dev/rel068-log/cs-probe.txt`). Every
  warning raised while the sampling route builds its frames is
  repeated the same way. Found at the 0.68.0 consolidation.
- **`cs()` on a cumulative component of an ordinal mixture often stops
  in the optimizer.** On 20 seeds of a two-class mixture with no `cs()`
  effect (n = 300; the release review's `dev/relrev-csfail.R`):
  `mixture(cumulative, sratio)` with `cs(z)` gave a NaN-gradient error
  on 7, no convergence on 11, 1 clean fit and 1 warned;
  `mixture(cumulative, cumulative)` with `cs(z)` gave 4 errors, 16 not
  converged, 0 clean. Controls: `mixture(sratio, sratio)` with `cs(z)`
  0 errors (6 not converged, 3 clean, 11 warned), and
  `mixture(cumulative, sratio)` without `cs()` 0 errors (5 not
  converged, 14 clean, 1 warned). Plain `cumulative()` and the probit
  with `cs()` fit 20 of 20. The NaN gradient comes from rows crossing
  during the search; brms's warning ("may have convergence issues")
  covers it. The release review, m5.
- **`mo()` point estimates are not always the maximum.** On 52 of 200
  seeds of `brms_monotonic`'s `ls ~ mo(income) * age` the log-likelihood
  is more than 0.01 below the exact profile maximum, up to 1.95 below;
  on 20 the fit stops with code 0 on a softmax plateau that a step
  toward a simplex vertex improves by up to 0.70 (prototype
  `dev/nanse-mo-escape.R` in lane nanse's worktree), and on 38 it is at
  another local maximum after that escape. An escape, or a simplex
  parameterization with a reachable boundary, changes every saturated
  `mo()` fit and the sampler's parameterization: a lane of its own.
  Lane nanse, `dev/nanse-findings.md`, "Defect found, not fixed"; it
  corrects the vigport review's "valid optimum" (corrected on that page).

### Open - low (found by the review of the 0.68.0 consolidation)

- **Unseeded data in the ported brms suite** moves the recorded
  messages of `standata:310`, `:315` and `:316` at every regeneration,
  so `dev/rel068-msgdiff.R` lists them each time. Mask numbers in its
  `norm()` (for example `gsub("-?[0-9.]+(e-?[0-9]+)?", "#", s)`) or
  seed each block in `brms_port()`. The verdicts are not affected. The
  release review, m9.
- **`(cs(1) | g)` with brms attached after frmtmb dies with an internal
  error**, "variable lengths differ (found for '.frm_cs(1)')", on every
  family (`dev/relrev-cs1g.R`); without brms attached it gets frmtmb's
  refusal. Pre-existing, the same on rellib-r5. And
  `conditional_effects(method = "predict")` on `categorical()` says it
  "has no meaning on an ordinal family". The release review, m10.
- **A sum-to-zero component of an `order = "none"` mixture lists no
  class `"Intercept"` rows**, where brms lists 4, as for one sum-to-zero
  family (the threshold vector brms declares before centering is not a
  frmtmb parameter). The NEWS mixture bullet says so. The release
  review, m11.
- **`test-perf.R` "fit time grows with n and stays within a linear
  envelope" is a wall-clock bound** (`t_large < 100 * max(t_small,
  0.01)`), which the rules forbid. Under load it failed once in the
  release's first gated run after lane nanse's merge (1.29 s against
  1.00; 22 test processes beside two study scripts), and passed 3 of 3
  alone (`dev/release/perf1..3.log`). nanse's fit-time check adds 8 to
  24 percent to a small fit, which narrows the margin. Replace it with
  a count (objective evaluations, AD nodes) or a ratio to a control
  timed in the same process. Found at the 0.68.0 nanse merge.

### Open - low

- **vigport defect 3: `plot()` of a fit suggests `x` for brms's `N`**
  ("plot() has no argument `N`. Did you mean `x`?"). Repro
  `dev/vigport-repros.R` R2. Lane vigport.
- **vigport defect 5: `fixef()` has no `frm_multiple()` method**, and
  the pooled table names the intercept `(Intercept)` and sigma
  `sigma_(Intercept)`, where brms's `fixef()` of a `brm_multiple()`
  fit says `Intercept`. Repro `dev/vigport-repros.R` R4. Lane vigport.
- **ordmix: the degenerate-component check misses a component that
  stopped using an interior category** (8 misses of the review's
  sets, accepted at its re-check), and it costs 11 to 14 percent of
  the fit at n = 50000 (1.45 s of 12.9 s; 2.55 s of 17.8 s with
  `cs()`; `dev/ordmix-rev3-misc.R`). The ordmix final review, c2, c3.
- **ordmix: an `order = "mu"` mixture's `insight` parameter scale**
  stays internal, since its block list holds `tau_raw` twice; and a
  sum-to-zero component has no class `Intercept` row, which brms
  lists. `dev/ordmix-findings.md`, "Not done, and why".
- **fixes: F11 seed 8** (`bf(y ~ s(x1), sigma ~ s(x2))`) stops 0.142
  below rellib-r5's optimum, converged, with finite standard errors:
  sigma's smoothing SD sits at a local optimum (theta -7.99, -0.80)
  where 0.67.0 runs it to the boundary (-28.88, -9.87). `restarts = 3`
  does not move it, and dropping the smooth units brings back the
  NaN-SE runaway elsewhere (`dev/fixes-p2-f11.R`). A known local
  optimum of the new optimizer path. `dev/fixes-findings.md`, m11.
- **gpby: brms's other `gp()` kernels** (`cov = "matern32"`,
  `"matern52"`, `"exponential"`) are refused by name. The exact form
  needs the kernel in `gp_corr()` and `gp_cross_cov()`; the
  Hilbert-space form needs brms's spectral densities.
  `dev/gpby-findings.md` section 7.
- **gpby: `zgp` is not a draw name.** frmtmb samples the field and
  brms its standardized `zgp`; the map is the change of variables
  `L^-1 b`, which a column cannot be renamed into. Section 7.
- **gpby: the draws' column order** puts each sub-GP's `sdgp` beside
  its `lscale`; brms puts every `sdgp` first. Names match, order does
  not. Section 7.
- **gpby: core `predict()` holds a `gp()` curve at its mode** and draws
  neither its coefficients nor the kriging residual, so its interval
  past the positions is narrower than brms's. Section 7.
- **gpby: a nonlinear body refuses a contributing exact `gp()` at an
  unseen position** (`lp_basis_nl()`). Section 7.
- **gpby r4: `test-gp-by-draws.R`'s formula-route fixture chain is
  poorly mixed** (on the m1b construction, `iter = 600`, one chain: 215
  transitions over the maximum tree depth, R-hat 1.56,
  `dev/gpby-p1-m1b-lane.txt`). Its assertions read the kriging law at
  one fixed draw, so they are not at risk; a later test must not read
  posterior summaries off it. The gpby re-check, r4.

### Upstream, listed in `dev/upstream-bugs.md` for the user to file

- brms 2.23.0: brms-20 (an ordinal mixture's NaN gradient at a
  finite density) and brms-22 (crossing `cs()` thresholds), lane
  ordmix; brms-21 (`incl_thres` on `hurdle_cumulative()`), lane fixes;
  brms-23 (`threshold =` passed to `brm()` and dropped), lane vigport.
  Lane ordmix met brms-1 (on the `cs()` path) and brms-17 (in
  mixtures) again.

## Filed by lane ordmix after punch round 1, 2026-10-05

### Open - medium

- **The probit's log-odds form underflows past `|eta| = 38.2`**
  (`R/links.R`, `logit_eta` of `probit`): `log(pnorm(eta)) -
  log(pnorm(-eta))` is `-Inf` or `Inf` there, and the ordinal densities
  then give `NaN`. A plain cumulative probit fit at slope 20 has a
  `NaN` objective on both 0.67.0 and the lane build
  (`dev/ordmix-rev-probit-sat.R`); in an ordinal mixture a component
  can reach that region on its own, and the review's
  `mixture(cumulative("probit"), sratio("cloglog"), acat())` ended
  there with an infinite gradient (now a warning in
  `check_convergence()`). The fix, `pnorm(eta, log.p = TRUE) -
  pnorm(-eta, log.p = TRUE)`, holds to `|eta| = 60` and beyond, but it
  is not the current form bit for bit where that is finite: 483 of
  2001 grid values on `[-38, 38]` differ (at most 4.3e-12 relative),
  and 1903 of 1977 derivatives (at most 1.4e-13), through RTMB's tape
  (`dev/ordmix-p1-probit.R`, log `dev/ordmix-p1-log-probit.txt`). So
  it moves every plain probit fit, and the round's rule kept it out. A
  fix that keeps the plain fits needs a branch at `|eta| = 38` on the
  tape, or an accepted change to the probit fits.

### Open - low

- **A multi-start option for ordinal mixtures** (`frmtmb_control(starts
  = n)`, best of `n` jittered starts). Not built: in the review's 160
  two-component fits, 22 ended at a degenerate boundary, and in 16 of
  the 22 that point had the higher log-likelihood, so "best of n"
  picks the degenerate point more often than it escapes it. A useful
  option has to rank starts by something other than the likelihood
  (the degenerate check of `mixture_ord_degeneracy()` is a candidate)
  and report every start, which is more than a small change. Priors on
  the components (as brms has them) are the remedy the warning names.

## Filed at 0.68.1

Found by lane cifix (`dev/cifix-findings.md`) and its review
(`dev/reviews/2026-10-06-cifix.md`); none blocks the CI patch. "The
emulator" is `dev/cifix-openblas.sh`, R 4.6.1 with OpenBLAS 0.3.26.

### Open - high

- **The SE check's tier 1 accepts an inverse built on a Hessian row of
  pure noise** (review item 4). Instrumented, all 353 test files: 117
  fits in 58 files (core 102, coupling 4, learn 2, sample 6, spline 3)
  keep `solve(H)` with a row whose maximum is 2.0e-12 to 3.1e-05; 125
  of the 128 parameters are variance components at a boundary, whose
  reported SE (optimizer scale) is 884 to 750,000, median 4,960, with
  no warning. 5 of 104 change verdict with the BLAS, which is why
  `test-aliased-grouping.R:138` and `test-open-issues.R:52` now wrap
  the warning. Repro: `(1 | Subject/a)` on sleepstudy with
  `a = factor(Days %% 3)` (`dev/cifix-diagnest.R`: Hessian row
  +7.1e-12 with the reference BLAS, -7.1e-12 with OpenBLAS). Remedy the
  review recommends: apply tier 3's empty-row rule before tiers 1 and
  2, in a planned round. Budget 48 test files (most need
  `allow_warnings()`; `test-diagnostics-ux.R:99` and
  `test-car-spde.R:324` need a decision), and consider a shorter
  message for a variance component at its boundary (lme4: "boundary
  (singular) fit").
- **A nonlinear ridge on a fit with random effects reports a huge SE
  with no SE warning** (review item 4). `y ~ a + b, a ~ 0 + f,
  b ~ 1 + (1 | g), nl = TRUE` (`dev/cifixrev-refits.R`, R2, seed 1,
  k = 10 levels of f, 6 of g): tier 1 accepts the finite-difference
  Hessian and a prediction of `a` has SE 822,571; only the lane-fixes
  identification warning fires. Without `(1 | g)` the exact Hessian
  gives NaN. The empty-row rule does not catch it; it needs tier 2's
  eigenvalue test on the finite-difference Hessian.
- **Separation is not named at the default budget** (review m5).
  `set.seed(514); d <- data.frame(x = rnorm(240) * 1e-4, z =
  rnorm(240)); d$yb <- as.integer(d$z > 0); frm(yb ~ z + x, family =
  bernoulli(), data = d)`: with OpenBLAS nlminb stops at code 9 after
  1997 evaluations with "function evaluation limit reached", and the
  word separation never reaches the user; the reference BLAS stops at
  code 0 after 1939 and names it. `glm()` reports fitted probabilities
  of 0 or 1 after 25 iterations. The convergence check should name
  separation when a binomial-type mean coefficient has run past |10|
  with fitted probabilities at 0 or 1, whatever code nlminb returns.
  `test-se-check.R` now sets `eval.max = iter.max = 4000` for that
  test (`dev/cifix-diagsep.R`).

### Open - medium

- **`ranef(condVar = TRUE)` still reads sdreport()'s own random-effect
  covariance** (review m6), which is the indefinite inverse when an SE
  is lost. On the `gr(g, by = f)` fixture (`dev/cifixrev-r1.R`, seed
  11) condsd / sqrt(diag of the repaired V_bb) runs from 0 to 99, so
  `ranef()` and `frm_linpred()` describe the same `b` with different
  covariances. Route it through `get_joint_cov()` when the repair ran.
- **`test-cumulative-cs.R:132` and `test-ordinal-mixture.R:751` are
  fragile to rounding** (review item 5). With the emulator, on 0.68.0
  and 0.68.1 alike, the first stops with "NA/NaN gradient evaluation"
  (`mixture(cumulative, sratio)`) and the second lets "singular
  convergence (7)" escape; both pass on the Ubuntu runner. The
  `ubuntu-latest` label moves to Ubuntu 26 from 2026-10-19, and a new
  OpenBLAS there may surface them. Repro: run each file with the
  emulator, `OPENBLAS_NUM_THREADS=4` (`dev/cifix-par.sh`).
- **frmtmb.sample lets two warnings escape on the Ubuntu runner**
  (review item 5; check-frmtmb.sample on d24f7b86: FAIL 0, WARN 2,
  SKIP 18, PASS 2274). `test-brms-shapes-draws.R:106`: "Method
  'posterior_samples' is deprecated", from the runner's brms, more
  than once, past the test's `expect_warning()`; use
  `allow_warnings()`. `test-compat-preflight.R:95`: "no DISPLAY
  variable so Tk is not available"; find which call loads tcltk at
  preflight.

### Open - low

- **Core's gp() position key is 15 significant digits, not exact**
  (review m2). `pos_rowkey()` pastes the coordinates, so 1/3 and
  1/3 * (1 + 2^-52) share a key and are kriged as one position; the
  spline's `sp_row_key()` compares doubles exactly, so the two layers
  disagree on what one row is. Harmless numerically (the rows differ
  by 9.1e-11 at most there, the per-row rounding). Remedy: key on
  exact doubles, for example `sprintf("%a")` per coordinate
  (`dev/cifixrev-dedup.R`).
- **`gp_krig_cov()` forms `outer(w, w)` whole under `gp(x, by =
  <numeric>)`** (review m9): two n x n beside the result, which undoes
  part of the blockwise memory discipline of the gpby review's m5.
  Peak memory on a 2000-row grid was not distinguishable from load
  noise (`dev/cifixrev-krigmem.R`). Remedy: scale each column block
  after it is symmetrized.

## Filed by lane optima after punch round 1, 2026-10-07

`dev/optima-findings.md` and the review `dev/reviews/2026-10-07-optima.md`.

### Open - medium

- **cox()'s baseline simplex (`sbhaz_raw`) is a softmax and meets the
  plateau `mo()` met**: weights of 1.1e-8, 8.3e-9 and 5.9e-8 on seeds
  2, 5 and 6 of the review's `dev/optima-rev-cox1.R`. The argument of
  `mo_simplex()` applies, but the change is not local: the baseline is
  read by `cox_baseline()`, the predictions and the sampler, which
  would need the softmax-and-Dirichlet treatment frmtmb.sample now gives
  a `mo()` simplex (brms puts `dirichlet(1)` on `sbhaz`). Review m9.
- **`escape_stationary()`'s sibling paths**: the objective's state is
  now settled after an escape and after `mo_search()`
  (`obj_settle()`); any other post-fit refit that keeps a run other
  than the last must call it too.

### Open - low

- **`mo()` seed 194** of `dev/optima-mo-study.R` stays 0.0489 below the
  exact maximum, at the other sign's cone, which 0.68.1 found: the
  search's vertex start (chosen from the gradient at `b = 0` with the
  other parameters held) did not lead into that cone's optimum there.
  A start at each vertex would reach it at D times the cost. Review m8.
- **A `simo` prior row is refused by `frm()` and `frm_sample()`**: brms
  accepts `prior(dirichlet(c(2, 1, 1)), class = simo, coef =
  moincome1)`. `set_prior()` has no Dirichlet density; the sampler's
  `mo_simplex_nlp()` is where a concentration vector would go.
- **`as_tmbstan(fit)` on every `mo()` fit**, ML or MAP, samples the
  fit's own tape, whose simplex coordinates have no density: on a
  weakly identified simplex the draws collapse to the barycenter (sd
  1e-17, coordinates to 7e17; the review's re-check). By contract the
  route adds nothing; since punch round 2 it says so in a message and
  `?as_tmbstan` names `frm_sample()`. (`check_laplace()` retapes like
  `frm_sample()` since punch round 2, a MAP fit with its own prior.)
- **`mo(x) * f` with a factor and `mo()` in group-level terms**, both of
  which brms fits, are refused (pre-existing). Review m11.

## Reference

Full agent report with per-item repro sketches and issue links:
sourced from glmmTMB/lme4/brms GitHub test suites and NEWS, 2026-08-31.
Key files to mirror: glmmTMB test-predict.R, test-formulas.R,
test-varstruc.R, test-offset.R, test-weight.R, test-NAhandling.R (lme4),
tests.brmsformula.R (brms).
