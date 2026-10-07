## Filed at the 0.69.0 release (2026-10-07)

What the four lanes of the round of 2026-10-06 (`ciharden`, `surface`,
`setier`, `optima`) found and did not fix, and what their reviews and
the consolidation listed. Each item names its lane and the record with
the measurement: `dev/<lane>-findings.md` and
`dev/reviews/2026-10-07-<lane>.md`. `dev/round-20261007.md` is the
round's record. This section replaces the four sections the lanes
added ("Filed by lane ciharden", "Filed by lane surface", "Filed by
lane setier", "Filed by lane optima after punch round 1"), item for
item. Each item this round closed is marked in place above with
"CLOSED at 0.69.0".

### Closed by this round (were open at 0.68.1)

- The SE check's tier 1 accepted an inverse built on a Hessian row of
  noise; a nonlinear ridge with random effects kept a huge SE;
  separation was not named at the default budget; `ranef(condVar =
  TRUE)` read sdreport()'s indefinite covariance (setier).
- `test-cumulative-cs.R:132` and `test-ordinal-mixture.R:751` fragile
  to rounding; frmtmb.sample's two escaped warnings on Ubuntu;
  `gp_krig_cov()`'s `outer(w, w)`; the ported suite's unseeded data;
  `test-perf.R`'s wall-clock bound (ciharden). The `gp()` position key
  stays brms's 15-digit key, the user's decision of this round, since
  brms 2.23.0 groups new-data rows the same way (ciharden review, B2).
- vigport defects 2, 3, 4 and 5 (`stancode()`, `standata()` and
  `pp_mixture()` on a fit; `plot()` and `N`; `update()` and a stale
  prior; `fixef()` of `frm_multiple()`); `fitted()` with no fixed
  `disc` column; `frm_sample(fit)` on an exact `gp()`; `(cs(1) | g)`
  with brms attached and the categorical message (surface). The
  repeated frame-build warning of `frm_sample()` was a misreading: the
  six warnings were six different ones, and the frame-build warning is
  given once (surface, item 5, with a guard test). The sum-to-zero
  component's missing `Intercept` rows are a documented divergence:
  brms's rows name a vector its Stan program declares before
  centering (surface, item 7).
- `mo()` point estimates below the maximum; `cs()` on a cumulative
  component of an ordinal mixture stopping on a NaN gradient; the
  probit log-odds underflow past `|eta| = 38.2` (and the same in
  cloglog and softit) (optima). frmtmb.sample's `mo()` draws are named
  `simo_<term>1[k]` as brms names them, which closes surface's filing
  of the `zeta1_1` names in the same round.
- vigport defects 8 and 9, fixes' consolidation item 1 and the NaN
  `fixef()` standard errors were closed at 0.68.0 by lane nanse;
  setier measured them closed on rellib-r6 (`dev/setier-sx.R`, 140
  fits, 0 non-finite; `dev/nanse-item3.R`).

### Open - high

- **A gaussian random-intercept fit stops short of the maximum near
  sd = 0, at code 0.** `y ~ x + (1 | g)` (20 groups of 5, sigma 1) with
  the response multiplied by 1e-3, 1 and 1e3, 40 seeds each: the fit
  ends more than 1e-4 below lme4's log-likelihood on 4, 0 and 20 of
  the 40, by up to 0.825, 23 of the 24 at optimizer code 0 (seed 36 at
  1e-3: "false convergence (8)"). The sd sits near zero (2.5e-22 times
  sigma on seed 19 at 1e-3, where lme4 puts it at 0.31), so the
  log-likelihood is flat toward zero, and `escape_stationary()` does
  not reach it. rellib-r6 gives the same fits, and so does lane
  optima's build. Since lane setier the SE warning says "the fit
  stopped short of the maximum" (or the convergence warning stays);
  the answer itself is still wrong. Under `quadrature = TRUE` the
  short fits keep the "flat" SE warning or nothing, as on 0.68.1.
  Repro: `dev/setier-rev2-scale2.R <lib>`, `dev/setier-rev2-trap.R
  <lib>`, `dev/setier-rev2-short.R <lib>`, `dev/setier-rev3-trapq.R`
  (setier review, RB2 and RQ1). Remedy to try: start a group sd scaled
  to the response, or an escape that tries the sd at a fraction of the
  residual sd as `se_sd_gain_up()` does. Lane setier.
- **A prediction along a lost direction is flagged when it is
  estimable** (the consolidation review, B1 and m4). On the ported
  fixture 1 (`data-helpers:7`), where the fit converges under OpenBLAS
  LAPACK, `fitted(fit1, newdata)` warns that all 5 predictions move
  along a direction the fit does not determine. They do not: the
  `s(Age)` fixed column is exactly `-0.046394458 + 0.190041455 * Age`
  (residual 1.3e-16), so the alias is (Intercept, Age, `s(Age).fx1`),
  `max |X d|` along it is 1.1e-16, and every prediction is constant
  along it. The lost basis the prediction test reads
  (`get_joint_cov()$null`) is restricted to the named parameters (lane
  nanse's m1 follow-up): it keeps Age 0.293 and `s(Age).fx1` -1.545
  and zeroes the intercept's share (0.0455 of unit norm) because the
  intercept keeps its SE. That truncated direction is not flat
  (objective +1.0e-6, +9.9e-5, +9.8e-3 at steps 1e-3, 1e-2, 1e-1), so
  `jc_nonest()`'s cosine against `se_pred_tol = 1e-6` flags every row.
  Estimability must be judged against the full null direction, before
  the loadings of kept parameters are zeroed, or `se_pred_tol` scaled
  to the basis's own precision. The SE losses themselves (Age and
  `s(Age).fx1` flat, the correlation at its boundary) are right.
  Repro: `dev/relrev069-dh7diag.R` under
  `dev/ciharden-openblas.sh 0.3.26 lapack`; `dev/rel069-dh7.R`. The
  test must not be loosened to allow the warning. Consolidation.
- **A restart that raises discards a finite first run** (the
  consolidation review, item 3a and m2). In `optimize_obj()`, when the
  restart from a first run's optimum raises "NA/NaN gradient
  evaluation", `optimizer_from_best()` restarts from `last.par.best`,
  where the objective can be NaN, and the error propagates although
  the first run's result was finite (fuzz seed 20379118: first run at
  "false convergence (8)", objective -49.18, finite gradient). The
  minimal fix the review measured keeps the run before a restart that
  raises:

      opt2 <- tryCatch(run(opt$par), error = function(e) NULL)
      if (is.null(opt2)) { obj$fn(opt$par); break }

  (`obj$fn()` settles the objective's state at the kept point, as
  `obj_settle()` does). Over 300 seeds per cell of the spec's family
  (`dev/relrev069-fuzzsweep.R`) it removes 9 to 12 errors and adds
  none; on seed 20379118 under OpenBLAS 0.3.26 it returns code 1,
  objective -54.0241049, as 0.68.1 returns a fit. 0.69.0's error rate
  is not worse than 0.68.1's overall (spec family: 37, 42, 42 of 300
  against 31, 52, 52 under the reference BLAS, 0.3.26 and 0.3.32). It
  needs a test seen to fail and a run with OpenBLAS LAPACK.
  Consolidation.
- **A boundary message on an sd the likelihood does not read** (the
  consolidation review, item 6 and m3). On `test-diagnostics-ux.R` fit
  58 ("a theta in the flat set ...") the message says
  "sd_id__pk_Intercept is at zero" while `diagnose()` gives that sd as
  0.2199, and the objective is exactly unchanged moving its log by
  -10, -3, -log 7, +log 7 and +1: the sd is not at an edge, it does
  not enter the likelihood. 0.68.1 listed it as flat with the other
  three in one SE warning (`dev/relrev069-flatsd.R` on both builds).
  `se_at_edge()` and `se_boundary_names()` should require the sd to be
  near zero on the scale `se_sd_gain_up()` uses (below its smallest
  probe, 0.01 of the residual sd, or of 1 on the link scale) before a
  boundary verdict, and report it as "flat" otherwise. Lane setier.

### Open - medium

- **One fit in `test-id-kron.R` ends at one of two optima from run to
  run.** The merged-dpar fit ("merging across dpars of one response
  takes the same path") stops at code 7 with the log-likelihood flat in
  the sd-sigma correlation. Ten runs of the file at once: 8 of 10 end
  at objective 128.07435688 (largest gradient 1.6e-4), 2 of 10 at
  128.074356883 (gradient 5.6e-4); on rellib-r6, 6 and 2 of 8. Run
  alone it is always the first. Since lane setier the first is a
  boundary stop (the boundary message) and the second keeps the code 7
  warning, so the file's condition count changes between runs; the
  tests pass either way. Likely `dev/rtmb-pitfalls.md` item 21 (P and
  E cores). Repro: `dev/setier-idk-test.R`, ten at once;
  `dev/setier-idk.R <lib>`. Lane setier.
- **cox()'s baseline simplex (`sbhaz_raw`) is a softmax and meets the
  plateau `mo()` met**: weights of 1.1e-8, 8.3e-9 and 5.9e-8 on seeds
  2, 5 and 6 of `dev/optima-rev-cox1.R`. The argument of
  `mo_simplex()` applies, but the change is not local: the baseline is
  read by `cox_baseline()`, the predictions and the sampler, which
  would need the softmax-and-Dirichlet treatment frmtmb.sample now
  gives a `mo()` simplex (brms puts `dirichlet(1)` on `sbhaz`). Lane
  optima, review m9.
- **`escape_stationary()`'s sibling paths**: the objective's state is
  settled after an escape and after `mo_search()` (`obj_settle()`);
  any other post-fit refit that keeps a run other than the last must
  call it too. Lane optima.
- **A dense block with several scales still has a `+Inf` or `NaN` log
  density where a log sd underflows.** Lane surface floored the
  variance of the six one-scale dense structures (`sd2_floored()`), but
  at log sd -1137.64 with a nonzero field `us` with two or more
  coefficients and `cs` give `+Inf`, and `ar1` gives `NaN`, on 0.68.1
  and 0.69.0 alike (`us` with one coefficient and `diag` give `-Inf`;
  `dev/surface-rev-otherblocks.R`). On `frm_sample(fit)` every prior
  is flat on the fit route and frmtmb.sample keeps these blocks
  centered (`ncp_plan()`), so a long first leapfrog step could stick a
  chain there as it stuck the exact `gp()` chain. Latent: 16 sampler
  seeds of `(1 + x | g)` never went below theta -15.3
  (`dev/surface-rev-ussweep.R`). The remedy is the same floor per
  scale, or a log-scale density. Lane surface, review m4.
- **A degree-7 raw polynomial with `(1 | g)` keeps standard errors far
  from lmer()'s.** Seed 7 (`dev/setier-rev-deg7.R`) loses x3 to x5 and
  keeps five SEs up to 0.946 relative from lmer()'s (before lane
  setier it lost x2 to x6 and kept `(Intercept)` at 37.7 against
  lmer()'s 2774). The tier-3 threshold `tau = 10 * ||E|| / gap`
  explodes when the eigenvalue gap is 1.4e-11, so only dominant
  loadings are named. Pre-existing. Setier review, m7 and the final
  check.
- **Curved ridges with an eigenvalue above 1e-9 keep silent finite
  standard errors**, on 0.68.1 and 0.69.0 alike: `a + exp(b)` on 9 of
  20 fits, the same ridge with `(1 | g)` on 20 of 20, c0k on 16 of 18
  code-0 fits, `a * b` on 3 of 20. The objective along the exact ridge
  moves by at most 1e-8 (`dev/setier-rev-aexpb.R`,
  `dev/setier-rev-c0k-sweep.R`). The curvature probe never sees them,
  because their eigenvalues are above the threshold that sends a
  direction to it. Pre-existing. Setier review, m8.

### Open - low

- **Two results of the merged tree depend on OpenBLAS's LAPACK**,
  found at the 0.69.0 consolidation and corrected by its review
  (`dev/reviews/2026-10-07-release.md`, B1). Both fail under OpenBLAS
  0.3.26 and 0.3.32 when OpenBLAS also provides LAPACK (the
  ubuntu-24.04 and 26.04 runners' configuration) and pass with either
  version as BLAS only; neither file runs on CI. (1) `test-fuzz.R`:
  the spec of seed 20379118 (`Gamma(link = "log")`, `REML = TRUE`,
  `weights`, `ar1(tt + 0 | g)`, `mo(xo) * z`, `shape ~ 1 + z`) stops
  with "NA/NaN gradient evaluation" after the restart from the best
  point; with the reference BLAS it ends at code 1 (objective
  -4074.54), with 0.3.26 as BLAS only at code 1 (-354.88), and
  rellib-r6 at code 1 in each configuration measured; lane optima's build
  errors the same way. The trigger is optima's chart on a REML
  likelihood without a bound; the mechanism is 0.68.1's restart (the
  Open - high item above). It is not a rate regression: over 300 seeds
  of the generator per cell (`dev/relrev069-fuzzsweep.R`), the spec's
  family errors on 37, 42 and 42 under 0.69.0 (reference, 0.3.26,
  0.3.32) against 31, 52 and 52 under 0.68.1, and with `mo(xo)` alone
  on 35, 35 and 35 against 58, 48 and 48 (`dev/rel069-fuzz1.R`,
  `dev/rel069-log/fuzz1/`). (2) The ported row `data-helpers:7`
  (`expect_silent(fitted(fit1, newdata))` on fixture 1): with OpenBLAS
  LAPACK the merged build's fixture converges (code 0; base, setier
  and surface stop at code 1, optima alone at code 0 without the
  warning), and the prediction warns that it moves along a lost
  direction. That is a false alarm: the alias is (Intercept, Age,
  `s(Age).fx1`), the five predictions are estimable, and the lost
  basis the prediction test reads drops the intercept's share because
  the intercept keeps its SE (the Open - high item above;
  `dev/rel069-dh7.R`, `dev/rel069-log/dh7/`). Fix to try: test a
  prediction against the lost directions before the loadings of kept
  parameters are zeroed, or scale `se_pred_tol` to the basis's own
  precision.
- **A mixture precision with an infinite MLE fits with no warning.**
  `test-bcm-latent-mixtures.R`'s Malingering_2 fit,
  `mixture(beta_binomial, beta_binomial)` on `bcm_malingering_data()`:
  component 2's phi runs to the binomial limit, and nlminb stops where
  the flat ridge lets it, phi2 7.19e8 with the reference BLAS and
  2.13e9 with OpenBLAS 0.3.32. With the reference BLAS the fit gives
  no warning at all; with OpenBLAS only "Optimizer did not report
  convergence: false convergence (8)" (`dev/ciharden-rev-bcmmal.R`).
  The flat-direction warning does not fire on a parameter that runs to
  infinity. Lane setier's boundary rule covers random-effect blocks,
  not a family's precision. Lane ciharden, review m5.
- **`update(prior =)` replaces the stored prior; brms merges** the new
  rows with the stored user rows (a new row wins its slot) and drops
  what matches nothing. A stored `normal(0, 0.05)` on `b_x` is
  discarded by `update(fit, prior = <a sigma prior>)` (b_x 0.236 to
  0.766; merged 0.236). frmtmb could merge and still check the new
  rows strictly. Lane surface, review m1.
- **`mo()` seed 194** of `dev/optima-mo-study.R` stays 0.0489 below
  the exact maximum, at the other sign's cone, which 0.68.1 found: the
  search's vertex start (chosen from the gradient at `b = 0` with the
  other parameters held) did not lead into that cone's optimum there.
  A start at each vertex would reach it at D times the cost. Lane
  optima, review m8.
- **A `simo` prior row is refused by `frm()` and `frm_sample()`**:
  brms accepts `prior(dirichlet(c(2, 1, 1)), class = simo, coef =
  moincome1)`. `set_prior()` has no Dirichlet density; the sampler's
  `mo_simplex_nlp()` is where a concentration vector would go. Lane
  optima.
- **`as_tmbstan(fit)` on every `mo()` fit**, ML or MAP, samples the
  fit's own tape, whose simplex coordinates have no density: on a
  weakly identified simplex the draws collapse to the barycenter (sd
  1e-17, coordinates to 7e17). By contract the route adds nothing; it
  says so in a message and `?as_tmbstan` names `frm_sample()`. Lane
  optima.
- **`mo(x) * f` with a factor and `mo()` in group-level terms**, both
  of which brms fits, are refused (pre-existing). Lane optima, review
  m11.
- **The SE check misses an exactly flat simplex direction that is not
  along a coordinate**: with four categories and the third unobserved
  (`dev/optima-segap.R`), the flat direction warned on 6 of 9 fits and
  was silent on 3 on lane optima's build. Lane setier's tier-3-first
  rule may cover it; rerun the script on 0.69.0 before starting. Lane
  optima.
- **The `cs()` ordinal-mixture wall has no message naming its cause**:
  most `mixture(cumulative(), sratio())` fits with `cs(z)` stop with
  "false convergence (8)" where a row's category closes in one
  component; the convergence warning is right but does not say why.
  Lane optima.
- **vigport defect 7 and fixes F11 seed 8, measured again**: `fit_loss2`
  stops 0.0319 below the best of 40 jittered starts (-358.5132) and
  says so ("false convergence"); `restarts` does not move it
  (`dev/optima-loss2.R`). F11 seed 8 stays 0.1418 below the optimum a
  start at the boundary finds; a one-evaluation screen does not see it
  (`dev/optima-f11.R`, `dev/optima-f11-screen.R`). Both stay open
  above. Lane optima, item 4.
- **A multi-start for ordinal mixtures, measured again and not built**:
  best of 10 jittered starts beats the default on 20 of 40 fits, but
  on 12 of those 20 the better point is degenerate
  (`dev/optima-ordmix-starts.R`). Lane optima, item 4c.
- **`parnames.frmtmb_fit()` and `parnames.frmtmb_draws()` decide
  whether to warn from which generic the `parnames` binding resolves
  to**, not which generic dispatched; a call through a generic the
  binding does not name can still warn twice or not at all. Not
  measured to happen; `posterior_samples()` asks the dispatching frame
  instead, and the two could share that test. Lane ciharden.
- **`logitnormal_sd()` at extreme inputs** (`pp_mixture()`'s
  `Est.Error`): at p = 1e-30 with logit SD 10 the single `integrate()`
  misses the logistic's transition, 2.21e-6 against 1.41e-6 (relative
  0.57; absolute size negligible); at logit SD 0 the result is
  rounding noise (3e-16), not 0. Surface review, n2 and n3.
- **frmtmb.sample's exact `gp()` sampling mixes poorly from either
  route** (R-hat 1.56 to 1.89 at 600 iterations, 197 to 215
  transitions over the tree depth): the field is sampled centered, a
  funnel, where brms samples the standardized `zgp`. Confirms gpby r4
  (above). Lane surface.

### Upstream, listed in `dev/upstream-bugs.md` for the user to file

- RTMB-3: RTMB's `pnorm(log.p = TRUE)` derivative is 1.2 percent off at
  9.4e6 and NaN near 4.46e9 (`dev/optima-pnorm-scan.R`); frmtmb's
  probit continues as the tail's leading term past `|eta| = 1000`.
  Lane optima.

