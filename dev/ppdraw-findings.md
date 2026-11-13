# Lane ppdraw: predictive draws at a boundary (round of 2026-11-11)

Base e8abed0a (frmtmb 0.81.0), branch wt-r82-ppdraw, library
`C:/Users/adf44/source/r/wt-ppdraw-lib`. Brief
`round-briefs/r082-ppdraw.md`. Backlog items 3269 and 3277 (lines at
e8abed0a). Logs named below are local (gitignored); every number this
record rests on is pasted here.

## The cause of both defects

`predict(propagate_error = TRUE)` draws the outer parameters from
`N(theta_hat, vcov(full = TRUE))` (`predict_par_drawer()`). Both fits
have one parameter whose estimate sits at the END of its range on its
internal scale, where the likelihood is flat toward that end. The
standard error is finite and huge: it measures only how little the
likelihood curves there. The normal law then puts half its mass on the
other side, on values the data rule out. No condition was raised: the
SE check reports no lost parameter on either fit
(`sdr_of(fit)$se_lost` is NULL, `dev/ppdraw-probe1-base.txt`).

- **3269** (`zero_inflated_binomial`, seed 7, 200 rows,
  `dev/ppdraw-probe1.R`): `zi_(Intercept)` = -22.197, SE 4926.7. Half the
  draws put zi near 1. `predict()` mean 3.042962 against `fitted()`
  5.735; `propagate_error = FALSE` gives 5.735362, so the simulator is
  right and the draw law is not. Lane: 5.7257.
- **3269, second case** (`se(s, sigma = TRUE)` with sigma estimated at
  zero, the construction of `dev/ndaterm-probe.R`, seed 7;
  `dev/ppdraw-probe4.R`): `sigma_(Intercept)` = -10.656, SE 3598.17. On
  rellib-r20 the predict column sds are 1.60289e+141 and 4.999856e+141
  against simulate()'s 0.04713429 and 5.145934, with 16 "NAs produced"
  warnings and the masked-draw warning (159 of 400 per row). Lane:
  0.05699568 and 5.272561, the message naming `sigma_(Intercept)`.
- **3277** (`mixture(cumulative(), cumulative())`, seed 7, 400 rows):
  the second component's log threshold gap
  `log(mu2_Intercept[2] - mu2_Intercept[1])` = -13.104, SE 3840.47. A
  draw on the positive side makes the gap `exp(+thousands) = Inf`, the
  thresholds after it infinite, the category probabilities NaN, and the
  inverse-CDF draw of `ord_sim()` NA. NA share 0.4275 on newdata
  (`ndraws = 4000`, seed 13), 0.4175 in-sample, 0 under
  `propagate_error = FALSE`. Thresholds are not distributional
  parameters, so `dpv_row_finite()` never saw them, and nothing counted
  the NA. Lane: NA share 0 on newdata and in-sample.

The triage named "the boundary message and `frmtmb_boundary_fit`" as
the seam. Neither fires on these fits (the SE check's boundary rule
reads random-effect parameters only; lane boundfam extends it this
round). So the rule here belongs to predict() and reads the likelihood
directly. It composes with boundfam: a parameter that the SE check
holds has variance 0 in the draw law, and `draw_edge_scan()` skips it.

## What changed

`R/predict-brms.R` only, plus tests and NEWS.

1. `draw_edge_held()` / `draw_edge_scan()`: once per fit (the verdict
   is cached in `fit$cache$draw_edge_held`), the objective with each
   outer parameter moved by 2 of its own draw-law standard deviations,
   the others at the estimates. The normal law says the log-likelihood
   falls by at least 2 there. A change of at most `draw_edge_flat` =
   1e-4 on either side holds the parameter at its estimate in every
   draw, as a variance component at its boundary and a separated
   coefficient already are (`D[held, ] <- v0[held]`, beside
   `latres_held_pos()`). `predict()` says so with a message of class
   `frmtmb_draw_held` naming the parameter. A message and not a
   warning: a zero-inflation model fit to data without excess zeros is
   an ordinary, correct fit, as a variance component at zero is
   (`se_boundary_act()`).
2. `predict_report_sim_na()`: a cell that the family's simulator left NA
   at finite parameters is counted per row and warned about (class
   `frmtmb_predict_na`), so a draw that cannot be made is never a silent
   NA. A `cs()` response under `summary = TRUE` is exempt because
   `predict_category_props()` already warns with its own count (the
   user's decision of 2026-10-06). `summary = FALSE` on such a fit now
   warns where 0.81.0 returned the NA cells in silence (BREAKING, in
   NEWS; `test-cumulative-cs.R` updated to require the warning).

### The threshold, measured (`dev/ppdraw-probe6.R`, log `dev/ppdraw-probe6.txt`)

The first build used 0.05. The one-file-per-process run of every suite
on it held one parameter outside this lane's own tests: `theta_2` (a
random slope's log sd, -1.963, SE 2.757) of `test-re-formula-partial.R`'s
fit. That fit is interior: lme4 calls it not singular (`isSingular`
FALSE, slope sd 0.14038), and the slices change by 0.024 toward zero
and 55.13 away (`dev/ppdraw-probe5.R`). A check that fires on a correct
model is worse than none, so the smaller change of the two sides was
measured on every parameter of the sweep's fits, the backlog fits and
15 random-slope fits with true slope sds 0.05, 0.1 and 0.2 (seeds 1 to
5), printing those below 0.5:

- 14 estimates at the end of their range (zi, zoi or hu intercepts on
  data with no excess zeros in 11 families, the 3269 fit, the 3277 gap,
  the se(sigma = TRUE) sigma): 1.29e-08 to 1.08e-07.
- 18 interior estimates near an end: 0.0219 and up. Two were below
  0.05: slope sd 0.20 seed 2 (`theta_2`, 0.0219) and the
  re-formula-partial fit (0.024).

1e-4 sits three orders above the first group and two below the second.
`test-predict-draw-edge.R` pins the re-formula-partial fit as the
absent case.

CORRECTED in the fix round (review m3): the margin to interior
estimates is smaller than this says. The review's gaussian `(1 | g)`
fits, 6 groups of 20 (`dev/reviews/ppdraw-review-g.R`, seeds 1 to 40),
have interior group log sds that lme4 does not call singular with flat
sides 0.00298, 0.00727, 0.0129, 0.0132 and 0.0142: the smallest is 30
times 1e-4, not 200. The flat side of an approach to a limit is the
likelihood's gain over the limit, which is continuous down to 0, so
the rule is "a gain over the limit of at most about 1e-4", not "at the
end of the range". Below 1e-4, not drawing into that side is the
honest reading of the gain.

### Seen to fail on rellib-r20 (`dev/ppdraw-newbase.log`)

`test-predict-draw-edge.R` pass=6 fail=9 on rellib-r20, pass=15 fail=0
on the lane. The behavioral failures, not missing symbols:
`max(abs(z))` 18.8 (the zib boundary predict mean, in Monte Carlo
standard errors); the se(sigma = TRUE) sd ratios 3.40e+142 and
9.72e+140; `anyNA(p)` TRUE and categories outside 1:4 on the ordinal
mixture. `test-cumulative-cs.R` fail=2 on rellib-r20 (the NA warning
absent on `summary = FALSE`), pass=54 on the lane.

### Robust to rounding (`dev/ppdraw-ulp.R`, log `dev/ppdraw-ulp.txt`)

x perturbed by `1 + 2e-15 k`, k = 0..20, run on the 0.05 build. zib
zi = 0: 21 of 21 held (estimate -22.1970 every time); zi = 0.3 and 0.02:
0 of 21 held each. Ordinal mixture: the gap estimate moved between
-13.6958 and -13.0469, held in 21 of 21, NA in 0 of 21. Rerun on the
final (1e-4) build, log `dev/ppdraw-ulp-final.txt`: the same counts,
zi = 0 held 21 of 21, zi = 0.3 and 0.02 held 0 of 21 each, ordinal
mixture held 21 of 21 with 0 NA and 0 categories outside 1:4, gap
between -13.6958 and -13.0469.

The paired reference in the test: predict()'s mean is compared with the
mean of the expected response at the same parameter draws (same seed;
the drawer is the first thing predict() takes from the stream), which
differs from predict() by observation noise alone. Against `fitted()`
the zi = 0.02 fit sits at z = -3.73 on every k: that is the draw law's
spread around a skewed likelihood, not a defect, and the test does not
ask it.

### Cost (`dev/ppdraw-cost.R`, log `dev/ppdraw-cost.txt`, 0.05 build)

Two objective evaluations per parameter, once per fit. Poisson GLMM,
3000 rows, `(1 + x | g) + (1 | h)`, np = 7: scan 0.068 s against
0.197 s for `predict(ndraws = 200)` with the verdict cached. zinb with
`(1 | g)`, np = 7: scan 0.410 s against 0.463 s. Far points make the
inner solve slower than at the estimates (one objective there: 0.003
s). With the default 1000 draws the scan is under a fifth of the first
call. The machine was shared by twelve lanes; these are single runs.

## brms 2.23.0 on the same fits (`dev/ppdraw-brms.R`, log `dev/ppdraw-brms.txt`)

- 3269 fit, default priors (zi `beta(1, 1)`), 2 chains of 2000, seed
  11: zi posterior quantiles 0.5, 0.9, 0.99, 1 = 0.003811, 0.012190,
  0.023090, 0.035087. `posterior_predict` mean 5.70256,
  `posterior_epred` mean 5.701684, NA share 0. brms keeps zi near 0, as
  holding does; its mean is lower than frmtmb's 5.7257 because the
  posterior keeps zi off 0.
- 3277 fit: brms sets `order = 'none'`, and with flat threshold priors
  the chains run away (thresholds near -2.2e9 and 3.5e12, Rhat 1.83 to
  2.33, 3 divergences, 137 max-treedepth). `posterior_predict` on the 3
  new rows still gives no NA (2986, 781, 750, 1483 of 6000). brms never
  returns NA here, which is what the lane does.

## The sweep (`dev/ppdraw-sweep.R`, data `dev/ppdraw-sweep-gen.R`)

Every core family with a zero-inflation, hurdle or mixture parameter,
fitted on interior data (zi or hu 0.3, zoi 0.2, mixture weight 0.4) and
on boundary data (no excess zeros; one component), seed 101, 300 rows,
1000 draws, seed 5. z is predict()'s mean minus fitted()'s in Monte
Carlo standard errors of the draws' mean; for a category-valued
response, the largest |z| over the rows and category shares (a max over
up to 1200 cells, so 3.5 to 4 is ordinary). "No propagation" is
`propagate_error = FALSE`, which tests the simulator alone. Generated by
`dev/ppdraw-table.R` from the 0.05 build's TSVs (`dev/ppdraw-sweep-base.tsv`,
`dev/ppdraw-sweep-lane.tsv`). The lane arm was rerun on the final
(1e-4) build and the table regenerated (`dev/ppdraw-table.txt`): every
z and every held entry is the same.

| family | fit | z base | z lane | z, no propagation | NA base | NA lane | held (lane) |
|---|---|---|---|---|---|---|---|
| zero_inflated_poisson | interior | 0.27 | 0.27 | 1.26 | 0.00 | 0.00 |  |
| zero_inflated_poisson | boundary | -30.79 | 2.00 | -0.77 | 0.00 | 0.00 | zi_(Intercept) |
| zero_inflated_negbinomial | interior | -0.77 | -0.77 | 2.36 | 0.00 | 0.00 |  |
| zero_inflated_negbinomial | boundary | -31.50 | 0.86 | 0.70 | 0.00 | 0.00 | zi_(Intercept) |
| zero_inflated_binomial | interior | -1.37 | -1.37 | 1.83 | 0.00 | 0.00 |  |
| zero_inflated_binomial | boundary | -31.07 | 0.95 | -1.62 | 0.00 | 0.00 | zi_(Intercept) |
| zero_inflated_beta_binomial | interior | -0.65 | -0.65 | 2.00 | 0.00 | 0.00 |  |
| zero_inflated_beta_binomial | boundary | -31.53 | 1.11 | 2.16 | 0.00 | 0.00 | zi_(Intercept) |
| zero_inflated_beta | interior | -0.46 | -0.46 | 2.46 | 0.00 | 0.00 |  |
| zero_inflated_beta | boundary | -31.72 | 1.74 | 0.05 | 0.00 | 0.00 | zi_(Intercept) |
| zero_one_inflated_beta | interior | -0.74 | -0.74 | -0.18 | 0.00 | 0.00 |  |
| zero_one_inflated_beta | boundary | 0.56 | 0.56 | 0.56 | 0.00 | 0.00 |  |
| zero_inflated_asym_laplace | interior | 1.60 | 1.60 | 1.60 | 0.00 | 0.00 |  |
| zero_inflated_asym_laplace | boundary | 0.30 | 0.30 | 0.30 | 0.00 | 0.00 |  |
| hurdle_poisson | interior | 3.63 | 3.63 | 0.34 | 0.00 | 0.00 |  |
| hurdle_poisson | boundary | -30.93 | 3.37 | -2.79 | 0.00 | 0.00 | hu_(Intercept) |
| hurdle_negbinomial | interior | 5.09 | 5.09 | -0.28 | 0.00 | 0.00 |  |
| hurdle_negbinomial | boundary | -31.39 | 1.87 | -2.41 | 0.00 | 0.00 | hu_(Intercept) |
| hurdle_gamma | interior | -0.20 | -0.20 | 1.03 | 0.00 | 0.00 |  |
| hurdle_gamma | boundary | -31.63 | 0.71 | -0.84 | 0.00 | 0.00 | hu_(Intercept) |
| hurdle_lognormal | interior | 0.43 | 0.43 | 1.33 | 0.00 | 0.00 |  |
| hurdle_lognormal | boundary | -31.61 | 1.62 | -1.61 | 0.00 | 0.00 | hu_(Intercept) |
| hurdle_cumulative | interior | 3.59 | 3.59 | 3.41 | 0.00 | 0.00 |  |
| hurdle_cumulative | boundary | 1.39e+06 | 3.29 | 3.80 | 0.00 | 0.00 | hu_(Intercept) |
| mixture_gaussian | interior | -1.62 | -1.62 | -2.81 | 0.00 | 0.00 |  |
| mixture_gaussian | boundary | -21.86 | -21.86 | 0.15 | 0.00 | 0.00 |  |
| mixture_poisson | interior | 1.38 | 1.38 | -3.04 | 0.00 | 0.00 |  |
| mixture_poisson | boundary | -7.88 | -7.88 | -1.15 | 0.00 | 0.00 |  |
| mixture_cumulative | interior | 9.67 | 9.67 | 3.62 | 0.00 | 0.00 |  |
| mixture_cumulative | boundary | 3.90 | 3.90 | 3.90 | 0.00 | 0.00 |  |

fits with |z| > 4.5: base 14, lane 4 of 30. Fits with a parameter held
(lane): 10.

Reading it:

- Same cause, fixed: all ten zero-inflated and hurdle boundary fits
  that draw from a covariance (base z about -31, the draws half zero;
  hurdle_cumulative's 1.39e+06 is category 0, of fitted probability near
  0, drawn half the time).
- zero_one_inflated_beta, zero_inflated_asym_laplace and
  mixture_cumulative at the boundary: z equals the no-propagation z
  because predict() falls back to drawing at the estimates with its
  "covariance could not be recovered" warning (the family's own verdict
  or the fit's convergence warning made a parameter lost). Said out
  loud already; not this lane's.
- Different cause, filed: the four lane rows with |z| > 4.5
  (`dev/ppdraw-probe3.R`, 400 draws, seed 5). Against the mean of the
  expected response at each drawn parameter vector, predict() agrees:
  z 0.45 (mixture_gaussian boundary), 0.12 (mixture_poisson boundary),
  largest category z 3.21 (mixture_cumulative interior), -0.43
  (hurdle_negbinomial interior; hurdle_poisson interior 0.77). So the
  simulator is right. On the hurdle fits the gap to fitted() is the
  plug-in gap of an identified fit (the drop in log-likelihood at the
  draws has median 1.68 against chi2_p / 2 = 1.68). On the three mixtures
  the normal law is wrong in the interior: median drop 5.68, 11 and 10.2
  (99th percentile 85.84, 1182, 156.2) against 3.17, 2.18 and 4.17. No
  parameter is at an end there, so holding is not the remedy. Filed.
- The 3277 test does not compare shares with fitted() for that reason:
  the fit is a weakly identified mixture (largest category z 7.03
  against fitted, `dev/ppdraw-probe3.R` on the sweep's interior
  ordinal mixture; on the 3277 fit itself, 20000 draws: P(Y = 1) 0.5653
  against 0.3831 on row 1, `dev/ppdraw-probe2.R`, while
  `propagate_error = FALSE` gives 0.3814).

## Suites

On the 0.05 build, one file per process, 3 at a
time (`dev/ppdraw-run-par.sh`, runner `dev/ppdraw-run-tests.R`, which
also logs every hold as a HELD line): ALL ran 511 of 511 (core 318,
coupling 11, eam 33, latent 12, learn 16, ode 11, sample 87, spatial 8,
survival 15), pass 32694 fail 0 err 0 skip 289 warn 0, every log's
`lib:` line on wt-ppdraw-lib. HELD lines: this lane's own test file
(zi_(Intercept) 3, the ordinal gap 1) and test-re-formula-partial.R
(theta_2 1), which led to the 1e-4 threshold above. The final build
holds fewer parameters and changes the message text; the three files
that can see the difference were rerun on it (`dev/ppdraw-new.log`):
test-predict-draw-edge.R 15/0, test-cumulative-cs.R 54/0,
test-re-formula-partial.R 41/0, no HELD line for theta_2.

Gated tier on the final build: the 47 files whose plain run skipped on
a gate (`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS`), log
`dev/ppdraw-gated.log`: GATED ran 47 of 47, pass 4602 fail 0 err 0
skip 0 warn 0, no HELD line. The 15 scale-tier skips and 1 fuzz-tier
skip of the plain run were not run (other gates, not predict's).

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
OPENBLAS_NUM_THREADS 4; Rblas 74 symbols and Rlapack 616 to OpenBLAS),
final build: test-predict-draw-edge.R 15/0, test-cumulative-cs.R 54/0,
test-re-formula-partial.R 41/0 (`dev/ppdraw-ob-out/`). 0.3.32 not run:
no assertion sits near a rounding-decided branch (the flat sides are
below 1.1e-7 or above 0.02, against 1e-4). `dev/ciharden-out/` deleted.

One-session runs of core on the final build
(`dev/release/run-sessions.sh ppdraw 3 frmtmb`, FRMTMB_RELLIB = the
lane library, log `dev/release/ppdraw.log`):

    SESSION frmtmb forward files=305/318 pass=23504 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb reverse files=305/318 pass=23669 fail=0 err=0 skip=239 warn=0
    SESSION frmtmb preload files=305/318 pass=23772 fail=0 err=0 skip=241 warn=0

The 13 files with no result are gated whole (10 `test-brms-suite-*`,
`test-latres-newrow-brms.R`, `test-linpred-layers-brms.R`,
`test-me-newdata-brms.R`); 20 LEAK lines, every one tagged
`[with load of ...]`, none untagged. The eight extension suites were
NOT run in one session: no extension changed, this change registers
nothing and sets no option, and their one-file-per-process runs were
clean. That is 24 more sessions of hours each on a machine shared by
twelve lanes; the consolidation's session runs cover them.

`R CMD check --as-cran` of core on the final build (built and checked in
`dev/ppdraw-check/`, then deleted; output kept as
`dev/ppdraw-check-status.txt`): Status: 2 NOTEs, the examples timing
(`frm_hazard_reads` 5.25 s elapsed, machine load, not this change) and
the V8 math-rendering NOTE. Tests inside the check: OK. Extensions are
unchanged and were not checked.

Bump: minor (two BREAKING behavior changes in predict()). Floor: none
new; no extension calls the changed functions.

## Fix round (review `dev/reviews/2026-11-11-ppdraw.md`, verdict FIX)

Everything above this section is the first round's record; where the fix
round changed a rule, this section wins. Reviewer scripts and logs:
`dev/reviews/ppdraw-review-*`. Coordinator's decisions 1 to 4.

### Findings reproduced on this build before the fix

From the reviewer's logs (same scripts, same lane build of 09:12) and
from `dev/ppdraw-fix-cases.R` on rellib-r20 (log
`dev/ppdraw-fix-cases-base.txt`):

- M1: zero_inflated_poisson null data, seeds 9, 13, 29 (zi -5.07,
  -4.98, -4.77, not held at 1e-4): z -6.23, -7.07, -7.57; seed 3 (zi
  -7.011, flat side 0.0017): mean 1.9756 against fitted 3.0799, P(y = 0)
  0.3984, z -49.59 at 4000 draws. G3 seed 5 new level: sd 2.215e+07
  against 0.9953 at the estimates.
- M2: `zi ~ x` (A1) z -25.91 on base, -57.73 on the first-round lane;
  `hu ~ x` (A4) -26.16 and -118.06.
- B1: `ppdraw-review-c-boundfam.txt`, and on boundfam's library
  (`dev/ppdraw-fix-bf.R`, log `dev/ppdraw-fix-bf-boundfam.txt`): a
  "limit" parameter arrives with variance 0 in `cov_fixed_prop` and its
  row intact in `cov_fixed_raw` (zi intercept 2.427e+07), NaN in
  `cov.fixed`. With `zi ~ x` both coefficients are "limit"; their raw
  correlation is 6.605e5 / sqrt(2.068e7 * 2.218e7) = 0.031.

### What changed

1. **The side law** (decision 1, `draw_edge_law()` and helpers in
   `R/predict-brms.R`, replacing `draw_edge_held()`/`draw_edge_scan()`).
   Each side of each outer parameter, the others at their estimates:
   - flat (|change| <= 1e-4 at 2 draw sds): scale 0, as before;
   - steep (the log-likelihood falls by 2 within `draw_edge_c` = 0.5 of
     2 conditional sds, i.e. by 1 csd): the side is narrowed;
   - otherwise the normal law, and if no side of any parameter is flat
     or steep the draws are bit for bit 0.81.0's (the same `rnorm()`
     matrix through the same factor; test "keeps its draws bit for
     bit").
   A parameter with a flat side is drawn on its own (its other side
   mapped), the rest from the Schur complement given it, and the rest
   are scanned again under that conditional law (decision 1c).
   Said once per fit (m2, cache flag `draw_edge_said`).
2. **B1** (with the coordinator's two addenda): boundfam's "limit"
   covers two kinds, (a) a family parameter at its end, the family's
   limiting case, and (b) a mean or covariance parameter a fit with no
   maximum leaves undetermined (raw variances 2.4e137, 5575, 6e8).
   - `fit_draw_space()` (`R/brms-shapes.R`): `boundary_only =
     all(lost == "boundary" | (lost == "limit" & draw_limit_quiet()))`,
     the proposed wiring, and it returns `raw` (`cov_fixed_raw`) and
     `lost`.
   - `predict_par_drawer()`: any "limit" that is not quiet sends the
     draw to the existing fallback ("covariance could not be
     recovered", draws at the estimates); its raw row is never drawn.
   - `draw_limit_restore()` puts back the raw row only where
     `draw_limit_held()` is TRUE (kind (a), at its end), if finite; a
     family parameter the limit leaves "undefined" (`zi_x` beside a
     zi intercept at -22) keeps variance 0, as the check holds it.
   - Wiring choice: `draw_limit_quiet()` and `draw_limit_held()` call
     boundfam's `se_limit_quiet()` and `se_limit_draw_held()` through
     `get0(..., envir = asNamespace("frmtmb"))`, guarded: on this tree,
     which has neither, nothing is called "limit" and both return FALSE.
     So the merged build needs no edit for the wiring. The test
     "a 'limit' parameter is drawn only when it is the family's limit"
     simulates all three cases (kind (a) at its end, kind (a) undefined,
     kind (b)) by mocking `fit_draw_space()` and the two predicates.
3. **csd from the Hessian** (found by the core suite on the first fix
   build: `test-se-reparam.R`, "a degree-8 raw polynomial gets lm()'s
   covariance", sd ratio off by 0.95). `draw_edge_csd()` takes `1 /
   sqrt(H[j, j])` from the fit's cached outer Hessian
   (`fit$cache$hessian_fixed`, at the estimates), which is exactly what
   a one-parameter slice sees, under REML too. Without it, `solve(V)`;
   and where V is held as a factor (`se_reparam_cov()`) and no Hessian
   is cached, no law is built, because `solve()` of that matrix has no
   digits (it narrowed the polynomial's draws about 20-fold).
4. **m5**: the NA warning counts Inf too (`!is.finite`), and its text is
   family-neutral ("returned no finite response ... the family has no
   finite draw at them").
5. **m6**: NEWS rewritten and wrapped to 80 columns.

### Where this departs from the decision, and why (measured)

- **csd, not sd, is the reference** (1a says `s_side = min(sd, d / 2)`).
  On an uncentered covariate (`x ~ N(100, 1)`, poisson and logit, seeds
  1 to 5) the intercept's marginal sd is 95 to 110 times its conditional
  sd (`dev/ppdraw-fix-c.tsv`, column `sd_over_csd`). Its slice falls by 2
  at 2 csd, so `d / 2` against the marginal sd is about 0.01: the literal
  rule narrows every such interior fit a hundredfold. Against csd its
  ratio is 0.99 or 1 and nothing changes.
- **The steep side's map** (1a says one scale per side). One scale lets
  the normal tail past z = 2 run through a wall: se(sigma = TRUE)'s
  first row drew with sd 0.61 (brms 2.23.0: 0.1077). Matching every draw
  to the point where the log-likelihood falls by z^2 / 2 (knots at z = 1
  to 4) fixed that (0.080) but followed 3269's exponential wall so
  closely that `zi ~ x` and `hu ~ x` moved to z -8.12 and -7.96. The
  build takes the smaller of the two at each z (`draw_edge_map()`): the
  decision's scale inside 2, the likelihood's limit past it. The three
  variants, each from a run of `dev/ppdraw-fix-cases.R` on its build
  (the log keeps the last):

  | variant | 3269 mean / P(y=0) / 98% width | A1 z | A4 z | se(sigma) row 1 sd |
  |---|---|---|---|---|
  | one scale per side (1a as written) | 5.6982 / 0.0106 / 7.1953 | -3.37 | -3.30 | 0.61 (ratio 13.04 to simulate()) |
  | likelihood-matched knots | 5.7275 / 0.0055 / 6.7102 | -8.12 | -7.96 | 0.0799 |
  | smaller of the two (built) | 5.7317 / 0.0048 / 6.6553 | -3.11 | -3.23 | 0.06826 |
  | brms 2.23.0 | 5.7022 / 0.0097 / 7.0005 | | | 0.1077 |

### Measuring `draw_edge_c` (decision 1b; `dev/ppdraw-fix-c.R`, log `.txt`, `.tsv`)

The unscreened ratio r = d / (2 csd) of every side (d: where the fall
is 2), on 124 fits: the 30-fit sweep, the review's zip null seeds 1 to
30, gaussian `(1 | g)` 6 groups seeds 1 to 40, random slopes (sd 0.05,
0.1, 0.2; seeds 1 to 5), 10 uncentered-covariate fits and sleepstudy.
There is no gap: smallest r per fit runs continuously from 0.001 (at-end
zi) through 0.111 (zip seed 3), 0.16 (G3 seed 5), 0.25 to 0.5 (interior
group and slope log sds near 0, lme4 not singular) to 1. c = 0.5 means a
side is narrowed only when the log-likelihood falls at least four times
faster than the normal law says. Fits whose draws change at c = 0.5:

- sweep: 12 of 28 with a covariance (the 10 zi/hu boundary fits, the
  zoib boundary fit and the ordinal mixture boundary fit, the last two
  of which predict() still draws at the estimates by its fallback);
- zip null seeds: 22 of 30 (18 at the boundary, flat side; seeds 3, 9,
  13 and 29 near it, steep side only, r 0.111 to 0.416);
- gaussian `(1 | g)` 6 groups: 7 of 40 (seeds 1, 5, 11, 17, 23, 28, 40;
  r 0.160 to 0.444), all on the group log sd's upper side;
- random slopes: 7 of 15 (sd 0.05 seeds 3, 4, 5; 0.1 seeds 2, 4; 0.2
  seeds 2, 4; r 0.247 to 0.495), the slope log sd's upper side;
- uncentered covariates and sleepstudy: 0 of 11.

Final build (Hessian csd), log `dev/ppdraw-fix-c.txt`. The fits that
keep their draws have r of 0.505 or more on every side; the largest
smallest-r that changes is 0.495. (On the first fix build, csd from
`solve(V)`, slope sd 0.05 seed 2 read 0.495 and changed; it reads 0.505
now.) A fit whose two sides both fall by about 2 at 2 sds keeps its
draws bit for bit.

### Acceptance (decision 1d; `dev/ppdraw-fix-cases.R` on base and lane, `dev/ppdraw-fix-brms.R`)

| case | rellib-r20 | lane (fix round) | reference |
|---|---|---|---|
| 3269 zib, 4000 draws: mean / P(y=0) / 98% width | 2.9481 / 0.4879 / 8.6551 | 5.7317 / 0.0048 / 6.6553 | brms 5.7022 / 0.0097 / 7.0005; plug-in 5.7344 / 0.0042 / 6.5952 |
| A1 zib `zi ~ x`, z vs fitted (400) | -25.91 | -3.11 | within 4.5 |
| A4 hurdle_poisson `hu ~ x`, z | -26.16 | -3.23 | within 4.5 |
| zip null seed 3, 4000 draws: mean / P(y=0) | 1.9756 / 0.3984 | 3.0726 / 0.0710 | brms 3.0632 / 0.0802 |
| zip null seeds 9, 13, 29: z | -6.23, -7.07, -7.57 | 0.33, -0.03, -0.59 | |
| zip null boundary seeds (18): z | -18.62 to -18.70 | 0.41 to 1.11 | |
| G3 seed 5 new level sd | 2.215e+07 | 1.01 | 0.9953 at the estimates |
| se(sigma = TRUE) column sds | 1.6e+141, 5.0e+141 | 0.06826, 4.96529 | brms 0.1077, 5.0106 |

- 3269 against brms: the zero share is half brms's and the 98% width
  6.66 against 7.00. Not met as stated, and I think brms is not the
  arbiter here: brms's `beta(1, 1)` prior on zi is, on the logit scale,
  a density proportional to `exp(x) / (1 + exp(x))^2`, about e^17 times
  higher at logit -5 than at -22, so its posterior sits where the prior
  pushes it (median 0.0037 with no excess zeros in the data). The lane's
  draws stop where the likelihood falls by z^2 / 2. Reported, not
  widened.
- zip seeds 1, 2, 4, 6, 10, 15, 18, 22 (zi -3 to -4.3) are interior,
  keep their draws bit for bit (same means as base), and their z of
  -1.1 to -3.8 is the plug-in gap filed in the first round.
- Sweep (`dev/ppdraw-table.txt`, regenerated): the ten zi/hu boundary
  z are now 0.06 to 3.05; hurdle_cumulative boundary reads 8630.97
  because fitted() puts P(Y = 0) at 1.2e-10 and the z divides by its
  root: the draws have 85 zeros in 300000 cells, against 144900 on base
  (`dev/ppdraw-fix-hc.R`). Interior rows unchanged. |z| > 4.5: base 14,
  lane 5 (the four filed in the first round, and that artifact).
- Rounding (`dev/ppdraw-fix-ulp.R`, k = 0 to 20): A1 z -3.11, A4 -3.23,
  zip seed 3 narrowed and not flat, G3 sd ratio 1.015, on all 21; the
  perturbation moves none of these estimates.
- `zi ~ x` and `hu ~ x` at 4000 draws (`dev/ppdraw-fix-a14.R`): A1 mean
  5.6798, P(y = 0) 0.0126, 90% / 98% widths 5.1100 / 7.4553 (base 2.9082
  / 0.4941; brms 5.7036 / 0.0098 / 4.9707 / 7.1353); A4 12.3313 /
  0.0069 / 11.8605 / 17.7244 (base 6.3242 / 0.4907; brms 12.3835 /
  0.0033 / 11.6305 / 16.7241). The zero share is above brms's: the
  joint move the one-parameter slices do not see, filed (LOW).
- Merged build (lane + boundfam's current `R/se-check.R`, `R/confint.R`,
  `R/autocor-latent.R`, with `se_limit_draw_held()`; installed into the
  lane library from a scratch copy outside the worktree, then the lane
  build reinstalled and checked to lack `se_limit_draw_held`):
  `test-predict-draw-edge.R` pass 41, fail 0, warn 0 (two expectations
  fewer: the slope branch that boundfam's verdict holds). Kinds through
  predict() (`dev/ppdraw-fix-merged-b.R`, log
  `dev/ppdraw-fix-merged-b.txt`): 3269, quiet TRUE, held TRUE, the side
  law's message, mean 5.722; `zi ~ x`, intercept held TRUE and `zi_x`
  FALSE (undefined), mean 5.74; all-zero poisson, `(Intercept)` and `x`
  limit, quiet FALSE, the fallback warning, draws at the estimates;
  all-right-censored gaussian seed 3, three limit, quiet FALSE, the
  fallback warning. `dev/ppdraw-fix-cases.R` on the merged build: every
  case as on the lane alone except A1 z 0.63 and A4 z 0.10 (the slope
  held at its estimate by boundfam's verdict instead of narrowed).

### Merge note for the consolidation (B1)

Rerun `test-predict-draw-edge.R`, `dev/ppdraw-fix-cases.R` (A1, A4)
and `dev/ppdraw-fix-merged-b.R` on the merged build; my run of all three
is above. No edit is needed at the merge: the predicates are found by
name. Boundfam's fitted() on these fits warns "200 of 200 predictions
move along a direction the fit does not determine" (its NaN standard
error for the "limit" parameter). WITHDRAWN in the short check: the
test helper no longer allows that warning (coordinator's item 4; it is
boundfam's defect, sent to boundfam). On the merged build the lane's
test therefore fails on it until boundfam's fix lands: the run without
the allowance read pass 35, fail 0, warn 3, the three escaped warnings
being that fitted() warning on the 3269, `zi ~ x` and `hu ~ x` fits
(`dev/ppdraw-merged-out/`, first merged run). WAITING ON BOUNDFAM.

### Cost (m1; `dev/ppdraw-fix-cost.R`, log `dev/ppdraw-fix-cost.txt`)

| fit | np | fit | law | first predict (100 draws) | second | review's scan |
|---|---|---|---|---|---|---|
| gaussian `y ~ f(200) + x`, 20000 rows | 202 | 6.0 s | 0.13 s | 0.34 s | 0.16 s | 0.17 s |
| poisson `y ~ f(60) + x + (1 | g)`, 2000 levels | 62 | 13.0 s | 3.27 s | 3.58 s | 0.20 s | 4.17 s |
| categorical, 8 categories x 7 | 49 | 2.7 s | 0.61 s | 1.25 s | 0.57 s | 0.78 s |
| zinb `(1 + x | g)`, 300 levels, `zi ~ x` | 8 | 9.0 s | 1.22 s | 1.73 s | 0.22 s | 1.30 s |

(Final build. The first fix build measured 0.16, 3.78, 0.73 and 1.34 s
for the law.) "held messages" 0 on all four: nothing changes their
draws.

Single runs on a shared machine. One evaluation per side where nothing
changes, as before; up to about 50 on a steep side, which none of these
has. No screen built: the standard error is in the parameter's units,
and the curvature at the estimate is what the normal law already reads,
so nothing at the estimate says a side is steep. Stated in the
`draw_edge_law()` comment and filed (LOW).

### Suites (fix round, final build)

- Seen to fail on rellib-r20 (`dev/ppdraw-finalbase.log`):
  `test-predict-draw-edge.R` pass 14, fail 21, err 1 (the err is the
  wiring test's mock of a function base lacks, the weak form; the 21
  failures are behavioral: z 18.8, 25.9, 26.2 and 14.9 against fitted(),
  zero shares 0.47 and 0.39, sd ratios 2.2e7 and 3.4e142, NA draws);
  `test-cumulative-cs.R` pass 52, fail 2 (no NA warning).
- One file per process, 3 at a time, every suite (`dev/ppdraw-final.log`):
  FINAL ran 511 of 511 (core 318, coupling 11, eam 33, latent 12, learn
  16, ode 11, sample 87, spatial 8, survival 15), pass 32725, fail 0,
  err 0, skip 289, warn 0, every `lib:` line the lane library. The side
  law spoke in one file outside this lane's own:
  `test-re-formula-partial.R` (`theta_2` narrowed on its upper side).
  An intermediate fix build (csd from `solve(V)`) ran the same 511 with
  one failure, `test-se-reparam.R`, which led to change 3 above.
- Gated tier, the 47 files whose plain run skips on a gate
  (`dev/ppdraw-gatedfix.log`): 47 of 47, pass 4602, fail 0, err 0,
  skip 0, warn 0.
- OpenBLAS 0.3.26 with LAPACK (Rblas 74 symbols and Rlapack 616 to
  OpenBLAS, 4 threads; `dev/ppdraw-ob-out/`): test-predict-draw-edge.R
  43/0, test-se-reparam.R 79/0, test-cumulative-cs.R 54/0,
  test-re-formula-partial.R 41/0. 0.3.32 not run: the test cases that
  `dev/ppdraw-fix-c.R` measures sit far from `draw_edge_c` (r 0.111 for
  zip seed 3, 0.160 for G3 seed 5, 0.001 to 0.003 at the end) and from
  1e-4 (flat sides 1e-8 to 1e-7, or 0.0017 and up); A1 and A4 were not
  measured there, and their verdicts held on all 21 ulp perturbations.
  `dev/ciharden-out/` deleted.
- One-session core runs (`dev/release/run-sessions.sh ppdrawfix 3
  frmtmb`, log `dev/release/ppdrawfix.log`):

      SESSION frmtmb forward files=305/318 pass=23532 fail=0 err=0 skip=239 warn=0
      SESSION frmtmb reverse files=305/318 pass=23697 fail=0 err=0 skip=239 warn=0
      SESSION frmtmb preload files=305/318 pass=23800 fail=0 err=0 skip=241 warn=0

  20 LEAK lines, all tagged `[with load of ...]`; the 13 files without a
  result are gated whole, as in the first round.
- `R CMD check --as-cran` of core on the final build (built and checked
  in `dev/ppdraw-check/`, then deleted; output kept as
  `dev/ppdraw-check-status.txt`): Status: 1 NOTE, the V8
  math-rendering NOTE on the HTML manual. Tests inside the check OK
  (21 min). Extensions unchanged, not checked.
- Bump: minor (two BREAKING behavior changes in predict()). Floor: none
  new; no extension calls the changed functions, and the merge with
  boundfam needs no edit (the predicates are found by name).

### Short check (coordinator's FIX items after the coverage check)

The short check confirmed the side law against an exact coverage
reference (`dev/reviews/ppdraw-short-cov-sum.txt`: changed fits cover
0.890 to 0.920 at 90%, where 0.81.0 was two to seven times too wide).
Its items:

1. MAJOR, fixed. `draw_edge_csd()` took `1 / sqrt(diag(H))`; on a
   `(1 + x | g)` fit with the correlation held at its limit that entry
   reads about +-1e-9, and a negative one leaked base R's unclassed
   "NaNs produced" from predict() on 6 of 16 fits
   (`dev/reviews/ppdraw-short-nan2.R`, log `-nan2-lane.txt`: lambda 0.5,
   0.6, 0.7 among them). Now `1 / sqrt(pmax(diag(H), 0))` (and the same
   on the `solve(V)` path): the entry is Inf, the `is.finite()` filter
   drops it, and the marginal sd stands. Rerun on this build
   (`dev/ppdraw-short-nan2-final.txt`): "warnings from predict():
   (none)" on 16 of 16 lambdas, and no "NaNs produced" in the log. New
   test "a curvature that rounds below zero raises no warning" (lambda
   0.5, 0.6, 0.7, 1; testthat counts any escaped warning). The before
   is the reviewer's log on the previous lane build; base has no
   `draw_edge_csd()` and cannot fail this way.
2. MINOR, docs. `draw_edge_map()` reads the slice in conditional sds
   and scales by the draw sd, so "no draw goes past where the
   log-likelihood falls by z^2 / 2" holds exactly only where the two
   agree. The comments of `draw_edge_map()` and `draw_edge_law()`, and
   the NEWS bullet, now say what holds: a draw's distance in its own
   sds is at most the slice's in conditional sds; where the parameter
   is correlated (`s > cs`) it goes `s / cs` times as far along the
   axis, as the normal law's draws do.
3. MINOR, filed (LOW, "cliff at `draw_edge_c`"): on the reviewer's
   smooth family (`dev/reviews/ppdraw-short-cliff.R`) the 98% width is
   5.17 to 5.19 up to lambda 1.4 (fall at 1 csd 2.08) and 5.56 at 1.5
   (1.94), against the reviewer's exact reference of about 5.70. Fix
   route filed: a continuous scale, `min(1, t2 / (2 csd))`.
4. Done. The fitted() warning allowance is out of the test helper; on
   the merged build the test waits on boundfam's fix (merge note
   above).
5. Reruns on this build:
   - `test-predict-draw-edge.R` (lane): pass 47, fail 0, err 0, warn 0
     (`dev/ppdraw-one.log`).
   - `dev/reviews/ppdraw-short-nan2.R`: as in item 1.
   - OpenBLAS 0.3.26 with LAPACK (Rblas 74 symbols and Rlapack 616 to
     OpenBLAS, 4 threads): `test-predict-draw-edge.R` pass 47, fail 0,
     warn 0 (`dev/ppdraw-ob-out/edge-final.txt`); `dev/ciharden-out/`
     deleted.
   - `R CMD check --as-cran` of core (built and checked in
     `dev/ppdraw-check/`, deleted; output in
     `dev/ppdraw-check-status.txt`): Status: 1 NOTE (the V8
     math-rendering NOTE); tests OK (21 min).

The other suites were not rerun for this change: it touches only the
two `sqrt()` lines (a negative or zero curvature now reads as missing
rather than NaN, which the old filter also dropped, so the only change
is the warning) and comments, NEWS and the one test file.

## Not done, and why

- Weakly identified mixtures (above): filed, needs a user decision on
  whether predict() may reweight its draws.
- The plug-in gap between predict() and fitted() on identified fits:
  filed (LOW, documentation).
- `summary()` still prints the zi intercept's SE of 4926.7 with no
  word: that is the SE check's boundary rule, lane boundfam's.
- frmtmb.sample's `posterior_predict()` was not examined for silent NA:
  its draws come from a posterior, not from this law, and the file is
  lane sampfront's.
- (First round; RESOLVED in the fix round, change 2.) For boundfam's
  merge, `fit_draw_space()` treated only reason "boundary" as a hold.
- (Fix round) The side law reads one parameter at a time; a joint move
  (`zi ~ x`) is not seen: filed LOW with the numbers. Its first-call
  cost is an inner solve per side on a fit with random effects: filed
  LOW, no screen found.
- (Fix round) 3269's zero share and 98% width stay below brms's
  (0.0048 / 6.66 against 0.0097 / 7.00), which I read as brms's prior on
  zi; reported, not widened.
- (Fix round) The extensions' suites were run one file per process, not
  in one session (unchanged packages; the change registers nothing).

## Files

`R/predict-brms.R`, `R/brms-shapes.R` (`fit_draw_space()`, a shared
file owned by neither lane: the `boundary_only` line and two returned
fields), `NEWS.md`, `tests/testthat/test-predict-draw-edge.R` (new),
`tests/testthat/test-cumulative-cs.R`, `dev/test-backlog.md` (CLOSED
lines under 3269 and 3277; five items filed), and the lane's dev scripts
`dev/ppdraw-*.R`, `dev/ppdraw-run-par.sh`, `dev/ppdraw-run-tests.R`.
