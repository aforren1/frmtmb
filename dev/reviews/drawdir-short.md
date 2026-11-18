# Short check of lane drawdir, fix round 1 with option (b) (2026-11-17)

Review `dev/reviews/2026-11-17-drawdir.md`, the lane's section 12 of
`dev/drawdir-findings.md`, user decision (b). Lane library
`wt-drawdir-lib`, rellib-r22 behind it; reference BLAS. Scripts
`dev/reviews/drawdir-short-1.R` (parts A to E) and
`dev/reviews/drawdir-short-bit.R`; logs `drawdir-short-1ab.txt`,
`drawdir-short-1cde.txt`, `drawdir-short-bit-{lane,base}.txt`,
`drawdir-short-t-{pdc,grby}.txt`, all in `dev/reviews/`. Every `lib:`
line is the one claimed. Code read against `git show
b255350e:R/predict-brms.R` with a plain `diff` of copies. No lane file
edited, no git writes, at most two R processes at once.

## Verdict: SHIP

All six checks hold. Three LOW items to file (end of this file); none
blocks.

## 1. FIX 1: non-finite ends

`draw_check_dirs()` reads a non-finite end again at t = 0.75, then 0.5,
through the side law with deviates `t Z`. `draw_logratio()`'s excess is
`drop - |Zt|^2 / 2`, so at t it is measured against the law's
`t^2 k^2 / 2`; the lane divides it by t^2 and compares with 6, which is
the same as comparing the excess at t with `6 t^2`. The slice is read at
`t k` conditional sds and scaled the same way. Exact for a quadratic
misfit.

slpw:12 (part A, independent reading of the ends):

    +theta_3 (theta)  e(1) Inf | e(0.75) -0.6404 [thr 6t^2 3.375; /t^2 -1.138] | lane at 0.75 stat -1.138
    -theta_2 (theta)  e(1) Inf | e(0.75) -0.8129 [thr 6t^2 3.375; /t^2 -1.445] | lane at 0.75 stat -1.445
    any stat infinite: FALSE
    cached: fire TRUE dir 202.6839 dir_par theta_3 dir_n 18 warn_only FALSE
    any condition text contains 'Inf': FALSE

The review's -0.64 at t = 0.75 is reproduced. The fire is kept by
finite directions (-theta_3 202.7, +theta_1 48.4, -sigma_(Intercept)
48.0, +theta_2 26.0).

Constructed case (mine, not the lane's test): 2-D, correlation 0.9, a
quadratic misfit `a/2 (u'p)^2` along u = (1, 1)/sqrt(2), a = 4, the
objective NaN (or +Inf) where |u'p| > 1.6, so every end fails at t = 1
and 0.75 and reads at 0.5. Analytic stat `a/2 (u'Delta)^2 - a/2 (u'
2 csd e_1)^2` = 13.68:

    constructed (end NaN): stat 13.680000 want 13.680000 (unscaled e(0.5) 3.4200) n 16 fires TRUE
    constructed (end Inf): stat 13.680000 want 13.680000 (unscaled e(0.5) 3.4200) n 16 fires TRUE

Read unscaled against 6 it would be 3.42 and silent; the lane scales it.
n = 16 = 4 ends + 4 at 0.75 + 4 at 0.5 + 4 slices.

Holds. One consequence, not a defect of the fix: an objective that is
non-finite at every read t (a hard wall within half the distance) is
left unread, and the reading returns NULL (`constructed hard wall ...
NULL`). 0.83.0 had no reading there, so nothing is lost; LOW item 2.

## 2. Option (b): the set and "all" against "any"

Where theta-class parameters are created (`R/frame.R` 5580 to 5618):

- `theta` (5591): the parameters of every random-effect block's
  covariance structure. Core structures (`R/covstruct.R`): `us`,
  `diag`, `homdiag`, `us_t`, `diag_t` (log sds and correlation
  parameters; the t df is fixed in the block, not a parameter), `cs`,
  `homcs`, `ar1`, `hetar1`, `toep`, `homtoep` (sds and correlations),
  `ou` (sd, decay rate), `exp`, `gau`, `mat` (sd, range, Matern
  smoothness), `gr_cov`, `gr_prec` (scale), `gp`, `hsgp` (sdgp,
  length scales), `smooth` and `ps()` terms (homdiag log sd: the
  penalty as a scale), `rr` (factor loadings), `equalto` (0
  parameters), `latres` (latent residual sd and its autocorrelation,
  `R/autocor-latent.R:560`), and extension structures
  (frmtmb.spatial `car`, `spde`).
- `thetaac` (5599): autocorrelation parameters (ar, ma, cosy, unstr,
  sar).
- `thetar` (5603): residual correlations.
- `logsdme`, `thetame` (5608 to 5610): `me()` latent sds and
  correlations.
- Not in the set: `beta`, `betad` (every linear predictor's
  coefficients: sigma, mixture weights, zi, hu ...), `meanme`, and the
  `extras` components: `zeta*` (mo() simplexes, 5075), `bcs*` (cs()
  coefficients, 5191), `tau_raw` (ordinal thresholds, family
  `extra_pars`), `sbhaz_raw` (frmtmb.survival baseline hazard).

Every parameter in the set is a covariance parameter; none is a
mixture weight, a cs() or ordinal parameter, or a regression
coefficient. The set belongs. Two remarks, neither a defect: ranges,
Matern smoothness, ar correlations and rr loadings are covariance
parameters but not "a variance toward 0", so the stated reason (flat
toward a variance of 0) fits them only loosely (their unconstrained
transforms have flat tails, which is the same problem); and
`?predict.frmtmb_fit`'s list of variance components does not name
spatial ranges, rr loadings or the latent-residual sd.

Firing directions on the 12 fires and test-gr-by.R's fit (part B):

    slpw:12   +theta_1(theta) 48.4, +theta_2(theta) 26.0, -sigma_(Intercept)(betad) 48.0, -theta_3(theta) 202.7 | all FALSE any TRUE | warn_only FALSE | frmtmb_draw_pareto
    slpw:14   +theta_2 9.0, -theta_3 11.3          | all TRUE | warn_only TRUE | frmtmb_draw_wide_tail
    slpw:15   +theta_2 27.1                        | all TRUE | warn_only TRUE | wide_tail
    slpw:18   +theta_2 16.9, -theta_3 16.0         | all TRUE | warn_only TRUE | wide_tail
    slp10:10  +theta_1 8.5, -theta_3 11.6          | all TRUE | warn_only TRUE | wide_tail
    prsp:1    +theta_2 12.1                        | all TRUE | warn_only TRUE | wide_tail
    prsp:10   +theta_2 51.3, +theta_3 57.5         | all TRUE | warn_only TRUE | wide_tail
    ordrs:9   +theta_3 6.2                         | all TRUE | warn_only TRUE | wide_tail
    S0.1:14   +theta_2 42.2, -theta_3 46.1         | all TRUE | warn_only TRUE | wide_tail
    S0.1:16   +theta_2 25.7, +theta_3 29.0         | all TRUE | warn_only TRUE | wide_tail
    S0.1:24   +theta_2 12.7, -theta_3 11.6         | all TRUE | warn_only TRUE | wide_tail
    S0.2:8    +theta_2 7.0                         | all TRUE | warn_only TRUE | wide_tail
    gr-by     +theta_2 9.1, -theta_3 11.8          | warn_only TRUE | wide_tail

"All" against "any" differs on slpw:12 only (1 of 13). There "all"
reweights, which is right on the evidence: the normal law's new-level
w90 is 19471 against the reweighted 114.78 and a posterior's 120.1
(review section 6); "any" would keep the 19471. So "all" is the better
reading of the decision. But the reweighting message on slpw:12 names
`theta_3` (the largest excess) as the direction "so the draws are
reweighted", while `?predict.frmtmb_fit` says a variance component's
direction is not reweighted; the direction that caused the
reweighting, `sigma_(Intercept)`, is not named. LOW item 1.

## 3. frmtmb_draw_wide_tail per call

S0.1:14 (part C). The verdict, `warn_only` included, is cached in
`fit$cache$draw_check`; the warning is raised by `draw_reweight()`,
once per `predict_par_drawer()` call, and `predict_simulate()` calls
that once per `predict()` call:

    predict() first call (check computed)      wide_tail 1
    predict() second call                      wide_tail 1
    predict(summary = FALSE)                   wide_tail 1
    predict() in a loop of 3                   wide_tail 3
    predict(propagate_error = FALSE)           wide_tail 0
    fitted(); fitted(new level)                wide_tail 0
    conditional_effects(); method = 'posterior_predict'   wide_tail 0
    pp_check(); simulate()                     wide_tail 0
    predict() on a fresh cache                 wide_tail 1

One warning per call, never per draw batch. Only `predict()` draws the
parameters from the normal law: `fitted()` and `conditional_effects()`
(both methods; its predict band simulates at the estimates,
`R/conditional-effects.R` 3100 to 3125), `pp_check()` and `simulate()`
do not, and stay silent, which is correct. A frmtmb_fit has no
`posterior_predict()` method; `predict(summary = FALSE)` is that path
and warns once. A user calling `predict()` in a loop gets one warning
per call, as with 0.83.0's `frmtmb_draw_pareto`. Acceptable: the
warning is about the draws each call returns.

## 4. The new-level caveat

S0.1:14 (part D):

    predict(new level, allow_new_levels = TRUE)                    wide_tail 1 caveat 1
    predict(old level, allow_new_levels = TRUE)                    wide_tail 1 caveat 0
    predict(old level)                                             wide_tail 1 caveat 0
    predict() in-sample                                            wide_tail 1 caveat 0
    predict(new level, re_formula = NA, allow_new_levels = TRUE)   wide_tail 1 caveat 0
    predict(new level, sample_new_levels = 'gaussian')             wide_tail 1 caveat 1

The caveat appears only where the call predicts a level the fit did not
see. `allow_new_levels = TRUE` alone does not add it, and `re_formula =
NA` (no group effect read) does not either. Holds.

## 5. Bitwise

`drawdir-short-bit.R`, both arms pinned to cores 0 and 1 (`start
/affinity 3`), seeded, 200 in-sample and 100 new-data draws:

    grp    fire FALSE dir 4.376639   par identical TRUE | in-sample draws identical TRUE | newdata identical TRUE
    pois   fire FALSE dir 0.06221108 par identical TRUE | ... TRUE | ... TRUE
    sleep  fire FALSE dir 0.322998   par identical TRUE | ... TRUE | ... TRUE
    slope  fire FALSE dir 0.5825509  par identical TRUE | ... TRUE | ... TRUE
    mix    fire TRUE excess 8.985 (summed; frmtmb_draw_pareto on both)
                                     par identical TRUE | ... TRUE | ... TRUE

Unfired fits keep 0.83.0's draws bit for bit, and 0.83.0's summed path
reweights exactly as before.

A non-variance per-direction fire against the first build: the first
build is installed nowhere, so no bitwise comparison is possible. What
was checked instead: the reweighting path (`draw_logratio()`,
`draw_resample()`, the side law) is unchanged from b255350e in the
diff, and FIX 1 changes only which ends are read, not the draws. On
slpw:12 (fired by a `betad` direction) the lane's reweighted arm at
20000 draws, seed 1 (part E), equals the review's first-build line on
all five printed figures:

    lane    mae 7.4062 | frmtmb_draw_pareto | cov90 0.839 w90 114.78 cov98 0.957 w98 165.38 | warn_only FALSE
    review  mae 7.4062 | paretok 0.8        | cov90 0.839 w90 114.78 cov98 0.957 w98 165.38

## 6. Tests and the test-gr-by.R hunk

    RESULT test-predict-draw-check.R pass=62 fail=0 err=0 skip=0 warn=0
    RESULT test-gr-by.R pass=71 fail=0 err=0 skip=0 warn=0

The test-gr-by.R hunk wraps the refusal `expect_error(predict(fit,
newdata = <level "new">), "New levels in grouping factor `g`: new")` in
`allow_warnings(..., "far tail")`, not required. It does not weaken the
test: the error and its message are still asserted, and only a warning
whose text contains "far tail" is muffled; any other warning still
escapes and fails. The pattern is not unique to the new warning: a
`me()` warning (`R/me.R:1088`, "diverges in the draw's far tail") also
matches, but this fit has no `me()` term. "normal law's far tail" would
be exact; optional. The fit does fire (part B: theta_2 9.1, theta_3
11.8, `warn_only` TRUE). Not requiring it is right: it is incidental to
the test. A side effect it shows: `predict()` warns about its draw law
and then refuses the call (the refusal is raised in `R/predict.R`
after the drawer). LOW item 3.

## LOW items to file

1. On a fire with variance and non-variance directions (slpw:12), the
   reweighting message and warning name the largest direction
   (`theta_3`, a variance component), not the non-variance direction
   that made the draws reweighted (`sigma_(Intercept)`), which reads as
   a contradiction of `?predict.frmtmb_fit`. Name the non-variance
   direction where `vc_only` is FALSE, or both.
2. A direction whose objective is non-finite at t = 1, 0.75 and 0.5 is
   left unread: a hard support wall within half the distance is not
   reported by the reading (no regression against 0.83.0).
3. `predict()` on a call it then refuses (an unrouted new level on a
   `gr(by =)` fit) first raises `frmtmb_draw_wide_tail`; the draw law
   is read before the data are routed.
