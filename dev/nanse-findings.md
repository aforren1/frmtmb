# Lane nanse: fits that lose their uncertainty silently, and frm_allfit()

Base: 9e902909 (frmtmb 0.67.0), branch wt-nanse. Library
`C:/Users/adf44/source/r/wt-nanse-lib` (core only; every extension loads
from rellib-r5 against this core). "Base" below is rellib-r5. Scripts
are `dev/nanse-*`; logs are `dev/nanse-log/` and the suite directories
`dev/nanse-suite*/`, `dev/nanse-*-log/` (local, gitignored). Every count
below was printed by the script named beside it.

## Punch round 2b

Final check, RB4: on `y ~ x + f + (1 + x | g2)` (480 rows, 40-level
`f`, no g2 variation; lme4 singular) the check took the intercept and
f2..f40 on seed 11 and the intercept on seed 17
(`dev/nanse-rev3-falseloss.R`, `-falseloss-sweep.R`). The cause: after
theta_3's empty row was removed the remainder was positive definite,
but the flat threshold was 10 times the Frobenius norm of the whole
44 x 44 block's noise (0.033 and 0.113), above the identified
intercept-against-contrasts eigenvalue (0.022 and 0.012), whose own
noise is 0.0012 to 0.0014.

The change (`se_tier3()`, R/se-check.R):
- after the empty, bound and non-finite rows are removed, a remainder
  whose smallest unit-diagonal eigenvalue is above `se_flat_tol` of
  the largest is inverted as it is, with no flat test;
- otherwise each eigenvalue is compared with its own direction's noise,
  `||Es v_k||` (a bound on that eigenvalue's error), times
  `se_dir_mult = 3`, or `se_flat_tol` of the largest, whichever is
  larger. The negative-direction test uses the same per-direction
  threshold. `tau`, the projection threshold for naming, is unchanged.

Test "an identified remainder keeps its SEs after a flat row goes"
(seeds 11 and 17: only theta parameters lost; the Intercept and f2 SEs
equal to lme4's within 2 percent of the SE): 4 failures on the round-2
build (`p2b-se-check-prefix.txt`), passes now. `test-se-check.R` 114
of 114, warn 0 (in the suite run).

Reruns on this build (`dev/nanse-log/p2b-*`):
- `dev/nanse-rev3-falseloss-sweep.R lane` (`p2b-sweep.txt`): 0 of 20
  f40 fits and 0 of 20 plain fits lose an identified fixed effect.
  Intercept SE against lme4: seed 11 0.23770 vs 0.23735, seed 17
  0.28228 vs 0.28214. The fits that lose something lose no fixed effect
  (f40 seeds 11, 17, 19; plain seeds 8, 13).
- B1 attacks: spread k = 10, 40, 60, 120 lose every a_k and b
  (`p2b-cases.txt`); crossed k = 20 and 60 lose 40 of 40 and 120 of 120
  (`p2b-spread2.txt`); gr-by x SE 0.1028, theta_2 and theta_3 NaN
  (`p2b-grby.txt`); near-collinear 0.999 and 0.99999 keep the
  reference SE to 3.2e-7 / 1.6e-6 and 2.6e-7 / 1.4e-7; the 2-D flat
  subspace loses 34 of 36 with the slope at 2.9e-7 (`p2b-b1.txt`).
  `dev/nanse-rev2-b1-rho.R`: the finite-difference path now keeps the
  pair at cor 1 - 1e-7 (SE 126.8, the reference) instead of losing
  it; at 1 - 1e-9 both paths lose it, as designed (`p2b-b1-rho.txt`).
- B2 (`p2b-pred.txt`): fitted(dpar = "a"|"b") NaN with one warning,
  as before.
- 200-seed mo() study (`mo-lane-p2b.tsv`, `mo-lane-p2b-sum.txt`):
  every column identical to round 2's `mo-lane-p2.tsv`. 26 SE warnings
  plus 17 code-7 fits, 0 silent, 0 of 357 false alarms, SEs at the
  reference to 2.87e-5, print failures 0.
- Suites, every gate on, core, frmtmb.sample and frmtmb.spline
  (`dev/nanse-suite9`): 264 of 264 files, pass 19708, fail 0, error 0,
  warning 0, skip 2 (test-scale.R). `test-diagnostics-ux.R` pass 129.
  Firings per file identical to round 2's run except the 2 of the new
  test.

Files touched: `R/se-check.R`, `tests/testthat/test-se-check.R`, this
file.

## Punch round 2

Re-check "Re-check after punch round 1" in
`dev/reviews/2026-10-06-nanse.md` (NOT MERGEABLE on RB1 and RB2).
Numbers are from the lane build in `wt-nanse-lib` (core and, new this
round, frmtmb.spline) unless a line names another build. Logs
`dev/nanse-log/p2-*`. "Seen to fail": the new tests in
`test-se-check.R` on the punch-round-1 build gave pass 95, fail 9,
error 1 (`p2-se-check-prefix.txt`; the 10-failure cap hid the rest),
and the new spline test errored (`p2-spline-curve-prefix.txt`, the
"reading the seam wrongly" stop). On this build: `test-se-check.R`
pass 110, fail 0, warn 0 after the summary() fix below; spline
`test-curve.R` pass 75, fail 0.

### RB1: a deterministic decision, and a report that cannot be lost

The clock is gone from the decision. `optimize_obj()` counts the
objective and gradient evaluations the optimizer asks for
(`fit$opt$evals`, summed over restarts and recovery starts), and
`se_check_at_fit()` builds the finite-difference Hessian of a fit with
random effects at fit time when the fit has at most 10 outer
parameters (`se_check_np_free`) or its `2 * np` gradients are at most a
quarter of the counted evaluations (`se_check_share = 0.25`).
Otherwise the check waits. `fd_hessian()` has no budget and no clock;
`fit_assembled()` no longer reads `proc.time()`. A model without
random effects reads the exact Hessian as before.

    dev/nanse-rev2-determinism.R   ridge_re_s1   ridge_re_s11
      idle (p2-det-idle.txt)       20 of 20 W    20 of 20 W
      30 busy R processes          20 of 20 W    20 of 20 W
        (p2-det-loaded.txt, dev/nanse-p2-burn.R x 30)
      reviewer, round 1 build      18 W 2 D      17 W 3 D (loaded)

`test-diagnostics-ux.R`, 8 concurrent runs started while the 30
busy-loop processes ran: 8 of 8 pass 129, fail 0, warn 0
(`p2-dux-1.txt` to `p2-dux-8.txt`). The two test workarounds are
removed: `test-case-studies.R` no longer calls `vcov()` inside the
block, and `test-diagnostics-ux.R` no longer allows the warning around
`diagnose()`; both pass (`p2-cs.txt`, the dux runs).

The deferred report:
- `sdr_of()` delivers a pending report on every call, not only the one
  that computes the sdreport.
- `se_deferred_report()` keeps the report pending when it runs inside
  a `suppressWarnings()` that a frmtmb or frmtmb.* function wrote
  (`se_muffled_inside()`: the calling function's environment is that
  namespace). A user's own `suppressWarnings()` counts as delivered.
- The readers that suppress call `se_flush_deferred()` first, so their
  caller gets the report: `fixef()`, the summary fixed and spec tables,
  `fit_draw_space()`, `me_hyper_table()`, `interop_vcov()`, and
  `fitted_point_se()`, whose handler muffles every warning.

`dev/nanse-p2-deferred.R` (the reviewer's script with the instrument
moved to `se_check_at_fit()`, since the `budget` argument it patched is
gone; `p2-deferred.txt`), c0 + exp(a)^k losing 4 of 6 SEs, every fit
waiting:

    first use summary() 1, fixef() 1, vcov() 1, confint() 1, VarCorr() 1
    each followed by vcov() 0 and summary() 0

The reviewer's unchanged `dev/nanse-rev2-deferred.R` now breaks its own
instrument (it calls `fd_hessian()` with the removed argument) and
reports nothing meaningful (`p2-rev2-deferred.txt`).
`dev/nanse-rev2-c0k.R lane 6` (`p2-c0k.txt`): 6 of 6 warn at fit time,
none deferred (round 1: 1 of 6 deferred and silent).

Tests: "whether the check waits is decided by counted work" (the rule
on synthetic counts, and two identical frm() calls give identical
`evals`) and "a waiting check reaches the caller, even through
fixef()" (fixef() delivers; an internal suppression keeps it pending; a
user's suppression consumes it). The second failed 2 expectations on
the round-1 build, the first errored.

### RB2: VarCorr() and the other products of hyp_par_cov()

`hyp_par_cov()` now also returns `Vp` (the propagation covariance
without the lost directions) and `jc` (the null basis), and
`hyp_prop_var(pc, G)` gives `diag(G Vp G')` with NaN for the rows of
`G` that `jc_nonest()` flags. Routed through it: `VarCorr()` (one
warning, "k of n VarCorr() entries move along ..."), the thresholds of
`brms_fixef_extra_vcov()` (NaN rows and columns for a flagged row),
the ordinal delta rows of the summary spec table, the `frm_multiple()`
pooling of `hypothesis()` (one warning), and `hypothesis()` itself.
`interop_vcov()` stays on the shown covariance: it hands a covariance
to marginaleffects, which multiplies it itself, so NaN is the honest
value there.

    dev/nanse-rev2-varcorr.R lane (p2-varcorr.txt), seed 1, g2 at 0
      VarCorr g1 sd Est.Error 0.233082   (base 0.233, round 1 NaN)
      residual sd Est.Error 0.07940131   (base 0.0794, round 1 NaN)
      hypothesis(sd_g1__Intercept)       0.233082
      warnings 1 (the g2 entries)

On the same fixture (`dev/nanse-p2-nanprod.R`, base and lane):
- `summary()`: g1 and residual errors finite. Base and round 1 also
  raised two "NaNs produced" warnings from TMB's `summary.sdreport()`
  (`sqrt()` of `diag.cov.random`, 7 of 24 entries at -1.2e-8 for the g2
  effects). `par_est_se()` now muffles exactly that message when the
  fit has lost SEs, since the SE warning already named the block.
- `ranef(condVar = TRUE)`: identical to base (g1 0.431 to 0.456, g2
  0 to 7.8e-5).
- `conditional_effects(fit, "x", re_formula = NULL)`: 100 of 100 bands
  finite, base and lane.

Test "a lost group sd does not take the other variance components":
VarCorr g1 equals hypothesis() to 1e-6, residual finite, g2 NaN with
one warning, `summary()` g1 and sigma finite, the re_formula = NULL
bands finite. 4 failures on the round-1 build.
`dev/nanse-rev2-port.R` (`p2-port.txt`): the flat `(1 + x | g)` block
gives NaN for its 8 entries and one warning; the residual sd keeps its
error.

### m1 follow-up: the null basis and the plot

The null basis of a removed direction keeps only the coordinates of
the parameters it named (`se_tier3()`); a kept parameter's small
loading (zeta1_2 at 0.022 on seed 12) no longer flags every prediction
through it. `dev/nanse-rev2-seed12.R` (`p2-seed12.txt`): zeta2_2 still
lost as concave; `conditional_effects()` "income" 4 of 4, "age" 100 of
100 and "income:age" 12 of 12 finite, no warnings, and the interaction
prints. (`frm_lp_basis()` carries the `beta` columns only, with the
`mo()` design at the fitted simplex, so these bands do not include the
simplex's uncertainty; that is how they were on base too,
`dev/nanse-p2-seed12.R`.) The plot draws the curve alone when the band
is NaN: `ce_plot_one()` includes the estimate in `ylim`. Test "a kept
parameter's small loading does not take its predictions" (3 failures on
round 1, including "need finite 'ylim' values").

### The 200-seed mo() study

`dev/nanse-mo-lane.R`, 4 chunks, `mo-lane-p2.tsv`, summary
`mo-lane-p2-sum.txt`: every column identical to `mo-lane-p1b.tsv`
except `ce_ok`, which changed on 18 seeds. So: 26 SE warnings plus 17
code-7 fits, 0 silent, 0 of 157 + 0 of 200 false alarms, SEs equal to
the held-coordinate reference to 2.87e-5. "conditional_effects print
fails" fell from 18 to 0: seed 12 (the null basis) and the 17 code-7
fits, whose all-NaN bands now plot as a curve.

### frmtmb.spline

`sp_one_basis()` reads `frm_lp_basis()$se_nonest`, sets those rows'
`.se` to NaN, and `sp_cov_check()` leaves them out of the comparison
(muffling core's own prediction warning there, so the user gets one).
`sp_sim_crit()` already left non-finite rows out of the maximum; with
every row NaN the simultaneous critical value is NaN instead of the
"exactly zero" stop. `frm_curve_deriv()` takes a row as lost when any
stencil row is, `frm_curve_feature()` when any of its five is, and a
difference when either grid's row is (conservative). One warning per
call, class `frmtmb_se_lost_prediction`. `dev/nanse-rev2-curve.R lane`
(`p2-curve.txt`): `dpar = "a"` curve and derivative, pointwise and
simultaneous, 0 of 9 finite with 1 warning each; the mu curve and
derivative 9 of 9 finite with no warning. Spline NEWS has a
development-version bullet.

### Recorded, not changed: the ridge with random effects

`dev/nanse-rev2-ridge-re.R lane` (`p2-ridge-re.txt`): the exact `c + dd`
ridge keeps sdreport()'s finite SEs of 3.3e5 to 2.0e6 with no warning
on seeds 3, 4, 9 and 10, exactly as base does (the reviewer's
`ridge-re-base.txt`); the other 8 seeds warn. Tier 1 is kept whenever
solve() succeeds with a positive diagonal, because changing a fit whose
plain inverse works would break the bit-identity of every healthy fit.
These SEs are huge, not confident.

### Other rev2 scripts on this build

- `dev/nanse-rev2-b1.R lane` (`p2-b1.txt`): near-collinear pairs keep
  the reference SE (exact 3.2e-7 and 2.6e-7, finite-difference 1.6e-6
  and 1.4e-7); the 2-dimensional flat subspace loses 34 of 36 and the
  identified slope keeps its SE to 2.9e-7. Unchanged from round 1.
- `dev/nanse-rev2-b1-rho.R lane` (`p2-b1-rho.txt`): as the reviewer
  recorded.
- `dev/nanse-rev2-hyp.R lane` (`p2-hyp.txt`): `b_a_f1 + b_b_Intercept`
  0.1890631, `b_a_f1 - b_a_f2` 0.2673756, `b_a_f1` NaN with a warning.

### Cost

`dev/nanse-rev-cost.R 3` and `dev/nanse-cost2.R`, idle machine
(`p2-revcost.txt`, `p2-cost2.txt`), and the new `dev/nanse-p2-mid.R`
(`p2-mid.txt`) for mixed models near the boundary:

    model           pars  evals  2np/evals  decision  ratio  control
    glm_f300         301                    exact     1.230  0.967
    glmm_f100        102                    waits     0.991  1.009
    lmm_f400         403                    waits     1.005  0.990
    pois_p12          14    237    0.12     fit time  1.186  1.015
    pois_p20          22    292    0.15     fit time  1.233  1.028
    pois_p40          42    319    0.26     waits     1.013  0.980
    gaus_p20          23    294    0.16     fit time  1.188  0.993
    gaus_p40          43    242    0.36     waits     0.975  0.981
    seven small models (cost2)              fit time  1.117 to 1.227,
                                                      control 0.965 to 1.024

The glm pays the exact AD Hessian at fit time and sdreport()'s own
finite-difference Hessian at `summary()` (fit + summary 1.136); the
NEWS now says `summary()` reuses the finite-difference Hessian, which
is the case for models with random effects.

### Suites

`dev/nanse-suite.sh dev/nanse-suite8 16`, all eight suites, one file
per process, every gate on (`NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`, and
new in this round's runner `FRMTMB_DRMTMB_FIT_TESTS` and
`FRMTMB_FUZZ`), the lane core and the lane frmtmb.spline (every `lib:`
line): 340 of 340 files, pass 23399, fail 0, error 0, warning 0, skip
15 (the seven `test-scale.R` files, `FRMTMB_SCALE_TESTS`).

Firings (`dev/nanse-diag.sh` on the 28 files where the warning fired,
`dev/nanse-fire-sum.R dev/nanse-p2fire-log`, `p2-fire-sum.txt`): 118,
all at optimizer code 0 on fits where base's sdreport() has a
non-finite SE; by reason bound 10, flat 108, concave 0. Of these, 20
are `test-se-check.R`'s own, 52 are `test-fuzz.R` (the fuzz gate was
off in round 1), and 1 is the new spline test. The rest, 45, compare
with round 1's 41 and 43: the difference is the fits that the clock
used to defer (test-id-kron.R now 2 at fit time, and so on), and the
count is now the same on every run.

`sh dev/nanse-check.sh frmtmb` and `sh dev/nanse-check.sh
frmtmb.spline` on the final tree: Status: 1 NOTE each, V8 unavailable
for the HTML manual (`dev/nanse-check/`).

### Files touched in this round

- Core R: `R/se-check.R`, `R/fit.R`, `R/confint.R`, `R/methods-fit.R`,
  `R/brms-shapes.R`, `R/multiple.R`, `R/me.R`, `R/interop.R`,
  `R/predict.R`, `R/conditional-effects.R`; `man/frmtmb_control.Rd`
  regenerated.
- Core docs: `NEWS.md`, `vignettes/diagnostics.Rmd`,
  `vignettes/inputs.Rmd`.
- Core tests: `test-se-check.R` (four new tests, one rewritten),
  `test-case-studies.R` and `test-diagnostics-ux.R` (round-1
  workarounds removed).
- frmtmb.spline: `R/curve-cov.R`, `R/curve.R`, `R/curve-deriv.R`,
  `R/curve-feature.R`, `tests/testthat/test-curve.R`, `NEWS.md`.
- `dev/nanse-suite.sh` and `dev/nanse-diag.sh` (the two extra gates),
  new `dev/nanse-p2-*.R`, this file.

### Proposed change to lane fixes' `nl_flat_message()` (updated, not made)

Replaces the round-1 proposal (RB3). Flatness is decided by
`nl_flat_tol * big` alone; the noise term only enters `tau`. The
round-1 version put `10 * ||E||` into the flat threshold too, and on
`test-nl-rtmb-scope.R`'s identified `a * squash(b * x)` (smallest exact
unit-diagonal eigenvalue 9.0e-6, `dev/nanse-rev2-squash.R`) that
function's own finite-difference noise crossed it. The reviewer
validated this variant (`dev/nanse-rev2-log/*-variant.txt`): squash
silent, the 5 ridges of the fixes table warn, 16 identified models
silent, crossed design 120 of 120 at k = 60 and 40 of 40 at k = 20 on 4
seeds.

```r
  Hr <- H                       # the loop's matrix, before (H + t(H)) / 2
  H <- (Hr + t(Hr)) / 2
  ...
  S <- H[ok, ok, drop = FALSE] / outer(dg[ok], dg[ok])
  E <- abs(Hr - t(Hr))[ok, ok, drop = FALSE] / 2 / outer(dg[ok], dg[ok])
  ev <- eigen(S, symmetric = TRUE)
  big <- max(abs(ev$values))
  flat <- abs(ev$values) < nl_flat_tol * big
  if (!any(flat) || all(flat)) return(NULL)
  gap <- min(abs(ev$values[!flat]))
  tau <- max(sqrt(nl_flat_tol), 10 * sqrt(sum(E^2)) / gap)
  proj <- sqrt(rowSums(ev$vectors[, flat, drop = FALSE]^2))
  load <- labels[ok][proj > tau]
```

On the merge, in `fit_assembled()` after the `frm_warning()` of a
non-NULL `flat_msg`, set `fit$cache$se_explained <- "nl_flat"` (the
family form: it explains every lost parameter except a bound-held
one).

## Punch round 1

Review `dev/reviews/2026-10-06-nanse.md` (NOT MERGEABLE: B1, B2,
m1 to m10). Every number below is from the current lane build in
`wt-nanse-lib` unless the line names another build. "Seen to fail":
`test-se-check.R` on the round-0 build (`nanse-rev-lib`, the
reviewer's copy of round 0) gave pass 50, fail 24, error 1
(`dev/nanse-log/p1-test-se-check-round0.txt`), on base pass 24, fail
48, error 2 (`p1-test-se-check-base.txt`); the hypothesis/emmeans
test added last failed 6 expectations on the build before its fix
(`p1b-se-check-prefix.txt`).

### B1: the flat-direction rule

The 0.1 loading threshold is gone. Tier 3 of `cov_from_hessian()`
(`se_tier3()`) now works on the unit-diagonal free block and its
noise `E`, the asymmetry of the unsymmetrized Hessian (`fd_hessian()`
for a finite-difference Hessian, the AD Hessian's own asymmetry when it
is exact):

- A direction is flat when its eigenvalue is at most
  `flat_thr = max(se_flat_tol * largest, 10 * ||E||_F)`.
- A parameter loses its SE when its projection onto the flat subspace,
  `sqrt(sum_k v_ik^2)` over the flat eigenvectors, exceeds
  `tau = max(sqrt(se_flat_tol), 10 * ||E||_F / gap)`, where `gap` is
  the smallest kept eigenvalue. `sqrt(1e-9) = 3.2e-5`, so a direction
  shared by up to about 1e9 parameters still reaches every one of
  them; the second term is the Davis-Kahan bound on how far noise can
  turn an eigenvector, so noise cannot name a parameter either.
- An empty Hessian row (largest entry at most 10 times its noise, and
  on a finite-difference Hessian at most `se_empty_row = 1e-6`) is
  removed first. The 1e-6 floor no longer applies to an exact AD
  Hessian, whose rows are trusted to rounding.

Results (`dev/nanse-rev-cases.R lane`, `dev/nanse-rev-spread2.R lane`,
`dev/nanse-rev-grby.R lane`, logs `dev/nanse-log/p1b-*.txt`):

    spread y ~ a + b, a ~ 0 + f   k = 10, 40, 60, 120: every a_k and b
                                  lost (0 finite); round 0 kept 60 and
                                  120 a_k at SE 0.197 and 0.214
    crossed a ~ 0 + f, b ~ 0 + g  k = 20: 40 of 40 lost; k = 60: 120 of
                                  120 lost (round 0: 97 of 120 finite)
    gr-by fixture                 x SE 0.1028 (base 0.1023, lme4
                                  0.1009); theta_2, theta_3 flat
    y in 1e5 units, a + b ridge   a, b lost; c keeps SE 9249

Tests: "a flat direction shared by many coefficients takes all of
them" (k = 60 spread and k = 20 crossed; 5 failures on round 0) and
"a small loading on a downward direction keeps its SE" (gr-by; 3
failures on round 0).

`mo()` study (`dev/nanse-mo-lane.R`, seeds 1 to 200 in 4 processes,
`dev/nanse-log/mo-lane-p1b.tsv`, summary
`dev/nanse-mo-lane-sum.R dev/nanse-log/mo-base.tsv
dev/nanse-log/mo-lane-p1b.tsv` in `mo-lane-p1b-sum.txt`): unchanged.
26 SE warnings plus 17 code-7 fits on the 43 interaction seeds that
lose an SE, 0 silent, 0 of 157 interaction and 0 of 200 main-effect
fits with every SE finite warned (0 of 357). Every-SE-finite seeds
match the held-coordinate reference to 2.87e-5 (median 4.5e-7). 47 of
base's 72 all-NaN seeds get every SE back, 8 lose some, 17 are the
code-7 fits. The lost sets moved on 13 seeds, each to a subset of its
round-0 set: "concave" labels fell from 12 seeds to 2 (m1), seed 50
(6 parameters "concave" in round 0, the coefficients included) now
loses `moincome:age`, `zeta1_1` and `zeta1_2` as flat, and 7 seeds that
lost two simplex coordinates now lose one.

### B2: predictions along a lost direction

`sdr_rescue()` keeps the null basis of the removed directions
(`se_null`, optimizer units, with the unit vectors of removed rows).
`jc_nonest()` is predict()'s `alias_null` test for it: a prediction
whose gradient has a cosine above `se_pred_tol = 1e-6` with that basis
gets a NaN standard error, and the call gives one warning, "k of n
predictions move along a direction the fit does not determine",
class `frmtmb_se_lost_prediction`. Covered: `frm_linpred()`,
`predict()` (mean and ordinal probabilities), `fitted()`,
`conditional_effects()` (both routes), conditional smooths,
`emmeans()` (grid route by the cosine test; design route by the NaN
rows of `vcov()`, with the same warning), and `hypothesis()`, which
now reads the covariance without the lost directions and the same
test, so `a_f1 + b_Intercept`, which the data determine, gets the
ML SE (equal to sigma_ML / sqrt(5) to 1e-5) instead of NaN.
`frm_joint_cov()` shows lost rows as NaN.

    dev/nanse-rev-pred.R lane, k = 10 and 60:
      fitted(dpar = "a"), "b"   Est.Error NaN NaN, 1 warning each
      (round 0: 0.1818 and 0.0299, no warning)
    dev/nanse-p1-emm.R: emmeans(~ f) mu SE 0.189, no warning;
      emmeans(~ f, dpar = "a") NaN; hypothesis(a_f1) NaN
    mo() seed 71 conditional_effects("income:age"): every band
      finite, no warning (the gradient along a saturated coordinate
      is exactly 0)

Tests: "a prediction along a lost direction gets NaN and one warning"
(6 failures on round 0) and "hypothesis() and emmeans() treat a lost
direction the same way" (6 failures before its fix).

Not covered: `frm_curve()` in frmtmb.spline reads `lb$V` from
`frm_lp_basis()`, which now carries `se_nonest`, but frm_curve()
does not read it; the change belongs in frmtmb.spline. It reaches
only fits without random effects that lost an SE. `fitted()` on the
`mu` of that nonlinear model gives `Est.Error` NA on base and lane
alike (not this lane's).

### m1: concave is confirmed

A direction with eigenvalue below `-flat_thr` is probed:
`se_line_probe()` steps `sqrt(4 * grad_tol / |lambda|)` both ways in
the optimizer's units and calls it concave only when the objective
improves by more than `grad_tol`. Otherwise its dominant parameters
(loading at least half the largest) are flat. Probe and state are
restored (`obj_state_save()`). `dev/nanse-rev-concave.R lane`: the
spline `s(x, by = fac) + gp(x)` fixture (eigenvalue -0.308, the
objective moves at most 8e-9 along it) and the gr-by fixture (real
ascent worth 1.4e-4, below `grad_tol`) are both flat now. On the
`mo()` study, 2 seeds keep a confirmed concave `zeta2_2`.

### m2: cost

The check reads the exact AD Hessian on models without random effects
(and only computes the finite-difference one if that shows a
problem). With random effects it builds the finite-difference Hessian,
step for step `optimHess()` (identical to sdreport()'s, verified with
`identical()` on sleepstudy), only while the projected time is at most
30 percent of the fit or 0.05 s, measured after the first gradients;
past that it stops and the check is deferred to the first SE use
(`se_deferred_report()`, one warning, cached). The cached Hessian is
reused by `summary()`.

`dev/nanse-rev-cost.R 3` (the reviewer's script; arms interleaved,
blocks past 1.2 s, minimum of 3 rounds, control = check against
check) and `dev/nanse-cost2.R`, on an idle machine, final build
(`dev/nanse-log/p1d-revcost.txt`, `p1d-cost2.txt`; the earlier run of
this round in `p1-revcost.txt`, `p1-cost2.txt`):

    model           pars  round 0  this round      control  fit+summary
    glm_f300         301   1.908   1.189 (1.236)  0.980    1.074
    glmm_f100        102   1.505   1.033 (1.023)  1.005    1.031
    lmm_f400         403   1.716   1.030 (0.980)  1.005    1.015
    lm                      1.044   1.090 (1.087)  0.999
    lmm_slope               1.125   1.096 (1.114)  1.000
    pois_glmm               1.168   1.148 (1.144)  0.993
    pois_glmm_2000          1.110   1.109 (1.145)  1.008
    dist_sigma              1.061   1.099 (1.081)  1.020
    nl                      1.126   1.119 (1.129)  0.965
    negbin                  1.234   1.217 (1.226)  1.024

The two mixed models are deferred (their Hessian build is projected
past 30 percent of the fit), so their fit-time cost is the gradient
probe that measures that. The glm pays for the exact Hessian (301
reverse sweeps); that is the dependence on the parameter count, and it
is below the 1.3 the punch list set.

Kept on by default: the 400-parameter mixed model costs nothing at fit
time, because the check is deferred. The 301-coefficient glm pays
1.19 to 1.24, below the 1.3 the punch list set. The dependence on the
parameter count is stated in NEWS, `?frmtmb_control` (`check_se`) and
the `se_check()` comment. The dense Hessian is cached only when the
finite-difference build ran or a repair was needed, so none of the
three large models keeps one ("cached Hessian 0 Mb" in
`p1d-revcost.txt`; round 0 kept 2.1 MB at 301 and 3.8 MB at 403).

### m3: precedence per parameter

`check_re_structure()` returns what it warned about (the theta of a
one-level factor; the theta and the `sigma` coefficient of an OLRE),
and only when its action was "warning", so `"ignore"` explains
nothing. A family's `fit_check` hook that warns explains every lost
parameter except bound-held ones. A convergence warning explains the
whole fit. `vcov()` on a fit with `se_lost` uses the SE check's clauses,
so it cannot give the old "overparameterized" reason, and it stays
silent once the fit warned. `dev/nanse-rev-prec.R lane`: one-level
`g1` plus `x` held by `ub = 0.1` gives the single-level warning and
"Standard errors are not available for 1 of 4 parameters. x: a bound
holds it"; OLRE with `check_olre = "ignore"` names
`sigma_(Intercept), theta_1`. Test "an earlier warning explains only
the parameters it names" (6 failures on round 0). Two existing tests
expected silence under `check_olre = "ignore"` and one expected the
deferred check at fit time; they now allow the SE warning
(`test-diagnostics-ux.R`, `test-case-studies.R`, which calls `vcov()`
inside the block because a loaded machine can defer the check).

### m4: `mo()` fits below the maximum

Recorded: of the 17 interaction seeds more than 0.5 log-likelihood
units below the profile maximum, 12 are silent with every SE finite on
the lane (6 of those 12 were all-NaN on base). Their SEs belong to a
point 0.5 to 1.9 units short of the maximum. A cheap check was tried
(`dev/nanse-p1-m4.R`): a 1 percent step toward each simplex vertex
(the first round of `dev/nanse-mo-escape.R`). It catches 2 of the 17
and fires on 20 seeds in all, most of them within 1e-3 of the maximum.
Not added.

### m5: REML wording

Under REML or `profile = TRUE` the warning says "The other outer
parameters keep theirs; the coefficients' standard errors come from
the joint precision, which this does not repair". Test "REML does not
promise the coefficients' standard errors" (3 failures on round 0).
`dev/nanse-rev-reml.R lane`: 60 ML and REML fits with a variance at
0, none lost an SE or warned.

### m7: separation

Flat mean coefficients of a binomial or ordinal family at |estimate|
above 10 are labeled "the data separate the outcomes, so the estimates
run off toward infinity and the optimizer stopped where its tolerances
did" (`se_relabel()`, not under REML). `dev/nanse-rev-sep.R lane`:
predfix seed 514, the code-0 fit, is "3 of 3 parameters ... the data
separate the outcomes"; the other separated fits stop at code 1 or 9
and get the convergence warning. Test "separated data are named as
separation" (1 failure on round 0).

### m6, m8, m9, m10

m6 (the "converged elsewhere" split) is a suggestion; not done. m8
and m9 are for the merge: `extensions/frmtmb.sample/R/sample.R:69`
says "thirteen" names, and frmtmb.learn's floor must become the
release that carries this warning (its `test-engine.R` requires it).
m10 is lane fixes'.

### Suites

`dev/nanse-suite.sh dev/nanse-suite7 18` (all eight suites, one file
per process, every gate on, 18 at a time) on the build before the
hypothesis()/emmeans() change: 340 of 340 files, pass 23224, fail 0,
error 0, warning 0, skip 29. After that change,
`dev/nanse-diag.sh` reran the 70 files that call `hypothesis()` or
`emmeans()` or where the SE warning fired (`dev/nanse-p1c-log`): all
pass with no escaped warning except `test-portability.R`, whose
hypotheses on a flat `sd` block now warn by design; that test now
allows the warning (rerun clean, `dev/nanse-log/p1d-portability.txt`).

Firings (`dev/nanse-fire-sum.R dev/nanse-p1c-log`,
`dev/nanse-log/p1c-fire-sum.txt`): 57, of which 16 are
`test-se-check.R`'s own. All 57 are at optimizer code 0 on fits where
base's sdreport() has a non-finite SE. By reason: bound 10, flat 47,
concave 0 (round 0: 5 concave, all finite-difference noise or an
ascent below `grad_tol`, m1). Outside test-se-check.R the count is 41
in this run and 43 in the full suite run: the check defers on a loaded
machine, and a test that never reads an SE then never sees the
warning (test-id-kron.R, one brms-suite-methods fit, one
stan-identity fit).

`sh dev/nanse-check.sh frmtmb` on the final tree (installed and
roxygenised first): Status: 1 NOTE, V8 unavailable for the HTML
manual, as before. No extension's R code changed, so no extension was
rechecked.

### Proposed change to lane fixes' `nl_flat_message()` (not made)

It names a coefficient when `abs(v) > 0.1` on a flat eigenvector, the
rule B1 removed here, so on the crossed k = 60 design it names 23 of
120 and the SE check, explained by it, stays silent. Proposed: keep
the raw differenced matrix before it is symmetrized and use the same
projection rule.

```r
  Hr <- H                       # the loop's matrix, before (H + t(H)) / 2
  H <- (Hr + t(Hr)) / 2
  ...
  S <- H[ok, ok, drop = FALSE] / outer(dg[ok], dg[ok])
  E <- abs(Hr - t(Hr))[ok, ok, drop = FALSE] / 2 / outer(dg[ok], dg[ok])
  ev <- eigen(S, symmetric = TRUE)
  big <- max(abs(ev$values))
  flat <- abs(ev$values) < max(nl_flat_tol * big, 10 * sqrt(sum(E^2)))
  if (!any(flat) || all(flat)) return(NULL)
  gap <- min(abs(ev$values[!flat]))
  tau <- max(sqrt(nl_flat_tol), 10 * sqrt(sum(E^2)) / gap)
  proj <- sqrt(rowSums(ev$vectors[, flat, drop = FALSE]^2))
  load <- labels[ok][proj > tau]
```

On the merge, when it fires, set `fit$cache$se_explained <- "nl_flat"`
(the family form, which explains every lost parameter except a
bound-held one), not the whole-fit mark, so a bound-held coefficient
in the same fit is still named by the SE check.

## Summary

1. **Defect 8, cause.** `ls ~ mo(income) * age` loses every standard
   error because the interaction's simplex has a weight at 0. The
   softmax puts that weight at a coordinate of minus infinity. The
   optimizer stops at finite softmax coordinates (magnitudes of 15 to
   1927 on the seeds whose plain inverse fails) where the likelihood is
   flat along them to 1e-20 or less, so the outer Hessian has decoupled
   near-zero rows. sdreport() inverts it with solve(), which
   refuses when the reciprocal condition number is below machine
   epsilon and then returns NaN for every parameter. The rest of the
   Hessian is well conditioned. So the matrix is genuinely singular in
   the saturated coordinates, and the loss of the OTHER standard errors
   is numerical.
2. **Fix.** Where sdreport()'s covariance is not finite, frmtmb now
   inverts the same Hessian scaled to unit diagonal (with the exact AD
   Hessian on a model without random effects). Where that fails too,
   it drops only the directions the likelihood does not determine and
   names those parameters. `frm()` now builds the Hessian at fit time
   and warns, naming each parameter without a standard error and why.
   `summary()` prints the reason instead of a bare NaN.
3. **Defect 9.** `frm_allfit()` now starts every optimizer where lme4's
   `allFit()` does (the original fit's estimates by default; lme4
   2.0.6's `start_from_mle = TRUE`), or at the original fit's own start.
   The table marks a success code below the best log-likelihood as
   "converged elsewhere".
4. **Found, not fixed: the `mo()` point estimates are not always the
   maximum.** On 52 of 200 data sets frmtmb's log-likelihood is more
   than 0.01 below the exact profile maximum, up to 1.95 below. On 20
   of the 200 the fit sits on a softmax plateau: a step toward a
   simplex vertex raises the likelihood, by up to 0.70. After that
   escape 38 are still more than 0.01 below, at another local maximum.
   See "Defect found, not fixed".

## Item 1: defect 8

### The cause, measured

`dev/nanse-mo-sweep.R` (seeds 1 to 200 of brms_monotonic's data code,
rellib-r5), summarized by `dev/nanse-mo-sweep-sum.R`:

    form int (ls ~ mo(income) * age): 200 seeds
    all SEs non-finite (sdreport): 72; some: 18
      of the all-NaN: code 0 55, no warning at fit 55, both 55
      pdHess TRUE on all-NaN: 59
    seeds with a simplex component < 1e-8: 110; of all-NaN: 72
    finite-SE seeds with a saturated component: 38
    scaled inverse: all finite on 59 of the all-NaN seeds
    rcond(H) on all-NaN seeds: max 1.54e-16 ; on finite seeds: min 2.68e-16
    max |gradient| on all-NaN seeds: 0.00517
    form main (ls ~ mo(income)): 0 of 200 with any non-finite SE

- Every all-NaN seed has a simplex weight below 1e-8; 38 seeds with a
  weight that small keep finite SEs. What separates the two groups is
  the reciprocal condition number of the outer Hessian, exactly at
  solve()'s threshold (2.2e-16): at most 1.54e-16 on every all-NaN
  seed, at least 2.68e-16 on every finite one. That is a knife edge,
  not a property of the model.
- The null direction is one coordinate. `dev/nanse-mo-probe.R`, seed 7:
  smallest eigenvalue of the outer Hessian 1.76e-20, eigenvector
  `zeta2_1 = -1.000`, `zeta2 = (-46.04, -17.06)`, so the interaction's
  simplex is (1, 1.0e-20, 3.9e-8). Seed 8: three near-zero eigenvalues
  (2.1e-8, 1.7e-8, 5.3e-15) on `zeta1_2` and both `zeta2`.
  `optimHess()` and the exact `obj$he()` agree to three digits on each.
- The model is identified at these points in the sense that matters for
  the coefficients: holding the saturated coordinates fixed and
  inverting the exact Hessian over the rest gives finite SEs for every
  coefficient and sigma (the reference of `dev/nanse-mo-lane.R`).
- On some seeds the saturated coordinate curves slightly downward
  (seed 12: diagonal -4.1e-4, gradient -4.2e-4, `dev/nanse-mo-hess.R`).
  There the fit is not a maximum along that coordinate; the plain
  inverse then gives a negative variance (base's 18 "some" seeds).
- Why `mo(income)` alone never fails: its simplex is identified through
  a large coefficient. The interaction's coefficient is near 0 (the
  data have no interaction), so its simplex is barely identified and
  runs to the boundary.

### The change

`R/se-check.R` (new), with `autoscale_sdreport()` calling
`sdr_rescue()`:

1. sdreport()'s own inverse is kept whenever every variance is a
   positive number, bit for bit (the guard test "a healthy fit is
   untouched" asserts `identical()` with `RTMB::sdreport()`).
2. Otherwise, the same Hessian scaled to unit diagonal is inverted (the
   same matrix, so it only differs where solve() refused). On a model
   without random effects the exact AD Hessian replaces optimHess()'s
   first. This tier requires the scaled matrix to be positive definite
   with its smallest eigenvalue above 1e-9 of the largest.
3. Otherwise, the parameters whose standard error does not exist are
   named, each with a reason (a bound holds it; its Hessian row is not
   finite; its row is empty or it loads on a direction of eigenvalue
   below 1e-9 of the largest, "flat"; it loads on a direction below
   -1e-3 of the largest, "concave"), and the covariance of the others
   is the pseudo-inverse over the determined directions. `vcov()`,
   `confint()` and `summary()` read it with the lost rows NaN.
   Predictions read the finite pseudo-inverse, which gives the right
   variance to any prediction that does not move along the lost
   directions.
4. A fit that did not converge (optimizer code not 0, or a convergence
   warning) keeps sdreport()'s own covariance, unrepaired.

The fit-time check (`se_check()`) builds the outer Hessian exactly as
`autoscale_sdreport()` would, at the same point, with the objective's
state restored afterwards, and caches it, so a later `summary()` does
not build it again. `frmtmb_control(check_se = "ignore")` skips it.

### Numbers on the 200 seeds

`dev/nanse-mo-lane.R` on rellib-r5 and on the lane build, summarized by
`dev/nanse-mo-lane-sum.R dev/nanse-log/mo-base.tsv
dev/nanse-log/mo-lane-final.tsv`:

    == int: 200 seeds
    base: all SEs non-finite 72, some 18, all finite 110; SE warning at fit 0;
          any fit warning 17; conditional_effects print fails 72
    lane: all SEs non-finite 17, some 26, all finite 157; SE warning at fit 26;
          any fit warning 43; conditional_effects print fails 17
    lane: SE warning on seeds with a lost SE: 26 of 43; on seeds with every
          SE finite: 0 of 157
    lane: lost seeds whose fit already warned otherwise: 17
    lane: lost seeds with no warning at all: 0
    lane: max |SE/ref - 1| over the coefficients and sigma, seeds with every
          SE finite: 2.87e-05 (median 4.51e-07)
    lane: seeds whose base SEs were all finite and now differ: 0
    base->lane:   lane allNaN finite some
         allNaN        17     47    8
         finite         0    110    0
         some           0      0   18
    == main: 200 seeds: base and lane all finite on 200, no warning on 200

- Coverage: every fit with a lost standard error warns, 26 with the new
  warning and 17 with "Optimizer did not report convergence: singular
  convergence (7)". 0 of 43 are silent (base: 55 of 72 all-NaN fits
  were silent).
- False alarms: the new warning fires on 0 of the 357 fits whose
  standard errors are all finite (157 interaction, 200 main effect).
- The 17 that stay all-NaN are exactly base's 17 all-NaN seeds with
  optimizer code 7; precedence leaves them to the convergence warning.
- The recovered standard errors equal the reference that holds the
  saturated coordinates fixed to 2.87e-5 relative.
- Lost sets on the 26 (the table in the log): one or two simplex
  coordinates on 25 seeds, flat on 14 and concave on 11. On seed 50 the
  concave direction loads on the four coefficients too (eigenvalue
  -1.43 of the unit-diagonal Hessian), and the warning names them.

### Against brms 2.23.0

`dev/nanse-mo-brms.R` on the final build (brms default priors,
Dirichlet(1) simplexes, 4 chains of 2000, Stan programs from
`dev/stan-cache`; max Rhat 1.005; `dev/nanse-log/mo-brms-final.txt`).
frmtmb SE over brms posterior sd:

    seed  Intercept   age  moincome  moincome:age  sigma   base SEs
       1      0.869 0.850     0.905         0.860  0.909   finite
       7      1.112 1.128     0.865         0.838  0.919   all NaN
       8      0.534 0.501     0.224         0.205  0.923   all NaN
      11      0.718 0.693     0.922         0.900  0.906   all NaN
      13      1.061 1.077     0.909         0.916  0.915   all NaN
      71      1.192 1.218     1.042         1.053  0.917   all NaN (1)

(1) on the lane seed 71 loses zeta2_2 only.

Seed 1 (finite on base) sets the ordinary ML-to-posterior ratio, 0.85
to 0.91. Seeds 7, 13 and 71 are in that range or above it, seed 11 is
0.69 to 0.92 (below it for the intercept and age), and seed 8 is
0.2 to 0.5: frmtmb's interaction simplex is a vertex, (0, 0, 1), and its
SE is conditional on that vertex, while brms's posterior averages over
simplexes with sd 0.21 to 0.24 per weight. That gap is the difference
between a boundary ML estimate and a posterior, not a defect of the
repair.

## Item 2: defect 9

- lme4 2.0.6 (`lme4::allFit`, read in the installed package):
  `start_from_mle = TRUE` is the default and refits with
  `update(object, start = pars)` from the MLE; `FALSE` refits with
  `update(object)`, which reuses the original call, its `start` too.
  lme4's NEWS: "allFit(..., start_from_mle = TRUE) now specifies that
  models should be refitted using starting parameter values from the
  original model fit, which should increase efficiency". The brief's
  "lme4 starts every optimizer from the same start values as the
  original fit" describes `FALSE`; frmtmb now offers both with lme4's
  default.
- Why the MLE as default: allFit exists to separate optimizer trouble
  from the model at the reported optimum. From the MLE every optimizer
  either stays (agreement) or finds a better point nearby (the original
  stopped short). From the original start a nonlinear model's refits
  can wander far from the reported optimum, and a disagreement then
  mixes each optimizer's path with the question asked. `FALSE` stays
  available for the stronger test.
- The ported brms-suite ledger: no gated brms-suite file failed in the
  suites below, so this lane flips no row.
- `dev/nanse-allfit.R`, `fit_loss` (brms_nonlinear, `start = list(beta
  = c(5000, 1, 45))`, logLik -361.7056667):
  - rellib-r5: nlminb, optim, nloptr_lbfgs NULL; bobyqa -517.0224 code
    0 (`dev/nanse-log/allfit-base.txt`).
  - lane, default: all four -361.7057, spread 8.75e-12, "agrees".
  - lane, `start_from_mle = FALSE`: nlminb and nloptr_lbfgs agree; optim
    -362.5268 code 0 is "converged elsewhere" (0.821 below); bobyqa
    -362.0274 code 1 "did not converge".
  - controls: `b1 * exp(b2 * x)` and the epilepsy Poisson GLMM, all
    four agree under both settings (spreads 1.1e-13 to 1.7e-8).
- The table now has `vs best` (each refit's logLik minus the best any
  fit reached, the original included) and a verdict. An optimizer
  started at the optimum can fail its own line search there (optim code
  52 on the nonlinear test fixture) while standing on the best point;
  that row reads "agrees; code 52".

## Item 3

`dev/nanse-item3.R`, rellib-r5 against the lane
(`dev/nanse-log/item3-base.txt`, `item3-lane.txt`):

Each case: base SEs and warning, then lane SEs and warning.

- `a * exp(b * x)`, x 0 to 1e5: NaN NaN NaN, none; 0.00963 9.68e-08
  0.0408, none.
- control, x / 1e5: 0.00963 0.00968 0.0408, none; the same.
- `b` at its bound `lb = 0`: NaN NaN NaN, only R's "NaNs produced";
  0.0274 NaN 0.05, "b_(Intercept): a bound holds it".
- `a + b + log(c0)`, c0 = 5.9e-5: 6 NaN, none; 5 NaN and sigma 0.05,
  names a_*, b_* and c0 (flat).
- the same ridge with c0 near 1: 6 NaN, none; 5 NaN and sigma 0.05,
  names a_*, b_* and c0 (flat).
- `a + log(c0) * x`, c0 = 5e-5 (identified): NaN NaN NaN, none; 0.0215
  9.97e-07 0.05, none.
- `a + log(c0)` with `a ~ 1 + x` (a curved ridge): 4 NaN, none; NaN
  0.0192 NaN 0.05, names a_(Intercept) and c0 (flat).

- `a * exp(b * x)`: the cause is optimHess()'s absolute step of 1e-3 in
  `b = -2e-5`, which moves `b * x` by 100. Its `b` diagonal is 1.31e96
  against the exact 2.45e14; the exact Hessian's unit-diagonal
  eigenvalues are 1.75, 1 and 0.248, so the model is well identified.
  The repaired SE of `b` equals the rescaled control's (0.00968 / 1e5).
- `a + b + log(c0)`: optimHess() steps into `c0 < 0` and its `c0` row
  is NaN; the exact Hessian is finite (diagonal 6.86e11). The ridge is
  real (unit-diagonal eigenvalues 2.5e-14, 1.0e-16, -6.7e-16), so the
  five ridge parameters are named and sigma keeps its SE.

### dev/test-backlog.md entries

- **"The covariance machinery is not bound-aware"** (Open, minor):
  closed. The `ub = 0.1` construction (gradcheck-01, seed 101) now
  names `x` ("a bound holds it") and the free parameters keep the
  inverse of the free block of the Hessian (test "a bound-held
  parameter loses its SE and the others keep theirs").
- **"`y ~ mo(m) + m` is not identified and is fitted in silence"**
  (Open, medium): the ordered-factor case, whose 6 of 6 SEs were
  non-finite with only R's "NaNs produced", now warns at fit time. The
  unordered case, whose SEs "happen to stay finite", is not covered:
  this lane acts on non-finite SEs only.
- **skew_normal "diagnose() says nothing about a dpar whose information
  is singular"**: not covered (the SE there is finite, 5.09).

## Precedence with the warnings of lanes ordmix and fixes

One rule: the standard-error warning is the LAST verdict and stays
silent on any fit that an earlier, more specific verdict has already
explained. It reads two flags (`se_explained()` in R/se-check.R):

In order of precedence, the fit, what warns, and the SE warning:

1. Optimizer code not 0: check_convergence(), "Optimizer did not
   report convergence". SE warning silent, and the covariance is left
   as sdreport() gives it.
2. Large gradient (grad_verdict()): check_convergence(). Silent.
3. Non-finite gradient (lane ordmix): its new branch in
   check_convergence(). Silent, with no merge edit, because it is one
   of check_convergence()'s messages.
4. A grouping factor with one level, or a gaussian observation-level
   effect: check_re_structure() before the fit, also under "ignore".
   Silent.
5. A family's fit-end check warns: thres() categories no row takes,
   an ordinal `disc` intercept, zero_one_inflated_beta's `coi`, xbeta's
   `kappa`, and lane ordmix's degenerate mixture component (raised in
   `mixture_ord_fit_check()`). Silent, with no merge edit, because
   fit_end_checks() records any warning a family's `post$fit_check`
   raises.
6. A flat nonlinear direction (lane fixes' nl_flat_message()): fixes'
   warning. Silent after ONE line on the merged tree (below).
7. `se = TRUE`, a Hessian not positive definite, every SE finite:
   check_convergence(), "Hessian is not positive definite", unchanged.
   Not raised, because no SE is lost.
8. None of the above, some SE not finite: the SE warning, naming the
   parameters.

Integration onto the merged tree:

- **fixes**: in `fit_assembled()`, where its block reads
  `if (!is.null(flat_msg)) frm_warning(flat_msg, call. = FALSE)`, add
  `fit$cache$se_explained <- "nl_flat"` inside that `if`. Both run only
  on the user's own call. Its `nl_flat_message()` returns NULL when its
  Hessian is not finite (the `c0 = 5.9e-5` ridge), and the SE warning
  then covers that fit (item 3). Test file overlap:
  `extensions/frmtmb.ode/tests/testthat/test-ode-nlf.R`. Take fixes'
  version; with the line above the SE warning does not fire there.
- **ordmix**: nothing to add. Its non-finite-gradient warning is a
  check_convergence() message, and its degenerate-component warning is
  raised inside a family `fit_check` hook. Textual overlap in
  `check_convergence()`: keep its new `else if` branch after
  `grad_verdict`, then this lane's se = TRUE pdHess block.

Functions edited in R/fit.R: `frm()` (one line: the frame attribute
`se_explained` from `check_re_structure()`; and the `se` and RE3.0
documentation), `fit_assembled()` (new argument `check_se`, the
autoscale-choose branch runs the check once on the kept fit, `start`
stored on the fit, the convergence flag, the `se_check()` call),
`frmtmb_control()` (`check_se`), `check_re_structure()` (returns what
it found), `check_convergence()` (the non-finite branch moved to
se_check(); the pdHess branch now only when no SE is lost). Elsewhere:
`autoscale_sdreport()` (R/autoscale.R), `fit_end_checks()`
(R/fit-end.R), `get_joint_cov()` (R/predict.R), `vcov_estimated()`,
`summary.frmtmb_fit()`, `print.summary.frmtmb_fit()`
(R/methods-fit.R), `diagnose()` (R/confint.R, the smallest covariance
eigenvalue and its doc), `frm_allfit()`, `print.frmtmb_allfit()`
(R/allfit.R); new R/se-check.R.

## Tests

New: `tests/testthat/test-se-check.R` (8 tests),
`tests/testthat/test-allfit-start.R` (4 tests). Seen to fail on
rellib-r5 (`dev/nanse-testfile.R`, logs
`dev/nanse-log/test-se-check-base.txt`,
`test-allfit-start-base.txt`):

    lib: C:/Users/adf44/source/r/rellib-r5/frmtmb
    RESULT test-se-check.R pass=12 fail=29 err=1 skip=0 warn=1
      a mo() simplex weight at 0 no longer takes every SE    fail=4
      a lost simplex SE warns at fit time and keeps the rest fail=10
      a covariate spanning 1e5 in a nonlinear body ...       fail=2
      an unidentified ridge names its parameters ...         fail=8
      a bound-held parameter loses its SE ...                fail=4
      check_se = 'ignore' and 'stop' do what they say        err (no check_se)
      se = TRUE says it once, and a fit that stopped short   fail=1
    lib: C:/Users/adf44/source/r/rellib-r5/frmtmb
    RESULT test-allfit-start.R pass=2 fail=4 err=3 skip=0 warn=0
      a nonlinear fit's refits start at its estimates        fail=1, err
      start_from_mle = FALSE starts where the original ...   err (no argument)
      a success code at a worse optimum is not ...           fail=3
      a refit that errors keeps its message                  err (no $errors)

Every failure but the two missing-argument errors is behavioral: NaN
standard errors, no warning, NULL refits, no verdict in the table. The
one test that passes on base by design is the guard "a healthy fit is
untouched and never warns" (the absent case).

Changed existing tests, each because the new warning fires on a fit
whose SEs were already not finite on base (every one checked, see
"Firings"): core `test-brms-names.R`, `test-brms-parity-defects.R`,
`test-case-studies.R`, `test-id-kron.R`, `test-me.R`, `test-v15.R`;
frmtmb.eam `test-defects.R`; frmtmb.latent `test-hmm.R`; frmtmb.learn
`test-engine.R`, `test-reference.R`, `test-stan-identity.R`;
frmtmb.ode `test-ode-nlf.R`; frmtmb.sample `test-gr-by-draws.R`.
`test-diagnostics-ux.R`: three tests read the old vcov()-time warning;
they now read the fit-time one, and the flat-`mo()` test moved to a
fixture that still has an exactly flat row (seed 71 of the vignette's
data code), because its old fixture (seed 6 of `y ~ mo(mo) + x`) now
gets every SE, which the test asserts. `helper-fuzz.R`: the fuzz tier
counts the new warning as "the fit warned".

## The suites

All eight, one file per process, every gate of the ordinary tier on
(`NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`), `bash dev/nanse-suite.sh
dev/nanse-suite5 26` on the final build, summarized by
`dev/nanse-suite-sum.R dev/nanse-suite5`:

    logs: 340  with a RESULT line: 340  loading the lane core: 340

    frmtmb           files 204 pass 16426 fail 0 err 0 skip 14 warn 0
    frmtmb.coupling  files  11 pass   542 fail 0 err 0 skip 5 warn 0
    frmtmb.eam       files  29 pass  1743 fail 0 err 0 skip 3 warn 0
    frmtmb.latent    files  10 pass   360 fail 0 err 0 skip 2 warn 0
    frmtmb.learn     files  15 pass   501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode       files  11 pass   549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample    files  45 pass  2515 fail 0 err 0 skip 1 warn 0
    frmtmb.spline    files  15 pass   553 fail 0 err 0 skip 1 warn 0

The skips are the three gated files of other tiers
(`test-drmtmb-agreement.R` 13, `test-fuzz.R` 1, each package's
`test-scale.R`). Those ran separately with their gates
(`FRMTMB_DRMTMB_FIT_TESTS`, `FRMTMB_FUZZ`, `FRMTMB_SCALE_TESTS`;
`dev/nanse-gated-log/`): `test-drmtmb-agreement.R` pass 131,
`test-fuzz.R` pass 2, core `test-scale-contract.R` pass 37, and
`test-scale.R` of coupling 6, eam 18, latent 6, learn 14, sample 3,
spline 2, all with fail 0 err 0 skip 0 warn 0. frmtmb.ode's
`test-scale.R` was stopped there after 40 minutes, because the lane
library was reinstalled under it, and rerun alone against both builds
at once (`dev/nanse-log/ode-scale-{base,lane}.txt`): pass 2 fail 0
err 0 skip 0 warn 0 on both, 34 min 37 s on rellib-r5 and 35 min 1 s on
the lane, `fit_s` 1754.47 against 1774.24. That fit ends with optimizer
code 1, so the check does not run on it. The fuzz tier fired the new warning 53
times, every time on a fit with optimizer code 0 whose sdreport() SEs
were not finite.

## Firings in the suites

`dev/nanse-diag.sh` reran the 26 files where the warning fired, and
for every firing also ran base's plain `RTMB::sdreport()` on the same
fit (`dev/nanse-fire-sum.R dev/nanse-firediag2-log`):

    firings: 49  base sdreport() non-finite on: 49
    by reason: bound 9; concave 5; flat 35

Every firing is on a fit with optimizer code 0 where base had at least
one non-finite SE. Per firing (file, parameters lost of all, base's
non-finite count):

    core brms-names               3 of 9 flat         base 2 of 9
    core brms-parity-defects  x2  1 of 5 bound        base 2 of 5
    core brms-port                1 of 11 flat        base 1 of 11
    core brms-suite-emmeans       5 of 8 flat         base 1 of 8
    core brms-suite-methods   x4  5, 3, 3, 5 of 8 flat  base 1 of 8
    core case-studies             5 of 7 flat         base 5 of 7
    core diagnostics-ux       x4  3 of 6, 1 of 9, 3 of 6, 4 of 7 flat
                                  base 6 of 6, 9 of 9, 6 of 6, 7 of 7
    core grad-verdict         x4  1 of 3 bound        base 2 of 3
    core id-kron              x2  3 of 7 flat         base 3 of 7
    core me                       3 of 5 flat         base 2 of 5
    core portability              3 of 6 flat         base 1 of 6
    core predfix              x4  1 of 3 bound (2), 3 of 3 flat (2)
                                  base 3 of 3
    core se-check             x5  1 of 9, 5 of 6, 1 of 3 (bound),
                                  1 of 9, 1 of 9; base all non-finite
    core thres                    1 of 8 flat         base 8 of 8
    core v15                      3 of 4 flat         base 3 of 4
    eam defects               x2  4 of 6, 1 of 7 flat base 2 of 6, 7 of 7
    latent hmm-starts             4 of 6 flat         base 6 of 6
    latent hmm                x2  2 of 9 flat, 6 of 6 concave
                                  base 9 of 9, 3 of 6
    learn engine                  2 of 2 flat         base 2 of 2
    learn reference               1 of 6 concave      base 1 of 6
    learn stan-identity       x2  1 of 8 flat, 2 of 9 concave
                                  base 1 of 8, 1 of 9
    ode ode-nlf               x2  2 of 4 flat         base 2 of 4
    sample brms-suite-methods x2  5 of 8, 3 of 8 flat base 1 of 8
    sample gr-by-draws            3 of 9 concave      base 2 of 9
    sample laplace-draws          3 of 6 flat         base 3 of 6
    spline frailty                5 of 6 concave      base 3 of 6

`dev/nanse-fire-count.R`: the lane names more parameters than base left
non-finite on 15 firings, the same number on 13, fewer on 21. "More"
is a behavior change and deliberate: on those fits base's solve()
succeeded on a matrix with an empty or flat row and gave the parameters
along it a finite but meaningless SE (for example an `a_(Intercept)`
with a Hessian diagonal of 3.5e-9, which on a decoupled row is an SE
near 1/sqrt(3.5e-9) = 1.7e4; inferred, not printed). Those are now NaN
with the reason.

Inspected one by one. The concave cases are real: the latent `hmm()`
start on the label-symmetry axis (eigenvalue -22.1, the collapse its
own warning predicts), frmtmb.learn's `bandit4arm2_kalman_filter`
fixtures (-62.7 and an exact-Hessian diagonal of -4.2e-5),
`gr(g, by = f)` (-1.18); the frailty fixture is fitted with
`se = TRUE`, so its spectrum was not logged. The flat
cases are variance components at 0, ridges by construction (the ODE
`ke / V` ridge, the OLRE-like `gr(id, cov = I)` with one row per id,
`homcs` against sigma), a nonlinear peak off its support, a `me()`
latent sd, and an `x` the data do not see (`test-predfix`'s 3 of 3,
a Hessian of 1e-242). The bound cases are deliberate test bounds.

## Cost

`dev/nanse-cost2.R`: `frm()` with the check against
`check_se = "ignore"` on the same build, interleaved, each arm grown
past 1.2 s, minimum of 5 rounds; the control re-times the check arm:

    lm              check 0.0043s  ignore 0.0042s  ratio 1.044  control 1.060
    lmm_slope       check 0.0262s  ignore 0.0233s  ratio 1.125  control 1.022
    pois_glmm       check 0.0318s  ignore 0.0273s  ratio 1.168  control 0.982
    pois_glmm_2000  check 0.1025s  ignore 0.0923s  ratio 1.110  control 1.024
    dist_sigma      check 0.0053s  ignore 0.0050s  ratio 1.061  control 1.018
    nl              check 0.0074s  ignore 0.0066s  ratio 1.126  control 0.941
    negbin          check 0.0879s  ignore 0.0712s  ratio 1.234  control 0.949

So 4 to 23 percent, against a noise floor of 0.94 to 1.06. The Hessian
is reused by the later sdreport() (`autoscale_sdreport()` passes it as
`hessian.fixed` when it was taken at the same point), so a fit followed
by `summary()` pays for it once. The earlier estimate
(`dev/nanse-cost.R`, the Hessian alone against a fit, 2 to 20 percent)
agrees.

## R CMD check --as-cran

Once per changed package, built and checked in `dev/nanse-check/<pkg>/`
(`sh dev/nanse-check.sh <pkg>`; `R_LIBS` = lane, rellib-r5, user):

    frmtmb         Status: 1 NOTE  (V8 unavailable for the HTML manual)
    frmtmb.eam     Status: 1 NOTE  (the same)
    frmtmb.latent  Status: 1 NOTE  (the same)
    frmtmb.learn   Status: OK
    frmtmb.ode     Status: 1 NOTE  (the same)
    frmtmb.sample  Status: OK

The core tarball was built before the last edit to `R/se-check.R`,
which changed the comment above `se_empty_row` only (the measured
numbers of `dev/nanse-rowmax-sum.R`); no code changed after the build.

## Defect found, not fixed: `mo()` fits below the maximum

`dev/nanse-mo-profile.R`: given the two simplexes,
`ls ~ mo(income) * age` is a linear model, so its profile
log-likelihood over the simplexes is exact. Maximized over
stick-breaking coordinates in [0, 1]^4 with box bounds (so the simplex
boundary is a finite point) from 12 starts, against frmtmb's logLik
(rellib-r5, seeds 1 to 200):

    gap = profile max - frmtmb logLik: median 9.1e-08, 90% 0.40, max 1.95
    seeds with gap > 1e-3: 60 ; > 1e-2: 52 ; > 0.1: 39
    profile optimum on the simplex boundary (a weight < 1e-6): 192

`dev/nanse-mo-escape.R` (summary `dev/nanse-escape-sum.R`): at a
saturated simplex, move 1 percent toward each vertex in weight space
and reoptimize while that raises the likelihood.

    gap_fit    : > 1e-3 60  > 1e-2 52  > 0.1 39  max 1.95
    gap_escape : > 1e-3 46  > 1e-2 38  > 0.1 32  max 1.95
    seeds where a vertex direction improved: 20; gained up to 0.699

So on 20 of 200 seeds frmtmb reports code 0 on a softmax plateau: the
likelihood still rises toward a vertex, and the softmax gradient there
is of order the saturated weight (1e-8 or smaller), below every
stopping rule. On 38 seeds the fit is a different local maximum. This
corrects the vigport review's "the point estimates are at a valid
optimum": they are at least as good as optim's, but optim stops on the
same plateaus. Not fixed here: an escape (the prototype above) or a
simplex parameterization with a reachable boundary changes the
estimates of every saturated `mo()` fit and the sampler's
parameterization, which is a decision for its own lane. The lane's
warning does say "concave ... other starting values may reach a higher
one" on the seeds where the plateau curves downward.

## Not done, and why

- Random-effect standard errors (`diag.cov.random`) on a repaired fit
  are left as sdreport() computed them, from the failed inverse.
  Recomputing them needs sdreport()'s internals.
- The REML / `profile = TRUE` fixed effects come from the joint
  precision, whose failure keeps its existing once-per-fit warning
  (`solve_joint_precision()`); the fit-time check covers the outer
  parameters.
- A Hessian that is not positive definite while every SE is finite is
  not repaired or warned about under the default `se = FALSE` (as on
  base); under `se = TRUE` it keeps base's warning.
- A finite but enormous SE (a saturated simplex coordinate on the
  scaled inverse: 1.5e10 on seed 7) is reported, not flagged. That is
  what base reports on the 38 seeds whose rcond fell on the other side
  of machine epsilon.
- `extensions/frmtmb.sample/R/sample.R` says frmtmb_control() has
  "thirteen" names; it has fourteen now. A comment only; the test that
  reads the names (`test-stan-control.R`) checks disjointness and
  passes.

## Files touched

- R: `R/se-check.R` (new), `R/fit.R`, `R/autoscale.R`, `R/allfit.R`,
  `R/fit-end.R`, `R/methods-fit.R`, `R/predict.R`, `R/confint.R`,
  and in punch round 1 also `R/utils.R`, `R/interop.R`,
  `R/conditional-effects.R`, `R/conditional-smooths.R`.
- Rd (roxygenised): `man/frm.Rd`, `man/frmtmb_control.Rd`,
  `man/frm_allfit.Rd`, `man/diagnose.Rd`.
- Vignettes: `vignettes/diagnostics.Rmd` ("NaN standard errors"),
  `vignettes/inputs.Rmd` (the runtime table).
- `NEWS.md`.
- Core tests: new `test-se-check.R`, `test-allfit-start.R`; changed
  `helper-fuzz.R`, `test-diagnostics-ux.R`, `test-brms-names.R`,
  `test-brms-parity-defects.R`, `test-case-studies.R`,
  `test-id-kron.R`, `test-me.R`, `test-v15.R`, `test-portability.R`.
- Extension tests (no R code changed in any extension):
  frmtmb.eam `test-defects.R`; frmtmb.latent `test-hmm.R`;
  frmtmb.learn `test-engine.R`, `test-reference.R`,
  `test-stan-identity.R`; frmtmb.ode `test-ode-nlf.R`; frmtmb.sample
  `test-gr-by-draws.R`.
- `dev/nanse-*` scripts and this file.

Likely overlaps with the other lanes: `R/fit.R` (fixes, ordmix:
`check_convergence()`, `fit_assembled()`, `make_start()` is not
touched here), `R/confint.R` (fixes, ordmix), `R/methods-fit.R` (fixes,
ordmix), `R/predict.R` (fixes, ordmix), `R/autoscale.R` (fixes),
`NEWS.md`, `extensions/frmtmb.ode/tests/testthat/test-ode-nlf.R`
(fixes), `tests/testthat/helper-fuzz.R`.

Version: a minor bump (new control option, new `frm_allfit()`
argument, a new fit-time warning, breaking behavior). Floor: no
extension's R code needs this core, but two extension test files
REQUIRE the new warning (frmtmb.learn `test-engine.R`, frmtmb.ode
`test-ode-nlf.R`), so frmtmb.learn and frmtmb.ode need their frmtmb
floor raised to this release; the other changed extension tests only
allow the warning and pass on either core.
