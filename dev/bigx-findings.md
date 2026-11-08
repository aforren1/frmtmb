# Lane bigx: a covariate far from zero under the Laplace route

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-bigx`, library `wt-bigx-lib`. "Before" is `rellib-r18`.
Status: DONE, after a fix round. Sections 1 to 14 are the first round;
section 15 is the fix round after the review
(`dev/reviews/2026-11-06-bigx.md`, verdict FIX) and supersedes them
where they differ.

## 1. The defects, as measured on the base build

`dev/bigx-grid.R` fits `y ~ w` and `y ~ w + (1 | g)` with `w = z + m`,
`z ~ N(0, 1)`, 150 rows, 15 groups, for `m` from 1e3 to 1e8, families
poisson, bernoulli, gaussian and Gamma (log link), data seeds 1 to 5,
and compares each fit with the same model on `z`, an exact
reparametrization. `dev/bigx-summary.R` counts a cell "right" when the
slope is there, within 0.01 SE of the centered fit, the log-likelihood
within 1e-4 and the SE ratio within 1e-3 of 1. Logs:
`dev/bigx-out/grid-before.{log,tsv}`, `summary-1.txt`.

Base build, 240 cells: right 142, slope dropped 72, warned 12, silent
wrong 88 (a dropped slope comes with a message only, so it counts as
silent). By cell:

- With `(1 | g)`, poisson and Gamma at 1e6: 0 of 5 right each; slope
  off by up to 9.41 SE (poisson) and 9.96 SE (Gamma), up to 46.7 and
  38.4 log-likelihood units short; 4 and 2 silent.
- With `(1 | g)`, bernoulli from 1e5: 1 of 5 right at 1e5 and 0 of 5 at
  1e6; up to 2.58 SE and 3.69 units short; 4 of 5 silent at each.
- Without a group term, gaussian at 1e6: seed 2 is 12.0 SE off, 50.4
  units short, SE ratio 8.19, silent.
- Gamma with `(1 | g)` at 1e4, seed 1: 0.041 SE, 8.6e-4 short, silent
  (below grad_tol, recorded only).
- At 1e7 and 1e8: the frame drops the slope as rank deficient on 72 of
  80 cells and fits `y ~ 1`, 0.43 to 50.4 units below the model asked
  for, with a message, not a warning.

## 2. The cause, measured

`dev/bigx-cause.R` (poisson with `(1 | g)`, data seed 2), logs
`dev/bigx-out/cause-before.txt` and `cause-after.txt`:

| m | stop, code, iterations | short by | gradient at stop | outer Hessian condition at the max (unit diagonal) | design kappa |
|---|---|---|---|---|---|
| 1e4 | at the max, 0, 29 | -3.1e-10 | 1.8e-5 | 4.5e7 | 1.8e4 |
| 1e5 | at the max, 0, 32 | -2.1e-10 | 1.70 | 4.5e9 | 1.8e5 |
| 1e6 | slope -2.1e-7, 0, 2 | 30.88 | 1.6e-4, 5.6e-5, -3.7e-5 | 4.5e11 | 1.8e6 |
| 1e7 | slope -2.1e-8, 0, 2 | 30.88 | 8.5e-6, -71.8, -8.6e-5 | 4.0e13 | 1.8e7 |
| 1e8 | slope dropped by the frame | | | | |

- The start is the intercept-only point (log mean y, 0, 0). At 1e6 the
  least-squares start is skipped (design kappa above 1e5,
  `start-ls.R`), as it is by design.
- The tape is right: the raw tape at the centered maximum (intercept
  `a - b m`) reads 262.583113928 against the centered fit's
  262.583113929 at 1e6. So the inner solve and the Laplace objective are
  not the cause.
- The optimizer stops at the start's basin: the gradient there along
  the direction that moves intercept and slope together is the
  log-likelihood's per unit of intercept, about 1e-12 of the
  curvature per unit of slope (the Hessian's condition number at the
  maximum is 4.5e11). nlminb's relative function test passes after 2
  iterations, code 0. The gradient (1.6e-4) is below `grad_tol`, so the
  convergence check's trip-wire does not fire.
- newton_polish() is asked (`ill`, design kappa above 1e5) but takes no
  step: its natural Hessian with random effects is `optimHess()` with
  an absolute step of 1e-3, which on the slope moves the predictor by
  1e3, and 5 of its 9 entries are not finite at the stop.
- Autoscale does not engage: its default triggers on a column's spread
  (sd below 1e-3, `autoscale_small_sd`), and this column's sd is 1.12.
  A shift is not a scale problem; centering is what it needs, and the
  default plan never looks at the center. `autoscale = TRUE` centers
  and reaches the maximum, but on the base build it warns "false
  convergence (8)" at every offset from 1e3 to 1e6 (`probe-before.txt`).
- At 1e7 the frame's rank test (`qr()` at lm()'s 1e-7) drops the
  column: `x + 1e7` with sd 1 leaves 1e-7 of its norm once the
  intercept is out.

## 3. References (`dev/bigx-refs.R`, `dev/bigx-out/refs.txt`)

- brms 2.23.0: `standata()` keeps both columns (K = 2) at 1e5 to 1e8,
  and its Stan code centers X (`Xc`, `means_X`).
- lme4 2.1-0 `glmer()`: errors at 1e5 ("Downdated VtV is not positive
  definite"), 1e6 ("PIRLS step-halvings failed"), 1e7 and 1e8
  ("subscript out of bounds", after dropping the column). `lm()` drops
  the slope at 1e7 and 1e8. In the grid, from 1e5 up with `(1 | g)`,
  lme4 returned a slope on 0 of 20 poisson cells, 2 of 20 bernoulli
  (within 0.039 SE of the centered fit), 11 of 20 gaussian (`lmer`, up
  to 1.19 SE off) and 0 of 20 Gamma; without a group term `glm()`
  returned all 80, within 1.1e-4 SE (`grid-before.tsv`, column
  `lme4_slope_z`).
- The reference for the grid is therefore the same model on `z`, which
  is exact.

## 4. The change

(see section 6 for the files)

1. The optimizer steps in the coordinates of each ill-conditioned fixed
   design (`opt_basis()`, `opt_coords()` in `R/fit.R`): `p = M q`, with
   `M` the basis the standard errors already use (`reparam_basis()`),
   which `design_whitener()` gives a design only when its condition
   number is at least 1e4. A design below that keeps its natural
   coordinates and its fit bit for bit. The basis rides on the
   objective's environment, so every `optimize_obj()` call on that
   objective (recovery starts, the stationary-point escapes, the sd
   escape, `mo_search()`, quadrature's and importance's Laplace stage)
   uses it. A call whose bounds reach a coordinate of the basis, and a
   caller's optimizer function, keep the natural coordinates.
2. `reparam_basis()` carries `M`'s inverse (attribute `inv`), block by
   block from each design's own factor; `solve(M)` refuses `M` at 1e8.
3. The rank test of every fixed design (the frame's drop,
   `se_design_null()`, `design_whitener()`, `lp_design_kappa()`) is at
   glm.fit()'s tolerance, 1e-11 (`design_alias_tol`), not lm()'s 1e-7.
4. A restart in the design's coordinates is judged by its gain alone
   (`restart_result(any_move = TRUE)`): there the objective carries the
   rounding of the natural intercept, and a restart from the maximum
   wandered 5e-6 within the noise, gained nothing and stopped at
   "false convergence (8)" (poisson, no group term, 1e6, seed 1:
   `dev/bigx-trace2.R`, `dev/bigx-out/trace2-pois-1e6.txt`).
5. `grad_headroom()`'s last-resort Hessian takes `fit_outer_hessian()`'s
   steps, not `optimHess()`'s absolute 1e-3: with `autoscale = TRUE` at
   1e5 the lane's fit converged at code 0 and the headroom's Hessian was
   not finite, so every SE was NaN with a warning
   (`dev/bigx-auto.R`, `auto-before.txt`, `auto-after.txt`). With the
   change, `autoscale = TRUE` is silent with SE ratio 1.0000 at 1e4 to
   1e7 on data seeds 1 to 3 (`probe-after2.txt`); on the base build it
   warned "false convergence (8)" on 8 of the 9 fits at 1e4 to 1e6, with
   the right SEs (`probe-before.txt`; at 1e7 the frame dropped the
   column). At 1e3, where no basis engages, seed 1 still warns on both
   builds.
6. `natural_finish()`: a run in the design's coordinates that ends at a
   nonzero code is followed by one natural-coordinate run from its stop,
   kept only if it reports convergence and ends no higher. On the
   degree-8 raw polynomial of `test-curve-illcond.R`,
   `test-se-reparam.R`, `test-se-reparam-inv.R` and `test-se-flags.R`
   the run in the design's coordinates stopped at the maximum at "false
   convergence (8)" and warned (suite run `all1`); with this step those
   files pass with no escaped warning (run `six1`).
7. `options(frmtmb.opt_basis = FALSE)` keeps the natural coordinates;
   for the tests that pin `newton_polish()` and `sd_up_escape()`, whose
   premise (nlminb stopping short) the change removes.
8. `?frmtmb_control`, `autoscale`: one sentence that a shifted column
   needs no autoscale.

Why not autoscale's trigger: autoscale's plan centers and scales in a
pre-fit and then runs nlminb again in natural coordinates from the
mapped-back point. That second run is where the base build's
`autoscale = TRUE` warned "false convergence (8)" (item 5). By
construction its plan cannot remove an ill-conditioning that is not one
column against the intercept: it centers an interaction column `w:f` on
its overall mean, which leaves it collinear with the dummy of `f`, and
it only scales a design without an intercept, such as an ordinal one.
The design's own coordinates cover both (section 5, `extra-*.txt`), and
no pre-fit is run.

Why not a tape in the new coordinates: the objective at natural
parameters carries the rounding of the natural intercept (`b0 = a - b
m`). At 1e6 (poisson, no group term, seed 1) the objective along a
coordinate of `q` moves by 5e-10 to 1.4e-9 for steps of 1e-9 to 1e-6,
against a curvature near 300 (`trace2-pois-1e6.txt`). A tape taking `q`
would compute `X (M q)` and carry the same rounding; only a tape on
`X M`, with every reader of `fit$obj` and `fit$opt$par` mapped back,
would not. That is the residual "false convergence (8)" at 1e8 (section
5), and it was judged not worth the blast radius.

## 5. After

Grid, final build (`dev/bigx-out/grid-after5.{log,tsv}`, data seeds 1
to 5, `dev/bigx-summary.R`): 240 cells, right 240, dropped 0, silent
wrong 0, warned 11. All 11 warnings are at 1e8 and read "Optimizer did
not report convergence: false convergence (8). A Newton step ..."
(poisson with `(1 | g)` seeds 1 and 5, without seed 1; bernoulli
without seed 4; gaussian with seeds 4 and 5, without seeds 3 and 4;
Gamma with seeds 3 and 5, without seed 5). Largest |slope shift| over
all 240 cells 2.6e-4 SE (poisson without a group term, 1e8, seed 1),
largest |SE ratio - 1| 1.1e-5, log-likelihood difference -5.4e-8 to
+1.8e-7. Below 1e8 no cell warns.
Base on the same grid: section 1.

Shapes beyond the grid (`dev/bigx-extra.R`, data seeds 1 to 3,
`extra-before.txt`, `extra-after.txt`), each against the same model on
the centered column:

| shape | base build | lane build |
|---|---|---|
| cumulative `(1 | g)`, `5 z + 1e5` | seed 1: 2.5e-3 short, warned, SE NaN | all right, silent |
| cumulative `(1 | g)`, `5 z + 1e6` | 6.88, 0.11 and 4.50 short; 3.57 and 2.90 SE off silent (seeds 1, 3) | all right (\|dll\| <= 1.5e-10), silent |
| gaussian `sigma ~ w`, 1e5 and 1e6 | "false convergence (8)" on 6 of 6, up to 5.2e-6 short | right, silent |
| poisson `w * f + (1 | g)`, 1e6 | 18.0 to 30.0 short, 2.8 to 4.5 SE off, warned | right, silent |
| same at 1e7 | `w` dropped, silent, 18.0 to 30.0 short | right, silent |
| poisson `w1 + w2` at 1e6 and -5e5 | up to 2.2e-7 short | within 1.3e-9 |
| gaussian `s(t) + (x + 1e6)` | up to 1.8e-7 short | within 1.3e-9 |
| `year` with `(1 | g)`, poisson and bernoulli | right | right, identical numbers |
| poisson `(1 | g)`, 1e5, `b` with `lb = 0` | right | right (natural coordinates) |
| gaussian raw degree 6, `(1 | g)`, vs lmer | +7.1e-8 to +1.7e-6 | the same |

(The `w1 + w2` row's slope z of 1.8 to 4.0 and SE ratio 0.333 are the
units of `w2 = 3 z2 - 5e5` against `z2`, not an error; its
log-likelihood is the check.)

The defect's own designs, `test-se-reparam.R`'s raw powers with
`(1 | g)` and no group effect (`dev/bigx-edge.R`, `edge-before.txt`,
`edge-after.txt`), against lme4's ML maximum:

| degree, seed | base | lane, design's coordinates | lane, natural |
|---|---|---|---|
| 5, 4 | -1.07e-4, log sd -17.4, SE warning | +5.1e-12, log sd -5.23 (lmer sd 0.00534), silent | as base |
| 6, 4 | -1.05e-4, log sd -95.3, SE warning | +6.1e-9, log sd -5.24, silent | as base |
| 7, 4 | +3.8e-3, by `sd_up_escape()` | +3.8e-3 without the escape | as base |

So two of the base build's "sd at zero" stops were 1e-4 short of a
small interior maximum, which the lane reaches.

Cost (`dev/bigx-cost.R`, `cost.txt`): `opt_basis()` on a design it
does not engage on is a QR of the design: 0.0075 s on n = 1e4, p = 20
(0.9 percent of a 0.86 s poisson fit) and 0.45 s on n = 1e5, p = 50
(1.6 percent of 27.6 s); minimum of 5 interleaved rounds, the control
read 0.84 to 0.88, so the noise is about 15 percent.

## 6. Tests

One file per process, `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`,
`dev/bigx-run-par.sh` (a copy of `dev/release/run-par.sh` whose runner,
`dev/bigx-run-tests.R`, puts `rellib-r18` behind the lane library).
Every log's `lib:` line reads frmtmb 0.79.1 from `wt-bigx-lib`; the
extensions load from `rellib-r18` (none of their code changed).

Run `all1` (first build, before `natural_finish()` and the test
edits), `dev/bigx-out/all1.log`: 500 of 500 files ran; failing
test-grad-verdict.R (4), test-optrap.R (4), test-se-flags.R (2),
test-se-reparam.R (6); escaped warnings in test-curve-illcond.R and
test-se-reparam-inv.R; every extension file clean. Causes: section 4,
items 6 and 7.

Run `all2` (final build), `dev/bigx-out/all2.log`, 500 of 500:

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 309 | 25594 | 0 | 0 | 15 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1854 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 12 | 683 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 545 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 85 | 4326 | 0 | 0 | 1 | 0 |
| frmtmb.spatial | 8 | 561 | 0 | 0 | 0 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |

The 29 skips are other gates: 15 `FRMTMB_SCALE_TESTS`, 13
`FRMTMB_DRMTMB_FIT_TESTS`, 1 `FRMTMB_FUZZ` (run `gates1` below).

New test file `tests/testthat/test-design-offset.R`: on rellib-r18
pass=8 fail=11 err=3 (`test-offset-before.txt`: the poisson slope 7.55
SE off, bernoulli 1.84, the 1e8 slope dropped, no basis, autoscale's
"false convergence (8)" warning, `any_move` unknown); on the lane
build pass=30 fail=0 err=0 skip=0 warn=0 (run `all2`).

Test files edited, each because the change removes its premise (nlminb
stopping short in natural coordinates):

- test-grad-verdict.R, "a flat ridge warns on curvature the gradient
  does not show": the short stop now comes from a caller's optimizer
  function wrapping nlminb, which keeps natural coordinates.
- test-optrap.R, three tests that assert `newton_polish()` took steps:
  `options(frmtmb.opt_basis = FALSE)`.
- test-se-reparam.R, "an sd near zero is judged by the check": the
  natural coordinates for the check, plus a default arm that asserts the
  interior maximum (silent, log sd 5 or more above the stop, higher
  log-likelihood); "raising the optimizer's caps keeps the sd escape":
  natural coordinates for `f2`.

The three edited files on rellib-r18 (`base-edited.log`):
test-grad-verdict.R pass=74 fail=0, test-optrap.R pass=147 fail=0 (the
option is unknown there, so they run the old path), test-se-reparam.R
pass=76 fail=3: the new default arm fails on the base build (warns,
log sd -17.4, no higher log-likelihood), as it should.

Other gates (`gates1.log`, all gates on, lane build): test-curve-scale.R
2/0, test-drmtmb-agreement.R 131/0, test-fuzz.R 2/0, and the scale
tiers of coupling 6/0, eam 18/0, latent 6/0, learn 14/0, ode 2/0,
sample 3/0 (pass/fail; 0 error, 0 skip, 0 warn each).

## 7. OpenBLAS 0.3.26 with LAPACK

`dev/ciharden-openblas.sh 0.3.26 lapack` (Rblas 74 symbols and Rlapack
616 to OpenBLAS; dgemm 2000 in 0.11 s). With that R copy:
test-design-offset.R pass=30 fail=0 err=0 warn=0, test-se-reparam.R
pass=79 fail=0 err=0 warn=0 (`ob-*.txt`). The grid (`grid-ob.tsv`):
240 of 240 right, 5 warned, all at 1e8, "false convergence (8)"
(gaussian with `(1 | g)` seed 2, without seeds 3 and 4; Gamma with
seeds 3 and 5). Not run under OpenBLAS: 0.3.32, and the rest of the
suite. The new tests' designs have condition numbers of 2e5 and above,
far from the basis's threshold of 1e4; whether an existing test's
design sits near that threshold was not checked.

## 8. One session per suite

`dev/bigx-run-sessions.sh sess1 3` (a copy of
`dev/release/run-sessions.sh` with the lane's runner and `rellib-r18`),
`NOT_CRAN=true`, `dev/bigx-out/sess1.log`: 27 of 27 sessions, fail 0
and error 0 in every one; every LEAK line is tagged `[with load of
...]` (extension or package loads), none is a test's own.

| package | files | forward pass | reverse pass | preload pass |
|---|---|---|---|---|
| frmtmb | 296/309 | 22544 | 22694 | 22787 |
| frmtmb.coupling | 11/11 | 551 | 551 | 551 |
| frmtmb.eam | 33/33 | 1848 | 1848 | 1848 |
| frmtmb.latent | 12/12 | 683 | 683 | 683 |
| frmtmb.learn | 16/16 | 474 | 474 | 474 |
| frmtmb.ode | 11/11 | 549 | 549 | 549 |
| frmtmb.sample | 80/85 | 3962 | 3962 | 3962 |
| frmtmb.spatial | 8/8 | 478 | 478 | 478 |
| frmtmb.survival | 15/15 | 552 | 552 | 552 |

## 9. The warning at 1e8

The full text on three of the grid's warned cells
(`dev/bigx-warn1e8.R`, `warn1e8.txt`): "Optimizer did not report
convergence: false convergence (8). A Newton step from the stop gains
7.54e-09 [9.11e-10, 5.64e-08] in log-likelihood, within
frmtmb_control(grad_tol = 0.001), so the estimates are at the maximum
and their standard errors are judged as for a converged fit". So these
fits are right and say so; the warning is the existing policy for a
nonzero code at the maximum.

## 10. R CMD check

`dev/bigx-check.ps1`: core built from the worktree and checked with
`--as-cran` in `dev/bigx-check/` (`R_LIBS` lane library, rellib-r18,
user library; pandoc and TinyTeX on PATH;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`): `Status: 1 NOTE`, the V8
math-rendering NOTE on the HTML manual. Tests 20 min OK, examples and
`--run-donttest` OK, vignettes rebuilt OK, PDF manual OK. No extension
changed, so none was checked.

## 11. The equivalence survey

Every one of the 172 pairs remeasured on the lane build
(`dev/equiv-measure.R <id> dev/bigx-equiv-out`, `FRMTMB_EQUIV_LIB` the
lane library then rellib-r18; logs `dev/bigx-equiv-out/<id>.log`), and
compared with `dev/equiv-data/results.tsv` by `dev/bigx-equiv-compare.R`
and `dev/bigx-equiv-metrics.R` (`dev/bigx-out/equiv-compare.txt`,
`equiv-metrics.txt`): the L1 constant and gradient differ by more than
1e-6 relative on 0 pairs, the L2 metric on 0, the L4 z text on 0,
`frm_ok` and `brms_ok` on 0. The measured verdicts are the record's
`verdict_measured` on every pair; the rows the compare script lists
differ only in the adjudicated columns (`verdict`, `agree`, p069's
`level`), which a remeasure does not carry, and in `agree` read as NA
where the record stores FALSE on "cannot measure". No verdict changes,
so no sampled pair reaches the change.

## 12. Bump, floor, compat, ledger

- Bump: minor (two BREAKING changes: the optimizer's coordinates move
  every fit on an ill-conditioned fixed design, and the alias tolerance
  keeps columns lm() drops).
- Floor: none. No extension calls a changed function in its `R/`, and
  every extension suite passed against this core, loading the
  extensions' rellib-r18 builds (sections 6 and 8).
- Compat rows: none changed.
- Ported brms ledger: not regenerated (consolidation does that), not
  run separately; the equivalence survey shows no change (section 11).

## 13. Not done, and why

- A tape in the design's coordinates, which would remove the residual
  "false convergence (8)" at 1e8 (section 4, last paragraph).
- The ordinal paragraph of `vignettes/diagnostics.Rmd` ("Center a column
  whose mean is more than about 1e4 times its standard deviation") is
  about `importance =`. The importance rounds run on their own
  objective, which gets no basis; only their Laplace stage does. Not
  remeasured, so the paragraph stays.
- The date/time message of `report_datetime_columns()` still says such
  a column "can stop the optimizer converging". Left: it is also about
  the intercept's meaning.
- `autoscale = TRUE` at `x + 1e3` (seed 1) still warns "false
  convergence (8)" on both builds: the design's condition number there
  (2.2e3) is below the basis's threshold.
- `newton_polish()`'s own Hessian with random effects is still
  `optimHess()` with absolute steps, not finite on a shifted column (5
  of 9 entries at 1e6). With the basis the optimizer reaches the maximum
  and the steps are not needed on the grid; the defect stays for a fit
  whose bounds keep the natural coordinates.
- `profile()` or a refit with a coefficient held by a map builds a new
  objective; on a shifted coefficient that was not measured.
- A saved fit carries `M` and its inverse on `obj$env` (two np x np
  matrices; 2.6 MB at the cap of 400 parameters).
- OpenBLAS 0.3.32 and the whole suite under OpenBLAS (section 7).

## 14. Files

R: `R/fit.R` (opt_basis(), opt_coords(), natural_finish(),
optimize_obj(), optimizer_from_best(), restart_result(),
grad_headroom(), fit_assembled(), `?frmtmb_control`), `R/frame.R`
(design_alias_tol, the drop, the sparse screen's comment),
`R/se-reparam.R` (reparam_basis()'s `inv`, design_whitener()'s
tolerance), `R/se-check.R` (se_design_null()'s tolerance),
`R/autoscale.R` (lp_design_kappa()'s tolerance). `man/frmtmb_control.Rd`
(roxygen). `DESCRIPTION`: roxygen 8.1.1 wrote
`Config/roxygen2/version: 8.1.1`. `NEWS.md`. Tests:
`tests/testthat/test-design-offset.R` (new), `test-grad-verdict.R`,
`test-optrap.R`, `test-se-reparam.R`. Lane scripts `dev/bigx-*`.

Status: DONE.

## 15. Fix round (review of 2026-11-06, verdict FIX)

Every number below is from the fix build in `wt-bigx-lib`; the
reviewer's scripts (`dev/review-bigx-*.R`, `REV_ARM=after`) were rerun
unchanged and their logs written to `dev/bigx-out/fx*-*.txt`.

### B1. A bound on a shifted slope

Cause: `opt_coords()` returned the natural coordinates for the whole
fit when any finite bound touched the basis.

Change: `opt_basis()` keeps the basis block by block
(`reparam_basis()` now lists its blocks), with each design block's
centering map (`opt_center()`: the column means, the intercept's
position, the centered design's condition number). `opt_coords()`
gives a block no bound reaches its design's coordinates; a design block
with a bound on a coefficient but not on its intercept the centering
map (every coefficient but the intercept stays natural, so the box stays
a box); and a block with a bound on its intercept, with no intercept, or
whose centered design still has a condition number of 1e3 or more, the
natural coordinates, noted on the objective. `opt_coords_warn()` then
warns once, class `frmtmb_ill_design`, naming the columns far from zero
(mean above 100 sd). An nl body's block is not noted.

Reviewer's bound designs (`fx-guard-bounds.txt`, poisson `(1 | g)`,
1e6, seeds 1 to 5): `lb = -5` and `ub = 5` within 9e-10 log-likelihood
and 0.00 SE of the centered fit with the same prior, silent, 10 of 10;
`lb = 0.45` (binds on seeds 1 to 3) fits, at the centered fit's
log-likelihood within 6e-10, the only warning that the bound holds the
slope's SE. Base: 6.3 to 9.4 SE off, 4 of 5 silent, and `lb = 0.45`
stopped with "NA/NaN gradient evaluation" (the review's table). Pinned
in `test-design-offset.R` ("a bound on a shifted slope fits the
centered model", `lb` -5 and 0.45; "only an ill-conditioned design gets
the design's coordinates" now asserts the centering map and the
`frmtmb_ill_design` warning for a bound on the intercept).

### M1. Shapes with no basis

Rerun of `dev/review-bigx-shapes.R` (`fx3-shapes.txt`, seeds 1 and 2):

| shape | first build | fix build |
|---|---|---|
| nl body `b0 + b1 * w` | 7.58 SE, 24.32 short, silent; 9.89 SE warned | 1.4e-5 and 6.9e-7 SE, within 2.6e-9, silent |
| `w + mo(o) + (1 | g)` | 7.39 SE, 27.10 short, warned | within 1.9e-9 and 3.7e-5 SE; still warns "curvature unusable", SEs NaN |
| `acat`, `cs(w)` | 0.78, 1.44 SE off, warned | unchanged (filed) |
| `w + (1 + w | g)` | 0.72, 0.77 SE off, warned | unchanged (filed) |

- nl body: `nl_basis_descent()` runs a second descent in coordinates
  that also whiten the body's Jacobian at the first stop (`opt_basis()`
  with `est`), kept only when it ends lower; a Jacobian under the
  threshold gets no block, so other nl fits are left bit for bit.
  Pinned: "a nonlinear body linear in a shifted column fits the
  centered one".
- `mo()`: `reparam_basis(mo = TRUE)` (the optimizer's call) gives the
  design's non-`mo()` columns a block, and the fit reaches the maximum.
  Its SEs stay NaN with a warning: `outer_fd_steps()` skips `mo()`
  designs, so the Hessian's step on the shifted column is an absolute
  1e-3. Scaling those steps would move the SEs of every `mo()` fit with
  a column wider than one unit (lane adhess kept 1e-3 there on
  purpose), so it is filed; the standard-error route keeps `mo = FALSE`.
- `cs()` and the random slope: filed (NEWS, Limits). A random slope on a
  shifted column makes theta ill-conditioned, which no fixed basis
  reaches; centering Z is a change to the covariance blocks.
- NEWS lists what is fixed and has a Limits section for the rest.

### M2. The frame's test and design_whitener()

Change: one rank test, `design_rank_qr()` (QR of the design scaled to
unit root mean square columns, at `design_alias_tol` 1e-11), in the
frame's drop and in `se_design_null()` (its null vectors mapped back
from the scaled columns), as in `design_whitener()`; then
`design_kappa_keep()`: when the kept columns' scaled condition number
is above `design_kappa_max` = 1e10, the frame keeps, in their order,
each column that leaves it under the bound, with the base message
naming the dropped ones. Their prediction null space is the scaled
design's smallest right singular vectors (`design_near_null()`), since
a near alias has no exact one.

Why 1e10 (`dev/bigx-poltol.R`, `poltol.txt`): at the scaled test's
tolerance of 1e-11, 1e-10 or 1e-9 alone, 8 of 12 raw polynomials on
[1000, 1010] (degree 5 to 8, seeds 501 to 503) kept designs with
condition numbers of 4.7e13 to 2.8e16 and stopped at false convergence
with every SE NaN; the designs kept at 4.6e8 to 3.7e11 reached glm()'s
maximum on the same columns to 1.6e-5 with finite SEs; a shift of 1e9
(kappa near 2e9) fits the centered model.

After (`dev/bigx-polyfar.R`, `polyfar-after2.txt`): 12 of 12 keep 4
columns (base kept 4 as well), basis on, code 0, finite SEs, no
warning, with the drop message. Pinned: "a raw polynomial far from zero
keeps only columns it can estimate" (seed 501, degree 6: message, kept
kappa at most 1e10, finite SEs, glm()'s log-likelihood on the kept
columns to 1e-6).

### M3. Profiles

Change: `profile_run()` hands TMB::tmbprofile() and TMB::tmbroot() a
view of the fit's objective in the optimizer's coordinates (`p = M q`,
the profiled combination `M' l`), with the fit's objective state put
back; tmbroot's range is formed from the fit's own SE, since it would
call sdreport() on the view. Used by `confint(method = "profile")`,
`confint(method = "uniroot")`, `profile()`, `hypothesis(method =
"profile")` and `conditional_effects()`'s profile bands. A fit with no
basis is profiled as before, bit for bit (same call).

Reviewer's profile part (`fx-follow.txt`, seeds 1 to 3): the slope's
profile interval equals the centered one to 6 digits (0.26219 0.49987,
0.302343 0.513781, 0.268649 0.475463), theta's too, no warning. Base and
the first build: NA, NA. Pinned: "the profile of a shifted slope is the
centered one" (profile and uniroot, within 1e-3 SE).

### M4. The threshold

`opt_basis_kappa` = 1e3 for the optimizer; the SEs keep 1e4.

- The reviewer's 36 well-conditioned designs (`review-bigx-noreg.R`,
  tag `bigxfix`, against `before1` and `before2`,
  `fx2-noreg-cmp.txt`): parameters and objective identical on 36 of 36;
  none of them gets a basis at 1e3 either (the four "ill" ones are
  below it).
- The reviewer's 80 cells with scaled kappa 2e3 to 6e6
  (`review-bigx-noreg2.R`, `fx-noreg2-after.txt`): basis on 80 of 80
  (70 on the first build), max |dll| 6.6e-9, max shift 8.4e-5 SE, 0
  warned.
- The band (`review-bigx-thresh.R`, `fx2-thresh.txt`, offsets 3e3 to
  2e4, four families, seeds 1 to 5): 260 cells, basis on 260, max |dll|
  5.4e-9, max shift 4.2e-5 SE, max |SE ratio - 1| 4e-6, 0 warned. The
  first build had 72 cells without a basis, among them bernoulli 0.101
  SE and Gamma 0.042 SE off, silent.
- Below the threshold (`dev/bigx-thresh-low.R`, offsets 100 to 2000,
  `thresh-low-after.txt`, `thresh-low-before.txt`): 200 cells, basis on
  116 on the fix build and 0 on base; max |dll| 6.7e-9 on both, max
  shift 4.2e-5 and 2.8e-5 SE, 0 warned on both.
- The 172 equivalence pairs: section 15, "Equivalence" below.

### Minors

1. `natural_finish()` runs only after "false convergence (8)". The
   reviewer's iteration caps (`fx-guard-itermax.txt`, seeds 1 to 5,
   `iter.max` 2, 3, 5, 8): every capped fit keeps code 1 and the
   "iteration limit reached" warning; none is relabeled. Pinned: "an
   iteration limit keeps its code and its warning".
2. Unit test "natural_finish() keeps the stop when the natural run does
   not converge": with `iter.max = 1` the natural run ends at code 1 and
   `opt` must come back identical. The reviewer's mutant mD (keep any
   result) returns the natural run's result there, so the test fails on
   it by construction (`identical()`); the mutant itself was not
   rebuilt.
3. `frm_allfit()`: its wrappers carry `frmtmb_builtin` and take the
   basis (`user_optimizer()`). `fx-follow.txt`, seed 2, 1e6, from
   `start_from_mle = FALSE`: bobyqa 3.45e-5 SE and NLopt 4.5e-7 SE from
   the centered fit (first build: 7.55 SE, 33.1 and 33.6 short).
4. The comment in `grad_headroom()` now names the differences, not
   optimHess().
5. Quadrature and importance tapes: filed (NEWS, Limits). On the fix
   build bernoulli `quadrature = TRUE`, seed 1, is still 1.2e-5 short,
   0.0049 SE, silent (`fx3-shapes.txt`).

### Merge note for consolidation (nestquad)

nestquad's `freeze()` builds each round's tape with its own
`MakeADFun()`; attaching `lap_obj$env$.frmtmb_opt_basis` to those tapes
alone does nothing for the rounds, because each round's box is finite
in every coordinate. With this round's `opt_coords()`, a finite box on
a design block's coefficients gives that block the centering map, not
the natural coordinates, but a box on its intercept still gives the
natural ones and notes a warning. Consolidation must do one of:

- pass the round's box in the new coordinates: give `optimize_obj()` a
  box already in `q` (`to_q(par) +- radius / q-units`) for the rounds;
  or
- state that nested rounds run in natural coordinates: do not attach
  the basis to the round tapes (so no `frmtmb_ill_design` warning is
  raised by a round's box), and attach it only to the final tape, for
  `escape_stationary()` and `mo_search()`.

The second is the smaller change and keeps every nested round as it is
on nestquad's build.

### Equivalence, autoscale, near aliases

- The 172 pairs on the fix build (dev/bigx-equiv-out2, compare and
  metrics scripts with BIGX_EQUIV_OUT, equiv-compare2.txt,
  equiv-metrics2.txt): L1 constant and gradient, L2 metric, L4 z text,
  rm_ok and rms_ok differ on 0 pairs; the rows the compare script
  lists are the same 21 adjudication-only rows as the first round
  (identical id list). No verdict moves with the threshold at 1e3.
- utoscale = TRUE (dev/bigx-probe.R, probe-after3.txt): silent
  with SE ratio 1.0000 from 1e3 to 1e7, seeds 1 to 3 (15 of 15); the
  first build still warned at 1e3, seed 1.
- Near aliases (dev/review-bigx-collin.R rerun, x-collin2.txt):
  x2 = x1 + 1e-9, 1e-8, 1e-6 noise and proportions summing to 1 to
  1e-9 are kept as glm() keeps them, SEs 6.5e4 to 1.6e8 within 1.1
  percent of glm's, log-likelihood within 1.2e-7, no warning, no
  message; to 1e-12 the proportion is dropped with the message. The raw
  polynomials on [1000, 1010] keep 4 columns with the message, finite
  SEs, no warning.
- Grid on the fix build (grid-after6.tsv): 240 of 240 right, 11 warned
  (all at 1e8), as in section 5.
- The five named test files on the fix build (ive1.log):
  test-design-offset.R 64/0/0/0/0, test-se-reparam.R 79/0/0/0/0,
  test-optrap.R 147/0/0/0/0, test-grad-verdict.R 74/0/0/0/0,
  test-autoscale.R 50/0/0/0/0 (pass/fail/error/skip/warn).

- R CMD check --as-cran on the fix build (dev/bigx-check.ps1, then
  the check tree deleted): `Status: 1 NOTE`, the V8 NOTE; tests 16 min
  OK, examples and `--run-donttest` OK, vignettes OK.

### Suite on the fix build

Run `all3` (one file per process, `NOT_CRAN`,
`FRMTMB_BRMS_FIT_TESTS`, `dev/bigx-out/all3.log`), 500 of 500, every
`lib:` line the lane core:

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 309 | 25628 | 0 | 0 | 15 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1854 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 12 | 683 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 545 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 8 |
| frmtmb.sample | 85 | 4326 | 0 | 0 | 1 | 0 |
| frmtmb.spatial | 8 | 561 | 0 | 0 | 0 | 0 |
| frmtmb.survival | 15 | 552 | 0 | 0 | 0 | 0 |

The 8 escaped warnings were deSolve's "Excessive precision requested"
and "Returning early" in test-ode-events.R and test-ode-nlf.R (4 each),
raised at trial points of `nl_basis_descent()`'s second descent. That
descent now muffles warnings (the fit's checks speak to the point it
keeps). After it (`ode2.log`): the 11 frmtmb.ode files pass with 0
warnings (test-ode-events.R 142/0/0/0/0, test-ode-nlf.R 6/0/0/0/0, the
scale file skipped by its gate), and test-design-offset.R 64/0/0/0/0.
The R CMD check above ran on the build before this muffle; the muffle
changes no value a test reads.

Not rerun on the fix build: the one-session runs (section 8),
OpenBLAS (round rule).

## 16. Fix round: files, bump, still open

Files touched in the fix round: `R/fit.R` (opt_basis(), opt_center(),
opt_coords(), opt_coords_warn(), user_optimizer(), opt_basis_kappa,
nl_basis_descent(), natural_finish()'s gate, fit_assembled(),
grad_headroom()'s comment), `R/frame.R` (design_rank_qr(),
design_kappa_max, design_kappa_keep(), design_near_null(), the drop),
`R/se-reparam.R` (reparam_basis()'s `blocks`, `min_kappa` and
`mo`; design_whitener()'s `min_kappa`), `R/se-check.R`
(se_design_null() on the scaled QR), `R/confint.R` (profile_run(), its
four callers), `R/conditional-effects.R` (the profile band's tmbroot
call), `R/allfit.R` (`frmtmb_builtin`), `NEWS.md`,
`tests/testthat/test-design-offset.R`.

Bump: minor. Floor: none (no extension calls a changed function; every
extension suite passed on this core). Compat rows: none changed.

Still open (filed in NEWS, Limits): `mo()` SEs on a shifted column
(`outer_fd_steps()` skips `mo()` designs); `cs()` and a random
slope on a shifted column; the quadrature and importance tapes; a
bounded nonlinear body. Above 400 outer parameters no basis (review
minor 5, unchanged, warned on both builds).

Status: DONE.
## 17. Short check (the review's last item, M2)

Cause: `design_kappa_keep()` judged the 1e10 bound on the raw scaled
design, so it dropped `w1:w2` and `I(w^2)` of shifted columns, which
the basis fits, and the frame called them "rank deficient".

Change: the bound is judged on the design with every column but the
intercept centered, as `opt_center()` centers it (the raw design where
there is no intercept). A column dropped by the bound gets a warning of
its own, class `frmtmb_design_capped` (`design_cap_warn()`): the design's
condition number at unit scale, centered, is above 1e10, the design has
full rank, and the model fitted is not the one written; it names the
columns. "Rank deficient" (the message) is kept for columns aliased at
the rank tolerance.

`dev/bigx-prod.R` (`prod-after.txt`, `prod-before.txt`), seeds 1 to 3:

| design | m | base | lane |
|---|---|---|---|
| gaussian `w1 * w2` | 1e4 | `w1:w2` dropped, 40.0 to 46.1 short | kept, within 2.5e-7; 1 of 3 warns "false convergence (8)" at the maximum |
| | 1e5 | dropped, same | kept, 1.0e-5 to 1.9e-5 above the centered fit (rounding); 3 of 3 warn false convergence at the maximum |
| | 1e6 | dropped | dropped as rank deficient (at 1e-11 scaled; glm drops it too) |
| poisson `w1 + I(w1^2)` | 1e4 | `I(w1^2)` dropped, 0.60 to 5.66 short | kept, within 4.4e-8, silent |
| | 1e5 | dropped | kept, within 7.9e-6; 3 of 3 warn false convergence at the maximum |
| | 1e6 | dropped | dropped as rank deficient |
| raw `poly(x, 4/6/8)` on [1000, 1010] | | | 5 columns kept, finite SEs, glm on the kept columns to 1.5e-5; the bound warns `frmtmb_design_capped` on degree 8 (3 of 3) and degree 6 seed 503 |

Reviewer's `review-bigx-fx.R` part `rank` (`sc-fx-rank.txt`): near
aliases in every order kept as glm keeps them; products and powers at
100 to 1e5 kept (at 1e4 and 1e5 within 9.4e-6 of the centered fit, with
false convergence at the maximum on 3 of 4); `y ~ w` kept up to 1e10.
Parts `bounds`, `nl`, `profile` (`sc-fx-bounds-nl-profile.txt`): as in
the review.

Pinned in `test-design-offset.R`: "products and powers of shifted
columns are kept and fitted" (`w1 * w2` and `w1 + I(w1^2)` at 1e5); the
polynomial test now uses degree 8 and requires the capped warning, the
rank message, finite SEs and a centered kappa at most 1e10.

NEWS: bullet 2 names products and powers and the new warning class;
Limits adds the curved body `b0 * exp(b1 * (z + m))` and a body with an
internal alias. `dev/test-backlog.md`, new section "Filed in the round
opened 2026-11-06": the bound on class `Intercept` (natural against
brms's centered intercept), not fixed.

Runs on this build: test-design-offset.R pass=71 fail=0 err=0 skip=0
warn=0; grid (`grid-after7.tsv`) 240 of 240 right, 11 warned, all at
1e8.
Core suite on this build (`core4.log`, one file per process, gated):
309 of 309 files, pass 25635, fail 0, error 0, skip 15 (other gates),
warn 0; every `lib:` line the lane core. R CMD check --as-cran on this
build: `Status: 1 NOTE` (V8); tests, examples, `--run-donttest` and
vignettes OK; the check tree deleted.

## 18. Confirmation (the review's last item)

Cause: the rank test ran on the uncentered scaled design, so at an
offset of 1e6 it called `w1:w2` and `I(w1^2)` aliased (rank 3
uncentered, 4 centered) and dropped them with "rank deficient", 32 to
47 units short, though the data identify them.

Change (`R/frame.R`): `design_kappa_keep()` runs the rank test on the
centered design (`design_center()`, which sets a column constant to
rounding to exactly zero, so a constant column is still an exact alias).
Those columns are `aliased` and get the message. Of the rest, a column
that would leave a design the optimizer cannot fit (centered condition
number above 1e10, or uncentered columns dependent at 1e-11, where
`design_whitener()` refuses) is `capped` and gets the warning
`frmtmb_design_capped`, which now says the data identify the column and
the model fitted is not the one written.

Reviewer's `dev/review-bigx-fx2.R` (`sc2-fx2.txt`), seeds 1 to 3, and
`dev/bigx-prod.R` (`prod-after2.txt`):

- 1e4 and 1e5: `w1 * w2` and `w1 + I(w1^2)` kept and right (within
  1.3e-5 of the centered fit, and glm's), some with "false convergence
  (8)" at the maximum.
- 1e6: `w1:w2` and `I(w1^2)` dropped with `frmtmb_design_capped`, no
  "rank deficient" message (6 of 6 in `prod-after2.txt`, seeds 1 to 3;
  10 of 10 in the fx2 script, seeds 1 to 5). That script's code-8
  summary on this build: 29 of 90 at code 8, gains 1.95e-9 to 9.92e-6,
  none short by more than 1e-6, none with another warning.
- Raw polynomials on [1000, 1010]: 5 columns kept, finite SEs; the
  columns aliased at 1e-11 after centering get the message, the others
  the capped warning.

Pinned in `test-design-offset.R`: `w1 * w2` and `w1 + I(w1^2)` at 1e6,
kept and right or capped with the warning, never a message; the
polynomial test requires the capped warning's text.

`dev/test-backlog.md`: the review's recommendation on code-8 fits at
the maximum (a classed message when the Hessian is usable, every SE
finite and the Newton gain below 1e-5; 29 of 90 battery fits, all at
the maximum, gains 2.0e-9 to 9.9e-6), filed as its own change.

NEWS bullet 2 says the rank test is centered and that products and
powers shifted by 1e6 meet the capped warning.

Runs: test-design-offset.R pass=75 fail=0 err=0 skip=0 warn=0; grid
(`grid-after8.tsv`) 240 of 240 right, 11 warned, all at 1e8.

Core suite `core5` (one file per process, gated): 309 of 309 files,
pass 25639, fail 0, error 0, skip 15 (other gates), warn 0; every `lib:`
line the lane core. R CMD check --as-cran on this build: see below.

`Status: 1 NOTE` (V8); tests 16 min OK, examples, `--run-donttest` and
vignettes OK; the check tree deleted.

Status: DONE.
