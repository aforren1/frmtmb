# Lane mfxse: marginaleffects standard errors on models with group-level terms

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1), worktree
`frmtmb-wt-mfxse`, library `wt-mfxse-lib`. Status: done, uncommitted. This
file is the record; every number below is pasted from the log named beside
it.

Bump: minor (a BREAKING change to what `get_coef()`, `get_vcov()` and the
marginaleffects standard errors return). Floor: none; no extension changed
or calls the changed functions. Files touched: `R/interop.R`, `R/compat.R`
(4 rows), `man/frmtmb-marginaleffects.Rd` (roxygen), `NEWS.md`,
`tests/testthat/test-marginaleffects-re.R` (new), `test-me-newdata.R`,
`test-v07.R`, `helper-warnings.R`; scratch `dev/mfxse-*` (scripts tracked,
logs local). `DESCRIPTION` left as it was: roxygen2 8.1.1 rewrote
`Config/roxygen2/version` to 8.1.1 and the lane set it back to 8.1.0.

## The defect, and where it lived

marginaleffects forms every standard error as `J V J'`, with `J` a forward
difference over `get_coef()` (it perturbs one entry, calls `set_coef()` and
`get_predict()`) and `V` from `get_vcov()` (marginaleffects 1.0.0,
`get_delta_jacobian()`, `std_error_from_jacobian()`; dumped with
`dev/mfxse-dump-me.R`). frmtmb's methods are in `R/interop.R`:

- `get_coef.frmtmb_fit()` returned `interop_coef_vector()`: `beta`, the
  estimated `betad`, the ordinal extras and the `mo()` simplexes.
- `get_vcov.frmtmb_fit()` returned `interop_vcov()`, their covariance, and
  warned with class `frmtmb_interop_se_part` (`interop_warn_se_part()`).

So a prediction that reads the group-level effects `b` (and, for an exact
`gp()`, the `theta` its field `sdgp L z` reads) held them at their modes.
The variance that `fitted()` and `frm_linpred(se.fit = TRUE)` hold, the joint
covariance of the fixed effects, the variance parameters and the predicted
effects (`get_joint_cov()`: `[Qrr^-1 + J Vf J', J Vf; Vf J', Vf]`, or
`Vrest + root root'` where it is held as a factor), was missing.

## What each quantity must hold, and which call maps to which

- Group-level terms kept (`re_formula = NULL` or a formula keeping a term):
  the joint covariance of the coefficients, the variance parameters and the
  predicted effects. A prediction at a level the fit saw reads `b`; the
  variance parameters reach it through the `J Vf J'` part of `b`'s block.
  `predictions()` (a level), `avg_predictions()` (a mean of levels),
  `slopes()`/`avg_slopes()` (a difference of levels in `x` over `eps`),
  `comparisons()`/`avg_comparisons()` (a difference of levels), and
  `hypotheses()` on any of these all read the same `J V J'`.
- `re_formula = NA`: the coefficients alone, plus a smooth's coefficients,
  which `re_formula` does not drop.
- An exact `gp()`: `theta` (sdgp, lscale) as well as `z`, since the field is
  `sdgp L(lscale) z`; `rr` blocks and registered structures whose expansion
  reads `theta` the same way (`block_reads_theta()`).
- `hypotheses()` on the model's own coefficients: those coefficients, as
  `hypotheses()` on a brmsfit reads its population-level parameters.

## The change (R/interop.R)

- `interop_re_part()`: the parameters appended after `interop_coef_vector()`:
  every `b` of a block a prediction reads (all but a latent residual's; under
  `re_formula = NA` only the smooth, `gp()` and `hsgp()` blocks), then the
  `theta` of each block whose coefficients read it. Names are
  `frm_lp_basis()`'s `coef_names`, `b.<term>.<level>.<coef>` and
  `theta.<name>`. `NULL` when the objective holds no `b` (quadrature,
  importance) or the fit has no block.
- `get_coef()`: the vector plus that part; the vector alone inside
  `marginaleffects::hypotheses()` called on the fit (`interop_in_hypotheses()`,
  fix round; the first build keyed on `hypothesis_is_formula`, which the joint
  test never passes).
- `set_coef()`: writes the appended part into `estimates$b` and
  `estimates$theta`, by name.
- `get_vcov()`: `interop_vcov_re()` keeps `interop_vcov()` as the coefficient
  block (so `vcov()`, `summary()` and `hypotheses()` errors are unchanged) and
  adds the cross block `M V[, a]` and the effects' block `V[a, a]` from
  `get_joint_cov()`, with `M` (`interop_joint_map()`) the map from the joint
  rows to the coefficient rows (a selection, or the cumulative-weight Jacobian
  of a `mo()` gap simplex, `mo_gap_cum_jacobian()`, factored out of
  `interop_vcov()`).
- User `vcov =`: a covariance over the outer vector (`vcov(full = TRUE)`,
  `vcov_cluster(full = TRUE)`, `~ g`; `"HC0"` and `"HC1"` are refused on any
  mixed fit, so they never reach this) is carried into the
  effects by `interop_joint_from_outer()`, `Qrr^-1 + J Vf J'`, the formula
  `joint_cov_repair()` uses for the model's own `Vf`. A matrix over
  `get_vcov()`'s rows is read as given. A matrix over the coefficients alone
  has nothing to carry them from: a warning of class `frmtmb_interop_se_part`
  (`interop_user_re()`).
- `interop_warn_se_part()` now names only what no parameter holds: a new row's
  `me()` draw, a new row's latent residual, and the effects of a fit whose
  objective integrates them out (or whose joint covariance failed).
- `interop_rows_arm()` / `interop_rows_check()`: `get_vcov()` arms and keeps
  the estimates; the first `get_predict()` at other estimates (a Jacobian
  call, fix round) checks once, and warns (same class) with a count when rows
  read a grouping level the fit did not see or an exact `gp()` at an unseen
  position with a kriging variance above `1e-4` of the field's
  (`interop_krig_share`).

No shared SE builder changed. `R/predict.R`, `R/autocor-latent.R` are
untouched; the lane reads `get_joint_cov()`, `joint_precision()`,
`b_coef_labels()`, `lp_eta_design()` and `lp_extra_var()` as they are.

## The off-rate on the 168 fits, before and after

`dev/mfxse-rate.R` is the review's `dev/review-predread-r2-mfxrate.R` on this
lane's library setup (the review's copy sourced a worktree that is gone), with
`fitted()`'s `Est.Error` and the warnings kept beside each row. Same designs:
bernoulli, gaussian, lognormal, poisson x 14 designs (8 or 30 groups, 10 or 40
rows, group SD 0.2 or 0.8, `(1 | g)` or `(1 + x | g)`, without 8 x 40 at 0.2)
x data seeds 7501 to 7503; reference = SD of 150 draws (seed = data seed +
1000) of the whole joint covariance through the same quantity; "off" is a
ratio outside 0.8 to 1.25. Summaries by `dev/mfxse-ratesum.R <arm>`, logs
`dev/mfxse-rate-out/sum-before.log` and `sum-after.log` (local).

Before (rellib-r18, 0.79.1); this reproduces the review's table to the digit:

    By quantity:
          qty ratio.n ratio.off ratio.low ratio.high ratio.med ratio.min ratio.max
     avg_pred     168      0.77      0.00       0.77      2.03      0.92     11.86
        cmp_f     168      0.22      0.00       0.22      1.05      0.87      2.92
        cmp_w     168      0.02      0.01       0.01      1.01      0.56      1.80
         pred    2016      0.67      0.41       0.26      0.91      0.18      9.35
      slope_x     168      0.52      0.00       0.52      1.30      0.86      9.86
    A fit with any material quantity, share: any_avg 0.798, any_pred 0.905
    Levels against fitted()'s Est.Error: n 2016 min 0.264269 median 0.934660
      max 9.960300 outside 1e-3: 1799
    Warnings: frmtmb_interop_se_part on 168 of 168 fits

After (wt-mfxse-lib):

    By quantity:
          qty ratio.n ratio.off ratio.low ratio.high ratio.med ratio.min ratio.max
     avg_pred     168      0.00      0.00          0      1.00      0.84      1.19
        cmp_f     168      0.00      0.00          0      1.00      0.85      1.18
        cmp_w     168      0.02      0.02          0      0.97      0.66      1.16
         pred    2016      0.01      0.01          0      1.00      0.25      1.55
      slope_x     168      0.00      0.00          0      1.00      0.86      1.18
    A fit with any material quantity, share: any_avg 0.024, any_pred 0.077
    Levels against fitted()'s Est.Error: n 2016 min 0.999998 median 1.000000
      max 1.000004 outside 1e-3: 0
    Warnings: none on 168 of 168 fits
    Convergence codes (both arms): 0 on 149, 1 on 19

The material rows left after (`dev/mfxse-rateoff.R after`, `off-after.log`):
24 rows on 16 fits, 4 of them averages (all `cmp_w`, `avg_comparisons()` of
the group-level covariate `w + 1`, at 8 groups: bernoulli 7502 with a slope,
lognormal 7501 with and without, poisson 7501), 20 of them levels. On every
material level `fitted()`'s ratio to the draws is the same number
(`mfx_fit` = 1): bernoulli levels near p = 0.93, poisson levels at sd 0.8.
So what is left is the delta method against simulated draws of a nonlinear
response, which `fitted()` shares, not a gap in the interop path. The 4 fits
are why the help says "on 164 of the fits".

## brms 2.23.0 on a subset

`dev/mfxse-brms.R`, log `dev/mfxse-rate-out/brms.log`: 5 of the 168 designs
at seed 7501 (poisson 30 x 10 at sd 0.8 with and without a slope, gaussian
30 x 10, bernoulli 30 x 40, lognormal 30 x 10, all `(1 | g)` but the one
slope fit), brms default priors, 2 chains x 2000, max Rhat 1.006 to 1.032.
brms's posterior SD of each quantity from `posterior_epred()` (re_formula
NULL), against each arm's marginaleffects SE:

| family, design | quantity | brms SD | before | after |
|---|---|---|---|---|
| poisson 30x10 | avg_pred | 0.1010 | 0.4320 | 0.1003 |
| poisson 30x10 | slope_x | 0.1043 | 0.2436 | 0.1076 |
| poisson 30x10, slope | avg_pred | 0.1102 | 0.5001 | 0.1082 |
| poisson 30x10, slope | slope_x | 0.1336 | 0.5574 | 0.1293 |
| gaussian 30x10 | avg_pred | 0.04459 | 0.14982 | 0.04479 |
| bernoulli 30x40 | avg_pred | 0.01254 | 0.02931 | 0.01297 |
| lognormal 30x10 | avg_pred | 0.1533 | 0.4962 | 0.1466 |
| lognormal 30x10 | slope_x | 0.1391 | 0.2643 | 0.1346 |

All 20 non-level quantities: after / brms median 0.977, range 0.872 to 1.104;
before / brms median 1.246, range 0.936 to 4.537. Levels (12 per fit), after /
brms medians 1.002, 1.000, 0.985, 1.019, 0.977 (ranges within 0.949 to 1.045);
before medians 0.839, 0.727, 0.874, 0.684, 1.214 (ranges 0.399 to 1.817).

## Structures, readers and the warned cases

`dev/mfxse-cmp.R` (seeds 7411, 8150; `mfxse-cmp-{before,after}.log`):
`predictions()` SE / `fitted()` Est.Error on 4 rows each.

| fit | before | after |
|---|---|---|
| poisson `(1 | g)` | 1.5294 2.9437 0.8681 0.8844 | 1.0000 x 4 |
| poisson `(1 + x | g)` | 1.2639 2.8162 0.8534 0.6518 | 1.0000 x 4 |
| gaussian `(1 | g)` | 2.0146 1.7129 1.8798 1.7071 | 1.0000 x 4 |
| gaussian REML | 2.0984 1.7868 1.9590 1.7722 | 1.0000 x 4 |
| lognormal `(1 | g)` | 2.2357 1.9152 2.0892 1.8824 | 1.0000 x 4 |
| bernoulli `(1 | g)` | 0.7738 0.4656 0.6856 0.8159 | 1.0000 x 4 |
| `s(x)` | 6.8613 0.5111 6.5458 6.8303 | 1.0000 x 4 |
| `gp(x, k = 15)` | 7.5893 9.2184 7.8929 2.4689 | 1.0000 x 4 |
| poisson, no group | 1.0000 x 4 | 1.0000 x 4 |
| exact `gp()`, 3 rows | 6.3337 8.2422 1.0273 | 1.0000 1.0000 1.0010 |

`re_formula = NA` reads 1.0000 on both arms for the group-level fits; `s(x)`
read 6.86 0.51 6.55 6.83 before (a smooth is not dropped by `re_formula`).
`hypotheses(fit)` on poisson `(1 + x | g)` reads 0.22110271 and 0.04996868 on
both arms, `sqrt(diag(vcov(fit)))`.

The exact `gp()` row at x = 7 (1.0010): `dev/mfxse-probe5.R` takes a central
difference with a fresh cache as referee. Central difference at h = 1e-3,
1e-4, 1e-5 reads 0.533873 every time; `frm_lp_basis()`'s `A V A'` 0.533872;
marginaleffects' default forward difference 0.534407. At x = 12, 0.550415,
0.550414, 0.550973. So the 1e-3 is marginaleffects' forward difference, and
`numderiv = "fdcenter"` is the remedy; the help says so.

`dev/mfxse-structs.R` (seed 9101; `mfxse-structs-{before,after}.log`), 5 rows:
after, crossed `(1 | g) + (1 | h)`, nested `(1 | h / g)`, `(x || g)`, `mm()`,
`t2()`, `gp(x, k = 10) + (1 | g)`, `(1 | p | g)` in `mu` and `sigma`,
zero-inflated poisson with `zi ~ (1 | g)`, `mvbf` with `(1 | g)` in each,
`cumulative()` with `(1 | g)` (20 category cells): every ratio 1.0000. Before:
1.65 to 1.84, 0.69 to 0.84, 1.06 to 1.22, 1.15 to 1.30, `t2()` 1.0000 (fully
unpenalized at k = 4), 7.70 to 10.61, 1.06 to 1.35, 0.76 to 1.00, 1.06 to
1.22, 0.56 to 1.00. The nonlinear body `a ~ 1 + (1 | g)` reads NA because
`fitted()`'s Est.Error is NA on a body (filed at 0.79.0, not this lane); the
marginaleffects SE there is finite. `gr(g, by = f)` did not fit on this
design (f varies within g); not measured.

`dev/mfxse-ext.R` (seeds 9101, 31, 7411; `mfxse-ext-{before,after}.log`),
with frmtmb.spatial from rellib-r18 on the lane's core:

| fit | before | after | warning after |
|---|---|---|---|
| `car(type = "esicar")` | 0.3360 0.3356 0.4230 0.4452 | 1.0000 x 4 | none |
| `car(type = "icar")` | 0.3361 0.3358 0.4232 0.4453 | 1.0000 x 4 | none |
| `car(type = "escar")` | 1.8602 1.8420 1.8388 1.7728 | 1.0000 x 4 | none |
| latent residual + `(1 | s)`, new rows | 0.4883 0.4896 | 0.4408 0.4548 | get_vcov(), latent residual |
| `me(xo, sx) + (1 | g)`, new rows | 0.7347 0.7572 | 0.7573 0.7479 | get_vcov(), me() draw |
| `(1 | g)`, `quadrature = TRUE` | 1.0000 x 3 | 1.0000 x 3 | get_vcov(), integrates out |
| `(1 | g)`, Laplace | 0.9953 1.3738 1.0362 | 1.0000 x 3 | none |
| a new level, `allow_new_levels` | 0.9920 0.2859 | 1.0000 0.2859 | get_predict(), 1 row |

Under quadrature `fitted()` is itself conditional on the modes and says so
(`frmtmb_modes_conditional_se`), so the two agree and both warn.

## The gaussian-field avg_slopes() (MINOR)

`dev/mfxse-gpslope.R` (data seed 8150, draws 8151;
`mfxse-gpslope-{before,after}.log`); "basis" is the delta method by hand
through `frm_lp_basis()` at x +- h/2:

| fit | before se | after se | basis | 150 draws |
|---|---|---|---|---|
| `y ~ z + gp(x)` | 3.08e-07 | 0.0308 | 0.03079 | 0.06699 |
| `y ~ gp(x, gr = FALSE)` | 1.27e-06 | 0.0516 | 0.05161 | 0.06578 |
| `y ~ gp(x, k = 15)` | NA | 0.0516 | 0.05162 | 0.04997 |

After, marginaleffects is frmtmb's own delta method. WITHDRAWN at the fix
round: the "150 draws" column is not a reference on an exact `gp()`. The
delta method is invariant to centering the field (the `Qrr^-1` and `J`
blocks transform with the map); a Gaussian draw in the non-centered
`(theta, z)` coordinates is not, so 0.067 measures the coordinates, not the
uncertainty (review m6). The draws column stays as measured; the claim that
it is a reference does not.

brms 2.23.0 on the same data (`dev/mfxse-brms-gp.R`, `mfxse-brms-gp.log`,
2 chains x 2000, default priors, Rhat 1.008): posterior SD of the average
slope 0.15409 on `y ~ z + gp(x)` and 0.16135 on `gp(x, gr = FALSE)` (means
-0.14949, -0.15303 against frmtmb's -0.15657, -0.15661). So on an exact
`gp()` the delta method, frmtmb's own as much as marginaleffects', reads a
fifth of brms's posterior SD. The reviewer's profile of the conditional mode
over log `lscale` (inner problem re-solved) explains it: the average slope
moves by -0.0051, -0.0005, -0.0055, -0.0086 at -1, -0.5, +0.5, +1 SE, then
-0.0558 at +2 and -0.1546 at +3 SE; flat near the mode and steep in the upper
tail, which a posterior with mass there (brms's priors on `lscale`) widens
and a linearization cannot see. It is not the interop gap, which was a
factor of 1e5. Filed below, not fixed, no warning (one on every field slope
would be noise).

The field's conditional (kriging) variance at the unseen positions x +- h/2
is not in the slope (`dev/mfxse-gpkrig.R`, `mfxse-gpkrig.log`, `w' C w` with
`C = gp_krig_cov()`): kriging SD 0.00001, 0.00013, 0.00128 at h = 1e-2, 1e-3,
1e-4 against the coefficient SE 0.03082, total 0.03082, 0.03082, 0.03085. It
grows as 1/h because of the 1e-12 relative nugget, and at marginaleffects'
step it moves the third significant digit at most. This is the measurement
behind `interop_krig_share`: a slope does not warn.

## Time

`dev/mfxse-time.R`, data seed 7411, poisson `(1 + x | g)`, 40 rows per
level, minimum of 3 rounds, with a fixed control (20 solves of 300 x 300)
that read 0.22 and 0.23 s before, 0.24 and 0.21 s after
(`mfxse-time-{before,after}.log`):

| levels, rows, get_coef() length | call | before s | after s |
|---|---|---|---|
| 30, 1200, 2 -> 62 | predictions(), 12 rows | 0.03 | 0.19 |
| | avg_predictions() | 0.02 | 0.14 |
| | avg_slopes() | 0.03 | 0.30 |
| | avg_slopes(vcov = FALSE) | 0.01 | 0.03 |
| 500, 20000, 2 -> 1002 | predictions(), 12 rows | 0.02 | 2.83 |
| | avg_predictions() | 0.05 | 8.28 |
| | avg_slopes() | 0.09 | 16.26 |
| | avg_slopes(vcov = FALSE) | 0.04 | 0.03 |

The cost is marginaleffects' Jacobian, one `get_predict()` per parameter:
8.8 ms per `get_predict()` on 20000 rows (`dev/mfxse-prof.R`), of which the
design rebuild (`lp_eta_design()`, `model.matrix()`) is about 40 percent of
the samples. A design memo across perturbations would cut that share at
most; the count of predictions is inherent. In the 168-fit run the five
calls took a median 0.42 s per fit at 8 groups and 1.11 s at 30 (0.28 s and
0.30 s before). Stated in the help and NEWS, with `vcov = FALSE` and
`re_formula = NA` as the ways around it.

## Tests

New `tests/testthat/test-marginaleffects-re.R` (10 tests, 39 expectations).
Seen to fail on rellib-r18 (`MFX_ARM=before dev/mfxse-run-test.R`, log
`dev/mfxse-tests/re-before.log`): `RESULT pass=21 fail=18 error=0 skip=0
warn=18`. The behavioral failures, among them: predictions vs fitted
`max(abs(p/ft - 1))` 5.694 against a bound of 0.006; the average and slope
"1.00 >= 0.01" (the old SE IS the old convention); `s(x)` SE 0.479 0.030
0.490 1.215 against fitted 0.076 0.063 0.076 0.200; `gp(x, k = 15)` 0.2570 on
every row against 0.0692 0.0559 0.0698 0.2536; the gaussian field's
avg_slopes SE 3e-7 against a lower bound of 0.0154. On the lane's build:
`pass=39 fail=0 error=0 skip=0 warn=0` (`re-after-2.log`).

Changed: `test-me-newdata.R` (the field and group-level test now asserts
agreement with `fitted()` and no warning; the "either direction" test reads a
`me()` fit, gains `skip_if_not_installed("marginaleffects")`, and asserts a
group-level fit no longer warns), `test-v07.R` (glmmTMB's SE is matched
through the coefficient block as `vcov =`, the convention glmmTMB uses; the
default now differs from it), `helper-warnings.R` (the pin `mfx_se_part`
follows the new message).

Whole suites, one file per process, lane library first then rellib-r18
(`dev/mfxse-run-par.sh`, a copy of `dev/release/run-par.sh` writing to
`dev/` and reading rellib-r18; jobs `dev/mfxse-full.jobs`, split into
`mfxse-full` and `mfxse-full2` when the driver was restarted from 2 to 3
processes; totals summed from every RESULT line, `dev/mfxse-full-all.txt`):
500 files (core 309, coupling 11, eam 33, latent 12, learn 16, ode 11,
sample 85, spatial 8, survival 15), pass 31652, fail 0, err 0, warn 0,
skip 289 in 55 files (the gated brms tiers; FRMTMB_BRMS_FIT_TESTS unset).
All 500 `lib:` lines read `frmtmb 0.79.1 .../wt-mfxse-lib`.

After the last edit to `R/compat.R` (a sentence in the `gp_pred()` note) the
library was reinstalled and the compat and changed files rerun
(`mfxse-final.log`): test-compat-extensions 2, test-compat-register 103,
test-compat 759, test-marginaleffects-re 39, test-me-newdata 181, test-v07
27, test-marginaleffects-args 57 passes, every one fail 0 err 0 skip 0
warn 0.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`, build
log `mfxse-openblas-build.log`: "Rblas: 74 to OpenBLAS", "Rlapack: 616 to
OpenBLAS"; runs `dev/mfxse-tests/openblas-*.log`): test-marginaleffects-re
39, test-marginaleffects-args 57, test-v07 27 passes, each fail 0 error 0
skip 0 warn 0. `dev/ciharden-out/` deleted afterwards.

One session per suite, three orders (`dev/mfxse-run-sessions.sh mfxse-sess
3`, a copy of `dev/release/run-sessions.sh` on rellib-r18 behind the lane's
library; log `dev/mfxse-sess.log`): 27 of 27 sessions, every one fail 0 err
0 warn 0. Core: forward pass 22555, reverse 22705, preload 22798, each
files=296/309 (the 13 others are the gated brms files, which skip at file
level and leave no row). coupling 551 x 3, eam 1848 x 3, latent 683 x 3,
learn 474 x 3, ode 549 x 3, sample 3962 x 3 (files=80/85, gated), survival
552 x 3, spatial 478 x 3. 46 LEAK lines, every one tagged `[with load of
...]`; no untagged LEAK line.

## Equivalence survey and brms-port ledger

`dev/equiv-data/results.tsv` holds no pair whose feature is
`marginaleffects` (grep), so no sampled pair reaches this change and none
was remeasured. `dev/brmsport-verdicts.tsv` holds no marginaleffects row;
no row flips.

## Compat rows changed (R/compat.R, `compat_hand_rules_tbl()`)

- `kind:covstruct` x `marginaleffects`: conditional -> works.
- `s()` x `marginaleffects`: conditional -> works.
- `t2()` x `marginaleffects`: conditional -> works.
- `gp_pred()` x `marginaleffects`: stays conditional, note rewritten (an
  unseen position's kriging variance; the delta method's distance from brms
  on an exact field).
- `me()` x `marginaleffects`: unchanged.

## What I did not do, and why

- A test of `car()` in frmtmb.spatial. Measured by `dev/mfxse-ext.R` (1.0000
  on esicar, icar, escar), but frmtmb.spatial does not suggest
  marginaleffects, and its CI skip gate fails a skipped
  `skip_if_not_installed()`. Adding marginaleffects to its Suggests is the
  user's call; no extension file changed, so no extension check ran.
- The kriging variance at an unseen exact `gp()` position, a new level's
  variance, a new row's `me()` draw and latent residual: no parameter holds
  them, so `J V J'` cannot. Each warns by name; counted in this run: the
  168 fits raised none; `dev/mfxse-ext.R` raised the latent, `me()`,
  quadrature and new-level warnings once each where they apply.
- The cost at many levels. marginaleffects fixes the Jacobian at one
  prediction per parameter; frmtmb cannot hand it a factor or a batched
  Jacobian through any exported hook (`get_jacobian_analytic()` is
  internal to marginaleffects).
- (Superseded at the fix round.) `hypotheses()` on the model keyed off
  marginaleffects passing `hypothesis_is_formula` to `get_coef()`; the joint
  test never passes it (review M1). It now keys off the frame of the
  exported `hypotheses()` (below).

## Defects found, not fixed

- On an exact `gp()` the delta method of an average slope is far below
  brms's posterior SD (0.0308 against 0.154 on `y ~ z + gp(x)`, 0.0516
  against 0.161 on `gr = FALSE`). This is frmtmb's delta method on a field
  whose slope is flat in the length scale near the mode and steep two SEs
  above it (review m6's profile), shared by `frm_lp_basis()` and
  `emtrends()`, not the interop path. MINOR, no warning.
- The nonlinear body `a ~ 1 + (1 | g)`: `fitted()`'s Est.Error is NA (filed
  at 0.79.0); marginaleffects now gives a finite SE there, which nothing in
  frmtmb checks.
- Far from its data an exact `gp()`'s finite-difference Jacobian parts from
  `fitted()`'s exact one (review m5, `dev/review-mfxse-out/follow.log`,
  `gp(x) + (1 | g)`, x = 7, 8, 9, 10, 12): forward 0.9996, 1.0217, 1.0749,
  1.0232, 0.6265 times `fitted()`; `fdcenter` 1.0006, 1.0081, 1.2061,
  1.1244, 0.6441, no closer. My first help said "third digit" and named
  `fdcenter` as the remedy; both corrected. Rows 8 to 12 warn by count.

## Fix round (review dev/reviews/2026-11-06-mfxse.md: FIX, 3 MAJOR, 9 MINOR)

### M1, `hypotheses(joint = TRUE)` tested the group-level effects
Signal searched in marginaleffects 1.0.0. `joint_test()` calls `get_coef(object)`
and `get_vcov(object, vcov = vcov)` with no other argument
(`dev/review-mfxse-out/mesrc/joint_test.R`). `internal_call` reaches
`get_coef()` only from `avg_slopes()`; `predictions()` and
`avg_comparisons()` pass none (`dev/mfxse-probe2.R`, my first-round trace:
`get_coef [ variables,numderiv ]` for `avg_comparisons()`, `[ ]` for
`predictions()`), so it cannot key the effects in. No sound key exists in the
arguments. What is sound: `hypotheses()` is exported, and `model` is its
documented argument. `interop_in_hypotheses()` walks the frames for the
function identical to `getExportedValue("marginaleffects", "hypotheses")`
whose `model` is a `frmtmb_fit`; there `get_coef()` and `get_vcov()` give the
coefficients alone, no warning, no arming. A test of the model's own
coefficients is the only thing `hypotheses()` on a raw fit computes.
Reviewer's `dev/review-mfxse-hyp.R` on the fix build
(`dev/mfxse-fix-out/hyp.log`): `joint = TRUE` F 15.4386, df1 3, df2 234 (the
rellib-r18 value; the first lane build read F 328.3 on 27); `joint = c(2,3)`
21.9391 df1 2; `joint = "x"` 21.8519 df1 1; `"b2 = 0"`, `"b2 = b3"`,
`~ reference`, `~ pairwise`, `vcov = ~ g` unchanged against `fixef()`.
Test: `test-marginaleffects-re.R`, the joint statistic against
`t(b) solve(vcov(fit), b) / 3` by hand, df1 3, and a non-name pattern
`"^b[.]"` refused by marginaleffects.

### M2, the new-level warning was muffled in slopes and contrasts
`interop_rows_arm()` keeps the estimates; `interop_rows_check()` runs at the
first `get_predict()` whose `estimates` are not identical to them, a Jacobian
call, never the estimate pass that `get_predict_error()` muffles. Reviewer's
`dev/review-mfxse-vcovwarn.R` on the fix build (`mfxse-fix-out/vcovwarn.log`):
`predictions`, `avg_predictions`, `avg_slopes(x)`, `avg_comparisons(f)` with
two new-level rows: 1 warning each reaching the caller (the first build: 1,
1, 0, 0). Seen levels, `re_formula = NA`, `vcov = FALSE`: 0. After
`get_vcov(fit)` the flag stays armed, but a following `predictions(vcov =
FALSE)` at the estimates gives 0 warnings (it never perturbs). Test: the
`avg_slopes()` and `avg_comparisons()` rows require the count warning.

The kriging count beside a new level (review m5): it was counted; the
reviewer's log truncated the message at 80 characters
(`dev/mfxse-krigcount3.R`, the whole message: "1 read a grouping level ...,
and 1 read an exact gp() at a position the fit did not see"; the design's
kriging share at x = 9 is 4.95e-02, `dev/mfxse-krigcount.R`). Now pinned by
a test that requires both counts.

### M3, tests that catch `oneblock` and `sandown`
New tests in `test-marginaleffects-re.R`: crossed `(1 | g) + (1 | h)`
poisson and `s(x) + (1 | g)`, levels against `fitted()` and
`avg_predictions()` against `mre_avg_se()`, seed 9102; and a cluster
sandwich against the joint covariance built by hand (`Qrr^-1 + J Vf J'`, Vf
from `vcov_cluster(~ g, "CR1", full = TRUE)`, read through
`frm_lp_basis()`), seed 9301 heteroskedastic. Run against the reviewer's
mutant builds (`dev/review-mfxse-run-test.R`, `dev/mfxse-tests/mut-*.log`):
`oneblock` pass 46 fail 8, of which the two-block test 3 (levels
max|p/ft - 1| 0.1315 against a bound of 0.0002 and 4.105 against 0.004;
average 1.14 against 0.03); `sandown` pass 48 fail 6, of which the sandwich
test 1 (max|ps/ref - 1| 0.060 against 0.001). The other failures on both
mutants are the M1 and M2 tests, which those builds (made from the first
lane source) predate. On the fix build: pass 55, fail 0.

### m1, the dense inverse in `interop_joint_from_outer()`
Measured (`dev/mfxse-cost-outer.R`, seed 7411, 500 levels, 1000 inner rows,
`mfxse-fix-out/cost-outer.log`): `Matrix::solve()` of the sparse block 0.00
to 0.02 s, a dense Cholesky inverse 0.14 to 0.16 s; `predictions(2 rows,
vcov = ~ g)` 19.47 s, of which `vcov_cluster(~ g, full = TRUE)` 15.84 s and
`get_vcov(vcov = ~ g)` 16.19 s; `predictions(2 rows)` default 4.14 s. So
the reviewer's 21 s is the cluster sandwich (core code older than this lane)
and the Jacobian, not the inverse. Kept `Matrix::solve()`, added the guard:
above `frmtmb.latres_dense_rows` inner rows the carry is not formed, and a
warning of class `frmtmb_interop_se_part` says the given covariance could
not be carried (test with the option at 5). latperf's `frm_joint_sparse`
covers the model's own joint covariance on latent fits, not a user `Vf`;
not merged. No `get_predict()` memo: it would save about a third per
parameter (review m1) at the cost of a cache keyed on which `b` moved, not
cheap and safe enough for this round.

### m2, m4, m5, m6, m9
- m2: help and NEWS no longer say `"HC0"`, `"HC1"` reach a mixed fit; they
  say only the outer part becomes robust and the conditional part stays
  model-based.
- m4: latent `b` stay out of `get_coef()`, with the reason corrected (one
  per fitted row; a new row under `cov = TRUE` does read them, through its
  conditional mean; the review measured 0.989 to 1.064) and the warning
  clause naming "the fitted latent values it continues from".
- m5: help paragraph replaced with the measured forward and central ratios
  and the `eps` sensitivity (0.0308 at the default and at 1e-3, 0.0328 at
  1e-4, 0.1124 at 1e-5); the row warning adds, when it counts a kriging row,
  that the finite-difference Jacobian is unreliable there.
- m6: the 0.067 draws claim withdrawn above; profile explanation recorded.
- m9: NEWS, `kind:covstruct` compat note (joint, slopes and contrasts warn,
  sandwich outer-only, HC0/HC1 refused), help.
- m3 left as is, a question for the user.

### m7, the merge with lane readers (not done here; for consolidation)
After the merge, `get_predict()` at a new level runs `fitted_me_moments()`'s
Gauss-Hermite integral over the level's effect (readers' `fitted()` on 3
new-level rows: 6.0 ms before, 189.0 ms after, review m7). Two consequences
for this lane's path:
1. The estimate at a new level then reads the block's `theta` (the effect's
   SD), which `get_coef()` holds only for `gp` and `rr`. The user's dots reach
   `get_coef()` in the Jacobian call (my first-round trace: `re_formula`
   arrives there), so consolidation should add, in `interop_re_part()`, the
   `theta` of every non-field block when `isTRUE(dots$allow_new_levels)`:
   a few parameters, read only then. The new-level warning stays true (the
   draw's own variance is still in no parameter) and must keep this round's
   Jacobian-call arming.
2. Cost: one integral per `get_predict()`, one call per parameter: at 62
   parameters about 12 s per call, at 1002 about 190 s (review's arithmetic,
   not run). If that is too slow, `interop_rows_check()` already knows the
   new-level rows; a cheaper course is to integrate only those rows, which
   readers' `read_mean()` would have to accept. Do not edit `predict.R` here.

### Reruns on the fix build (one file per process, `dev/mfxse-fix.log`)
Every core file that reads marginaleffects, 17 files, all `lib:` lines on
`wt-mfxse-lib`: test-arg-refusal 122, test-autocor-sigma 79,
test-brms-shapes-punch 55, test-compat 759, test-famgaps 64,
test-latres-newrow 37, test-latres 298, test-linpred-layers 83,
test-marginaleffects-args 57, test-marginaleffects-re 55, test-me-newdata
180, test-mo-coding 164, test-ordinal-fitted 129, test-resp-read 111,
test-review-v29 186, test-se-reparam 76, test-v07 27 passes; every one fail
0, err 0, skip 0, warn 0. `test-me-newdata.R` lost one expectation: a
coefficient test no longer warns, so its `require` of the warning went.

The 168-fit table was not rerun: the default standard-error path is
unchanged (the fix touches the `hypotheses()` context, the user `vcov =`
carry, the arming of the row check and messages), and none of the 168
designs has a new level or an unseen `gp()` position.

## R CMD check --as-cran (core)

`dev/mfxse-check.sh` (built and checked in `dev/mfxse-check/`, since deleted;
`00check.log` kept as `dev/mfxse-check-00check.log`, local): `Status: 2
NOTEs`. The examples-timing NOTE names `frm_hazard_reads` at 6.26 s elapsed
(2.86 s user), not an example this lane touched, on a machine running 12
lanes; the other is the V8 math-rendering NOTE. Tests under check: `[ FAIL 0
| WARN 0 | SKIP 270 | PASS 22308 ]`, 51 minutes. Vignettes rebuilt OK.
Extensions: none changed, so none checked.

## Short check of the fix round (FIX, one MAJOR)

### S1, `hypotheses(fit, hypothesis = <function>)` held the effects fixed
marginaleffects hands a function `hypothesis` the model and differentiates
what it returns, so it is not a test of the coefficients. `interop_in_hypotheses()`
now reads `hypothesis` from the `hypotheses()` frame whose `model` is a fit
and returns FALSE when it is a function. Reviewer's
`dev/review-mfxse-sc-hyp.R`, seed 7411 (`(1 + x | g)` poisson), before
(`dev/review-mfxse-out/sc/hyp.log`, the fix-round build) and after
(`dev/mfxse-fix-out/sc-hyp-after.log`):

| call | before | after | the call alone |
|---|---|---|---|
| `hypotheses(fit, function(x) avg_predictions(x))` | 1.77417 | 0.14299 | 0.14299 |
| `hypotheses(fit, function(x) predictions(x, nd))` | 1.76369 1.91840 17.14256 | 0.55222 0.80605 2.56104 | 0.55222 0.80605 2.56104 |

The joint test stays F 15.4386 df1 3 by all 8 routes of that script. New test
"hypotheses(fit, <function>) holds the group-level effects" (reviewer's design,
`mre_data()`): on the build with the old condition it fails, 11.41 against a
bound of 0.01 (`dev/mfxse-tests/re-sc-old.log`); on the fix build it passes.

### s1, a user `vcov =` never armed the row check
`get_vcov()` now arms the check on the user branch too (not inside a
coefficient test). Test: `predictions(vcov = ~ g)` and
`predictions(vcov = vcov(fit, full = TRUE))` at two new-level rows each
require the count warning; both failed on the old build
(`re-sc-old.log`, 2 failures "expected a warning containing '2 read a
grouping level ...'") and pass now. The reviewer's M2 rows on the fix build:
1 warning each from predictions, avg_predictions, avg_slopes,
avg_comparisons, slopes and comparisons at new levels; 0 at seen levels.

### s2, the option's wider meaning
`?frmtmb-autocor` (`R/autocor.R`, the paragraph on
`frmtmb.latres_dense_rows`) now says the option also bounds the dense inner
block a carried `vcov =` inverts on any fit with group-level effects, latent
or not, and that past it the effects are held with a warning. `R/autocor.R`
is outside this lane's file; one paragraph added, nothing else touched.

### Reruns (one file per process, `dev/mfxse-sc.log`, all `lib:` on wt-mfxse-lib)
test-marginaleffects-re 59, test-marginaleffects-args 57, test-me-newdata
180, test-compat 759 passes; each fail 0, err 0, skip 0, warn 0.
