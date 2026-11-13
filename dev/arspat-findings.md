# Lane arspat: residual correlation and spatial refusals

Round opened 2026-11-11, base e8abed0a (frmtmb 0.81.0), branch
wt-r82-arspat. Brief: `round-briefs/r082-arspat.md`. Status: DONE
(items 1 to 4 and 3897, 4110, 4115, 4120, 4320; not done: 3901, the
`ar(cov = FALSE)` and `newdata2` halves of 3903; one question for the
user, section 4).

## What brms 2.23.0 writes (run, not recalled)

`dev/arspat-brms-probe.R`, log `dev/arspat-out/brms-probe.txt`
(`make_stancode()` only, brms 2.23.0 from rellib-r20):

- `y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)`: fits,
  `normal_time_hom_se_lpdf()`, which adds `diag(se2)` to each group's
  `Cov = (sigma L)(sigma L)'`. The same for `ma()`, `cosy()`, `unstr()`
  (`_flex`), `sigma ~ x` (`_het_se`) and `student()`
  (`student_t_time_hom_se_lpdf()`).
- `y | se(s) ~ ... ar(cov = TRUE)` (sigma = FALSE): refused, "Please set
  argument 'sigma' of function 'se' to TRUE when modeling time-series
  covariance matrices."
- `ym | mi() ~ x + ar(t, g, cov = TRUE)`: fits, `Yl[Jmi] = Ymi` and
  `normal_time_hom_lpdf(Yl | ...)`. The same for `cosy()` and `student()`.
  `mi(s)` (measurement error) puts every row's `Yl` latent.
- `fcor(M)` with `mi()` fits; with `se()` it is refused ("Invalid addition
  arguments for this model").
- `sar()` with `mi()` fits (error and lag); with `ar(cov = FALSE)` it fits,
  `mu[n] += Err[n, 1:Kar] * ar` and then the SAR density of `Y` at that
  `mu`; with `ar(cov = TRUE)` it is refused ("Can only model one covariance
  matrix of natural residuals"); with `se()` it writes
  `sqrt(square(sigma) + se2)` into a `real sigma` slot, which stanc
  rejects (brms-38 in `dev/upstream-bugs.md`).
- `ar(t, g, cov = TRUE)` on groups with an interior gap: no `Jtime_tg`
  in the call (`normal_time_hom_se_lpdf(Y | ..., nobs_tg, begin_tg,
  end_tg)`), so brms indexes by position within the group. Only `unstr()`
  passes `Jtime_tg` (the `_flex` functions).

## What changed

### 1. `se()` beside the covariance form (backlog 4108 = 4172)

`check_autocor_response()` (R/autocor.R) no longer refuses `se()` for the
ar() family; it refuses `se()` without `sigma = TRUE` in brms's words, and
still refuses `se()` beside a registered term (`fcor()`: brms refuses;
`sar()`: brms's code does not compile). The objective
(`autocor_loglik_se()`, R/autocor.R, called from R/objective.R) takes
each group's covariance `D R D + diag(se^2)` and factorizes it per group,
since the added diagonal differs by row. `simulate()` and `predict()` add
an independent `N(0, se^2)` to the correlated draw
(`autocor_draw_resid_se()`; under `student()` one chi-square scales both,
as brms's `multi_student_t` scale). The pointwise log-density
(`autocor_group_loo_parts()`, R/autocor-loglik.R) inverts `D R D +
diag(se^2)` per group, as brms's `get_cov_matrix_ac()` adds `se2`.

### 2. `mi()` beside the covariance form and registered terms (4108)

- The ar() family: a missing value that no other density reads is taken
  out of the block's patterns (`autocor_drop_rows()`, applied in
  `build_objective()` beside `mi_gone`). The density of the rows left is
  the margin of the joint one, which is the same block over the time points
  the group keeps, because lags are read off the global time levels. This
  is exact on `student()` too (a margin of a multivariate t is a
  multivariate t with the same nu). A value another density reads stays
  latent and enters the joint draw through `mivals`.
- `fcor()` and registered terms (`sar()`): `mi_rows_joint()` (R/frame.R)
  now holds every missing row of such a response latent, since the map
  takes every row at once, and the seam density reads `mivals`. Laplace,
  exact on gaussian. (First pass: "approximate on student". WITHDRAWN in
  the fix round: the student error is the data-free constant `c(nu, n, m)`,
  which the objective now takes off; see "Fix round", MAJOR-1.)

### 3. Cluster sandwich over the series (4299)

`cluster_guard()` (R/sandwich.R) refuses `fcor()` and registered terms
(one density over all rows) and, for the ar() family under either form,
refuses only a clustering that splits a series (`autocor_row_groups()`).
The cluster weight enters the covariance form once per group
(`autocor_loglik(w =)`, `mvt_std_loglik(w =)`, `autocor_loglik_se(w =)`);
the `cov = FALSE` form already took it per row, which is right once a
series is inside one cluster. The HAC refusal text (4320) now says what is
missing (one time order over all rows) and names `~ g`.

### 4. Interior time gap (4104)

brms counts rows within a group under `ar()`, `ma()` and `arma()` (code
above). frmtmb keeps nlme's reading (lags from the time levels, as the gls
agreement tests pin) and now says so: `autocor_note_gaps()` gives a
message of class `frmtmb_autocor_gap` when a group skips a level that other
groups have. A message, not a warning, because the levels reading is a
correct model. Not for `cosy()` (one correlation for every pair) or
`unstr()` (brms reads the time index too), nor under `cov = FALSE` (counts
rows, as brms).

QUESTION FOR THE USER (recommendation: keep the levels reading and the
message). Matching brms would mean counting rows: a BREAKING change of the
likelihood on ragged data, which would break the nlme `gls()` agreement
the validation suite pins, and would treat a missing time point as no gap,
which is wrong for a process in time. The message makes the difference
visible to a ported brms script, and the help says how to get brms's model
(the within-group position as the time variable).

### 5. Same cluster

- 3897 (frmtmb.spatial): `sar()` takes a `W` with complex eigenvalues.
  `sar_build()` keeps the complex ones apart, `sar_logdet()` adds
  `0.5 log((1 - rho Re)^2 + (rho Im)^2)` per eigenvalue, and rho's interval
  comes from the real eigenvalues alone, spatialreg's rule. (First pass:
  "from the real parts of every eigenvalue, spatialreg's reading", which
  was not spatialreg's rule and could be narrower; CORRECTED in the fix
  round, MINOR-2.) brms refuses such a W
  (brms-37), so the reference is the dense normal with `determinant()`.
- 4110: a jacobian slot that gives `O K` (O orthogonal) is accepted when
  the Gram matrices of the two check directions agree; such a block is
  marked `jacobian_rotated`, and the pointwise log-density then reads
  `Q e` as `K'(K e)` with `e = y - unwhiten(0)`, since `K'u` is wrong there.
- 4115: `options(frmtmb.autocor_dense_max = )` (default 1000) replaces the
  constant; frmtmb.spatial's test sets it to 60 instead of reading
  `frmtmb:::autocor_loo_dense_max`, which also makes that test smaller.
- 4120: `fcor()`'s pointwise log-density reads `diag(M^-1)` from the
  block (computed once at build) and `Linv'u`, so no `n x n` matrix is
  formed per draw. Measurement below.
- 4320: `predictions(newdata =)` on an `ar(cov = FALSE)` fit failed because
  marginaleffects' `pad()` prepends one row per factor level, copied from
  the first row, so with the factor as the term's group the padded rows
  repeat a time of the rows asked for (`dev/arspat-4320.R`,
  `dev/arspat-out/4320.txt`, marginaleffects 1.0.0). `get_predict()` now
  predicts the padding rows apart from the rows asked for under a
  `cov = FALSE` term.

## L1 and L4 against brms 2.23.0

Scripts `dev/arspat-l1b.R` (logs `dev/arspat-out/l1b-all.txt`) and
`dev/arspat-l1c.R` (sar cases, `dev/arspat-out/l1c-sar.txt`), table by
`dev/arspat-l1-summary.R`. Data: `ac_data()` in the script, seed 29, 12
groups of 6 (30 groups and multivariate t residuals, df 4, on the student
cases; 24 groups of 4 with 4 rows dropped on `se_unstr`). Joint cases read
brms's `log_prob` (no Jacobian) at frmtmb's estimate and at two points
stepped 0.05 off it in every coordinate (`eq_l1()` of
`dev/equiv-lib.R`); "marginal over Ymi" cases read frmtmb's logLik
against brms's program integrated over `Ymi`: its Laplace value on
gaussian (exact, the program is quadratic in `Ymi`), and one-dimensional
`integrate()` per missing value on student. brms's SAR densities leave out
`0.5 N log(2 pi)` (66.16357 at N = 72), which is subtracted. L4 is
`rstan::optimizing()` of brms's program under flat priors from 6 starts
against frmtmb's maximized logLik (`eq_l4_opt()`).

<!-- generated by dev/arspat-l1-summary.R; do not edit -->
| case | identity | max abs const over 3 points | brms gradient at estimate | L4 brms max - frmtmb max |
|---|---|---|---|---|
| se_ar | joint | 1.4e-14 (3 points) | 3.0e-06 | 7.1e-14 |
| se_ar_het | joint | 1.4e-14 (3 points) | 2.6e-04 | 6.8e-10 |
| se_ma | joint | 2.8e-14 (3 points) | 9.7e-06 | 9.5e-13 |
| se_arma | joint | 2.8e-14 (3 points) | 1.7e-04 | 3.9e-10 |
| se_cosy | joint | 1.4e-14 (3 points) | 1.3e-05 | 1.3e-12 |
| se_unstr | joint | 1.4e-14 (3 points) | 2.3e-04 | 6.6e-10 |
| se_ar_re | joint | 4.3e-14 (3 points) | 6.7e-16 | - |
| se_ar_student | joint | 1.1e-13 (3 points) | 5.1e-04 | 4.9e-10 |
| se_cosy_student | joint | 1.1e-13 (3 points) | 5.4e-05 | 1.6e-11 |
| mi_ar | marginal over Ymi | 4.3e-14 (3 points) | - | - |
| mi_cosy | marginal over Ymi | 4.3e-14 (3 points) | - | - |
| mi_unstr | marginal over Ymi | 2.8e-14 (3 points) | - | - |
| mi_ar_student | marginal over Ymi | 1.7e-13 (3 points) | - | - |
| mi_fcor | marginal over Ymi | 1.4e-14 (3 points) | - | - |
| mi_sar_error | marginal over Ymi | 6.4e-11 (3 points) | - | - |
| mi_sar_lag | marginal over Ymi | 6.4e-11 (3 points) | - | - |
| mi_fcor_student | marginal over Ymi | 2.1e-10 (3 points) | - | - |
| mi_sar_error_student | marginal over Ymi | 1.2e-10 (3 points) | - | - |
| mi_sar_lag_student | marginal over Ymi | 2.2e-10 (3 points) | - | - |
| skew_se_ar | joint | 2.8e-14 (3 points) | 1.7e-14 | - |
| skew_se_ar_cond | ERROR: Please set cov = TRUE in ARMA structures when including known standard errors.  | - | - | - |
| skew_se_cosy | joint | 0.0e+00 (3 points) | 9.8e-13 | - |

Joint cases: 11; max |const| 1.1e-13; max |L4| 6.8e-10 over 8.
Marginal cases: 10; max |const| 2.2e-10.
<!-- end generated -->

The sar residual is the precision of the two printed values (`%.10f`).
`se_ar_re` has `(1 | h)`, so its gradient is read on brms's inner block
and it has no L4 (the posterior mode is not the marginal maximum).

A dissolved signal, kept: the first run (`dev/arspat-l1.R`, log
`dev/arspat-out/l1-all.txt`) used gaussian data on the student cases. nu
ran to about 1e9 there and the constant was 1.5e-5, 3.7e-6 and 1.3e-5 at
the three points of `se_ar_student`, and 1.8e-5 on `mi_ar_student`: not an
identity. With multivariate t data (df 4) it is 1.1e-13 and 1.7e-13. At a
nu of 1e9 the difference of two `lgamma()` values loses about 1e-5 nats
in one of the two programs; frmtmb's head is `lgamma_shift_diff()`, so the
loss is presumably Stan's `multi_student_t`. Not pursued: brms fits at
such a nu are not where the two would be compared.

## Tests seen failing on rellib-r20

`dev/arspat-test1.R` with `ARSPAT_LIB=""` (logs
`dev/arspat-out/t-se-mi-base.txt`, `t-sar-base.txt`):

- `tests/testthat/test-autocor-se-mi.R`: pass=4 fail=1 error=7 on the base;
  the seven errors are the refusals ("se() cannot be combined with a
  residual correlation term", the same for mi(), "vcov_cluster() does not
  support the residual correlation term"), the failure is the gap message
  that the base never gives. Lane build: pass=30 fail=0.
- `extensions/frmtmb.spatial/tests/testthat/test-sar.R`: two errors on the
  base, "the eigenvalues of 'M' are complex" and "sar(): mi() cannot be
  combined". Lane build: pass=121 fail=0 skip=4 (brms rows, gated).

## frm_compat() rows changed

`dev/arspat-compat-probe.R`, before (`ARSPAT_LIB=""`,
`dev/arspat-out/compat-before.txt`) and after
(`dev/arspat-out/compat-after.txt`), with frmtmb.spatial loaded:

| pair | before | after | where |
|---|---|---|---|
| `ar()` x `se()` | refused | conditional | R/compat.R, the `for (nm in c("ar()", "ma()", "arma()"))` rows |
| `ma()` x `se()` | refused | conditional | same row |
| `arma()` x `se()` | refused | conditional | same row |
| `fcor()` x `mi()` | refused (`fcor()` x kind:aterm) | conditional | R/compat.R, new row after `fcor()` x kind:aterm |
| `sar()` x `mi()` | refused (seam kind:aterm) | conditional | R/autocor-seam.R `autocor_seam_compat_rules()`, new row for every registered term |

Notes changed, status kept: `group:autocor` x kind:aterm (conditional;
it now says se() and mi() work on gaussian and student, so `cosy()` and
`unstr()` x `se()`/`mi()`, which resolve to it, read right), and `ar()`,
`ma()`, `arma()` x kind:aterm (conditional). No survey pair in
`dev/equiv-data/sample.tsv` is one of the five; the two that share a
changed note, p112 (`trunc()` x `ma()`) and p154 (`cens()` x `cosy()`),
keep their status and their code path (still refused under the
covariance form), so they were not remeasured.

## fcor() pointwise log-density cost (4120)

`dev/arspat-fcor-cost.R`, log `dev/arspat-out/fcor-cost.txt`, seed 1,
n = 1500, M = 0.5^|i - j|, one process, lane build; the old route is the
same build with `diag_inv` taken off the block:

- the two routes differ by 2.74e-16 relative;
- peak Vcells memory over one call: 0.2 Mb new, 34.5 Mb old (two
  1500 x 1500 doubles are 36 Mb);
- seconds per draw, minimum of 5 interleaved rounds of 128 calls: new
  0.0099, old 0.0398, control (new against itself) 0.0091; old / new
  4.01, control / new 0.91.

## OpenBLAS

`dev/ciharden-openblas.sh 0.3.26 lapack`, then `dev/arspat-openblas.sh
0.3.26` (OPENBLAS_NUM_THREADS = 4, the runners' count), logs
`dev/arspat-out/ob0.3.26-*.txt`:

| file | pass | fail | error | skip |
|---|---|---|---|---|
| test-autocor-se-mi.R | 34 | 0 | 0 | 0 |
| test-autocor-loglik.R | 78 | 0 | 0 | 0 |
| test-autocor.R | 195 | 0 | 0 | 0 |
| test-sandwich.R | 58 | 0 | 0 | 0 |
| test-autocor-seam.R | 80 | 0 | 0 | 0 |
| frmtmb.spatial test-sar.R | 121 | 0 | 0 | 4 (brms rows, gated) |

0.3.32 was not run: no change of this lane sits near a rounding-decided
branch (the Gram check's tolerance is 1e-6 relative, the identities are
1e-10 and tighter).

## What was not done, and why

- `sar()` beside `ar(cov = FALSE)` (half of 3903): a response holds one
  residual term in `frame$autocor`, and 33 readers in 11 files branch on
  `autocor_is_cond()` of that slot. Much larger than the item's size;
  filed in the backlog section of this round.
- `newdata2` for `sar()` (3903): not in the brief's list; left.
- 3901 (rho's interval past `1 / max eigenvalue`): an interval on the
  logit chart of the bounds needs `lo`/`hi` carried through five rbind
  and back-transform sites in four files shared with other lanes; filed.
- Counting rows instead of time levels under the covariance form (4104):
  a user decision; see the question in section 4.
- 0.3.32 OpenBLAS: see above.

## Suites, one file per process

`bash dev/release/run-par.sh arspat-all dev/arspat-out/jobs-all.txt 3`
with `FRMTMB_RELLIB` = the lane library and `NOT_CRAN=true`, log
`dev/release/arspat-all.log` ("ARSPAT-ALL ran 511 of 511"). Every file's
`lib:` line names the lane's core; frmtmb.spatial is the lane's build, the
other extensions rellib-r20's (none calls a function this lane changed).

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 318 | 23533 | 0 | 0 | 240 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 687 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 87 | 4050 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 493 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 556 | 0 | 0 | 0 | 0 |

The skips are the gated brms tiers. Gated, on the lane build
(`dev/release/arspat-gated1.log`, `FRMTMB_BRMS_FIT_TESTS=true`):
test-brms-likelihood.R pass=544 skip=0, frmtmb.spatial test-sar.R
pass=141 skip=0, test-autocor-se-mi.R pass=34. On rellib-r20 the new row
18f errors (the refusal), `dev/arspat-out/t-brmslik-base.txt`: pass=540
error=1.

## One session per suite

`bash dev/release/run-sessions.sh <tier> <N> <package>` with
`FRMTMB_RELLIB` = the lane library:

- frmtmb.spatial (`dev/release/arspat-sess-spatial.log`): forward,
  reverse and preload each files=8/8 pass=493 fail=0 err=0 skip=6 warn=0.
  One LEAK line each, `test-spde.R envvars: GDAL_DATA [with load of ...
  sf ...]`, tagged with a package load and from a file this lane did not
  touch. No option leak from test-sar.R's `frmtmb.autocor_dense_max`.
- frmtmb (`dev/release/arspat-sess-core.log`): forward files=305/318
  pass=23533, reverse pass=23698, preload pass=23801 skip=242, each fail=0
  err=0 warn=0. Every LEAK line is tagged `[with load of frmtmb.sample]`,
  `frmtmb.survival` or `frmtmb.eam` (test-arg-refusal.R,
  test-cox-home.R, test-families-extreme.R); none is untagged. The 13
  files with no row are the ones gated as a whole (test-brms-suite-*.R,
  test-latres-newrow-brms.R, test-linpred-layers-brms.R,
  test-me-newdata-brms.R), which read pass=0 skip=0 in the per-file run
  too.
- The other seven extensions ran one file per process only (above): no
  extension calls a function this lane changed, and every test of theirs
  passed against the lane's core.

## R CMD check --as-cran

`bash dev/arspat-check.sh <pkg>`, built with vignettes in
`dev/arspat-check/`, `R_LIBS` = lane library, rellib-r20, user library,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`:

- frmtmb: Status: 2 NOTEs. The V8 math-rendering NOTE, and the examples
  timing NOTE (`frm_hazard_reads` 3.55 s CPU, 8.06 s elapsed; `ps` 3.37
  and 5.36; `dharma_residuals` 1.28 and 6.88): elapsed far above CPU on
  examples this lane did not touch, the load signature lane-rules.md
  describes. Tests 69 min, OK.
- frmtmb.spatial: Status: 1 NOTE (V8).

## Bump and floor

Core: a minor bump (new combinations fit, a new option, a new message
class; the `se()` refusal is gone, which NEWS states). frmtmb.spatial: a
minor bump; its `mi()` test needs the core of this round, so its floor on
frmtmb is that new core version (the consolidator sets the number; the
DESCRIPTION is not edited here).

## Fix round (review dev/reviews/2026-11-11-arspat.md: FIX, 1 MAJOR, 8 MINOR)

Each finding was reproduced from the reviewer's scripts and logs in
`dev/reviews/arspat-review-*` before the fix.

### MAJOR-1: student() registered terms with mi()

Reproduced: `dev/reviews/arspat-review-laplace.log` (objective minus the
dense `dmvt()` margin equals `c(nu, n, m) = lgamma((nu + n)/2) -
lgamma((nu + n - m)/2) + (m/2) log(2/(nu + n))` to every printed digit, fcor
and sar, m = 3, 10, 20, nu 3 to 30) and `-laplace3.log` (prior gamma(2,
0.1): nu 28.11, logLik - exact -1.3614, closed form -1.3614, `tk` -0.927).

Fix: `autocor_seam_t_laplace()` (R/autocor-seam.R) is `c`, written with
`lgamma_shift_diff()` so it keeps its digits at a large nu, and
`build_objective()` adds it to the negative log-likelihood of a registered
term's student() response. `autocor_seam_t_missing()` decides where: the
maximum-likelihood route only (`mi_ml_frame()` has run; a sampler keeps the
joint density and gets no correction), no measurement-error `mi(sd)`, no
other density reading the response's missing values (`mi(y)` as a
predictor), every latent value in the model the response's own, and no
random-effect block or `me()` group. With another random effect the joint
Laplace error is not this constant (the Schur complement of the joint
Hessian carries the derivative of the inner log-determinant), so nothing is
corrected there and `mi_laplace_note()` still measures it.

Composition with `mi_laplace_exact()`/`mi_laplace_note()` (R/frame.R):
`mi_laplace_exact()` now also returns TRUE for a response that
`autocor_seam_t_missing()` corrects, so the note does not run on a case
that is now exact. Asserted in the test: `mi_laplace_exact(f$frame)` TRUE
and `f$mi_laplace` NULL on the reviewer's prior case; on the same model
with `(1 | g)` the correction is off and `mi_laplace_exact()` FALSE.

Tests (seen failing: on rellib-r20 each errors at the `mi()` refusal; on
the pre-fix lane build the reviewer's numbers above are the failure):
- `test-autocor-se-mi.R` "mi() beside fcor() on student() integrates out
  exactly": nu held at 3 and 10 by the optimizer, m = 3, 10, 20, logLik
  against `mvtnorm::dmvt()` on the observed rows, relative 1e-9; the prior
  case against `dmvt()` at the fitted nu; the absent case with `(1 | g)`.
- frmtmb.spatial `test-sar.R` "mi() beside sar() on student() integrates
  out exactly": lag and error, m = 3, 10, 20, nu 3 and 10, against a dense
  multivariate t written in the test (mvtnorm is not in spatial's
  Suggests), relative 1e-9.

L1 rows added (`dev/arspat-l1d.R`, logs `dev/arspat-out/l1d-student.txt`,
`l1d-fix.txt`): brms's program integrated over its `Ymi` by Laplace (Newton
to the mode, finite-difference Hessian, h = 1e-4) less `c`, against
frmtmb's logLik at nu held at 5, 10 and 30, one group of 72 rows, df 4
data; brms's student SAR densities leave out `0.5 N log(pi)` = 41.21027,
subtracted. The table at the top of this file is regenerated with them:
`mi_fcor_student` 2.1e-10, `mi_sar_error_student` 1.2e-10,
`mi_sar_lag_student` 2.2e-10 (the precision of the finite-difference
Hessian).

Not changed (reviewer's related note): an `ar()`-family value that another
density reads stays latent on student() under Laplace, and `summary()` does
not say so; `mi_laplace_note()` measures that case as before.

### MINOR-2: spatialreg's interval

`sar_build()` now takes rho's interval from the real eigenvalues alone,
`1 / range(Re(eig[real]))`, with the lane's `1e-8 * max Mod` test for real.
On a W built so that a complex real part (-0.5, a directed 3-cycle) is more
extreme than every real eigenvalue (+/- 0.2, an edge scaled by 0.2), the
interval is (-5, 1) where the real parts gave (-2, 1). `dev/arspat-sar-ends.R`
(log `dev/arspat-out/sar-ends.txt`): `sar_logdet()` against `determinant()`
      rho -4.999995: sign +1, log|det| -33.144340, |diff| 9.32e-10
      rho -2.500000: sign +1, log|det| 10.092902, |diff| 7.11e-15
      rho +0.000000: sign +1, log|det| 0.000000, |diff| 0.00e+00
      rho +0.500000: sign +1, log|det| -0.574327, |diff| 3.11e-15
      rho +0.999999: sign +1, log|det| -51.030885, |diff| 8.01e-10
    3-NN W: interval (-1.682443, 1.000000); 1/range(real eig) (-1.682443, 1.000000); 1/range(Re all) (-1.682443, 1.000000)
      rho -1.682441: sign +1, log|det| -27.688571, |diff| 1.06e-10
      rho -0.841221: sign +1, log|det| -2.949121, |diff| 2.66e-15
      rho +0.000000: sign +1, log|det| 0.000000, |diff| 0.00e+00
      rho +0.500000: sign +1, log|det| -1.552474, |diff| 1.33e-15
      rho +0.999999: sign +1, log|det| -46.886304, |diff| 9.63e-10

The sign is +1 at every point; the difference is at most 9.3e-10 at the
ends, where the log-determinant is -33 to -51. The test in `test-sar.R`
asserts the (-5, 1) interval and these values, which the real-parts rule
would fail. Docs (`?sar`, "Arguments of the term", "The range of rho"),
the spatial NEWS, the code comment and the backlog CLOSED line now say the
real eigenvalues. The refusal message now says "a negative and a positive
real eigenvalue" (test updated).

### MINOR-3: skew_normal() with se()

`dev/arspat-l1d.R`, log `dev/arspat-out/l1d-fix.txt` (12 groups of 6, seed
29): `y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)` on
`skew_normal()`, the latent residual, joint at the conditional modes:
const 2.8e-14 at three points, brms gradient 1.7e-14; the same with
`cosy(t, g)`: 0 at three points, gradient 9.8e-13. So under `cov = TRUE`
it agrees, and the compat notes (`group:autocor` x kind:aterm, `ar()`,
`ma()`, `arma()` x kind:aterm and x `se()`) now say so. Under
`cov = FALSE` brms refuses it ("Please set cov = TRUE in ARMA structures
when including known standard errors") and frmtmb fits it (logLik
-89.9534, the same on rellib-r20 per the reviewer). That is a parity
question, not an L1 disagreement: reported to the coordinator, filed in
the backlog, behavior not changed; the notes state it as an open parity
question. The vignette sentence now names the five terms that take
`se()` and says `fcor()` and registered terms such as `sar()` do not.

### MINOR-4: NEWS

Bullets added for `options(frmtmb.autocor_dense_max = )` (4115), the
`O K` jacobian slot (4110), `fcor()`'s pointwise cost (4120) and
`predictions(newdata =)` under `cov = FALSE` (4320). The `mi()` bullet
now states the student correction. No BREAKING marks: each change lifts a
refusal, adds a message or option, or fixes a failure.

### MINOR-5: the nu -> 1e9 gap

Closed in the backlog as not a frmtmb defect, with the reviewer's numbers
(`dev/reviews/arspat-review-se.log`): the naive lgamma difference errs by
-5.3e-7 (k = 3) and -1.1e-6 (k = 6) at nu = 1e9; `lgamma_shift_diff()`
by at most 2.1e-14. Nothing goes upstream.

### MINOR-6: dev/ paths

Removed from the `group:autocor` x kind:aterm note, from the comment in
`check_autocor_response()` and from the comment in `sar_build()`. The
pre-existing `(dev/arcov-validate.R)` in the `ar()`/`ma()`/`arma()` x
kind:aterm note is not this lane's and stays.

### MINOR-7

No action asked (merge notes). Spatial's floor is the new core, as stated.

### MINOR-8: gap message noise

`dev/arspat-gapnoise.R` names the test that raises each message (log
`dev/arspat-out/gapnoise-before.txt`): test-autocor.R 3 (ragged groups
against gls), test-autocor-cond.R 1 (`cond_data()` is ragged),
test-autocor-loglik.R 15 (`ll_ac_data()` skips a level on purpose),
test-autocor-se-mi.R 2 (the reference fit on `d[-miss, ]`). Each fixture
skips a level on purpose, so each place now asserts the message with
`expect_message(class = "frmtmb_autocor_gap")`; in test-autocor-loglik.R
through `ll_gap()`, which asserts it under `ar()`, `ma()` and `arma()` with
`cov = TRUE` and lets the other terms through. No `suppressMessages()`.

After: `dev/arspat-out/gapnoise-after.txt`, 0 uncaught in each of the four
files (the outer handler sees a message only when no `expect_message()`
caught it).

### Runs after the fix

- Changed test files on rellib-r20 (seen failing): test-autocor-se-mi.R
  pass=4 fail=1 error=9 (the refusals, the gap message, the new student
  fcor test); test-autocor-loglik.R pass=66 fail=13 error=2 (no gap
  message on the base, no `se()`); frmtmb.spatial test-sar.R pass=101
  fail=1 error=4 (complex W refused, `mi()` refused twice, the new
  refusal text, the option the base core lacks). Logs
  `dev/arspat-out/fixbase-*.txt`.
- Changed test files on the lane (`dev/release/arspat-fix.log`,
  `arspat-fix2.log`): test-autocor-se-mi.R 53/0/0/0, test-autocor-loglik.R
  93/0/0/0, test-autocor.R 198/0/0/0, test-autocor-cond.R 58/0/0/0,
  test-compat.R 868/0/0/0, test-mi-families.R 194/0/0/0, spatial
  test-sar.R 157/0/0/4 (pass/fail/error/skip), warn 0 each. The first run
  had 6 and 12 escaped "Large maximum absolute gradient" warnings from the
  nu-held fits and one error from the old refusal text; both fixed
  (`allow_warnings()` on the held fits, the text in the test).
- Full suites, one file per process (`dev/release/arspat-fixall.log`,
  "ran 326 of 326", every `lib:` line the lane's): frmtmb 318 files
  pass=23571 fail=0 err=0 skip=240 warn=0; frmtmb.spatial 8 files pass=529
  fail=0 err=0 skip=6 warn=0. The other seven extensions were not rerun:
  the fix touches the student seam density, `mi_laplace_exact()`, notes,
  docs and tests, and no extension test uses a registered term with
  `mi()` on student (their earlier run against this lane's core is above).
- frmtmb.spatial in one session (`dev/release/arspat-fix-sess-spatial.log`):
  forward, reverse and preload each 8/8 pass=529 fail=0 err=0 skip=6
  warn=0, one tagged GDAL_DATA LEAK from test-spde.R each, as before.
- OpenBLAS 0.3.26 with LAPACK (`dev/arspat-out/ob-0.3.26-fix-summary.txt`,
  4 threads): test-autocor-loglik.R 93, test-autocor-se-mi.R 53,
  test-autocor-seam.R 80, test-autocor.R 198, test-sandwich.R 58, spatial
  test-sar.R 157 + 4 skip, all fail 0 error 0. `dev/ciharden-out` deleted.
- R CMD check --as-cran (`dev/arspat-check.sh`): frmtmb Status: 1 NOTE
  (V8 math rendering; examples 56 s, tests 46 min, OK; the examples timing
  NOTE of the first pass did not recur, consistent with load); frmtmb.spatial
  Status: 1 NOTE (V8).
