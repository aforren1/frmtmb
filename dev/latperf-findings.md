# Lane latperf: latent residuals at scale

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-latperf`, library `wt-latperf-lib`, base library `rellib-r18`.
Status: DONE, fix round DONE (2026-11-06). This file is the record if the
report does not arrive. The review is `dev/reviews/2026-11-06-latperf.md`
(FIX: M1, m1 to m7); the section "Fix round" answers it and corrects the
claims below that it found wrong.

## Shared functions this lane changes (named early, for the merge)

- `R/predict.R`: `get_joint_cov()` (new argument `dense`, sparse branch for a
  latent residual), `joint_cov_repair()` (its non-`V` part moved, unchanged,
  into the new `joint_cov_repair_meta()`), `frm_linpred()` (one line: the
  design stays sparse when the joint covariance is), `linpred_layers()` (a
  skip above the row limit), `predict_mean_se()` (sparse gradient when the
  joint covariance is sparse), `fitted_point_se()` (the skip applies to the
  dense route only), `frm_joint_cov()` (forms `V` from the sparse form),
  `frm_lp_basis()` (forms `V` from the sparse form, through `lb_dense()`),
  `lp_basis_out()` and `lb_rowvar()` (sparse form), roxygen of
  `frm_lp_basis()`.
- `R/se-reparam.R`: `jc_rowvar()`, `jc_quad()` (one dispatch line each).
- `R/se-check.R`: `jc_nonest()` (a sparse gradient stays sparse).
- `R/confint.R` `smooth_edf()` and `R/methods-fit.R` `ranef()`: `diag(jc$V)`
  became `jc_diag(jc, bj)`, one line each.
- `R/autocor-latent.R`: `latres_skip_se()`.
- `R/autocor.R`: roxygen of `?frmtmb-autocor`.
- New: `R/joint-sparse.R`.

No `<<-` site was touched.

## Profile before (base build rellib-r18)

Script `dev/latperf-prof.R` (one process per design; seed 1; the generator of
`dev/review-latres-fittedprof.R`): `y ~ x + ar(t, g, cov = TRUE)`, one
poisson series ("one"), and the same plus `(1 | id)` with 50 levels ("two").
`dense` moves both row limits to Inf so the base build takes its dense route.
Logs `dev/latperf-out/base-*.txt`. Times in seconds (elapsed), heap is R's
"max used" since the step's `gc(reset = TRUE)` in MB, ws the process's peak
working set in MB, saved the serialized fit after the readers.

| n | design | route | fit | fitted() | residuals() | linpred se | heap | ws | saved MB |
|---|---|---|---|---|---|---|---|---|---|
| 2000 | one | dense (default) | 1.97 | 4.69 | 4.55 | 4.39 | 378 | 508 | 36.21 |
| 2000 | two | dense (default) | 3.30 | 5.94 | 5.20 | 5.00 | 385 | 554 | 38.57 |
| 4000 | one | dense (limit moved) | 3.80 | 40.49 | 41.14 | 39.40 | 838 | 953 | 129.92 |
| 8000 | one | dense (limit moved) | 8.39 | 332.22 | 376.61 | 512.91 | 2676 | 2778 | 500.43 |
| 8000 | one | default: skipped, NA | 1.50 | 0.00 | 0.00 | 0.00 | 176 | 415 | 9.56 |
| 32000 | two | default: skipped, NA | 10.94 | 0.00 | 0.03 | 0.00 | 202 | 617 | 35.59 |

At 32000 rows, design two: `fitted(newdata = 50 rows)` gave 50 of 50
`Est.Error` NA with the message "This standard error reads the whole joint
covariance of a latent residual as a dense matrix, about 16 GB at 32000
rows", and `frm_linpred(newdata, se.fit = TRUE)` stopped with that error
(item 2). At 8000 and 32000 rows the default skips every fitted-row standard
error (item 1). The 8000-row dense fit itself took 8.39 s against 1.50 s
because the Laplace check also reads `frmtmb.latres_se_rows`.

Where the time goes, base, 2000 rows, design two, `fitted()` alone
(`dev/latperf-prof-two.R base 2000 two`, Rprof at 0.01 s,
`dev/latperf-out/base-rprof-2000-two.txt`): sampled 5.95 s, of which `%*%`
5.61 s self (94.3%), under `rowSums` in `lb_rowvar()`: the dense design of the
fitted rows times the dense joint covariance. The same at design one
(`dev/latperf-out/base-2000-one.txt`, interval 0.02): `%*%` 2.94 of 3.04 s
under `lb_rowvar()`. The cache after `fitted()` at 2000 rows: `Vjoint` 32.27
MB of the 38.6 MB fit.

## The change

`R/joint-sparse.R`. On a fit with a latent residual, `get_joint_cov()` holds
the joint covariance sparse: the supernodal Cholesky factor of the inner block
`Q_rr` of the joint precision, `J = -Q_rr^-1 Q_rf` (n x |f|), and the outer
covariance `Vf` (the one the SE check kept for propagation, else the inverse
of the Schur complement, which is what `solve(Q)` holds). That is the block
inverse `joint_cov_repair()` already writes. A row then has
`a' V a = a_r' W a_r + g' Vf g` with `g = J' a_r + a_f`, and `a_r' W a_r`
reads `W = Q_rr^-1` only at pairs of positions one row loads together. Where
the row's likelihood couples a pair in `Q_rr`, it is in the pattern of the
factor, where TMB's selected inverse (`solveSubset`, the routine `sdreport()`
uses for `diag.cov.random`) is exact. A family that separates its dpars does
NOT couple every pair (review M1, "Fix round" below). The readers:

- fitted rows (`lb_rowvar()`, `jc_rowvar()`): per row, the pairs of its own
  entries looked up in the selected inverse; the pairs outside the pattern
  read columns of `W` solved once per distinct position (`jcs_W_pairs()`,
  since the fix round; before it, such a row was solved directly, which is
  quadratic); without TMB's routine, every row is solved directly in 32 MB
  chunks.
- a block `V[i, j]` (new rows, `frm_curve()`, `fd_joint_cov()`): the
  `[` method of class `frm_joint_sparse`, one solve per inner position on the
  shorter side, refused above `frmtmb.latres_dense_rows^2` numbers.
- `lp_basis_out()` drops the columns no row loads (a new row loads no fitted
  latent), and above 200 inner positions keeps the sparse form for
  `lb_rowvar()`; `frm_lp_basis()` and `frm_joint_cov()` form the dense block
  they return.
- the sparse route does not apply where the covariance is held as a factor
  (`cov_fixed_root`, `joint_root`) or where the inner factor or the Schur
  complement fails; `get_joint_cov()` then takes the dense route as before,
  with both limits.

`latres_skip_se()` now skips only where a dense block over every latent is
formed: the differenced category probabilities of an ordinal family's fitted
rows (`fit_fd_se()` with every latent), the stacked `cs()` layers, and a fit
the sparse route does not reach.

## Profile after (lane build)

Same script and seeds, `dev/latperf-out/lane-*.txt`:

| n | design | fit | fitted() | residuals() | linpred se | heap | ws | saved before | saved after |
|---|---|---|---|---|---|---|---|---|---|
| 2000 | one | 2.67 | 0.17 | 0.00 | 0.02 | 253 | 536 | 5.12 | 6.02 |
| 8000 | one | 2.24 | 0.42 | 0.03 | 0.02 | 269 | 423 | 9.56 | 13.43 |
| 32000 | one | 6.13 | 1.67 | 0.11 | 0.11 | 304 | 605 | 27.75 | 43.16 |
| 2000 | two | 4.38 | 0.14 | 0.02 | 0.01 | 262 | 524 | 5.66 | 7.51 |
| 8000 | two | 3.38 | 0.50 | 0.03 | 0.01 | 267 | 440 | 11.55 | 19.17 |
| 32000 | two | 15.31 | 2.80 | 0.14 | 0.15 | 354 | 656 | 35.59 | 65.98 |

New rows, design two: `fitted(newdata = 50 rows)` 0.17 s at 2000, 0.22 s at
8000, 0.44 s at 32000, every `Est.Error` finite (base at 32000: 50 of 50 NA).

Values against the base build where both run: at 2000 rows the first and
median `Est.Error` of `fitted()` print identically to 10 digits (one:
0.5218028311, 0.7217602799; two: 0.4932818549, 0.7561478484), and the 50
new-row values of design two at 2000 rows too (0.6880915469, 1.397609028). At
8000 rows, design one, base dense 0.528786047 and 0.7195509206, lane the same.
Within one process, ratio against the dense route on the same fit
(`dev/latperf-agree.R`, seed 11, 30 groups of 12 rows, 360 rows, log
`dev/latperf-out/agree.txt`), max |sparse/dense - 1| over every finite row:

<!-- generated from dev/latperf-out/agree.txt -->
```
ar1 poisson            sparse TRUE | fitted 4.44e-16 (360, mis 0) | linpred 4.44e-16 (360, mis 0) | resid 4.44e-16 (360, mis 0) | newdata 0.00e+00 (4, mis 0)
ar1 + (1|s) poisson    sparse TRUE | fitted 4.44e-16 (360, mis 0) | linpred 4.44e-16 (360, mis 0) | resid 4.44e-16 (360, mis 0) | newdata 2.22e-16 (4, mis 0)
ma1 poisson            sparse TRUE | fitted 1.78e-15 (360, mis 0) | linpred 1.78e-15 (360, mis 0) | resid 1.78e-15 (360, mis 0) | newdata 0.00e+00 (4, mis 0)
arma1 poisson          sparse TRUE | fitted 6.66e-16 (360, mis 0) | linpred 5.55e-16 (360, mis 0) | resid 6.66e-16 (360, mis 0) | newdata 0.00e+00 (4, mis 0)
cosy poisson           sparse TRUE | fitted 4.44e-16 (360, mis 0) | linpred 4.44e-16 (360, mis 0) | resid 4.44e-16 (360, mis 0) | newdata 0.00e+00 (4, mis 0)
unstr poisson          sparse TRUE | fitted 9.99e-15 (210, mis 0) | linpred 9.99e-15 (210, mis 0) | resid 9.99e-15 (210, mis 0) | newdata 0.00e+00 (4, mis 0)
condar poisson         sparse TRUE | fitted 1.11e-15 (360, mis 0) | linpred 1.11e-15 (360, mis 0) | resid 1.11e-15 (360, mis 0) | newdata NA (0, mis 0)
ar1 binomial trials    sparse TRUE | fitted 5.55e-16 (360, mis 0) | linpred 4.44e-16 (360, mis 0) | resid 5.55e-16 (360, mis 0) | newdata 0.00e+00 (4, mis 0)
ar1 zip                sparse TRUE | fitted 4.44e-16 (360, mis 0) | linpred 6.66e-16 (360, mis 0) | resid 4.44e-16 (360, mis 0) | newdata 0.00e+00 (4, mis 0)
ar1 cumulative         sparse TRUE | fitted 1.33e-15 (1440, mis 0) | linpred 1.22e-15 (360, mis 0) | resid NA (0, mis 0) | newdata 0.00e+00 (16, mis 0)
```

This is an IDENTITY measured as a numerical check: both routes compute the
same block inverse, one through a sparse solve of the whole `Q`, the other
through the factor of `Q_rr` and the Schur complement, so the residual is
rounding. `(n, mis k)` is the count of finite rows compared and of rows finite
on one route only. `unstr` reads 210 of 360: on this seed it loses a standard
error, and 150 rows are NaN with the lost-direction warning on both routes,
which exercises `joint_cov_repair_meta()` and `jc_nonest()` with a sparse
gradient. The `cumulative` fitted rows are the differenced category
probabilities (4 x 360, the dense block below the row limit). `condar`
(cov = FALSE) refuses new rows, as brms does. The newdata rows of `binomial`
read a dense block of 4 positions.

REML, where `beta` is an inner parameter (`dev/latperf-reml.R`, seed 11, the
same 360 rows, `pos ~ x + (1 | s) + ar(t, g, cov = TRUE)`, log
`dev/latperf-out/reml.txt`), `frm_linpred(se.fit = TRUE)` against the dense
route: lognormal 1.31e-13 (360 rows, 0 mismatched), Gamma(log) 2.05e-14.

How often a fitted row's pairs leave the pattern of the factor, so that it is
solved directly (`dev/latperf-pattern.R`, log `dev/latperf-out/pattern.txt`):
0 rows on every design above and on designs one and two at 8000 rows (seed 1).
Each call read the selected inverse once. WITHDRAWN as a general claim: none of
those families separates its dpars. The review found every row of a hurdle with
`hu ~ (1 | id)`, and the rows with y > 0 of a zero-inflated fit with
`zi ~ (1 | id)`, outside the pattern (M1); see "Fix round".

Where the time goes after, design two, `fitted()` alone
(`dev/latperf-prof-two.R lane <n> two`, `dev/latperf-out/lane-rprof-*.txt`):
2000 rows sampled 0.11 s; 8000 rows 0.35 s; 32000 rows 2.00 s, of which
`sdreport()` for the joint precision 1.27 s (63.5%), the selected inverse
and its lookup 0.34 s, the Cholesky factor 0.25 s.

What the fit kept (same logs), 32000 rows, design two: `Vjoint` 25.68 MB
(factor 9.73, `J` 1.22, the selected inverse's lookup 14.05), `Qjoint` 3.55,
`sdr` 1.13. At 2000 and 8000 rows `Vjoint` is 1.63 and 6.43 MB: linear in the
rows. Base at 2000 rows: `Vjoint` 32.27 MB. The fix round drops the lookup
from the fit (review m3): `Vjoint` 11.63 MB at 32000 rows, see "Fix round".

## Tests

New file `tests/testthat/test-latres-sparse.R` (7 tests). Seen to fail on the
base build (rellib-r18, `dev/latperf-out/t-sparse-base.txt`):
`pass=8 fail=7 err=5`. The behavioral failures, before any missing symbol:

- fitted rows above the limit: `said$n` 3 where 0 (the skip message), and
  `Est.Error` and `se.fit` not finite;
- new rows beside `(1 | s)` above the dense limit: the error
  `frmtmb_latres_dense` "about 0.0021 GB at 360 rows";
- the cache after `fitted()`: 1111494 bytes against the bound 267912
  (`0.25 * 8 nb^2`): the dense `V`;
- `frm_lp_basis(newdata)`: its `coef_pos` holds fitted latents;
- the ordinal route: the predictor's `se.fit` skipped too (message 1, not 0).

The missing-symbol errors (`jc_is_sparse`, `jcs_joint_cov`,
`jcs_rowvar_inner`) are the weak form, on the assertions against the dense
route and the unit test of the fallbacks.

Lane: `pass=53 fail=0 err=0 skip=0 warn=0`. The last test constructs the
guard's absent case: rows whose pairs leave the pattern (five far pairs on a
tridiagonal precision, confirmed absent from the lookup) and TMB's routine
mocked away; both give the dense values to `sqrt(eps)` relative. The test
"without the sparse route the skip and the refusal hold" mocks the sparse
route away and sees the old skip on the fitted rows and on new rows.

Changed tests, each because it pinned the behavior this lane removes:
`test-latres.R` "the SEs of many fitted rows are skipped, with a message"
(now: read from the sparse factor, equal to the value below the limit, no
message) and "a joint covariance past the limit is refused, by name" (now:
`frm_lp_basis()` and `frm_joint_cov()` on the fitted rows refuse;
`frm_curve()` on new rows reads); `test-latres-newrow.R` "a standard error
refused above the dense limit says so" (now: equal values, no message) and
the second half of "new rows above the row limit beside another block read
the joint covariance" (now: equal to the value below the limit).

Core suite, one file per process, 3 at a time (`dev/latperf-par.sh`, runner
`dev/latperf-run-tests.R`, `NOT_CRAN=true`, log
`dev/latperf-out/latperf-core.log`): ran 309 of 309; pass=22564 fail=1 err=0
skip=239 warn=0; every `lib:` line names wt-latperf-lib. The one failure,
`test-arg-refusal.R` "no registered S3 method has a `...` it never touches",
named the two new methods; they now call `frm_check_dots()`, and the rerun
after reinstall reads pass=122 fail=0 (`latperf-re1.log`, with
test-latres-sparse.R 53/0 and test-message-uniqueness.R 6/0). The 239 skips
are all environment gates, tallied from the DETAIL lines: 224 for
`FRMTMB_BRMS_FIT_TESTS` (brms fits and the Stan tier), 13 for
`FRMTMB_DRMTMB_FIT_TESTS`, 1 scale tier, 1 fuzz tier. The gated latent
files are run gated below.

Extension suites, one file per process (`dev/latperf-out/latperf-ext.log`,
extensions from rellib-r18, core from wt-latperf-lib on all 191 `lib:`
lines): ran 191 of 191, no failure, error or warning anywhere.

<!-- from dev/latperf-out/latperf-ext.log -->
| package | pass | fail | err | skip | warn |
|---|---|---|---|---|---|
| frmtmb.coupling | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 683 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 3962 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 552 | 0 | 0 | 0 | 0 |

Gated, `FRMTMB_BRMS_FIT_TESTS=true` (`dev/latperf-out/latperf-gated.log`):
test-latres-brms.R pass=281 fail=0 err=0 skip=0 warn=0;
test-latres-newrow-brms.R pass=2 fail=0 err=0 skip=0 warn=0.

All 50 extension skips are environment gates (`FRMTMB_BRMS_FIT_TESTS` 36,
`FRMTMB_SCALE_TESTS` 14, from the DETAIL lines).

OpenBLAS 0.3.26, BLAS and LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
then `dev/latperf-openblas-tests.sh`, `OPENBLAS_NUM_THREADS=4`, log
`dev/latperf-out/openblas-0.3.26/`): test-latres-sparse.R 53/0/0/0/0,
test-latres.R 298/0/0/0/0, test-latres-newrow.R 38/0/0/0/0,
test-predict-lp-basis.R 29/0/0/0/0, test-curve.R 91/0/0/0/0,
test-latres-brms.R 30/0/0/5 (ungated)/0 (pass/fail/err/skip/warn).
`dev/ciharden-out/` deleted afterwards. 0.3.32 not run: no test of this lane
sits at a rounding-decided branch; every comparison is a ratio at
`sqrt(eps)` against values that agree to 1e-13.

`R CMD check --as-cran` of core (`dev/latperf-check.sh`, built and checked in
`dev/latperf-check/`, deleted afterwards; log `dev/latperf-out/check.txt`):
`Status: 1 NOTE`, the expected "Skipping checking math rendering: package
'V8' unavailable" on the HTML manual. Examples 114 s, tests OK.

All suites in one session, three orders (`dev/latperf-sessions.sh latperf-sess 2`,
runner `dev/latperf-run-session.R`, log `dev/latperf-out/latperf-sess.log`):
LATPERF-SESS ran 27 of 27; no failure, error or escaped
warning in any session, and 0 LEAK lines without `[with load of ...]`.

<!-- generated from dev/latperf-out/latperf-sess-sessions/*.txt -->
```
SESSION frmtmb forward files=296/309 pass=22565 fail=0 err=0 skip=239 warn=0
SESSION frmtmb preload files=296/309 pass=22808 fail=0 err=0 skip=241 warn=0
SESSION frmtmb reverse files=296/309 pass=22715 fail=0 err=0 skip=239 warn=0
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam preload files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam reverse files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward files=12/12 pass=683 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=12/12 pass=683 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=12/12 pass=683 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn preload files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn reverse files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward files=80/85 pass=3962 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.sample preload files=80/85 pass=3962 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.sample reverse files=80/85 pass=3962 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.spatial reverse files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
```

`files=296/309` and `80/85` count files with at least one test record; the
13 and 5 others skip whole at their gate. After the last cosmetic edits (two
test titles, one wrapped line) the tree was reinstalled and
`test-latres-sparse.R` 53/0 and `test-latres.R` 298/0 rerun
(`dev/latperf-out/latperf-re2.log`).

## Equivalence record and ported brms ledger

The pairs of `dev/equiv-data/results.tsv` that name `ar()`, `ma()`, `arma()`,
`cosy()` or `unstr()` (`dev/latperf-equiv-pairs.R`): p028, p112, p128, p132,
p154, p156, p167. Their levels are L1 (log density and gradient) and L4
(draws), or a refusal; none reads a standard error of a prediction, and this
change moves no log density, fit or draw. So no pair is reached and none was
remeasured. `dev/brmsport-verdicts.tsv` has no row on a latent residual; no row
is known to flip. No compat rule changed.

## What I did not do, and why

- Fits without a latent residual keep the dense `get_joint_cov()`. The review's
  control, `(1 | id)` with one level per row, has the same cost (10.9 s at
  2000 rows, `dev/reviews/2026-11-02-latres.md`, M7). `R/joint-sparse.R` reads
  nothing latent-specific, so it could serve every fit with many group
  effects, but that moves every fit's standard errors by rounding across nine
  suites, and lanes `readers` and `mfxse` work in the same readers this round.
  It is a decision for the user, not a fix for this lane.
- The category probabilities of an ordinal family's fitted rows
  (`fit_fd_se()` with every latent differenced through `fd_eta_chain()`), and
  the stacked `cs()` layers, still form a dense block over every latent: above
  3000 rows they skip, with the message, as before. A sparse version needs the
  chain-rule columns of `fit_fd_se()` held sparse.
- A fit whose covariance of the estimates is held as a factor
  (`cov_fixed_root`, from an ill-conditioned design, or `joint_root` under
  REML) keeps the dense route and both limits. The factor route's rounding is
  the reason that path exists; mixing it with the Schur complement was not
  measured.
- `frmtmb.latres_dense_rows` is not removed. It now limits the dense blocks a
  function returns or must form (`frm_joint_cov()`, `frm_lp_basis()` on the
  fitted rows, the ordinal route below 3000 rows), counted as
  `limit^2` numbers. `frmtmb.latres_se_rows` stays for the ordinal route and
  for the Laplace check, whose first-order correction (`latres_tk()`) forms
  the latents' covariance `v` densely: O(n^2) by its formula, not by a reader.
- REML's `fit_draw_space()` inverts the joint precision densely through
  `joint_cov_whole()` for `predict()` draws and `fit_fd_se()`; a REML fit
  with a latent residual at many rows pays that. Not measured; pre-existing.
- The selected inverse is TMB's unexported `solveSubset()`, reached by
  `get0()` on TMB's namespace, as `R/zzz.R` reaches `isNullPointer`. If a TMB
  release drops it, every fitted row is solved directly (tested by mocking),
  which is O(n^2) time in 32 MB chunks: correct, and slow at 32000 rows.

## Defects found, not fixed

- None in the readers beyond the items above. The 8000-row dense base fit took
  8.39 s against 1.50 s with the default limits, because the Laplace check
  reads `frmtmb.latres_se_rows` too: moving that option to read the fitted rows
  densely also moves the check onto its dense correction. That coupling goes
  away with this change for the readers, and stays for the check.

## Bump and floor

Minor: two behavior changes are marked BREAKING in NEWS (standard errors of
fitted rows above 3000 rows, and of new rows beside another block above the
dense limit, were `NA` and are now numbers). No extension changes, and no
extension's code calls a changed function: every extension suite passed from
rellib-r18 against this core, per file and in session. So no extension floor
moves. Core's own dependencies are unchanged: TMB's `solveSubset()` exists in
the TMB 1.9.25 of the user library, and its absence falls back to direct
solves, so no TMB floor is added.

## Files touched

R: `R/joint-sparse.R` (new), `R/predict.R`, `R/autocor-latent.R`,
`R/autocor.R` (roxygen), `R/se-reparam.R`, `R/se-check.R`, `R/confint.R`,
`R/methods-fit.R`. Generated: `NAMESPACE` (two S3 methods),
`man/frm_joint_cov.Rd`, `man/frm_lp_basis.Rd`, `man/frmtmb-autocor.Rd`.
`NEWS.md`. Tests: `tests/testthat/test-latres-sparse.R` (new),
`tests/testthat/test-latres.R`, `tests/testthat/test-latres-newrow.R`.
`DESCRIPTION` is unchanged in content (roxygen 8.1.1 bumped its version line;
reverted). Scripts in `dev/`: `latperf-prof.R`, `latperf-prof-two.R`,
`latperf-agree.R`, `latperf-reml.R`, `latperf-pattern.R`,
`latperf-selinv-probe.R`, `latperf-ordprobe.R`, `latperf-equiv-pairs.R`,
`latperf-install.sh`, `latperf-check.sh`, `latperf-openblas-tests.sh`, and
lane copies of the release harness with rellib-r18 as base:
`latperf-run-tests.R`, `latperf-par.sh`, `latperf-run-session.R`,
`latperf-sessions.sh` (the originals in `dev/release/` still name rellib-r17;
not changed, that is the consolidation's to update).

## Fix round (review of 2026-11-06)

Changed in this round: `R/joint-sparse.R` (M1, m3, m5, m6 comment),
`R/predict.R` and `R/autocor.R` (roxygen only: m4, m5, M1), `NEWS.md`,
`tests/testthat/test-latres-sparse.R`, `test-latres.R`,
`test-latres-newrow.R`, and the new `tests/testthat/helper-latres-ref.R`.
New scripts: `dev/latperf-split.R` (from the reviewer's pattern script),
`dev/latperf-mutant.R`, `dev/latperf-unstr.R`.

### M1: hurdle and zero-inflated families are linear again

The pairs outside the factor's pattern now go to `jcs_W_pairs()`. Each pair
takes its column of `W` at the end that more pairs share (a level of `hu` or
`zi` over a latent). Each column is solved once, in 32 MB chunks, and the rows
read `W[latent, level]` from those columns. `dev/latperf-split.R <n>
<design>`, seed 1, the reviewer's designs (one poisson series,
`ar(cov = TRUE)` on `mu`, `(1 | id)` with 50 levels on the gate), logs
`dev/latperf-out/split-*.txt`:

| n | design | `fitted()` before (review) | after | rows solved one by one | columns solved |
|---|---|---|---|---|---|
| 1500 | `zi ~ (1 \| id)` | 0.28 s, 713 rows | 0.24 s | 0 | 50 |
| 1500 | `hu ~ (1 \| id)` | 0.41 s, 1500 rows | 0.19 s | 0 | 50 |
| 8000 | `zi ~ (1 \| id)` | 5.26 s, 3935 rows | 0.70 s | 0 | 50 |
| 8000 | `hu ~ (1 \| id)` | 5.06 s, 8000 rows | 0.78 s | 0 | 50 |
| 32000 | `zi ~ (1 \| id)` | not run | 3.45 s | 0 | 50 |
| 32000 | `hu ~ (1 \| id)` | 57.48 s, 32000 rows | 3.25 s | 0 | 50 |

Against the dense route on the same fit at 1500 rows, `fitted()` Est.Error
max |rel|: 6.66e-16 (hurdle) and 6.66e-16 (zi). Every `Est.Error` is finite
at every size. The "before" column is the reviewer's measurement
(`dev/review-latperf-out/pattern-*.txt`), under another load on the machine.

Pinned by two tests in `test-latres-sparse.R`:

- "a pair outside the factor's pattern costs one solve per column": a
  tridiagonal precision with five far pairs gives one call and at most 5
  columns.
- "a hurdle's or a zi's pairs cost one solve per level, not per row": 400
  rows, 12 levels, both families. `jcs_rowvar_solve()` is mocked to stop, the
  columns are at most the levels, and the values are checked against the
  dense route.

The mutant `perrow` (the route before the fix, `dev/latperf-mutant.R`) fails
both: `test-latres-sparse.R` pass 56, fail 6, err 1.

### m1 and m2: tests against an independent reference

`tests/testthat/helper-latres-ref.R` forms the joint covariance densely from
`joint_precision()` alone. It is `solve(Q)`, or, where the SE check kept a
propagation covariance, `joint_cov_repair()`'s block formula written out
again from `sdr$cov_fixed_prop`. The edited assertions compare with it:

- `test-latres.R` now uses 240 rows, above `jcs_block_max`, so the fitted rows
  take the selected inverse.
- `test-latres-newrow.R` compares with it too. For the response-scale new
  rows of its third edit, the reference is the dense route on the same fit.

New test in `test-latres-sparse.R`: "the standard-error check's propagation
covariance is the one read". It uses the reviewer's `boundary` design (seed
20261106, 420 rows; the sd of `(1 | id)` at 0, `cov_fixed_prop` kept), at a
tolerance of `1000 * .Machine$double.eps`. That is below the 2.7e-10 by which
reading `solve(Q)` moves it. Mutants, one process each
(`dev/latperf-mutant.R`, logs `dev/latperf-out/mut/`), pass / fail:

| mutant | `test-latres-sparse.R` | `test-latres.R` | `test-latres-newrow.R` |
|---|---|---|---|
| none | 66 / 0 | 299 / 0 | 39 / 0 |
| `novf`: drop `g' Vf g` | 55 / 11 | 297 / 2 | 38 / 1 |
| `norepair`: ignore the SE check's covariance | 64 / 2 | 299 / 0 | not run |
| `perrow`: one direct solve per row | 56 / 6, 1 err | not run | not run |

`norepair` cannot fail `test-latres.R`, because none of its fits keeps a
propagation covariance. Before the fix round, `novf` passed `test-latres.R`
and `test-latres-newrow.R` whole (the review's table).

### m3: the selected inverse's lookup is not kept

`jcs_selinv_keys()` is rebuilt by each reader (0.14 s at 32000 rows, the
review's measure), and the memo is gone from the object.
`dev/latperf-prof.R`, seed 1, logs `dev/latperf-out/fix-lane-*.txt`, saved
fit in MB:

| n | design | before `fitted()` | after, with the lookup | after, now |
|---|---|---|---|---|
| 2000 | one | 5.12 | 6.02 | 5.73 |
| 8000 | one | 9.56 | 13.43 | 12.27 |
| 32000 | one | 27.75 | 43.16 | 38.52 |
| 8000 | two | 11.55 | 19.17 | 15.67 |
| 32000 | two | 35.59 | 65.98 | 51.92 |

At 32000 rows, design two, `Vjoint` is 11.63 MB (factor 9.73, `J` 1.22;
`dev/latperf-out/fix-rprof-32000-two.txt`). `fitted()` takes 3.00 s there,
and every value equals the earlier run to 10 digits (0.5557698336,
0.7049856461). NEWS now says that the 500 MB at 8000 rows was with both
limits moved, and that with the default limits the base fit stayed at 9.6 MB
with every SE `NA`.

### m4 and m5

NEWS and `?frm_lp_basis` now say that the trim drops every column that no row
loads: a latent, an unnamed level, or a coefficient whose covariate is 0 on
every row.

The refusal of a dense block now states the peak: "about X GB, and about 4X
GB while it is formed". The factor 4 is the review's 4.0 to 4.1 times at 1204
and 3004 positions, and `?frm_joint_cov` says the same. I did not lower the
peak. The solved columns, `J Vf J'`, the block and its copy are each the
size of the block, and removing one copy does not justify a rewrite of `[`.

### m6: the Schur complement that base `solve()` refuses

I kept the fallback and say so in the code. Below base `solve()`'s rcond
limit, the fit takes the dense route, which inverts the whole `Q` as before.
On the reviewer's case (`unstr(t, g)` with `(1 | id)`, 95 lost SEs, rcond
1.6e-101; `dev/latperf-unstr.R`, logs `dev/latperf-out/unstr-*.txt`), the
lane leaves the sparse route and reads what the base build reads. Both give
420 of 420 finite, first 0.6352123012, median 0.4812177155, sum
202.3574618. I did not use `solve(S, tol = 0)`: there the inverse of an
indefinite Schur complement (eigenvalue -9.9e-9) is rounding, and the dense
route is the reference that the base build gives.

### m7: TMB floor

The signature `solveSubset(Q, L = Matrix::Cholesky(Q, super = TRUE, perm =
TRUE), diag = FALSE)`, with `stopifnot(is(L, "dCHMsuper"))` and the result in
the original order, is in `TMB/R/utils.R` of TMB 1.7.14, 1.7.20, 1.9.0 and
1.9.7. I read the CRAN archive sources and then deleted them. frmtmb needs
RTMB (>= 1.9), and RTMB 1.9 imports TMB (>= 1.9.7), so the effective floor
already has the routine, and no `TMB (>= x)` is added. `get0()` stays, and if
the routine is absent the code still falls back to direct solves.

### Item 8: how a caller reads the sparse form (for lane mfxse)

On a fit with a latent residual, `get_joint_cov(fit)` returns the usual list
with `root = NULL`, `Vrest = NULL` and `V` of class `frm_joint_sparse`.

- `jc$V[i, j]` (and `drop =`) returns a dense block with the dense route's
  dimnames. It costs one sparse solve per distinct inner position on the
  shorter side, about `nr * 8` bytes each, plus `J Vf J'` for the block.
- A block of `length(i) * length(j)` above
  `getOption("frmtmb.latres_dense_rows", 20000)^2` numbers stops with class
  `frmtmb_latres_dense`. The peak is about four times the block.
- `jc$V[, apos]` with `apos` a few levels is `p x levels`: cheap in time, and
  `p` rows of memory.
- `diag(jc$V)` does not work: use `jc_diag(jc, idx)`.
- `as.matrix(jc$V)` forms all of `V`, under the same limit.
- For `diag(G V[pos, pos] G')` use `jc_rowvar(jc, G, pos)`. For
  `G V[pos, pos] G'` use `jc_quad(jc, G, pos)`, which is for a few rows. Both
  read the sparse form at linear cost and take a dense or sparse `G`.
- `jc_is_sparse(jc)` tells the two forms apart.

Do not form an inverse of the joint covariance densely from `jc$V`. The joint
precision is `joint_precision(fit)`, and it is sparse. The reads of
`jc$V[, apos]` in `interop_vcov_re()` will go through `[`. They can now stop
with `frmtmb_latres_dense` outside its `tryCatch()`, at more than the square
of the limit.

### Tests of the fix round

One file per process, `NOT_CRAN=true`. Log
`dev/latperf-out/latperf-fix2.log`; every `lib:` line names wt-latperf-lib.
Counts are pass/fail/err/skip/warn:

| file | result |
|---|---|
| test-latres-sparse.R | 66/0/0/0/0 |
| test-latres.R | 299/0/0/0/0 |
| test-latres-newrow.R | 39/0/0/0/0 |
| test-predict-lp-basis.R | 29/0/0/0/0 |
| test-arg-refusal.R | 122/0/0/0/0 |
| test-bracket-access.R | 36/0/0/0/0 |
| test-message-uniqueness.R | 6/0/0/0/0 |

The code changed only in `R/joint-sparse.R`. `R/predict.R` and `R/autocor.R`
changed in roxygen only, so I did not rerun the core suite. The tree was
reinstalled after the last roxygen edit.

### After the short check (SHIP)

- s2: the repair-branch test now asserts `cov_fixed_prop` is present before
  the value check, so a fit that stops inside the boundary fails it.
- s1, not fixed, filed in `dev/test-backlog.md` ("Filed in the round of
  2026-11-06"): with many levels on `hu` or `zi` the cost is about
  levels x rows, one column of `W` per level. At 32000 rows, 3.25 s at 50
  levels, and (the check's measures) 3.5 s at 500, 6.8 s at 4000, 24.1 s at
  16000. The route: add the pairs the rows load to `Q_rr` as structural
  zeros before factoring.
