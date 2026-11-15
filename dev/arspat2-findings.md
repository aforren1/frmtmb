# Lane arspat2: residual terms and spatial

Round opened 2026-11-13, base ec8178d3 (frmtmb 0.82.0), branch
wt-r83-arspat2. Brief: `round-briefs/r083-arspat2.md`. Status: DONE
for 5289 (and two related silent defects found here), 5264, 4002; NOT
DONE, with reasons, for 4483 (blocked on lane priors2's class rename),
5256 (scoped and filed), 3979. Fix round after review
`dev/reviews/2026-11-13-arspat2.md` (FIX): done, section "Fix round"
below. The first pass's margin construction (item 1, "The fix") and
logit interval (item 2) are WITHDRAWN there; their sections stay as the
record of the first pass and are marked.

## Fix round (review FIX: F1 HIGH, F2 HIGH, F3 LOW, F4)

### F2: the margin, now analytic, with the missing values out of the fit

First pass (withdrawn): the missing values stayed latent and a standard
normal device in them made the joint Laplace exact. The device coupled
every missing value with every other (`d = G^-1 K_m'u`), so the inner
Hessian was dense `m x m`, and its tape grew past m^2: the reviewer's
n = 400, m = 300 fcor fit hit `std::bad_alloc` at 37 GB and surfaced as
an NA/NaN gradient.

Now (the reviewer's suggestion (a)): a missing value of a registered
term's response that no other density reads LEAVES the fit
(`mi_rows_held()` in R/frame.R: a seam response no longer holds its
rows; `mi_rows_joint()` still keeps every row's predictor read, since
the lag SAR's mean is `A^-1 mu`). The objective
(`autocor_seam_margin_loglik()`, R/autocor-seam.R) writes the exact
margin of the remaining rows with the missing rows filled by `mu` (any
value works): `Q_o = |u|^2 - v'G^-1 v`, log scale determinant
`logjac - 0.5 log|G|`, `G = K_m'K_m`, `v = K_m'u`. Nothing is latent, so
the inner Hessian holds the random effects only. One rule for gaussian
and student, with or without other integrated values: the 0.82.0
constant (`autocor_seam_t_laplace()`), `autocor_seam_t_missing()`,
`autocor_seam_scale_integrated()` and the `mi_inner_beta` flag in
R/fit.R are gone (R/fit.R is back to ec8178d3). A value another density
reads stays latent as before; the sampler holds every value
(`mi_sampling()`).

`G` and `v` come from a new optional seam slot,
`frmtmb_register_autocor(margin = function(rows, aux))`, run once per
objective and returning `function(u, sigma, theta) list(G, v)`:

- `fcor_margin()`: `G = [M^-1]_mm` scaled by `1/sigma` both sides
  (data, read once), `v = L^-1[, m]'u / sigma`.
- `sar_margin()` (frmtmb.spatial): `G = (I - rho (W_mm + W_mm') +
  rho^2 W_m'W_m) / sigma^2`, three products of `W` read once, and
  `v = (u_m - rho W_m'u) / sigma`.
- Without a slot: `K_m` from one `whiten()` call per missing value
  (`autocor_seam_margin_cols()`), on the tape. Filed LOW.
- A given slot is checked when the objective is built, at two columns
  and a generic point (sigma 1.37 or `exp(0.3 sin(i))`, theta moved by
  0.3), against those `whiten()` calls, and refused by the term's name
  when it disagrees (test: a slot that forgets sigma in `G` is refused;
  the right one fits and equals fcor()).

Tried and dropped (no gain, measured on sar n = 400, m = 300): only the
nonzeros of `G` on the tape (6.39 to 6.62 s, peak 0.62 to 0.64 GB) and
`solve(G) %*% v` in place of `solve(G, v)` (6.00 to 7.19 s, 0.61 GB).

Exactness, objective against the observed-rows fit or the dense
reference:

- `dev/arspat2-5289-probe.R` (log `dev/arspat2-out/5289-probe-fix.txt`):
  largest |diff| 4.3e-14 over the five cases (RE sd 0 and 1.5, no RE,
  REML, profile), three points each.
- `dev/arspat2-5289-sar.R` (`5289-sar-fix.txt`): fcor, sar error, sar
  lag against the dense t margin with RTMB's Laplace over b: largest
  |diff| 4.3e-14.
- The reviewer's cases, `dev/arspat2-rev-margin.R` and
  `dev/arspat2-rev-sar.R` (copies of `dev/reviews/arspat2-review-*.R`
  with `dev/arspat2-rev-shim.R` standing in for the removed
  `autocor_seam_t_missing()`; logs `rev-margin-fix.txt`,
  `rev-sar-fix.txt`): every margin case at or below 5.7e-14 absolute
  except the near-singular M (kappa 2.2e4; 1.2e-12) and gaussian
  sigma ~ (1|p|h) (1.5e-10 on -130.06, relative 1.1e-12); binary-W sar
  against the dense reference, lag relative 1.4e-11, 1.9e-12, 1.2e-13,
  error 4.2e-13, 0, 2.4e-13 (the reviewer's first-pass figures to the
  digit); AD gradient against finite differences relative 2.5e-8 (sar
  lag), 3.9e-9, 5.1e-10, 1.3e-9; multivariate sums -1.1e-8 and -9.2e-10;
  `sigma ~ s(z)` and `gp(z)` take the margin and are exact.
- Fits the change reaches without a model change
  (`dev/arspat2-rev-ident.R`, `rev-ident-fix.txt`, against the
  reviewer's `ident-base.rds`): the four `mi()` seam fits move by
  rounding only (logLik within 2.8e-14, parameters within 1.9e-12),
  since the margin replaces latent values that were exact there; the
  other five are bit-identical in parameters, logLik, SEs and messages;
  the sar rows of `confint_varcorr()`/`summary()` change type only.

Cost, `dev/arspat2-cost2.R` (the reviewer's designs: rook lattice,
M = solve(crossprod(I - 0.4 W)), `(1 | h)` 8 levels, student, gamma(2,
0.1) on nu, data seed 11, missing seed m; ONE fit per process, so the
peak working set is that fit's). Lane and base measured here except the
two fcor n = 1000 base rows, taken from the reviewer's
`cost-sar-base.txt` (process peak; a base fit at n = 1000, m = 200 needs
25 GB, above this round's 12 GB probe limit):

| n | m | term | base s | lane s | ratio | base peak GB | lane peak GB |
|---|---|---|---|---|---|---|---|
| 400 | 100 | sar | 1.01 | 1.06 | 1.05 | 0.37 | 0.37 |
| 400 | 200 | sar | 1.91 | 2.86 | 1.50 | 0.42 | 0.48 |
| 400 | 300 | sar | 3.31, 3.42, 3.45 | 5.72, 6.52, 6.69 | 1.73 to 1.95 | 0.48 to 0.49 | 0.60 to 0.61 |
| 1000 | 50 | sar | 2.09 | 1.85 | 0.89 | 0.44 | 0.35 |
| 1000 | 200 | sar | 4.54 | 4.59 | 1.01 | 0.66 | 0.52 |
| 400 | 100 | fcor | 16.63 | 4.14 | 0.25 | 2.88 | 0.66 |
| 400 | 200 | fcor | 26.20 | 6.56 | 0.25 | 5.34 | 0.80 |
| 400 | 300 | fcor | 48.97 | 12.04 | 0.25 | 7.06 | 1.11 |
| 1000 | 50 | fcor | 47.25 (reviewer) | 19.20 | 0.41 | 8.60 (reviewer) | 1.92 |
| 1000 | 200 | fcor | 196.61 (reviewer) | 25.31 | 0.13 | 25.14 (reviewer) | 2.24 |

Every lane fit: convergence 0, finite gradient at the optimum. The one
cell near the acceptance line is sar at n = 400, m = 300 (75 percent
missing): time 1.7 to 1.95 times base, peak 0.60 to 0.61 GB against
0.48 to 0.49 (+22 to +27 percent, +0.12 GB). The extra is the dense
300 x 300 solve and determinant of `G` on the tape; base held those
values in a sparse inner Hessian, but was not exact there (its logLik
-230.55 against the margin's -185.37 in this design). At m = 200 the
same design is 1.5 times and +14 percent.

bad_alloc: the margin route builds no inner Hessian in the missing
values, so its inner tape no longer grows with m (largest lane peak
2.24 GB). TMB's habit of printing `std::bad_alloc` and returning NaN
from `obj$fn` is generic to any model whose inner Hessian tape is too
large, and catching it is a fit-end change in R/fit.R (other lanes'
area); filed LOW with the reviewer's reproducer rather than built.

### F1: sar's rho interval: the raw Wald interval, clamped

First pass (withdrawn): a Wald interval on the logit of `(lo, hi)`. It
covered 0.665 at rho = 0.9 on a 5 x 5 lattice (the reviewer's
`arspat2-review-cover.R`).

Now: type `"clamped"`. `autocor_trans_rows()` keeps rho on its own scale
with the delta-method SE (bit-identical to 0.82.0's), carries `lo`/`hi`,
and `varcorr_untrans()` maps `"clamped"` to `pmin(pmax(v, lo), hi)`;
`varcorr_nat_deriv()` is 1 there (its `lo`/`hi` arguments are gone again).
`frm_multiple()` pools rho on its own scale (Rubin) and clamps the
pooled interval's ends: the same scale as the single-fit interval.

The reviewer's design rerun on the lane (200 seeds each; logs
`dev/arspat2-out/cover-0.9.txt`, `cover-0.5.txt`):

| lattice | true rho | coverage clamped | coverage 0.82.0 raw | upper end at hi (clamped) |
|---|---|---|---|---|
| 5 x 5 | 0.9 | 0.935 | 0.935 | 133 |
| 10 x 10 | 0.9 | 0.930 | 0.930 | 4 |
| 5 x 5 | 0.5 | 0.925 | 0.925 | 14 |
| 10 x 10 | 0.5 | 0.925 | 0.925 | 0 |

Identical coverage, as it must be: the clamp moves only ends past
`(lo, hi)`, where the true rho cannot be. Test (frmtmb.spatial
test-sar.R) "rho's interval is the raw Wald interval held inside its
range": seed 3 (raw upper end past hi, held at hi exactly) and seed 1
(inside, the raw interval), against the delta method from `cov.fixed`.

### F3: notes, comments, flag

- Both `mi()` compat notes (R/compat.R fcor x mi(), R/autocor-seam.R
  registered term) now say: a value no other density reads leaves the
  fit, which takes the exact margin; a value another density reads stays
  latent under Laplace.
- The `R/fit.R` comment and the `mi_inner_beta` flag are gone: the
  margin is exact whatever else is integrated, so no flag changes the
  answer. (The reviewer was right that `profile` maximizes beta with no
  Laplace term; the first pass's "profile integrates beta" was wrong.)
- Cost statements in NEWS, `?frmtmb_register_autocor` ("The margin"),
  `?frmtmb-autocor` (fcor) and `?sar` now quote the measured numbers.

### F4: frmtmb.spatial's floor

frmtmb.spatial now calls `frmtmb_register_autocor(margin =)`, which
0.82.0 does not have: installed on a 0.82.0 core it fails at load
(unused argument), and its new tests fail there. Consolidation MUST raise
`Depends: frmtmb (>= ...)` to the new core's version. Its NEWS says it
needs the frmtmb that takes `margin =`. The spatial tests of this fix
round ran against the lane's build of frmtmb.spatial (every `lib:`
line).

### Fix round 2 (short check `dev/reviews/2026-11-13-arspat2-short.md`)

- N1: the core NEWS bullet said the margin is "the fit over the observed
  rows with the term's matrix cut to them". True for `fcor()` only; for
  `sar()` the margin is the exact marginal of the observed rows under the
  full-lattice model, which is not `sar(W_oo)` (short check: lag -53.637
  against -55.204, rho 0.50 against 0.10). Reworded per term, core and
  spatial NEWS; my backlog CLOSED line for 5289 too.
- N2: the build-time slot check compared `G` at {1, m} x {1, m} and `v`
  at {1, m}; a sar slot right only for a symmetric W passed on a
  row-standardized W, 0.011 nats off. `autocor_seam_margin_setup()` now
  compares every entry of `G` and `v` against `whiten()` (m calls, on
  numbers, once per objective), relative tolerance 1e-8 of the largest
  entry. Cost (`dev/arspat2-slotcheck-cost.R`, minimum of 5): the whole
  setup with the full check is 0.010 s (one 10 ms tick) for sar error
  n = 400, m = 300 and for fcor n = 400, m = 200, against a fit of 3.49 s
  and 3.59 s. The roxygen, `?frmtmb_register_autocor` ("The margin") and
  NEWS say the check is in full at one value of the parameters and is
  not a proof for every value. New spatial test "a margin slot right
  only for a symmetric W is refused" (the short check's slot and data):
  seen failing on the pre-fix lane build (test-sar.R pass=178 fail=1,
  `dev/arspat2-out/fx2-sar-before.txt`), passing after; frmtmb.spatial's
  own slot is the absent case.
- N3: the NEWS error figures now give the measured range on 0.82.0,
  0.21 nats (gaussian, random effect in sigma) through 4.66 (student
  fcor, half of 60 rows missing) to 45.17 (student sar, 300 of 400
  missing).
- LOW: both spatial NEWS bullets are marked BREAKING (`estimates$miss`;
  interval type `"raw"` to `"clamped"`).
- The short check's cost table, recorded: sar n = 400 at m = 100, 200,
  300 took 1.15, 1.67 and 1.90 to 2.06 times base, with peaks 0.38, 0.50
  and 0.64 GB against 0.37, 0.43 and 0.49 (+31 percent, +0.15 GB at
  m = 300); fcor 0.25 times base and 4.5 times smaller. Accepted by the
  coordinator: base was 45.17 nats off at sar m = 300. The NEWS sar cost
  sentence now says 1.7 to 2.1 times and about 30 percent more memory.
- Runs on the lane library (`dev/arspat2-out/fx2-test-*.txt`):
  test-autocor-seam.R 90, test-autocor-se-mi.R 85, frmtmb.spatial
  test-sar.R 179 + 4 skip, test-car.R 132; fail 0, error 0, warn 0 each.
  `tools::checkRd()` is clean on the three changed Rd files. R CMD check
  was not rerun: R/ changed only in the check and its roxygen.

### Fix-round tests

Reference BLAS, lane core and lane spatial (`dev/arspat2-out/fx-*.txt`):

| file | pass | fail | error | skip | warn |
|---|---|---|---|---|---|
| test-autocor-se-mi.R | 85 | 0 | 0 | 0 | 0 |
| test-autocor-seam.R | 90 | 0 | 0 | 0 | 0 |
| test-autocor.R | 198 | 0 | 0 | 0 | 0 |
| test-mi-families.R | 194 | 0 | 0 | 0 | 0 |
| frmtmb.spatial test-sar.R | 177 | 0 | 0 | 4 | 0 |
| frmtmb.spatial test-car.R | 132 | 0 | 0 | 0 | 0 |

OpenBLAS with LAPACK (`dev/ciharden-openblas.sh <ver> lapack`, then
`dev/arspat2-openblas.sh`, OPENBLAS_NUM_THREADS = 4; logs
`dev/arspat2-out/ob0.3.26-*.txt`, `ob0.3.32-*.txt`,
`ob0.3.26-Haswell-*.txt`, the last with OPENBLAS_CORETYPE=Haswell; every
`lib:` line the lane's core and spatial): 0.3.26, 0.3.32 and 0.3.26
Haswell each read exactly the reference-BLAS table above (85, 90, 198,
194, 177 + 4 skip, 132; fail 0, error 0, warn 0). `dev/ciharden-out/`
deleted.

Seen failing on rellib-r21 (`fxbase-*.txt`): test-autocor-se-mi.R
pass=59 fail=26 error=0 (the margin assertions: values still latent,
and the objective points); test-autocor-seam.R pass=86 fail=2 error=1
(the two `unread` assertions; the slot test errors at the unknown
`margin =` argument, the weak form for a new argument); spatial
test-sar.R pass=162 fail=15 (the `unread` assertions, the student RE
margin points, and the clamp test's type).

Full suites after the fix, one file per process (`bash
dev/release/run-par.sh arspat2-fix dev/arspat2-out/jobs-fix.txt 3`,
log `dev/release/arspat2-fix.log`; core and spatial the lane's,
frmtmb.sample rellib-r21's on the lane's core):

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 330 | 24430 | 0 | 0 | 240 | 0 |
| frmtmb.spatial | 8 | 549 | 0 | 0 | 6 | 0 |
| frmtmb.sample | 91 | 4238 | 0 | 0 | 18 | 0 |

The machine's C: drive filled during this run (4.0 GB, then 2.9 GB
free, other lanes' output). Two core files were hit: test-negbinomial2.R
read warn=8 ("cannot create dir ... No space left on device") and
test-nl-body-vars.R died with "Fatal error: cannot create 'R_TempDir'"
(the log's "ran 428 of 429"). Both were rerun alone once space was back
(135 GB free): 69 and 23 passes, warn 0 (`dev/arspat2-out/rerun-*.txt`);
the table counts the reruns. Gated frmtmb.spatial test-sar.R:
pass=197 fail=0 skip=0 (`t-sar-gated-fix.txt`).

R CMD check --as-cran after the fix (`bash dev/arspat2-check.sh`, logs
`dev/arspat2-check/check-*.txt`, `.Rcheck` trees deleted): frmtmb
Status: 1 NOTE (V8 math rendering), tests 51 min OK; frmtmb.spatial
Status: 1 NOTE (V8), tests 62 s OK. The sessions (run-sessions.sh) were
not rerun in the fix round; the coordinator's list did not ask for them.

The removed test assertions: `autocor_seam_t_missing()` route flags
(the function is gone) and "the student() margin stays off where it is
wrong" (rewritten as "a registered term's missing values stay latent
where read": another density reads them, and the sampler's frame).

Process note: running the reviewer's `arspat2-review-cover.R` unmodified
wrote its `.rds` files into `dev/reviews/arspat2-review-out/`. It
overwrote the reviewer's `cover-5-0.5-lane.rds` and
`cover-10-0.5-lane.rds` (first-pass lane results; the reviewer's text
log `cover-lane-0.5.txt` still holds their numbers) and created two new
`cover-*-0.9-lane.rds`, which I moved to `dev/arspat2-out/`. My attempt
to stop those processes was denied by the permission system (the runs
had already finished); nothing was worked around.

## Files touched

- Core: R/autocor-seam.R (`autocor_seam_margin_setup()`,
  `autocor_seam_margin_cols()`, `autocor_seam_margin_loglik()`,
  `fcor_margin()`, the `margin` argument and its docs section, the mi()
  compat note; `autocor_seam_t_laplace()` and `autocor_seam_t_missing()`
  removed), R/objective.R (seam branch), R/frame.R (`mi_rows_held()`,
  `mi_rows_read()`'s peel test, `mi_laplace_exact()`), R/autocor.R
  (`autocor_types()`, `autocor_trans_rows()` `lo`/`hi` columns, fcor
  docs), R/confint.R (`varcorr_trans_rows()`'s `add()`,
  `varcorr_untrans()`, `confint_varcorr()`; three small hunks),
  R/methods-fit.R (summary autocor block, `summary_nat_frame()`,
  `varcorr_nat_deriv()`), R/multiple.R (pooled varcorr and a comment),
  R/autocor-latent.R (`ac_report_rows()` columns), R/compat.R (one note
  string), NEWS.md, man/frmtmb-autocor.Rd, man/frmtmb_register_autocor.Rd,
  tests/testthat/test-autocor-se-mi.R, tests/testthat/test-autocor-seam.R,
  dev/test-backlog.md. R/fit.R is unchanged from ec8178d3.
- frmtmb.spatial: R/sar.R (`sar_margin()`, its registration, docs),
  man/sar.Rd, NEWS.md, tests/testthat/test-sar.R.
- dev/arspat2-*: scripts and this record.

## frm_compat() rows

Statuses unchanged. Notes changed: `fcor()` x `mi()` (R/compat.R) and the
registered term's `mi()` row (R/autocor-seam.R
`autocor_seam_compat_rules()`). No equivalence-survey pair is an `mi()`
pair of these terms (`dev/equiv-data/sample.tsv`: p158, p159, p164
only, all refused), so none was remeasured. No ported-ledger row was
expected to flip; the ledger was not regenerated.

## Bump and floor

Core: minor (a new seam argument; BREAKING: `mi()` values of a
registered term leave `fit$estimates$miss`, and the bounded parameter's
interval type and ends). frmtmb.spatial: minor (it registers the new
`margin` slot); its floor on frmtmb MUST be the new core version (F4).

## 1. Backlog 5289: student() with fcor()/sar(), mi() and a random effect

### What was wrong

The Laplace value of a multivariate t integrated over `m` of its `n`
coordinates, AT FIXED VALUES of everything else, is the exact margin
times `exp(c(nu, n, m))`; 0.82.0 took `c` off when the missing values
were the whole inner problem. With another integrated value `b` (a
random effect, or beta under REML or `profile`) TMB's Laplace is joint
over `(y_m, b)`. Its Hessian in `y_m` is `(nu + n)/(nu + Q_o(b)) K_m'K_m`,
which moves with the observed rows' quadratic form `Q_o(b)`, so the
joint value is not the nested one (exact over `y_m`, Laplace over `b`)
plus a constant. 0.82.0 then applied no correction and said nothing.

A second defect of the same gate, found here: under `REML = TRUE` (and
`profile = TRUE`) beta is integrated too, the gate did not know it, and
the constant was applied and `mi_laplace_exact()` said TRUE. The fit
was wrong by 0.28 to 0.38 nats with no RE at all.

### The fix (first pass; WITHDRAWN in the fix round, see above)

`autocor_seam_t_margin()` (R/autocor-seam.R). With `K_m` the Jacobian
columns of `whiten` at the missing rows (`whiten(e_j) - whiten(0)`, one
call per missing value, on the tape), `P` the projection onto them,
`d = (K_m'K_m)^-1 K_m'u`:

    log t_{n-m}(|(I - P) u|^2; nu) + logjac - 0.5 log|K_m'K_m|
      - (m/2) log(2 pi) - 0.5 |d|^2

The first line is the exact margin of the observed rows (the Schur
identity `|Sigma| = |Sigma_oo| |S|`, `S^-1 = K_m'K_m`); `d` is
`y_m` minus its conditional mode, so the second line is a standard
normal in `y_m` with Hessian `I` whatever `b` is. A Laplace over
`(y_m, b)` of `A(b) - 0.5 |y_m - mu(b)|^2` equals the nested one exactly
(Schur complement of the joint Hessian is `A''`), so the fit is the one
over the observed rows.

`autocor_seam_t_missing()` now returns `list(m, rows, margin)` per
response; `margin = FALSE` (the 0.82.0 constant) only when the
response's missing values are the whole inner problem, `TRUE` when the
frame has random-effect blocks, a `b` template, `me()` groups, another
response's missing values, or `frame$mi_inner_beta`, which
`fit_assembled()` (R/fit.R, one line after `mi_ml_frame()`) sets under
REML or `profile`. A charted missing value (`mi_chart_maps()`) gets no
correction. The margin route is exact in every configuration, so a
frame that errs toward it is still right; only the cost differs.

### Measurements

`dev/arspat2-5289-probe.R` (data seed 4, n = 60, m = 20, M = 0.5^|i-j|,
`(1 | h)` 6 levels; logs `dev/arspat2-out/5289-probe-{base,lane}.txt`):
the mi() fit's `-obj$fn(p)` against the fit over the observed rows with
`fcor(M_oo)`, at three outer points with nu = 3, 6, 21.

| case | base diff (3 points) | lane diff (3 points) |
|---|---|---|
| (1 \| h), data RE sd 0 | -1.981, -1.878, -1.493 | -2.8e-14, -1.4e-14, -5.7e-14 |
| (1 \| h), data RE sd 1.5 | -3.074, -2.907, -2.288 | -1.4e-14, +1.4e-14, +2.8e-14 |
| no RE (constant route) | +1.4e-14, +7.1e-15, -2.8e-14 | the same |
| no RE, REML | -0.382, -0.361, -0.284 | -7.1e-15 (each) |
| (1 \| h), profile | -3.077, -2.909, -2.291 | -1.4e-14, +1.4e-14, +1.4e-14 |

Base `mi_laplace_exact()` FALSE with tk -2.9e-9 and -2.5e-9 on the RE
cases (no warning), TRUE on the REML case (wrong).

`dev/arspat2-5289-sar.R` (data seed 7, 7 x 7 rook lattice, 12 missing,
`(1 | h)` 7 levels; logs `dev/arspat2-out/5289-sar-{base,lane}.txt`):
frmtmb's objective against a reference written in the script, the dense
t margin of the observed rows (mean `(B (X beta + Z b))_o`, scale
`sigma^2 [(A'A)^-1]_oo`, `B = A^-1` lag, `I` error) with `b` integrated by
RTMB's Laplace; fcor is the control of the reference's parameter map.

| case | base diff (3 points) | lane diff (3 points) |
|---|---|---|
| fcor | -0.885, -0.832, -0.641 | -1.4e-14, +4.3e-14, -7.1e-14 |
| sar error | -0.885, -0.832, -0.641 | +1.4e-14, -1.4e-14, -2.8e-14 |
| sar lag | -0.885, -0.832, -0.641 | -1.4e-14, +4.3e-14, +1.4e-14 |

The base error is the same for the three terms at a point: the joint t
Laplace error depends on the data only through `Q_o`, and here on
nothing else that differs between them.

Cost (`dev/arspat2-cost.R`, seed 11, n = 400 on a 20 x 20 lattice,
`(1 | h)` 8 levels, prior gamma(2, 0.1) on nu, seconds per `frm()`,
minimum of 3 rounds; logs `dev/arspat2-out/cost-{base,lane}.txt`, the
two run at once in two processes):

| m | fcor mi base | fcor mi lane | sar error mi base | sar error mi lane | fcor on observed rows |
|---|---|---|---|---|---|
| 10 | 2.13 | 2.42 | 0.14 | 0.16 | 1.48 to 1.50 |
| 40 | 3.86 | 7.65 | 0.21 | 0.39 | 1.14 to 1.17 |

fcor's `whiten()` is a dense `n x n` product, so the route costs about
`m n^2` per evaluation; sar's is sparse. Only the newly exact cases pay
it: a fit whose missing values are the whole inner problem keeps the
constant.

### Tests (seen failing on rellib-r21)

- `tests/testthat/test-autocor-se-mi.R`: the end of "mi() beside fcor() on
  student() integrates out exactly" (its absent case with `(1 | g)` now
  asserts the margin route), and two new tests: "mi() beside fcor() on
  student() with more integrated is exact" (RE sd 0 and 1.5 with
  `profile`, REML; objective against the observed-rows fit at 3 points,
  relative 1e-10), "the margin of a registered term stays off where it
  is wrong" (another density reads the values; gaussian; the sampler's
  frame from `mi_sampling()`). Base: pass=54 fail=12 error=2; lane
  pass=73 fail=0 error=0 warn=0 (before the gaussian test below was
  added; the file's final counts are in that section).
- `extensions/frmtmb.spatial/tests/testthat/test-sar.R`: "mi() beside
  sar() on student() with a random effect is the observed rows' margin"
  (lag and error, 3 points each, reference as in the script). Base:
  pass=157 fail=6 (`t-sar-base.txt`, relative errors 0.024 to 0.034);
  lane pass=163 fail=0 skip=4 warn=0.

### The same defect on gaussian(), found here

`gaussian()` is exact over the missing values at fixed `sigma`, but its
Hessian in them is `K_m'K_m`, which moves with `sigma`; with
`sigma ~ 1 + (1 | h)` the joint Laplace is again not the nested one.
`dev/arspat2-sigre.R` (data seed 9, n = 60, m = 15, `fcor(M)`, 6
groups, sd(sigma) held at 0.3, 0.6 and 1.2; logs
`dev/arspat2-out/sigre-{base,lane}.txt`), the mi() fit's `-obj$fn`
against the observed-rows fit:

| family | base diff (sd 0.3, 0.6, 1.2) | lane diff |
|---|---|---|
| gaussian | -1.103, -1.735, -1.992 | -2.1e-13, -3.6e-14, 0 |
| student | -1.103, -1.735, -1.992 | -2.1e-13, -5.0e-14, +7.1e-15 |

The margin route now also takes a gaussian response when
`autocor_seam_scale_integrated()` finds a random-effect block, an `me()`
term or an `mi()` reader in a predictor of the response other than `mu`
(a block whose components carry no key falls back to its `dpar`). Its
absent case, gaussian with `(1 | h)` in `mu` only, keeps the joint
density, which is exact there (asserted against the observed-rows fit).
`sar()` takes one sigma (`sigma_varies = FALSE`), so only `fcor()` and a
registered term with a predicted sigma reach this.

New test "mi() beside fcor() with a random effect in sigma is exact"
(test-autocor-se-mi.R): base adds 8 failures (2 route flags, 6 objective
points); file totals base pass=58 fail=20 error=2, lane pass=85 fail=0
error=0 warn=0 (`dev/arspat2-out/t-semi-{base,lane}.txt`).

### Not done for item 1

- `fcor()`'s margin route costs one dense `whiten()` per missing value.
  A slot giving `K`'s columns directly (or `K` on the tape) would cut it
  to `m n`; filed LOW.
- No L1 row against brms: the change is in the maximum-likelihood
  objective only. The joint density a sampler reads, which is what brms's
  program is, is unchanged (`mi_sampling()` frames take no correction,
  asserted), and its L1 rows are the previous lane's.
- `me()` with the margin route was not measured: the observed-rows fit is
  not the reference there, because the dropped rows' `me()` observations
  still inform the `me()` hyperparameters in the mi() fit.

## 2. Backlog 5264 (= 3992): sar()'s rho interval (first pass; WITHDRAWN in the fix round, see above)

`autocor_types()` marks a registered term's bounded parameter
`"bounded"`; `autocor_trans_rows()` puts it on `qlogis((v - lo)/(hi -
lo))`, its internal parameter, with columns `lo` and `hi` (NA on every
other row, so `varcorr_trans_rows()`'s rows in R/confint.R and
`ac_report_rows()`'s latent rows in R/autocor-latent.R carry them too).
`varcorr_untrans()` (R/confint.R) and `varcorr_nat_deriv()`
(R/methods-fit.R) take `lo` and `hi`. The five sites pass them:
`confint_varcorr()` (R/confint.R), `summary()`'s autocor block and
`summary_nat_frame()` (R/methods-fit.R), `frm_multiple()`'s pooled table
(R/multiple.R), and the row builder in R/autocor-latent.R. `summary()`'s
`Est.Error` is the same delta method on rho as before,
`(hi - lo) p (1 - p) se_t`.

`dev/arspat2-5264.R` (seeds 1 to 20, 5 x 5 rook lattice,
row-standardized, rho = 0.9, lag; logs
`dev/arspat2-out/5264-{base,lane}.txt`): intervals past an end of rho's
interval, base 13 of 20 (largest upper end 1.0576, seed 18), lane 0 of
20. Seed 3: base (0.6459, 1.0398), lane (0.5011, 0.9572), estimate
0.8428.

Test "rho's interval is the Wald interval on the logit of its range"
(frmtmb.spatial test-sar.R, seed 3): the type, both ends inside, the
ends equal to `lo + (hi - lo) plogis(theta +/- z sqrt(V))` from
`cov.fixed` to 1e-8 of the interval's width, the summary row the same,
`Est.Error` the delta method on rho. Base: 4 failures; with the item 1
test the file reads base pass=159 fail=10, lane pass=169 fail=0 skip=4
warn=0 (`dev/arspat2-out/t-sar-{base,lane}.txt`).

## OpenBLAS (first pass)

`dev/ciharden-openblas.sh 0.3.26 lapack`, then `dev/arspat2-openblas.sh
0.3.26` (OPENBLAS_NUM_THREADS = 4), logs
`dev/arspat2-out/ob0.3.26-*.txt`, every `lib:` line the lane's core and
spatial:

| file | pass | fail | error | skip | warn |
|---|---|---|---|---|---|
| test-autocor-se-mi.R | 85 | 0 | 0 | 0 | 0 |
| test-autocor-seam.R | 80 | 0 | 0 | 0 | 0 |
| test-autocor.R | 198 | 0 | 0 | 0 | 0 |
| test-mi-families.R | 194 | 0 | 0 | 0 | 0 |
| frmtmb.spatial test-sar.R | 169 | 0 | 0 | 4 | 0 |
| frmtmb.spatial test-car.R | 132 | 0 | 0 | 0 | 0 |

0.3.32 was not run: no assertion of this lane sits near a
rounding-decided branch (the identities are 1e-10 relative and hold to
about 1e-15). `dev/ciharden-out/` deleted.

## Suites, one file per process (first pass)

`bash dev/release/run-par.sh arspat2-all dev/arspat2-out/jobs-all.txt 3`
with `FRMTMB_RELLIB` = the lane library and `NOT_CRAN=true`, log
`dev/release/arspat2-all.log` ("ARSPAT2-ALL ran 529 of 529", no "NO
RESULT" line). Every `lib:` line names the lane's core; frmtmb.spatial
is the lane's build, the other extensions rellib-r21's (none of their
code changed; they reach the changed core through `summary()`,
`confint_varcorr()` and the objective). Totals summed from the log:

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 330 | 24420 | 0 | 0 | 240 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1851 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 13 | 697 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 91 | 4238 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 541 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 16 | 575 | 0 | 0 | 0 | 0 |

The skips are the gated brms tiers. Gated on the lane build
(`FRMTMB_BRMS_FIT_TESTS=true`): frmtmb.spatial test-sar.R pass=189
fail=0 skip=0 (`dev/arspat2-out/t-sar-gated.txt`); test-brms-likelihood.R
pass=544 fail=0 skip=0 (`t-brmslik-gated.txt`).

## One session per suite (first pass; not rerun in the fix round)

`bash dev/release/run-sessions.sh arspat2-sess 3` (all nine packages,
forward, reverse and preload), `FRMTMB_RELLIB` = the lane library, log
`dev/release/arspat2-sess.log`, "ARSPAT2-SESS ran 27 of 27":

| package | files | forward pass | reverse pass | preload pass | fail, err, warn |
|---|---|---|---|---|---|
| frmtmb | 317/330 | 24420 | 24585 | 24688 | 0 each |
| frmtmb.coupling | 11/11 | 551 | 551 | 551 | 0 each |
| frmtmb.eam | 33/33 | 1851 | 1851 | 1851 | 0 each |
| frmtmb.latent | 13/13 | 697 | 697 | 697 | 0 each |
| frmtmb.learn | 16/16 | 474 | 474 | 474 | 0 each |
| frmtmb.ode | 11/11 | 549 | 549 | 549 | 0 each |
| frmtmb.sample | 86/91 | 4238 | 4238 | 4238 | 0 each |
| frmtmb.survival | 16/16 | 575 | 575 | 575 | 0 each |
| frmtmb.spatial | 8/8 | 541 | 541 | 541 | 0 each |

46 LEAK lines, every one tagged `[with load of ...]`; no untagged LEAK
and no DETAIL line. The files with no row (13 in core, 5 in
frmtmb.sample) are the ones gated as a whole, which read pass=0 skip=0
in the per-file run too.

## R CMD check --as-cran (first pass)

`bash dev/arspat2-check.sh frmtmb frmtmb.spatial`, built with vignettes
in `dev/arspat2-check/`, `R_LIBS` = lane library, rellib-r21, user
library, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE` (logs
`dev/arspat2-check/check-*.txt`; the `.Rcheck` trees and tarballs
deleted afterwards):

- frmtmb: Status: 2 NOTEs. The V8 math-rendering NOTE, and the examples
  timing NOTE (`frm_hazard_reads` 3.38 s CPU, 11.50 s elapsed;
  `dharma_residuals` 1.51 and 6.75): elapsed far above CPU on examples
  this lane did not touch, the load signature of lane-rules.md (14 lanes
  on the machine). Tests 92 min, OK.
- frmtmb.spatial: Status: 1 NOTE (V8). Tests 55 s, OK.

## 3. Backlog 4483: car()'s rho flat on its logit. NOT DONE

The density needed is `logistic(0, 1)` on car's `theta_2`, which is
exactly brms's `(flat)` on rho in (0, 1). Two findings:

- No new seam slot is needed to find the position: car's `natural` slot
  already returns type `"prop"` with a unit gradient row for `theta_2`
  (extensions/frmtmb.spatial/R/car.R, escar and bym2 alike), which says
  "this coordinate is the logit of a quantity on (0, 1)". A rule "a
  covariance position no class reaches, whose natural row is `prop`
  with a unit gradient, takes `logistic(0, 1)` among the sampler's
  defaults and under `prior = "flat"`" covers it. (Core's latent
  residual `cosy` rows are `prop` too, but class `cosy` reaches them.)
- Every line that would apply it names class `"theta"`:
  `default_prior()`'s `add("theta", coef = nm)` (R/priors.R;
  `dev/arspat2-4483-probe.R` shows `(flat) theta theta_2` on route
  `"sample"` for escar and bym2) and frmtmb.sample's
  `default_priors_for()` and `prior = "flat"` stack, beside
  `mix_theta_nlp()`. Lane priors2 owns those lines this round: it renames
  class `"theta"` everywhere and replaces `mix_theta_nlp()`'s fixed
  `dirichlet(1)`. A change here would collide with both, or name a class
  that no longer exists after the merge. Left for after that merge, and
  filed with the rule above.

## 4. Backlog 5256 (= 3994): sar() beside ar(cov = FALSE), and newdata2

Scoped, not built. What brms 2.23.0 takes (`dev/arspat2-brms-probe.R`,
`make_stancode()` only, log `dev/arspat2-out/brms-probe.txt`), gaussian
and student alike:

| second term beside `ar(t, g, cov = FALSE)` | brms 2.23.0 |
|---|---|
| `sar(W, type = "lag")`, `sar(W, type = "error")` | fits: `mu[n] += Err[n, 1:Kar] * ar`, then the SAR density of `Y` |
| `fcor(M)` | fits: the same shift, then `normal_fcor_hom_lpdf` |
| `cosy(t, g)`, `unstr(t, g)` | "Can only model one time-series term." |
| `ma(t, g, cov = TRUE)` | "Can only model one term per autocorrelation class." |
| `ma(t, g, cov = FALSE)` beside `sar(W)` | fits (`Kma`) |

So the second slot holds a registered term only, beside a `cov = FALSE`
ARMA. frmtmb holds one residual term per response in
`frame$autocor[[resp]]`. `autocor_is_cond()` is read at 33 places in 11
files (R/autocor.R 14, R/predict.R 5, R/priors.R 4, R/objective.R 2,
R/predict-brms.R 2, and R/autocor-loglik.R, R/curve-cov.R,
R/families.R, R/frame.R, R/gp-start.R, R/simulate-newdata.R one each),
plus every reader of `frame$autocor` that assigns `theta_idx`, names
`ar[k]` and `lagsar`, reports rows or refuses. The objective is the easy
part: R/objective.R's shift loop runs before the seam density and would
read a second slot. The rest is not: the fitted mean of the lag SAR
with the shift (`A^-1 (eta + shift)`, where the shift reads the observed
response); `simulate()`, where the shift reads earlier draws while the
SAR draw is joint over all rows, so the recursion is implicit;
`predict()` on `newdata`; `log_lik()`; the priors on `ar` beside
`lagsar`; frmtmb.sample's names; the compat rows; the cluster sandwich.
Much larger than a lane item; filed with this scope.

`newdata2`: frmtmb has no `newdata2` argument on any prediction path (in
R/ the word is in messages only). Adding it changes the API of
`fitted()`, `predict()`, `simulate()` and their frmtmb.sample
counterparts, in files other lanes own this round (R/predict*.R). Left
open.

## 5. Backlog 3979 and 4002

- 4002 (covr could not install frmtmb.spatial from inside R): settled,
  not reproduced. `dev/arspat2-4002.R` (log `dev/arspat2-out/4002.txt`,
  covr 3.6.5, R 4.6.1, `R_LIBS` = lane library, rellib-r21, user
  library): `covr::package_coverage(type = "none")` installed it 5 of 5
  times (8.6 to 9.9 s each), and `type = "tests"` ran the whole suite,
  coverage 81.88 percent in 51.4 s. Closed in the backlog.
- 3979 (rr()'s factor-space expansion, gp()/hsgp() as specials and the
  us_t / diag_t sharing not on the published entry shape): still so at
  ec8178d3 (R/covstruct.R branches on `bk[["covstruct"]] == "rr"` and
  `c("rr", "gp")`; `student_shared_lpdf()` in R/objective.R). Not done:
  each port changes the tape of every fit that uses the structure, gp()
  is lane gpfield's this round, and the item is M on its own. Left open.
