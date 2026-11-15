# Lane selimit2: residues of the standard-error check's limit rule

Round of 2026-11-13, base ec8178d3 (frmtmb 0.82.0), branch
`wt-r83-selimit2`. Brief: `round-briefs/r083-selimit2.md`. This file is
the record; logs named below are local (`dev/selimit2-out/`,
`dev/selimit2-suite-*/`, `dev/release/selimit2-*`, gitignored).

## Items and status

| item | what | status |
|---|---|---|
| 4800 | joint runaway called "flat"; `zi ~ f` silent at SE 5374 | fixed |
| 4814 | zi sent to 1 on a level called the family's limit | fixed |
| 5433 | `jc_nonest()`'s move test is local (new x = -10) | fixed on fitted()'s mean route for rows outside the fitted span; other routes filed; fitted rows filed (decision) |
| 4824 | negbinomial shape "nocurv" behind a noisy gradient | fixed |
| 5423 | asym_laplace NaN fitted values: warning words | fixed (a named separation keeps its words) |
| 2829 | "NaNs produced" from vcov()/summary() under `check_se = "ignore"` | does not reproduce on 0.82.0; closed |
| 2840 | warning side of a correlation past the clamp: SEs 4.7% low | does not reproduce on 0.82.0; closed |
| 3401 | nl body with random effects: SEs 2.05e-4 off | 5.3e-7 on 0.82.0; closed |
| 2888 | correlation verdict depends on where a flat-ridge stop lands | reproduces on 0.82.0 unchanged; open |
| 3370 | joint_cov_repair() inverts Qrr holding beta under REML | met (silent, up to 3.3%); fixed |
| 3020 | fd_hessian() keeps a PD step whatever its rise | already fixed in 0.81.0; closed |
| 2754 | c0k keeps one ridge SE | reference BLAS: 0 of 19; OpenBLAS seed 17 keeps 2 (warned); open |
| 2757 | `s(x) + s(x2)` named on one BLAS only | reproduces (verdict differs by seed and BLAS); open |
| 2759 | unconfirmed downward direction keeps zeroed loadings | ridge verdicts named on all seeds; loadings not inspected; open |
| 3404 | autoscale decides with optimHess() | read, filed (R/autoscale.R is slopescale's) |
| 2834, 5473 | dense us() density digits; Cholesky row held whole | not done: each changes every us() fit's rounding or the SE parametrization (M) |
| 3389, 3575 | asym_laplace: no Hessian SE, code-8 stops | words fixed (5423); refusal or sandwich filed for a user decision |

## What changed (R/se-check.R unless named)

### 4800: joint ends

- `se_limit_joint_dirs()`: per linear predictor, the block of the
  Hessian in natural units (se_limit_units()) over the candidate
  coefficients. A Cholesky factor of the block less `se_edge_screen`
  (0.1) times the identity rules it out with no evaluation. Otherwise the
  candidate is the estimate's own direction inside the eigenvectors at or
  below `se_edge_screen`, weighted by one over each eigenvalue (a Newton
  step inside that subspace). Where single coefficients are already at
  their end, the search runs over the others with those adjusting (the
  Schur complement of their block), so that a level's run through a
  coefficient that also ran alone is found (`zi ~ f` with group a at 1
  and group b at 0, below).
- `se_limit_joint_dir()`: the ray is asked as se_limit_dir() asks one
  coefficient, and two more conditions hold: both sides level at 2
  units, and the loss on the way back to zero more than
  `se_joint_quad` (10) times the quadratic prediction from the
  curvature along the direction and its slope. The first version lacked
  both and called test-se-tier.R's collinear poisson pair (correlation
  of 1 minus 5e-11, SEs 7135, the fit 0.07 short) a limit: 3 failures in
  `dev/release/selimit2-focus1-files/frmtmb--test-se-tier.txt`.
- A joint end is not held coefficient by coefficient: its direction goes
  to tier 3 (attribute `joint` of se_edge_limit()'s result,
  `se_tier3(joint =)`), which removes the eigenvector closest to it
  (cosine above `se_joint_cos`, 0.9) whatever its eigenvalue. So group
  b's mean in `y ~ g` keeps its standard error (test: glm()'s within
  1e-3). A joint end through a coefficient that ran alone loses its rows
  as a single end does.
- se_relabel() then calls the removed coefficients "limit" (kind "end",
  `joint` TRUE in `fit$cache$se_limit`), and the clauses say "as they
  move together further".

### 4814: a probability sent to 1

`se_limit_separating()` takes the ends with their forward direction. A
zi, hu or zoi predictor that the run moves up on some rows (by more than
`se_sep_move`, 1e-3, of the largest move) is a separation of those rows
(`sep1`), with the clause "the probability goes to 1 on some rows, whose
outcomes are all the inflated value"; the other coefficients keep their
own verdicts, as the bernoulli's fb does. Both directions at once keeps
boundfam's rule and its "sepother". Toward 0 alone stays the family's
limit. A single coefficient that separates both ways keeps the joint
search out of the rest of its predictor (they are "sepother" already).

### 5433: the move test away from the data

- `se_limit_ranges()`: per limit direction (each limit parameter, a
  joint end as one, at most `se_range_dirs` = 50), the level range found
  by doubling from 0.5 to 1024 natural units each way while the
  objective changes by at most 1e-6, and a grid of 41 points over it.
  Computed in get_joint_cov() (R/predict.R, both routes, and
  R/joint-sparse.R) through `se_limit_ranges_of()`, once per fit, only
  when a prediction asks.
- `jc_limit_extrap()`: a row is undetermined when the quantity, with the
  predictor moved along the direction anywhere on the grid and the rest
  held, changes by more than `se_limit_move` (0.01) of its standard
  error. Called in predict_mean_se() (fitted()'s expected response) with
  the family's response_mean(). Other routes keep the local test (filed).

### 4824: a limit's gradient

`se_boundary_stop(fit, g, control, gvec)`: with the gradient vector, a
gradient above grad_tol is accepted when every component above it lies
on a linear predictor's coefficient (else the analysis is not forced),
the analysis loses only boundary and limit parameters, and the gradient
over the parameters other than level ends ("end", not rising) is within
grad_tol. check_convergence() (R/fit.R, one hunk) asks it on a code-0
fit whose gradient exceeds grad_tol too, and reads the remaining
gradient in the gradient test.

### 5423: the words

jc_nonest() marks (attribute `nocurv`) rows whose only cause is a lost
parameter of reason "nocurv", read from `jc$lost_why` (get_joint_cov(),
joint_cov_repair_meta(), the sparse route). se_pred_warn() gives such
rows: "read a parameter whose variance the curvature where the fit
stopped does not give (a stop short of the maximum, or a log density
with a kink)". `|` drops attributes, so predict_mean_se() carries it
across its added rows. test-mi-families.R allows both texts.

### Second agent's changes (after the hand-over)

1. `jw` (se_edge_limit()'s cache): se_limit_ranges() read `li$jw`, the
   joint ends' directions, but nothing set it, so a joint end's
   parameters were dropped from the non-local move test altogether.
   The cache now carries the joint ends' forward directions, and
   se_limit_ranges() takes a column only when it runs through a
   parameter asked for.
2. The `el` hunk at the top of cov_from_hessian()'s limit step: same
   logic, laid out so that no line of ours sits between quadseam2's
   changed lines 387 to 389 (its `se_edge_probed()` gate); `el <- NULL`
   moved above `held`, `joint <- attr(el, "joint")` after the block. A
   three-way merge with quadseam2's tree then has an unchanged line
   between every pair of hunks (checked against
   `frmtmb-wt-quadseam2/R/se-check.R`, read-only). Note for the
   consolidator: quadseam2 lets se_edge_limit() run under one-block
   quadrature and swaps the probe's `fn` for se_probe_fn(); the joint
   search uses that `fn` too, but se_limit_ranges() calls `obj$fn`
   directly, which under quadrature is the frozen tape.
3. eam regression (found by the extension suite): frmtmb.eam
   test-defects.R "max_ndt above min(rt) ... mixture" failed (1 fail, 1
   warning): `residuals(fit, type = "response")` gave 232 of 300 NaN
   SEs. The wiener weight's logit is at its end (-23.4), mu1, bs1 and
   ndt1 are "undefined" and level over 1024 units each way
   (`dev/selimit2-eamprobe.R`, `dev/selimit2-out/eamprobe.txt`: ranges
   [-1024, 1024], [-256, 512], [-1024, 1024], theta1 [-1024, 4]); the
   mixture's mean is unbounded in the boundary separation, so the
   non-local test moved every row's mean. The first agent's premise
   ("the fitted rows stay where they are over the level range") holds
   for each row's probability of its outcome, not for its mean.
   Fix: se_limit_span() records, per direction, the range of the
   predictor's move over the fitted rows; jc_limit_extrap() asks only
   a row outside that span (relative tolerance 1e-8). An intercept
   moves every row alike, so no row is outside it; `zi ~ x` at new
   x = -10 is outside and still NaN. After: test-defects.R pass=63
   fail=0 warn=0, test-se-limit.R 362/0, test-mi-families.R 194/0
   (`dev/selimit2-suite-lane-f3/`). No core test reproduces the eam
   case; it is pinned by the extension's own test.
4. OpenBLAS (found by the 0.3.26/0.3.32 runs): test-se-limit.R "a
   fitted value that barely moves along a limit keeps its SE" failed on
   both versions (361/1). The separated bernoulli (seed 514) stops at
   code 1 there with every coefficient "nocurv"; the fit warning names
   the separation, but fitted() took item 5423's new words ("a stop
   short of the maximum, or a log density with a kink"), which
   misdescribe it. Reference BLAS: code 0, lost "separation", old words
   (`dev/selimit2-obsep.R`, `dev/selimit2-out/obsep-*.txt`; rellib-r21
   under 0.3.26 gives the old words). Fix: se_lost_why() reads a
   "nocurv" coefficient that the separation check named
   (`se_sep_named`) as "separation" for the prediction warning
   (R/predict.R, the two `lost_why` lines); `sdr$se_lost` is unchanged.

## Cluster 1 of part B: reruns on 0.82.0

- 2829 (`dev/selimit2-2829.R`, logs `2829-base.txt`, `2829-lane.txt`):
  rs20 seeds 7 and 12 under `check_se = "ignore"`. vcov() and summary()
  give 0 warnings on rellib-r21 and the lane; `vcov(full = TRUE)` gives
  its one designed warning ("Some standard errors are not finite").
  Seed 7: code 1, theta (-0.4263, -2.667, 53.61); seed 12: code 1,
  (-0.8143, -3.413, -190.9). No "NaNs produced". Closed.
- 2840 (`dev/selimit2-2840.R base 0 20`, `2840-base.txt`; lane
  `2840-lane.txt`): test-review-v29.R's `(1 | q | id)` fit on 21
  ulp-moved data sets (k = 0 to 20). 17 stop with a convergence warning
  (the warning side), yet every fit loses only theta_3 ("boundary") and
  reads the log-sd SEs 0.21767 or 0.21768 and 0.22244 or 0.22245,
  the rank-one model's 0.217708 and 0.222476 to 1.3e-4. The k = 0
  warning is the code-8 text that says the stop is at the maximum
  (Newton gain 2.6e-5). Base and lane lines are identical (`diff`).
  Closed: the 4.7 percent no longer reproduces.
- 2888 (`dev/selimit2-2888.R <lib> 4`, a copy of
  dev/relrev071-rssd.R on rellib-r21; `2888-base.txt`, `2888-lane.txt`):
  seed 4 rs_sd0.02 stops at t = 19.66 (code 1, singular convergence),
  logLik 5.68e-8 below lmer, the correlation parameter's SE 1.8639e+04,
  nothing lost, the code-7 "at the maximum" warning. Identical on both
  builds. lmer's rho is 0.99926 (t about 26), an interior maximum on a
  flat ridge, so neither verdict is wrong; the dependence on where the
  stop lands remains. Open; it needs the correlation-space pull of 2834.
- 3401 (`dev/selimit2-3401.R base`, `3401-base.txt`): lane adhess's
  nl_re (seed 11) against a Richardson Jacobian of the AD gradient
  (refself 5.5e-11): the fit's SEs are within 5.345e-7, the value
  adhess measured for its "curv" candidate, which fd_hessian()'s retake
  has applied unconditionally since 0.81.0. Closed.

## Tests and runs so far

- New tests in `tests/testthat/test-se-limit.R` (five `test_that()`
  blocks, items 4800, 4814/4800, 5433, 4824, 5423). On rellib-r21:
  pass=272 fail=87 err=1 (`dev/release/selimit2-sl-base-files/`), every
  new block failing behaviorally. Lane: pass=362 fail=0 err=0 warn=0.
- test-se-tier.R 112/0, test-mi-families.R 194/0 (its allowed warning
  text widened), 18 more SE and prediction files 0 fail 0 warn
  (`dev/release/selimit2-focus1.log`).
- Traced core suite on rellib-r21 (`dev/selimit2-suite-base/`, 330
  files): pass=24420 fail=87 err=1 skip=240 warn=0; the only failing
  file is test-se-limit.R (the new tests).

## SE comparison against rellib-r21 (core, traced)

Both arms one file per process, 330 files, through
`dev/selimit2-run-traced.R` (`dev/selimit2-par.sh base|lane
dev/selimit2-core.jobs 3`; every lane log's `lib:` line is
wt-selimit2-lib). The lane arm ran the build before the second agent's
changes (see "Second agent's changes" and item 3370). Those touch only
prediction standard errors (joint ends, rows inside the fitted span, a
separation's words, REML fits on a whitened design with a lost
parameter) and the place of two lines; none enters cov_from_hessian(),
so the COV records above are those of the final build. The final
build's runs are under "Final runs".

- Lane suite: pass=24510 fail=0 err=0 skip=240 warn=0 (base:
  pass=24420 fail=87 err=1, all in test-se-limit.R's new blocks).
- `Rscript dev/boundfam-cmp.R dev/selimit2-suite-base
  dev/selimit2-suite-lane` (`dev/selimit2-out/cmp-core.txt`):

      files compared 330
      COV bitwise=5724 ulps=1 moved=24 nan=10 length=18 count=4
      VCT bitwise=160 ulps=0 moved=0 nan=0 length=0 count=0

  Every non-bitwise record explained:
  - "count" 4 files, and with them every "length" and all "moved"
    records but one, are record misalignment, not changed SEs. A diff
    of the COV/LIM sequences with the record numbers stripped shows
    only inserted or removed records: test-ce-bands.R +3 and
    test-predfix.R +3 (the analysis that se_boundary_stop() now forces
    on a fit whose gradient is above grad_tol on linear-predictor
    coefficients, item 4824), test-se-reparam-inv.R +5 (the same, on
    degree-7 designs; nothing is lost, so the stop is not accepted),
    test-se-flags.R -1 (base computed record 28 twice with identical
    values, the lane once). Every other record of those files is
    identical, so no shown SE changed.
  - "moved" test-sugar.R record 2: `y ~ x` on `rnorm(50)` with no seed
    (line 40), different data each run, as boundfam found.
  - "nan" 10: 6 in the misaligned files (alignment only) and 4 real:
    test-se-limit.R records 83, 84, 87 (the new 4800/4814 tests:
    `zi ~ f`/`hu ~ f` SEs 5374, 5466, 5298 to NaN, the intended change)
    and test-ordinal-mixture.R records 28, 29 (below).
- `Rscript dev/boundfam-classify.R dev/selimit2-suite-lane
  dev/selimit2-suite-base` (`dev/selimit2-out/classify-core.txt`):

      TOTAL fired=166 already_nan=159 newly_nan=5 no_base=2

  Base fired 139 times, the lane 166. Per file the firing counts are
  equal except test-se-limit.R (43 to 68, the new tests) and
  test-ordinal-mixture.R (5 to 7). The 2 "no_base" are test-predfix.R
  LIM 32 and 33, the same firings as base LIM 29 and 31 (record numbers
  moved by 3; the LIM lines are equal but for the new fields). The 5
  "newly_nan":
  - test-ordinal-mixture.R #28 and #29: the mixing logit `theta1 ~ z`
    of `mixture(cumulative, sratio)` and `mixture(cumulative,
    cumulative)` ran off jointly to (15423, 14087) and (4015, 3745),
    ray gain 4.8e-9 and 1.6e-8, moving rows both ways (sep=1): the
    mixing weight separates. Base kept SEs of 1.44e6/1.31e6 and
    3.97e5/3.70e5 with no word; the other 3 SEs are kept at ratio 1.
    The file passes.
  - test-se-limit.R #92, #93, #96: the new `zi ~ f`/`hu ~ f` tests.
  Of the lane firings carrying the new `joint` or `sep1` flags, all
  but the 2 ordinal mixtures are test-se-limit.R's new tests (poisson
  `y ~ g`, the 20 ulp-moved `y ~ x` fits, the zi/hu factor fits). No
  firing outside those two files changed kind, estimate or SE.
- False positives on finite MLEs: none found. Every added firing has
  an estimate beyond 3700 units on its link or is a constructed limit
  of the new tests, with a ray gain below 2e-8. test-se-tier.R's
  collinear poisson pair (the joint rule's known trap, standard errors
  7135) passes 112/0.

## Cluster items taken by the second agent

- **3370, fixed.** Built where it acts (`dev/selimit2-3370.R <arm>
  <deg> <seeds>`, 256-bit reference with Rmpfr): gaussian REML,
  `y ~ poly + (1 | g)` raw, 200 rows, no group effect, so the sd runs to
  0 and is lost "boundary", `fit$beta_whiten` set. vcov() is within
  4e-11 (degree 6) and 1.1e-9 (degree 7) of the exact conditional
  covariance on both builds. fitted(re_formula = NA) at x = 1 to 2,
  seeds 1 to 5, max |SE / exact - 1|:
  - rellib-r21: degree 6 9.8e-5, 1.1e-4, 9.5e-5, 5.0e-5, 7.6e-5;
    degree 7 1.4e-2, 2.9e-2, 9.0e-3, 1.2e-2, 3.2e-2; no warning
    (`3370-base.txt`, `3370-base-d7.txt`). Qrr's condition number
    1e15 to 1.9e16, two with a negative eigenvalue in double.
  - inverting Qrr in gamma only: degree 6 3.7e-6 to 1.7e-5, degree 7
    2.6e-3 to 3.8e-3 (`3370-lane2*.txt`): the dense natural-coordinate
    `a' Wrr a` still cancels.
  - final (gamma inverse kept as a factor, `root`, Vrest NULL): degree
    6 3.2e-11 to 1.7e-10, degree 7 7.2e-10 to 1.8e-9
    (`3370-lane3.txt`, `3370-lane3-d7.txt`).
  Change: R/se-reparam.R sdr_unwhiten() keeps the precision it was
  handed (`sdr$joint_precision_white`, one line); R/predict.R
  joint_cov_white_inner() and a branch in joint_cov_repair(). New test
  in test-se-reparam-reml.R (degree 7, 20 ulp-moved responses):
  rellib-r21 pass=40 fail=20 (every draw, 0.68 to 3.3 percent;
  `dev/selimit2-suite-base-reml/`), lane 60/0, OpenBLAS 0.3.26 and
  0.3.32 60/0.
- **3020, closed.** fd_hessian() (R/se-check.R, unchanged here) retakes
  every column whose rise exceeds `tau_band` and keeps the retake
  whenever it is finite, singular or not; the condition "only where the
  first matrix did not factor" was removed in 0.81.0 (its comment, and
  3401's measured 5.3e-7 above).
- **2754, open.** `dev/selimit2-c0k.R <arm> 1:20`. Reference BLAS,
  rellib-r21 and lane alike: 19 fits with a flat ridge, 0 keep a
  finite ridge SE (seeds 1 and 7 lose all three as "flat"); seed 20
  stops 3.1 short with code 1 and warnings. OpenBLAS 0.3.26, lane:
  seed 17 (code 1, two warnings) keeps a_(Intercept) 48 and
  k_(Intercept) 0.017 with theta_1 NaN; seeds 1 and 7 are named
  (`c0k-*.txt`). The reference part no longer reproduces; the
  OpenBLAS seed does, under a convergence warning.
- **2757, open.** `dev/selimit2-smooth2.R`, lane: reference BLAS names
  the two log sds on seed 1 and keeps SEs 604 to 1040 on seeds 2 to 5;
  OpenBLAS 0.3.26 names seeds 1, 3, 4 and keeps 898 and 611 on 2 and 5.
  Every fit warns (the fx1 columns are lost). The trade of the two sds
  is a curved ridge (s1^2 + s2^2 fixed): ±0.5 in log raises the
  objective by 0.28 to 0.30, second order. Not fixed.
- **2759, open.** On the same c0k runs every ridge parameter is named
  on seeds 3 and 7 (both BLAS); the shown covariance's loadings were
  not inspected, so whether zeroed loadings remain is not measured.
- **3404, filed (not edited).** R/autoscale.R is lane slopescale's.
  Reading: autoscale_prefit_converged() and autoscale_fit_sound() take
  optimHess() (step 1e-3 in natural units). obj$he() exists only
  without random effects; a fit with them would need
  fit_outer_hessian()'s steps. The cutoffs (reciprocal condition 1e-4
  on the coefficient block, 1e-8 whole) were measured on optimHess()
  matrices (dev/predfix-p2-prefit.R) and would need remeasuring. Filed
  as a LOW item in the lane's backlog section.
- **3389, 3575, filed for a decision.** `dev/selimit2-asym.R lane`
  (`asym-lane.txt`): dbg14 (seed 43, quantile 0.3), gaussian data
  (seed 5), a fixed quantile 0.5, zero_inflated_asym_laplace (seed 9).
  Every fit loses every parameter as "nocurv"; codes 8, 8, 8, 0;
  largest gradients 1.38, 1.39, 0.127, 0.511; each warns (false
  convergence and/or the gradient); fitted() gives every SE NaN with
  item 5423's words. Nothing is silent. A refusal by name or a
  quantile-regression covariance is a user decision (filed MEDIUM).
  Recommendation: refuse by name at fit time, saying the log density is
  piecewise linear in the location, and point to frmtmb.sample for
  uncertainty; that is what brms users get (posterior SDs).
- **2834, 5473, not done.** 2834 evaluates every us() block from its
  Cholesky factor, which changes the rounding of every us() fit (the
  item itself says a minor); 5473 needs a Cholesky row in polar
  coordinates for the SE. Each is M and outside this lane's remaining
  time; 2888 waits on 2834.

## Final runs (the final build, wt-selimit2-lib)

- Focus files after the hand-over edits (`dev/selimit2-suite-lane-focus2/`,
  21 files): all fail=0 err=0 warn=0, among them test-se-limit.R 362,
  test-se-tier.R 112, test-mi-families.R 194, test-se-check.R 134,
  test-numerical-robustness.R 723, test-predfix.R 92.
- After the span and separation-words fixes
  (`dev/selimit2-suite-lane-f3/`): test-se-limit.R 362/0,
  test-mi-families.R 194/0, frmtmb.eam test-defects.R 63/0.
- OpenBLAS with LAPACK (`dev/selimit2-ob.sh`, OPENBLAS_NUM_THREADS=4,
  `dev/selimit2-out/ob-*`; builds checked: dgemm 2000 in 0.48 s
  single-threaded against 11.4 s for the reference): 0.3.26 and 0.3.32
  each test-se-limit.R 362/0, test-se-tier.R 112/0,
  test-mi-families.R 194/0, test-se-reparam-reml.R 60/0, frmtmb.eam
  test-defects.R 63/0. Before the separation-words fix both versions
  read test-se-limit.R 361/1 (above).
- Gated core files and every extension file, one per process,
  `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`
  (`dev/selimit2-suite-lane-final/`, 261 files, every `lib:` line the
  lane's core; skips are only the scale and drmTMB tiers):

      frmtmb (62 gated files) pass=6381 fail=0 err=0 skip=13 warn=0
      coupling  11 pass=551  fail=0 err=0 skip=5 warn=0
      eam       33 pass=1857 fail=0 err=0 skip=3 warn=0
      latent    13 pass=697  fail=0 err=0 skip=2 warn=0
      learn     16 pass=545  fail=0 err=0 skip=2 warn=0
      ode       11 pass=549  fail=0 err=0 skip=1 warn=0
      sample    91 pass=4603 fail=0 err=0 skip=1 warn=0
      spatial    8 pass=612  fail=0 err=0 skip=0 warn=0
      survival  16 pass=575  fail=0 err=0 skip=0 warn=0

- Earlier ungated extension run on the pre-span build
  (`dev/selimit2-suite-lane-ext/`): all fail=0 warn=0 except eam
  test-defects.R 62/1 warn 1 (the regression fixed above).
- Core in one session, three orders (`NOT_CRAN=true
  FRMTMB_RELLIB=wt-selimit2-lib bash dev/release/run-sessions.sh
  selimit2-core 3 frmtmb`, `dev/release/selimit2-core.log`):

      SESSION frmtmb forward files=317/330 pass=24551 fail=0 err=0 skip=240 warn=0
      SESSION frmtmb reverse files=317/330 pass=24716 fail=0 err=0 skip=240 warn=0
      SESSION frmtmb preload files=317/330 pass=24819 fail=0 err=0 skip=242 warn=0

  The 13 files without rows are the gated brms-suite and `-brms` files,
  which run (and pass) in the gated run above. Every LEAK line is
  tagged "[with load of ...]" (test-arg-refusal.R, test-cox-home.R,
  test-families-extreme.R loading an extension); none untagged.
- R CMD check --as-cran, core (built and checked in
  `dev/selimit2-check/`, `R_LIBS` lane lib, rellib-r21, user library,
  `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH;
  logs kept as `dev/selimit2-out/check-run*`, the tree deleted):
  - Run 1: `Status: 1 ERROR, 1 NOTE`. The tests stopped after 17 min
    with no testthat summary and no `proc.time()` line, so R ended
    abnormally (no Application Error event logged); before that,
    test-gp-noncentered.R lines 481, 527, 528 saved `_problems`.
    Those three pass one file per process (lane twice, rellib-r21,
    and under `OMP_THREAD_LIMIT=2 _R_CHECK_LIMIT_CORES_=true`) and in
    all three sessions. `dev/selimit2-gpnc.R`, 9 processes (lane,
    rellib-r21, the check's installed copy), reads bitwise the same
    values far inside the bounds: 481 1.989e-9 (bound 1e-3), 527
    1.472e-12 (1e-8), 528 5.046e-10 (1e-5)
    (`dev/selimit2-out/gpnc-repeat.txt`). Lines 527 and 528 call no
    code this lane changed (an RTMB tape and gp_tape_scope()). The
    first build attempt also failed, because the OpenBLAS tree was
    being deleted while R CMD build copied it.
  - Run 2, same tarball: `Status: 2 NOTEs` (examples over 5 s elapsed:
    frm_hazard_reads 6.59 s, dharma_residuals 6.03 s, CPU 3.7 and
    1.6 s, under 14 lanes' load; the V8 math-rendering NOTE). Tests
    OK in 78 min: `[ FAIL 0 | WARN 0 | SKIP 273 | PASS 24295 ]`.
  Run 1 is unexplained; it needs the user's eye only if it recurs.

## Files touched (both agents)

R/se-check.R, R/fit.R (one hunk, check_convergence()), R/predict.R,
R/joint-sparse.R, R/se-reparam.R (one line, sdr_unwhiten()), NEWS.md,
dev/test-backlog.md, tests/testthat/test-se-limit.R,
tests/testthat/test-mi-families.R, tests/testthat/test-se-reparam-reml.R;
new dev scripts `dev/selimit2-*`. No extension, no man/ change
(roxygenise wrote nothing: every new function is `@noRd`). R/predict.R
was CRLF in the checkout (core.autocrlf) and is LF after a sed edit;
git normalizes line ends on commit, so only content differs.

## Bump

Minor (BREAKING verdict changes: joint ends, one-sided separations,
NaN extrapolated fitted() SEs). No new floor.

## Not done, and what needs the user

- Decisions filed MEDIUM in the backlog: (1) whether a limit's
  undetermined parameters make the fitted rows' means NaN (the level
  set's answer) or keep them (the limit model's, today's behavior for
  rows inside the fitted span); (2) asym_laplace(): refuse the Hessian
  SE by name, or give a quantile-regression covariance.
- Open, with the numbers above: 2754 (OpenBLAS seed 17), 2757, 2759
  (loadings not inspected), 2888, 2834, 5473; 3404 filed (slopescale's
  file).
- The non-local move test covers fitted()'s mean only (filed LOW by the
  first agent).
- No brms-suite ledger row edited; every gated brms-suite file passes.

## Fix round (review `dev/reviews/2026-11-13-selimit2.md`, verdict FIX)

Before = the build the review read (installed 10:48); after = this
build. Logs `dev/selimit2-out/fix/`.

### F1 (HIGH): a joint end hid a stop short

se_boundary_stop() dropped the whole gradient of every "end"
coefficient. Now se_boundary_stop_rest() removes only the projection of
the gradient (optimizer units) on the span of the level ends'
directions: each single end not rising (unit vector) and each joint
end's `jw` column (divided by `par_units`), by a QR; the largest
remaining component is judged against grad_tol.
`dev/reviews/selimit2-review-grad2.R lane` (poisson `y ~ g`, group a
all zero, seed 1, start (-25, 24.5)):
- before: iter 1 code 1, gradient 94.2 94.2, group b's mean 85 percent
  off, 1 warning (the limit words only) (`grad2-lane.txt` of the
  review); base 2 warnings.
- after (`fix/grad2-lane-final.txt`): iter 1 the 2 base warnings
  ("would still gain 21.6", "Large maximum absolute gradient (94.2")
  and no limit words; iter 2 (gradient 2.57) the same 2; iter 3
  (gradient 0.00614, gain 1.7e-7, within grad_tol) the code-1 note and
  the limit words.
- New test "a stop short across a joint end keeps its convergence
  warnings" (test-se-limit.R): passes on base and after; fails on the
  reviewed build (its log above).

### F2 (MEDIUM-HIGH): joint verdict on a finite MLE

New se_limit_joint_cert() (R/se-check.R): a joint direction is kept
only when the rows it moves prove it. Per row, a move its outcome
rules out is `bad`: mean of poisson, negbinomial(2), geometric,
binomial, zero-inflated count: any row up, or down with y > 0;
bernoulli: down with y != 0 or up with y != 1; hu: up with y != 0 or
down with y = 0; zoi: up off {0, 1}, down on {0, 1}; coi: among y in
{0, 1}, up with y = 0, down with y = 1. zi: up with y != 0 is bad, a
zero moved down is `soft` (a bounded loss). A row counts as moved above
1e-8 of the largest move. Rows moved by at most `se_sep_move` (1e-3) of
the largest, the wrong way, are rounding of the run: the direction is
first projected (QR) to hold them; if nothing is left, no end. Any bad
row left: no end ("flat" as base). Soft rows, and every mixing weight
`theta`, are checked far out instead: the objective with every row
moved by more than 1e-3 of the largest put 40 units past its own
predictor must not rise by more than 1e-6 (each row's log density is
monotone and bounded along the ray). Other dpars are not certified.
The certified (projected) direction is the one stored and used.

`dev/reviews/selimit2-review-fp.R lane pois zi`
(`fix/fp-lane-final.txt`):
- delta 1e-3: no word, finite SEs (base the same); delta 1e-4, 1e-5:
  "flat" (base the same). Before: all three "the maximum likelihood
  estimate does not exist".
- one count at the top: still the joint limit; `zi ~ f` / `zi ~ 0 + f`
  / `hu ~ f` with level a all zero: still limits (sep1); one or two
  nonzero counts in level a: no sep1; `hu ~ f` with level b free of
  zeros: hu_fb a single end (quiet message, filed LOW).
- The 4800 test's new absent case runs these three deltas.
- First version of the projection held every row moved by less than
  1e-3: on test-ordinal-mixture.R it held the rows near the mixing
  boundary and lost #28/#29 (5 LIM lines instead of 7). Holding only
  rows moved the wrong way restored them (7 LIM lines, the same two
  joint ends, `dev/selimit2-suite-lane-f4/`).

### F3 (MEDIUM, silent wrong): the span rule

se_limit_span() is removed; jc_limit_extrap() tests every row. In
se_limit_ranges() an undetermined parameter's direction now carries
`base`, the ends moved to the far end of their level ranges toward
their limits (single ends with their sign, joint ends by `jw`; halved
until the combination is level), and its own level range is measured
from there. jc_limit_extrap() takes every linear predictor the mean
reads (`preds`) and moves them together, so an end in one dpar and an
undetermined parameter in another (the eam mixture) are moved jointly.
R/predict.R predict_mean_se() builds `preds` once.
- `dev/reviews/selimit2-review-span.R lane`
  (`fix/span-lane-final.txt`): rows (1, 1), (1, 0.8), (1.5, 0.9) NaN
  with the warning (before and base: 0.4857, 0.8433, 0.2625); (0, 0.5)
  and (3, 0) keep 0.3270 and 0.1434.
- frmtmb.eam test-defects.R: 63/0 (the fitted rows finite: at weight 0
  the wiener parameters do not enter the mean).
- test-se-limit.R 5433 block: x = -10 and 10 NaN, x = -1, 0, 3 and the
  fitted rows finite, unchanged.
- New test "a row inside the data that the limit moves has no SE".
- NEWS bullet rewritten. The lane's question (a) is settled this way.

### F4 (MEDIUM): 3370's cost

joint_cov_white_inner(): a sparse Cholesky of `Qg[r, r]`
(Matrix::Cholesky, perm), `Lg = P' L^-T`, `Wg` and `Jg` by solves, `T`
applied to the beta rows (and columns of `Wg`) only, returns NULL above
`joint_root_max`; joint_cov_repair() no longer rebuilds `V` from the
root. `dev/reviews/selimit2-review-3370.R lane cost 3000`
(`fix/3370-cost-lane-3000.txt`): fitted() 1.15 s, max used 598 MB
(base 0.73 s, 459 MB: 1.6 and 1.3 times; before 70.53 s, 701 MB); SEs
0.030084767 0.008550877 0.031072024 (reference 0.030084767
0.008550877 0.031072023). test-se-reparam-reml.R 60/0.

### LOW

- 4824 NEWS bullet marked BREAKING, and says the gradient across the
  end is still tested.
- hu ~ f wording asymmetry filed LOW (my backlog section); the MEDIUM
  item on fitted rows is replaced (settled by F3).
- se_limit_ranges() and se_probe_fn(): se_probe_fn() is in quadseam2's
  tree, not this one; filed LOW for the merge.
- The 5423 block's skip_if_not() is replaced: it asserts some "nocurv"
  is lost and accepts the general words where another reason is lost
  too.
- test-bracket-access.R caught `r$family$family` in the new
  certificate (a hazard read); now `[["family"]][["family"]]`.
- Memory notice: the first fix-round core run was stopped (my driver
  and one test) and resumed at 2 processes; it was then rerun whole
  after the F2 change.

### Fix-round runs (final build)

- The six files, reference BLAS (`dev/selimit2-suite-lane-f4/`,
  `-fix3/`): test-se-limit.R 378/0, test-se-tier.R 112/0,
  test-mi-families.R 194/0, test-se-reparam-reml.R 60/0,
  test-ordinal-mixture.R 99/0, frmtmb.eam test-defects.R 63/0, warn 0.
- OpenBLAS 0.3.32 with LAPACK, `OPENBLAS_CORETYPE` Haswell and
  Sandybridge (every log's load line reads "Core: Haswell" /
  "Core: Sandybridge", 6 of 6), `dev/selimit2-out/ob-0.3.32-*`: the same
  six files, the same counts, fail 0 warn 0 on both kernels.
  `dev/ciharden-out/` deleted.
- test-se-limit.R on rellib-r21 (`dev/selimit2-suite-base-sl2/`):
  pass=286 fail=89 err=1; the F3 block fails (2), the F1 block and the
  F2 absent case pass on base, as they should (both fail on the
  reviewed build, per the review's logs).
- Full core, traced, one file per process, 2 at a time
  (`dev/selimit2-suite-lane-fin/`, every `lib:` line the lane's):
  pass=24567 fail=0 err=0 skip=240 warn=0.
- `Rscript dev/boundfam-cmp.R dev/selimit2-suite-base
  dev/selimit2-suite-lane-fin` (`dev/selimit2-out/cmp-fin.txt`):

      COV bitwise=5704 ulps=7 moved=24 nan=18 length=24 count=6
      VCT bitwise=159 ulps=1 moved=0 nan=0 length=0 count=0

  New against the first comparison: test-se-limit.R and
  test-se-reparam-reml.R gain records (the new tests: count, and the
  nan of the shifted records); ulps 7 and 1 are cross-process
  rounding; test-id-kron.R record 10 reads theta_2 NaN against base
  0.9613 with no firing in either arm (0 LIM lines both): the
  hybrid-CPU flake boundfam recorded for that very record (three base
  runs gave NaN twice and 0.9613 once), not this lane's code.
- `Rscript dev/boundfam-classify.R dev/selimit2-suite-lane-fin
  dev/selimit2-suite-base` (`dev/selimit2-out/classify-fin.txt`):

      TOTAL fired=168 already_nan=155 newly_nan=5 no_base=8

  The 166 firings hold (every file's LIM count equal to the first lane
  run but test-se-limit.R, 68 to 70, the new tests); the 2 ordinal
  mixtures' joint ends hold with the same estimates and the same NaN
  (#28, #29); the 3 other newly NaN are test-se-limit.R's zi/hu factor
  tests (record numbers moved by the new tests, so their base pairing
  is shifted; the kept-SE ratios 0.95 and 1.26 are that pairing). The
  8 no_base are the same shift (test-predfix.R 2, test-se-limit.R 6).
- R CMD check --as-cran, core, the final build (`dev/selimit2-check/`,
  deleted; logs `dev/selimit2-out/check-fix-*`): `Status: 1 NOTE` (the
  V8 math-rendering NOTE). Tests OK in 29 min. Run 1's abort did not
  recur.
  `[ FAIL 0 | WARN 0 | SKIP 273 | PASS 24311 ]`.
- Question (a) of "Not done" is settled by F3; question (b)
  (asym_laplace) waits for the user, unchanged this round. Files
  touched in the fix round: R/se-check.R, R/predict.R, NEWS.md,
  dev/test-backlog.md (my section only), tests/testthat/test-se-limit.R;
  new dev scripts dev/selimit2-certprobe.R, dev/selimit2-lpprobe.R.

## Fix round 2 (short review `dev/reviews/2026-11-13-selimit2-short.md`)

### F4b (silent, pre-existing): whitened REML above joint_root_max

- Before (rellib-r21 and the fix-round build alike,
  `dev/reviews/selimit2-short-out/3370ref-*-4500-*`): 4510 joint
  parameters, sd lost, fitted(re_formula = NA) SEs 0, 0.006780297, 0
  (seed 2; conditional 0.026547773 0.006987372 0.024168042) and 7
  percent high (seed 3), no warning.
- Fix: joint_cov_white_inner() returns the beta positions as
  `too_big` above `joint_root_limit()` (new,
  `options(frmtmb.joint_root_max)`, default the constant 4000);
  joint_cov_repair() then sets those rows and columns of V to NaN,
  drops the root, and adds them to the lost positions with reason
  "size". jc_nan_lost() and jc_nonest() carry an attribute `size`, and
  se_pred_warn() words those rows: "read coefficients whose covariance
  this REML fit cannot give accurately above 4000 joint parameters, so
  their standard errors are NaN; raise options(frmtmb.joint_root_max =
  ), at memory that grows with its square, or draw from the fit with
  frmtmb.sample::frm_sample()". predict_mean_se() carries the
  attribute across its added rows.
- After, the same case (`dev/selimit2-out/fix/3370ref-lane-4500-s2.txt`):
  fitted SE NaN NaN NaN with that warning.
- New test in test-se-reparam-reml.R, the limit lowered to 20 with
  `withr::local_options()` (a 30-parameter fit): fails on rellib-r21
  (pass=43 fail=22: 2 here, 20 the 3370 test) and on the fix-round
  build (pass=63 fail=2), `dev/selimit2-suite-*-f4b0/`; passes after.
- The sparse `a' Qrr^-1 a` route is not built: filed HIGH with the
  numbers and the exact fix. sdr_unwhiten() still uses the constant
  (vcov()'s route; the whole-route seed 1, 4.5 percent low, is in the
  same item).
- NEWS: the 3370 bullet is bounded by size, with the NaN and warning
  above it marked BREAKING.

### Filed

- MEDIUM: a mixing weight's end is not certified (the short review's
  `selimit2-short-theta.R`, 0.082 at x1e4, 43 at x1e6, a finite MLE
  called separated on both builds).

### LOW

- The 5423 block's fit now runs under `allow_warnings()` with the stop's
  own texts in place of `suppressWarnings()`.

### Runs (reference BLAS, final build)

- test-se-reparam-reml.R 65/0, test-se-limit.R 378/0, test-se-tier.R
  112/0, warn 0 (`dev/selimit2-suite-lane-f5/`).
- R/se-check.R changed besides R/predict.R and NEWS, so R CMD check
  core was rerun: `Status: 1 NOTE` (V8), tests OK in 26 min,
  `[ FAIL 0 | WARN 0 | SKIP 273 | PASS 24316 ]`; logs
  `dev/selimit2-out/check-fix2-*`, the check tree deleted.
