# Short check of lane squant, fix round 1 (round of 2026-11-16)

Reviewer: short check, 2026-11-17. Worktree `frmtmb-wt-squant`, lane
library `wt-squant-lib` (installed `elf_cdf_log()` identical to the
worktree source, `squant-short-tape.txt` line 2), extensions from
rellib-r22. Every number is from a log in `dev/reviews/` written by the
script of the same stem (`squant-short-<name>.R` writes
`squant-short-<name>.txt`). No package code edited, nothing committed,
no `git diff`, at most two R processes of mine at once. One of my
processes (`squant-short-compat.R`, stuck in its gp() fit for 19
minutes) was stopped by its Windows PID after matching its command line
on this worktree's script path.

## Verdict: FIX (one item)

The tape fix, the left-censored refits, `elf_quantile()`, the docs and
the test file all hold. One "works" compat row warns and is very slow at
the default smoothing.

1. **(MEDIUM, must fix, compat) `gp_pred()` x `quantile_elf` says
   "works"; at the default `smoothing = 0.25` the fit warns and takes
   200 to 800 times longer.** `ys ~ gp(x) + x2`, quantile 0.5
   (`squant-short-gp.txt`):

   | rows | gaussian() | quantile_elf() (0.25) | quantile_elf(smoothing = 1) |
   |---|---|---|---|
   | 40 | 0.6 s | 129.4 s, warning | 0.6 s |
   | 80 | 0.3 s | 500.6 s | 3.1 s |
   | 150 | 1.2 s | not done after 19 min (stopped) | |

   The warning at 40 rows: "The Laplace approximation is not reliable
   at this fit's optimum for the exact gp() term gp(x): its first-order
   correction is 1.79 in log-likelihood". The profile
   (`squant-short-gpprof.txt`, 97.9 s run) puts 99 % of the time in
   `gp_second_descent()` -> `nlminb` -> `obj$fn` -> TMB's inner
   `newton()` (`iterate`, `updateCholesky`): 21 outer evaluations at
   about 3.7 s each, so the inner problem, not the outer search, is
   slow. Fix: an explicit `r("gp_pred()", "quantile_elf",
   "conditional", ...)` row that says to use `smoothing = 1` (and
   why), and one sentence in "Group-level terms". Recommended, same
   reason: the covariance-structure rows (`us`, `diag`, `homdiag`, `cs`,
   ...) and `s()`, `t2()` reach "works" from the generic
   covstruct-by-family rule, while the family page itself says that at
   a small smoothing a group-level fit "can stop on a local maximum,
   with code 0 and no warning". "conditional" with a pointer to that
   section is the honest status; the lane decides the row text.

Not defects of this lane, for the record:

- (LOW, pre-existing, file it) `cs(1 + x | g)` with the correlation at
  its limit gives "false convergence (8)" and a large-gradient warning,
  where `(1 + x | g)` on the same data gets the boundary message. The
  same under `gaussian()` (gradient warning) and `student()` (both
  warnings) (`squant-short-compat2.txt` (a)). On data with a slope sd
  of 1.5, `us` and `cs` at quantile 0.5 and 0.9 fit at code 0, no
  condition, 20 of 20 (5 seeds, (b)).
- (LOW, optional wording) "Group-level terms" opens "the Laplace
  approximation has a maximum in the group sd"; the earlier review
  noted the supremum can be the boundary sd = 0.

## 1. The tape derivative of `elf_cdf_log()`

`squant-short-tape.R`: tape `d/dmu` and `d/dsigma` of `lF` and `lS` at
mu 0, sigma 1, against `-f/F`, `f/S`, `-q f/F`, `q f/S`, with `f` from
an independently written density and `log F`, `log S` from numerical
integration in log space (geometric pieces, each tail from its own
side). 113 values of `t` in [-2000, 2000], dense at +-690 (steps 1e-6 to
10) and at 0; p 0.01, 0.1, 0.5, 0.9, 0.99; smoothing 0.05, 0.25, 1, 4;
`q` chosen so `q / lambda` is exactly +-690 and 0.

- SUMMARY A: worst value error 9.09e-13; worst relative error of `d/dmu`
  3.14e-12, of `d/dsigma` 3.14e-12 over 40 (p, smoothing, side), where
  the analytic value is above 1e-12 of the slope scale; 0 points above
  1e-8; 0 non-finite. Below that threshold the worst absolute error is
  1.4e-23 of the slope scale.
- The band between 1e-30 and 1e-12 of the slope scale
  (`squant-short-tape2.txt` (a), t in steps of 1 to +-689): worst
  relative error 2.1e-4 (p 0.1 and 0.9, smoothing 1, at t +-668, a
  derivative of 9.6e-31; 5 points above 1e-6 in each of the two settings);
  6.5e-07 at t +-690, smoothing 0.05 (a derivative of 1.45e-17);
  3.2e-09 and 1.8e-11 at t +-690 elsewhere; 1.9e-14 or less in the
  other settings.
- Single-point tapes at exactly t = -690, 0, 690 (`squant-short-tape.txt`
  C): the analytic value to 10 digits wherever it is above 1e-30
  (e.g. p 0.01, smoothing 0.05, t 690: lF mu -0.02346336045 against
  -0.02346336045). The review's -2.09e+74 at t = -690 is gone.

Jumps at t = +-690 and 0 (`squant-short-tape.txt` B, 11 points at
1e-6 spacing, 20 settings x 3 points x 2 sides):

- Values: no jump. The ratio of each value step to `h` times the tape
  slope is 0.999998 to 1.000001 everywhere.
- Derivatives: no jump where the derivative is above about 1e-30 of the
  slope scale. The step across the point equals the steps beside it
  (e.g. p 0.5, smoothing 1, t0 0: 3.2e-07 both; t0 -690 for lF: 1.1e-16
  both).
- The one discontinuity: the log of the LARGE side (lF at t = 690,
  lS at t = -690) where its derivative is below about 1e-30 of the
  slope. There the tape is exactly 0 inside, up to and including
  |t| = 690, and the right value past it (`squant-short-tape2.txt` (b):
  p 0.5, smoothing 0.25, lF: 0 at t <= 690, -8.9e-39 past it, analytic
  -8.9e-39 throughout; also 1.06e-31, 6.75e-77, 2.15e-299). A jump of
  at most about 1e-30 in absolute terms changes no fit. It is not
  `sign_pos()`'s: the zeros start inside the edge (t = 690 - 5e-6),
  where the weights are exactly 0 and the pair is the inside one
  (`log_ibeta_two()`) read at `t`; at exactly |t| = 690 the weights
  also select that pair, as the code comment says.

## 2. The review's left-censored fits (`squant-short-fit.txt`)

`squant-review-tape4.R`'s 16 fits (bounds 95, 90, 80, 60 in data near
100; p 0.1, 0.5; smoothing 0.25, 1), refitted. Finite differences:
central, one Richardson step, absolute step 1e-4 (a step scaled by
|intercept| = 100 spans the kernel width and was wrong in my first
run).

- Tape against finite differences away from the optimum (each
  parameter moved by 0.02 to 0.05): worst relative difference 1.55e-08.
  The review measured 0.0185 against 1.39e-05 at the optimum before.
- At the optimum: max |finite-difference gradient| / |objective| 5.29e-07
  at worst (L 95, p 0.5, smoothing 0.25: gradient 1.73e-04), 1.8e-08 or
  less in 13 of 16. `nlminb` from the fit with the tape gradient moves
  nothing (logLik change 0, parameter change 0, 16 of 16); BFGS on `fn`
  with the finite-difference gradient gains at most 2e-10 in logLik and
  moves the parameters at most 4.2e-06. Hessian smallest eigenvalue 0.459
  to 9.64, all positive. All code 0, no warning.

## 3. `elf_quantile()` (`squant-short-quant.txt`, `-quant2.txt`)

- The review's three cases and their mirrors: 229.34002832,
  68.31361329, 68.30968610 against the review's 229.34003, 68.313613,
  68.309686 (differences within the review's 8-digit rounding); no
  warning.
- Round trip pr -> q -> pr, pr from 1e-300 to 1 - 1e-15 (11 values), p
  0.001 to 0.999 (7), smoothing 0.001 to 4 (6), (mu, sigma) (0, 1) and
  (3, 2): 924 points, `t` up to 6.9e8 in size. Against integration of
  the independent density: worst |log tail(q) - log pr| 5.51e-13, 0
  above 1e-8, 0 non-finite (quant2). Against the package's own CDF:
  2.7e-13. `pbeta()` is -Inf past |t| = 745, so it is no reference
  there.

## 4. Compat rows (`squant-short-compat.txt`, `-compat2.txt`, `-gp.txt`)

`frm_compat("quantile_elf")`, core only: 84 rows (works 38, conditional
32, refused 13, untested 1). Fits on new data (seed 20261117, 400 rows,
25 groups):

| row | model | result |
|---|---|---|
| nl | `a * xp^b`, quantile 0.75 | code 0, no warning (message: starts at the prior locations) |
| mixture | two ELFs, quantiles 0.25 and 0.75 | code 0, no warning |
| us | `(1 + x | g)`, smoothing 1 | code 1, boundary message (correlation at its limit); same under gaussian() and student() |
| diag | `diag(1 + x | g)`, quantile 0.75, smoothing 1 | code 0, no warning |
| cs | `cs(1 + x | g)`, smoothing 1 | code 1, 2 warnings; same under student(), see above |
| s() | `s(x) + x2`, quantile 0.9 | code 0, no warning, 0.7 s |
| t2() | `t2(x, x2)`, quantile 0.25 | code 0, no warning, 2.1 s |
| gp_pred() | `gp(x) + x2`, quantile 0.5 | item 1 |

The us and cs results are the data (the same under gaussian() and
student()), and well-identified data fit cleanly (compat2 (b)).

## 5. Docs (`R/families-quantile.R`, `man/quantile_elf.Rd`, `squant-short-se.txt`)

What the review asked is there:

- "Standard errors" says the default `summary()`, `fixef()`, `confint()`
  and `hypothesis()` are the ELF model's Hessian standard errors and
  cover about 0.7 to 0.8 away from the median. It gives the steps (row
  index, `V <- vcov_cluster(fit, ~ row, "CR0", full = TRUE)`,
  `summary(fit, vcov = V)`, `confint(fit, vcov = V)`) and
  `vcov(fit, cluster = ~ row)`. The example runs all four. Measured:
  Hessian SEs 0.1013, 0.1902 against sandwich 0.1522, 0.3130; the three
  sandwich routes agree to 6 digits.
- "How to choose" step 2 gives the call.
- "Group-level terms" gives the SE per route: Laplace,
  `vcov_cluster(fit, ~ g, "CR0", full = TRUE)` (runs: 0.2251, 0.3539
  against Hessian 0.2320, 0.3943); quadrature, refused by
  `vcov_cluster()` (measured: "vcov_cluster() cannot use a fit made
  with quadrature = TRUE ..."), so only uncalibrated Hessian SEs. The
  sawtooth, the lower-bound counts and 3.42 / 1.76 are in.

## 6. Tests

`tests/testthat/test-families-quantile.R`, NOT_CRAN=true:

- lane build: pass=136 fail=0 err=0 skip=0 warn=0, 16 tests, 2.9 s
  (`squant-short-testrun-lane.txt`).
- rellib-r22: pass=0 err=16 of 16, every one "could not find function
  "quantile_elf"" (or `elf_quantile` NULL) (`-testrun-base.txt`): the
  missing-symbol form only. The behavioral "before" for the derivative
  pin is the lane's pre-fix log (`dev/squant-out/tapegrad-before.txt`)
  and the review's tape5.

## Files of this check

`dev/reviews/squant-short-{tape,tape2,fit,quant,quant2,compat,compat2,
gp,gpprof,se,testrun}.R`, their `.txt` logs
(`squant-short-testrun-lane.txt`, `-base.txt` for the test run;
`squant-short-compat.txt` ends at the stopped gp() fit), and this file.
