# Short check of lane famend's fix round (round of 2026-11-13, base ec8178d3 = 0.82.0)

Verification only: no edit to the lane's source, tests, NEWS or backlog.
Tree `frmtmb-wt-famend` (branch wt-r83-famend), library `wt-famend-lib`,
base `rellib-r21`. Scripts `dev/reviews/famend-short-*.R` and
`famend-short-check.sh`; logs in `dev/reviews/famend-short-out/`
(gitignored). Reference BLAS only: no OpenBLAS run in this check.

## Verdict: FIX

F1 is incomplete. The `NaN` gradient band of item 3067 is still there
at `xi` from about 0.44 to 0.84 (the lane and 0.82.0 alike), and NEWS
and the CLOSED line of the item (test-backlog.md line 3078) say it is
fixed at `xi` from 0.06 to 10. F2 to F6 and both filings SHIP.

## Step 0: library against tree

`famend-short-libcmp.R` sources each changed R/ file (families-extreme,
compat, families-cdf, families, fit, predict, priors) and compares every
function, source references removed, with the installed namespace:
0 differ, 0 missing in all seven files. Control on `rellib-r21`: 13
differ and 7 missing (`frechet_ridge_gain`, `gev_ratio_top`,
`lp_outer_names`, ...), so the compare discriminates. R/families-extreme.R
is dated 14:20 and the install 14:06; the edit after the install
changes no code. Roxygenise on a scratch copy of the tree gives `man/`
and NAMESPACE identical to the tree. No reinstall was needed.

## F1 (HIGH): FIX

What holds:

- Cap back at `exp(700)`: `exp_neg_capped()` is `exp(-ad_max(w, -700))`.
- `dev/famend-gevbit.R` rerun (`gevbit.txt`): inside rows 13257,
  bit-identical value, gradient and Hessian 13257, differ 0; outside,
  6000 rows per `xi` at 0.001, 0.01, 0.06, 0.3, 1, 3, 10, non-finite
  gradient 0 (0.82.0: 16 of 1200 at `xi` >= 0.06). The lane's numbers.
- `famend-review-gevquad.R lane 3`, reference BLAS (`gevquad-ref.txt`):
  code 0, objective 214.7801659, logLik -214.7801659, finite; with the
  "did not settle" and large-gradient warnings, as on 0.82.0.
- `gev_ratio_top` = 34.539575992340424; the ratio at the clamp
  34.539575992340914 (1.42e-14 relative above). At the clamp plus 1, 2,
  3, 4 ulps the ratio is 0.105, 0.201, 0.288, 0.368 below the constant,
  so no inside row reads it (`famend-short-gevclamp.R`).
- Rows whose tape `z` is 1 to 12 ulps above the clamp, `xi` 0.06 to 10
  and -0.5, -0.9, three slots: 288 of 288 bit-identical to 0.82.0
  (`famend-short-gevinside.R`).
- `xi` near 0 (-1e-4, -1e-8, -1e-12, 0, 1e-12, 1e-8, 1e-4), inside rows
  across the support, three slots: 4422 of 4422 bit-identical to 0.82.0
  (`famend-short-gevclamp.R`). 589 of them have non-finite entries, the
  same on 0.82.0: `lccdf` at `x` from 694.5 to 1e12 above the location
  (`xi` <= 0), value `-Inf` (588 rows) and gradient `NaN`
  (`famend-short-xi0.R`). Pre-existing, not this lane's; see "Not
  blocking".

What fails. `gev_parts()` forms `zc = ad_max(z + out, cl)` with
`cl = -1 + 1e-15` and `out = ad_pos(cl - z)`. Where `cl - z` lies in
[8, 16) its rounding grid is 2^-49 = 1.78e-15, and `z + out` lands on
`-1 + 2^-49`, 7.77e-16 (7 ulps) ABOVE the clamp. The second clamp does
not lift it, the ratio there (33.96) is below `gev_ratio_top`, so
`ad_min()` passes it with its live slope of about 5.6e14, and the
adjoint overflows as on 0.82.0. That is every outside row with `z` in
about (-17, -9]: at the band's `x` (-19.8 to -20.6, where `h` lies in
(-700, -673)), `xi` from 9 / 20.6 to 17 / 19.8.

`famend-short-gevscan.R` (`gevscan.txt`), the test's setting (mu 0,
sigma 1, `xi` through log1p), `x` from -19.6 to -20.8 by 0.01, lane
build:

```text
xi 0.42: non-finite gradient   0 of 121 rows
xi 0.44: non-finite gradient  15 of 121 rows, x from -20.46 to -20.60
xi 0.46: non-finite gradient  79 of 121 rows, x from -19.82 to -20.60
...      (every xi from 0.46 to 0.82 by 0.02: 79 of 121)
xi 0.84: non-finite gradient  42 of 121 rows, x from -19.82 to -20.23
xi 0.86: non-finite gradient   0 of 121 rows
total non-finite rows: 1558
xi 0.6, mu 3, sigma exp(0.7): non-finite rows 40 of 41
```

`famend-short-gevband.R` (`gevband-a.txt`, `gevband-b.txt`), `x` -19.4
to -20.3 by 0.001, lane against 0.82.0's `gev_parts()`:

```text
xi 0.5    rows 901 | non-finite gradient: lane 490, 0.82.0 490 | rows with zc above the clamp 901, of them non-finite (lane) 490; at the clamp, non-finite 0
   x -19.811 zc - clamp 7.77e-16 gradient NaN NaN NaN
xi 0.7    rows 901 | non-finite gradient: lane 490, 0.82.0 490 | rows with zc above the clamp 901 ...
xi 0.06, 0.13, 0.2, 1, 2, 5, 10: lane 0, 0.82.0 803; rows with zc above the clamp 0
```

So at `xi` 0.44 to 0.84 the lane changes nothing: the same 490 rows are
`NaN`. The grids of the lane's checks (0.001, 0.01, 0.06, 0.3, 1, 3, 10)
and of the new test (0.06, 1, 10) miss that range. NEWS ("at `xi` from
0.06 to 10. It is now finite") and the CLOSED line of 3067 are false
there. The rows at the clamp are fixed; the rows at the clamp's
rounded neighbor are not.

Fix direction (the lane's call): make every outside row read the
constant. For example, form the outside `zc` so that it is exactly `cl`
(the `z + out` sum is the rounding source; `ad_max(., cl)` alone is
exact), or hold the ratio at a constant below the ratio at
`-1 + 2^-49` (that also holds inside rows with `1 + z` below 1.8e-15,
which must then be shown harmless). Add `xi` 0.5 and 0.7 to the test
and to `famend-gevbit.R`, scan `xi` finely (0.4 to 0.9), and rerun
OpenBLAS on test-families-extreme.R and test-family-fitend.R.

Also LOW, wording: NEWS says "The density, its CDF and their gradients
at every other row are unchanged to the bit." Outside rows away from
the band change: `lpdf` relative change 4e-13 to 8.2e-12 at `xi` 0.3, 1,
3 (`famend-short-gevinside.R`), since `h` reads the constant 64 ulps
below the old ratio. Inside rows are unchanged to the bit. Say "every
row inside the support".

## F2 (MEDIUM): SHIP

`frechet_ridge_gain()` saves with `obj_state_save()` and restores with
`on.exit(obj_state_restore(), add = TRUE)` around both `fn()` calls.
`famend-review-frechet.R` rerun (`frechet.txt`): on 9031 fixed (se
FALSE and TRUE), 9001 fixed, 9031 / 9039 / 9001 with `(1 | g)`, 9031
weighted and 9031 censored, `last.par.best - opt$par` 0 and
`value.best - objective` 0 with the check on, and identical
`last.par.best`, `last.par`, vcov, ranef, logLik and fitted on and off
on every pair (round 1: 3 0 -3 and -1.37e-8).

## F3 (MEDIUM): SHIP

test-family-fitend.R lines 179-209: each spelling either (a) has the
`frmtmb_limit_fit` or `frmtmb_boundary_fit` message and NO warning
(`fe_only(r, character(0))`), or (b) has code 1, a warning whose text
contains "did not report convergence", and only warnings whose text
contains that or "Large maximum absolute gradient" (`fe_only(r, conv)`,
fixed-string match). Both then need logLik agreement below 1e-8
relative. No other ending passes: code 0 without the message fails, a
message plus any warning fails, and any third warning fails. Messages
other than the limit message are not checked in branch (b); harmless.

## F4 (LOW): SHIP

`famend-short-gevmsg.R` (`gevmsg.txt`): the call is read out of the
message and evaluated, on four fits, each code 0 and no warning:

```text
scalar xi -0.63             confint(fit, parm = "xi_(Intercept)", method = "profile")  [-1.3167, -0.9042]
xi ~ x (xi -0.81 to -0.60)  confint(fit, parm = "xi_(Intercept)", method = "profile")  [-1.6333, -1.0384]
xi held at -0.7             confint(fit, parm = "(Intercept)", method = "profile")     [-0.0914, 0.1525]
mu ~ x, xi held at -0.7     confint(fit, parm = "(Intercept)", method = "profile")     [-0.0991, 0.1328]
```

## F5 (LOW): SHIP

NEWS.md (development section) and the `gev_parts()` comment say `xi`
from 0.06 to 10 for 0.82.0's band, which is what the review measured.
The Rd shows the `parm =` call. (That the band is now gone at those
`xi` is F1's defect, not F5's.)

## F6 (LOW): SHIP

`fe_only()` fails on any caught warning that matches no allowed class
and no allowed fixed text. It is used on the Frechet edge fit
(`frmtmb_frechet_nu_edge` only), 9001 (none), both cs() fits and both
nb fits; the GEV fits assert `length(r$warnings) == 0`. The prior test
uses `allow_warnings()`.

## Filings: SHIP

"### From lane famend" in test-backlog.md has the MEDIUM GEV far-tail
tape density item (with the review's reproducer and owner) and the LOW
nb code-1 stop without the limit message (owner: the SE check's limit
rule), each with a severity.

## Step 2: core suite

`dev/release/run-par.sh famend-short dev/famend-core.jobs 3`,
`NOT_CRAN=true`, `FRMTMB_RELLIB=wt-famend-lib`, one file per process,
3 at a time, 06:42 to 06:58 (log `dev/release/famend-short.log` and
`dev/release/famend-short-files/`, deleted after the totals were taken
from them). Every one of the
331 `lib:` lines names `C:/Users/adf44/source/r/wt-famend-lib`.

```text
FAMEND-SHORT ran 331 of 331
pass 24447  fail 0  err 0  skip 240  warn 0
```

No failing file. The lane's files: test-family-fitend.R 45,
test-families-extreme.R 124, test-osa-cdf.R 71, test-cens-trunc.R 100,
test-linpred-layers.R 83, test-prior-table-brms.R 101 pass.

## Step 3: R CMD check --as-cran

`dev/reviews/famend-short-check.sh`: `R CMD build` of the tree and
`R CMD check --as-cran` of `frmtmb_0.82.0.tar.gz` in a scratch folder
outside the repo, `R_LIBS` = wt-famend-lib, rellib-r21, user library,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, `NOT_CRAN=false`, pandoc and
TinyTeX on PATH. 06:59 to 07:35.

```text
Status: 1 NOTE
* checking HTML version of manual ... NOTE: Skipping checking math rendering: package 'V8' unavailable
tests: Running 'testthat.R' [23m] OK  [ FAIL 0 | WARN 0 | SKIP 532 | PASS 19973 ]
```

The NOTE is the environment's (no V8). The check tree is deleted.

## Not blocking

- Pre-existing, the same on 0.82.0: the GEV log survivor at `xi` near 0
  or below, `x` past about 745 scale units above the location, is
  `-Inf` with a `NaN` gradient (and at `x` 694.5 the value is finite
  and the gradient `NaN`), where the exact value is about `-x`. A
  right-censored row there gives an infinite objective. Worth a LOW
  backlog item.
- The nested-quadrature GEV case on reference BLAS: the review had
  0.82.0 at 214.7792046 and OpenBLAS at 214.7801659; the fix-round
  findings have the lane at 214.7792046 under OpenBLAS and 214.7801659
  on reference, and this check reproduces 214.7801659 on reference.
  The values swapped sides, likely through the 1e-12 changes of the
  outside rows on an unsettled fit (not measured). Both are finite.

## Disk

Free on C: at start 261,120,942,080 bytes (243.2 GiB); at end
266,042,413,056 bytes (247.8 GiB). Deleted: the check tree, the suite's
per-file logs. Kept: `dev/reviews/famend-short-out/` (90 KB,
gitignored), the logs this record cites.

## Re-check after fix round 2

The change: in `gev_parts()`, `r = ad_max(z + out, cl)`,
`ex = r - cl`, and `zc = r - (ad_pos(ex) - ad_pos(ex - out))`. All runs
on reference BLAS, wt-famend-lib. Logs are `famend-short-out/r2-*`.
Free on C: at the start of the re-check 266,051,993,600 bytes; at the end
266,048,552,960.

### Verdict: SHIP

F1 is fixed. No outside row has a non-finite gradient on any grid I
ran. Inside rows of the family's three slots are bit-identical to
0.82.0. Fits whose optimizer tries outside points end where they did
on 0.82.0, or 8e-10 lower in one case.

### Library against tree

`famend-short-libcmp.R`: 0 differ, 0 missing in all 7 changed R/
files. The installed `gev_parts()` is the tree's. The tree file is
dated 15 s after the install, and the deparse compare shows no code
change. Of the files the lane changed, only R/families-extreme.R is
newer than round 1.

### (1) Outside rows

- `famend-short-gevscan.R` (`r2-gevscan.txt`): xi 0.40 to 0.90 by
  0.02: 0 non-finite rows in total (round 1: 1558). xi 0.6, mu 3,
  sigma exp(0.7): 0 of 41 (round 1: 40).
- `famend-short-gevband.R` (`r2-gevband-a.txt`, `-b.txt`), 901 rows
  each: lane 0 at xi 0.5, 0.7, 0.44, 0.84, 0.06, 0.13, 0.2, 1, 2, 5,
  10. 0.82.0 at the same xi: 490, 490, 803, 456, 803 (0.06), 803, 803,
  803, 803, 803, 803.
- `famend-short-gevdense.R`, my own scan: outside offsets from 1e-6 to
  60 scale units, 40 random offsets from 1e-10 to 1e3, the band at
  0.01 steps, and 1e4, 1e8, 1e12, 1e16.

```text
PART 1 xi 0.40 to 0.90 by 0.005, 4 (mu, sigma): outside rows 147460, lpdf non-finite gradient 0; NaN Hessian 0 of 21008
PART 2 28 xi from -0.9 to 10, with +-1e-8, +-1e-4, +-1e-3, 0.01, 0.03 (4 (mu, sigma)): outside rows 28748, lpdf non-finite gradient 0; NaN Hessian 0 of 4076
PART 3 9 xi in 0.41 to 0.86, 3 more (mu, sigma) incl. mu 1e4, sigma exp(4): outside rows 9855, lpdf non-finite gradient 0; NaN Hessian 0 of 1404
```

  Part 2 on 0.82.0 (`r2-gevdense-2-base.txt`): lpdf non-finite gradient
  4256, NaN Hessian 792 of 4076, so the scan discriminates. The
  outside lcdf and lccdf gradients are non-finite in 448 of 5112 on
  both builds (part 2 only; parts 1 and 3: 0 of 26664 and 0 of 1782).
  This is the same on 0.82.0, not this change.
- `famend-short-gevclamp.R` (`r2-gevclamp.txt`): of the outside rows
  whose 0.82.0 `zc` was above the clamp with `h` in (-700, -600), 400
  tape-tested, 0 non-finite (round 1: 53). In the band block, 300 of
  300 are finite.

### (2) Inside rows

- `famend-short-gevinside.R` (`r2-gevinside.txt`): 288 of 288 rows 1 to
  12 ulps inside the clamp (the tape's z), three slots, bit-identical
  value, gradient and Hessian.
- `famend-short-gevclamp.R`: xi near 0: 4422 of 4422 bit-identical. The
  589 non-finite rows are the same as on 0.82.0 (the LOW filing).
- The "JUST INSIDE" block (124 of 138, 14 differ) tests `xi * x` on
  doubles. `famend-short-just.R` (`r2-just.txt`) classifies each row by
  z as a tape forms it (`expm1(p[3]) * ((y - p[1]) / exp(p[2]))`):
  inside 120 identical and 0 differ; outside 4 identical and 14 differ.
  R's `expm1()` and the tape's differ by an ulp, which moves those rows
  across the clamp. This agrees with the lane's `famend-fix2-just.R`.
- `famend-short-blend.R` (`r2-blend.txt`): the region where
  `log1p_ratio()` blends (|z| < 0.016), 11 xi x 12 z x 3 slots = 396:
  value, gradient and Hessian identical 396. The same with xi on the
  identity (`famend-short-blend-id.R`): 396 of 396.

### (3) The inside adjoint through ad_pos()

`famend-short-subnorm.R` (`r2-subnorm.txt`) tapes `w1 * l + w2 * h` of
`gev_parts()` alone. The weights include odd subnormals (2^-1074 times
1, 3, 5, 7, 9, 1001; 2^-1022 + 2^-1074; 3e-310), and `zc` runs from 1
to 6 ulps inside the clamp to 40, at 10 xi. That is 1370 tapes, 774 of
them with a subnormal gradient or Hessian entry: 1366 bit-identical to
0.82.0. The value and the gradient are identical in all 1370, so the
first-order adjoint cancels exactly, subnormals included.

The 4 that differ (`famend-short-subnorm2.R`, `r2-subnorm2.txt`) are
the Hessian of `h` ALONE (w = 0, 1) at |z| 0.001 to 0.01, the blend
region: 2 to 4 of 9 entries differ, by 1.7e-16 to 2.2e-15 relative. The
second-order sweep through the two `ad_pos()` terms is not bit-exact
there. No family slot shows it (the 396 + 396 above, the lane's 19287,
my 288 + 4422). NEWS claims bit-identity for values and gradients, not
for Hessians. Not blocking.

### (4) Nested-quadrature GEV case

`famend-review-gevquad.R lane 3` (`r2-gevquad-ref.txt`): code 0,
objective 214.7792046, logLik -214.7792046, finite. The objective is
0.82.0's reference value from the review. The fit still gives the
"did not settle" and large-gradient warnings.

### (5) Fits with outside trial points

`famend-short-gevfits.R` traces `stats::nlminb()` and records every
objective value. I count a trial point as outside when its objective
is finite and above 1e6, which only the barrier gives here. The fits:
xi 0.2, 0.5, 0.7, 0.9, -0.3; 12 seeds; `y ~ 1`, `y ~ x` and
`y ~ 1 + (1 | g)`; n = 200. That is 180 fits on each build.
`famend-short-gevfits-cmp.R` (`r2-fits-cmp.txt`):

```text
fits 180; with an outside trial point (either build) 122
code equal 180 of 180
objective identical 179; lane lower 1; lane higher 0; NA 0
estimates identical to 8 digits 178; max relative change 5.22e-05
 xi0 0.7 seed 1 (1 | g): objective 460.5027821793 -> 460.5027821785 (-8.0e-10), est 5.2e-5 relative, 2 outside trial points
 xi0 0.9 seed 5 (1 | g): objective identical, est 1.3e-8 relative, no outside trial point
```

All 180 have code 0 and the same warning classes (none) on both builds.

### (6) NEWS, the 3067 CLOSED line, the LOW filing

- The NEWS bullet now says the band was about 19.5 to 20.6 scale units
  below the location, at xi 0.06 to 10, and is now finite. It says
  values and gradients are unchanged to the bit at every inside row.
  On outside rows that were above the clamp the barrier is below -1e13
  and changes by up to a factor of 1.2e5; on the other outside rows it
  changes by at most 1e-11. My measurements agree:
  - band: x -19.82 to -20.60 on round 1; 0.82.0's 19.47 to 20.26;
  - other outside rows: 4e-13 to 8.2e-12 relative;
  - the x = -15, xi = 1 row (0.82.0's zc above the clamp): 5.6e3
    relative (`r2-gevinside.txt`), with h near -509, where the barrier
    is about -1e221.
- The CLOSED line of 3067 matches: every outside row on the clamp; 0 of
  6000 at 12 xi; 0 of 75531 on the dense grid; 19287 inside rows x
  slots bit-identical. I did not rerun the lane's 75531-row grid;
  my own dense scan covers that range.
- The LOW filing (GEV log survivor -Inf with a NaN gradient, x 744.8
  to 1e12, 588 rows, and x 694.5, 1 row) matches `famend-short-xi0.R`.

### (7) Tests, reference BLAS

`dev/famend-runtest.R`, one file per process, `lib:` wt-famend-lib:

```text
RESULT test-family-fitend.R pass=49 fail=0 error=0 skip=0 warn=0
RESULT test-families-extreme.R pass=124 fail=0 error=0 skip=0 warn=0
```

The full suite and `R CMD check` of the first part were not rerun; the
lane reports its own check after fix round 2 (`Status: 1 NOTE`, the V8
NOTE, PASS 19977).
