# Lane ciharden: CI through the Ubuntu 26.04 move

The goal: CI stays green when `ubuntu-latest` becomes Ubuntu 26.04,
and the test suites stop depending on platform rounding.

Base: 4f5ea39f (frmtmb 0.68.1), built in
`C:/Users/adf44/source/r/rellib-r6`. Worktree `frmtmb-wt-ciharden`,
private library `wt-ciharden-lib`. Logs are in `dev/ciharden-log/`
(gitignored); every number below is copied from them, with the script
that made it.

## Punch round 1 (review of 2026-10-07, `dev/reviews/2026-10-07-ciharden.md`)

This section supersedes what Items 3 and 4 below say where they
disagree. Logs of this round are in `dev/ciharden-log/p1/` and
`dev/ciharden-log/p1-*`.

**B1, test-perf.R.** The review is right on both counts. The tape's
operations per observation depend on the response (CppAD records
nothing for a product by a constant 0 or 1), so "a + b n" with a
constant b is false and my comment about TMB sharing operations was
wrong. The reviewer's seed sweep, rerun here
(`dev/ciharden-rev-perfseeds.R`, `p1/perfseeds-*.txt`): GLM, 30 pairs,
node ratio 98.618 to 101.040, 13 of 30 over the old bound of 100;
GLMM, 10 pairs, 76.944 to 79.256. And the node count cannot see the
item-14 loop: `dev/ciharden-rev-perfmut.R` (`p1/perfmut.txt`), the
loop inside `build_objective()`'s objective, nodes 9.92 at tenfold n
and 99.43 at hundredfold, bytes 97.0 and 9676.5. Changed:

- Both counts are bounded by `2 * fold`. Against the sweep the node
  maximum is 101.04 of 200, the bytes 97.31 of 200 (GLM; the bytes
  ratio is 97.31 for every pair in the sweep's process order, 58.81 in
  the test's, because one-time allocations land in whichever call comes
  first: the headroom is 2.06x, not 3.4x, as the review's m6 says).
- The comment states the real reasons, and that the allocation count
  is the item-14 canary while the node count guards only the tape's
  size.
- New guard test "the census catches an elementwise loop in the
  objective": the reviewer's mutation, mocked into `build_objective()`
  with `local_mocked_bindings()`, at a tenfold n (2 s): bytes ratio
  must exceed 20 (it is 97.0), nodes stay within 20 (9.92).
- `test-perf.R`: pass=10 fail=0, reference BLAS and OpenBLAS 0.3.32
  (`p1-t-*.sum`).

**B2, the gp() key.** Reverted to brms's rule, as the tiebreaker says.
brms 2.23.0 groups NEW rows at 15 significant digits too
(`dev/ciharden-rev-brmsgpnew.R`: `Jgp_1` = 1 1 1 2 2 for rows `u,
u (1 + 2^-52), u, 1/3, 1/3 (1 + 2^-52)`), and my exact key also made
frmtmb inconsistent with itself (a new row equal to a fitted position
at 15 digits was that position, two new rows equal at 15 digits were
two).

- `pos_rowkey()` is brms's paste key again, and the one key for the
  frame, the fitted-position match and the unseen-row grouping;
  `gp_fit_poskey()` is gone and `R/frame.R` is back to the base
  commit. `gp_krig_cov()`'s one-allocation change stays.
- The NEWS bullet that called brms's grouping a bug is gone.
- Tests: "the kriging key is exact" and "two unseen positions one bit
  apart are two positions" are gone. New: "unseen rows make one
  position as brms's new-data rule does", the reviewer's construction:
  frmtmb's grouping, read from which rows of `gp_krig_cov()` are
  identical, must be 1 1 1 2 2 and equal brms's `standata(newdata =
  )$Jgp_1` (an empty brmsfit, no compile). Seen to fail on the punch-0
  lane build (`p1-gpm-oldlane-ref.sum`: fail=2, `actual: 1 2 1 3 4`);
  passes on rellib-r6 and on this build.
- frmtmb.spline's `sp_row_key()`: not changed, on a measurement.
  `dev/ciharden-splinekey.R` (fit `y ~ x + gp(x)`, rows `u` and
  `u (1 + 2^-52)` for four values of u): core gives the two rows
  identical gp columns (one position) but different `x` columns, so
  core's linear predictor differs between them by an ulp or more in
  all four cases with the reference BLAS and in three of four with
  OpenBLAS 0.3.32, where u = 7/3 rounds to the same value for both rows
  (for example `0x1.ddb22d49b5622p-2` and `0x1.ddb22d49b5626p-2`; the
  probe compares the values with `unname()` since the review's r1). A
  15-digit spline key would merge three of
  the four pairs and report row 1's value for row 2, disagreeing with
  core's own prediction at row 2. The spline's exact key merges only
  rows that are identical in every column, which core treats
  identically everywhere: its merges are a subset of core's, so the
  two are consistent, and the 0.68.1 property (identical rows
  evaluated once) is untouched. spline `test-difference.R` 71/0 and
  `test-curve.R` 82/0 against this core, both BLAS builds.

**Minors.**

- m1: `base_r_function()` decides a package that is not loaded by
  `exists()` alone and never reads the value, which is what loaded it
  (`tclVar` loaded tcltk). Its non-function objects then count as
  functions; none of the registry's 27 call-shaped names exists in any
  of those packages (`dev/ciharden-regexists.R`). Roxygen corrected.
  The preflight test asserts the package is still not loaded after a
  hit: on the punch-0 build FAIL (`p1-pf-oldlane-ref.sum`,
  `isNamespaceLoaded(unl[1])` TRUE), on this build pass=86.
- m2: eam `test-sampling.R` asserts 0 divergent transitions and
  `|posterior mean - mode| / posterior sd < 3` per parameter (measured
  at most 0.152 and 0.096 by the review): pass=100, both BLAS builds.
- m3: R-CMD-check.yaml's pin comment is the dated one of the others.
- m4: brms-likelihood.yaml's pin comment says the Stan cache prefix
  must gain the image in the commit that moves the pin.
- m5: filed in `dev/test-backlog.md`, "Filed by lane ciharden": the
  Malingering_2 precision with an infinite MLE fits with no warning
  under the reference BLAS; lane setier's boundary work may cover it.
- m6: recorded above (B1).

Verification of this round (lane library reinstalled with this
round's core and frmtmb.sample; every `lib:` line names it):

- The named files, reference BLAS and OpenBLAS 0.3.32 LAPACK 4
  threads alike (`p1-t-*.sum`): `test-perf.R` 10/0, `test-gp-multidim.R`
  67/0, `test-gp-by.R` 86/0, spline `test-difference.R` 71/0 and
  `test-curve.R` 82/0, eam `test-sampling.R` 100/0, sample
  `test-compat-preflight.R` 86/0; all with 0 error, skip and warn.
- The full scan, all 354 files, every gate on (`p1-ref.sum`,
  `p1-ob0.3.32-lapack-t4.sum`, "ran 354 of 354" each; diff
  `p1/p1-diff.out`): 0 configuration-dependent signatures; the one
  pass count that differs is `test-v11.R` (1 against 2), the
  coverage-only row. Per package, pass/fail/error/skip/warn, reference
  then OpenBLAS: frmtmb 215 files 17265 and 17266/0/0/0/0; coupling
  548; eam 1764; latent 366; learn 515; ode 551; spline 592, all
  0/0/0/0; frmtmb.sample 2618/0/2/0/0 in both.
- The 2 frmtmb.sample errors are mine, not the code's:
  `test-loo.R:688` and `:736` stopped in `compileCode()` ("sink ...
  invalid connection") in BOTH runs, in files that finished at
  22:44:34 and 22:44:40, about 110 s after they started, which is when
  I killed (`taskkill /T`) a shell whose process tree held a
  mis-launched R CMD check, and with it, by all appearances, the two
  compilers brms had started under it. `test-loo.R` alone afterwards,
  each configuration: 111/0/0/0/0 (`p1-loo-ref.sum`, `p1-loo-ob.sum`).
  Corrected, the frmtmb.sample total is 2626/0/0/0/0 in both.
- `R CMD check --as-cran`, only the packages whose code changed this
  round: frmtmb, Status: 1 NOTE (V8), tests FAIL 0 WARN 0 SKIP 198
  PASS 14710; frmtmb.sample, Status: OK, FAIL 0 WARN 0 SKIP 12 PASS
  2393. (The first core check of this round died at build: the same
  `taskkill` took its vignette install; rerun clean.)

## Summary

- All eight suites, every file in its own process, every gate on, ran
  under five BLAS and thread configurations plus a second reference
  run, on the base and on the lane. On the base, 7 test files depend
  on the configuration (6 by outcome, 1 by coverage only); the two
  reference runs agree exactly. On the lane, every configuration gives
  FAIL 0, ERROR 0, SKIP 0, WARN 0 in all 354 files, and the only
  difference left is the coverage-only one (`test-v11.R`).
- OpenBLAS 0.3.32, ubuntu-26.04's version, gives exactly the
  failures 0.3.26 as BLAS and LAPACK gives, with 1 or 4 threads. The
  move to 26.04 adds no rounding failure that 24.04 does not already
  show here.
- Two frmtmb.sample warnings on the runner were package defects, now
  fixed: `posterior_samples()` warned twice with brms loaded, and the
  `frm_sample()` pre-flight loaded tcltk.
- `test-perf.R` counts tape operations and allocated bytes instead of
  timing.
- Core's gp() position key stays brms's 15-digit rule everywhere (an
  exact key for unseen rows was reverted in punch round 1).
  `gp_krig_cov()` makes 1 n x n allocation where it made 3.
- The ported brms suite is seeded per block; two records agree on all
  865 rows.
- Every workflow is pinned to `ubuntu-24.04`; a weekly
  `ubuntu-next.yaml` checks core on 24.04 and 26.04 side by side.

## What changes on 2026-10-19

GitHub's changelog of 2026-09-17 ("Ubuntu 26 generally available and
latest migration"): `ubuntu-latest` moves from Ubuntu 24.04 to Ubuntu
26.04 between 2026-10-19 and 2026-11-19, rolled out gradually; the
image `ubuntu-26.04` is generally available now. The runner-images
issue 14226 lists the default changes: Clang 21 (was 14), GCC 13, 14
and 15 installed, Python 3.14, and some tools removed (none this
project uses). It names no R, BLAS or LAPACK versions, because the R
jobs install those themselves:

- `r-lib/actions/setup-r` installs OpenBLAS from apt as both BLAS and
  LAPACK (`update-alternatives` to `openblas-pthread/libblas.so.3` and
  `liblapack.so.3`, seen in the 0.68.1 check-frmtmb.sample log,
  `dev/ciharden-log/ci-sample-0681.log`). On 24.04 that is
  `libopenblas0-pthread 0.3.26+ds-1ubuntu0.1`; on 26.04 (resolute) it
  is `0.3.32+ds-5` (packages.ubuntu.com, 2026-10-06).
- R itself and the CRAN binaries: Posit's R builds list `ubuntu-2604`
  and Posit Public Package Manager serves binaries for `resolute`
  (api.r-hub.io/rversions/linux-distros and
  packagemanager.posit.co/__api__/status, both 2026-10-06), so the
  move should not fall back to source builds of rstan or TMB.

So the numerical change that reaches this package is OpenBLAS 0.3.26
to 0.3.32, for BLAS and LAPACK, plus whatever the new compiler does to
the packages built from source on the runner (none of the R packages
are; the CRAN binaries are built by Posit for each distribution).

## The instrument: OpenBLAS builds of R on Windows

`dev/ciharden-openblas.sh <version> <blas|lapack>` generalizes
`dev/cifix-openblas.sh` (fixed to 0.3.26, BLAS only, and to a worktree
that no longer exists). It copies R 4.6.1 into
`dev/ciharden-out/Rob<version>-<mode>` and replaces `Rblas.dll` (and,
in `lapack` mode, `Rlapack.dll`) by a forwarding DLL that sends every
exported symbol OpenBLAS has to `libopenblas.dll` and the rest to the
renamed original. In `lapack` mode LAPACK is OpenBLAS's, as on the
runner, where OpenBLAS replaces getrf, potrf, trtri, lauum and others
with its own blocked code. `dev/ciharden-blascheck.R`
(`dev/ciharden-ob*.log`):

| Build | forwarded | dgemm 2000 (s) | dpotrf 3000 (s) | dgesv 2000 (s) |
|---|---|---|---|---|
| 0.3.26, BLAS | Rblas 74 of 76 | 0.06 | 1.08 | 1.18 |
| 0.3.26, BLAS and LAPACK | Rblas 74, Rlapack 616 of 626 | 0.05 | 0.37 | 0.22 |
| 0.3.32, BLAS and LAPACK | Rblas 76, Rlapack 617 of 626 | 0.04 | 0.11 | 0.27 |

Configurations scanned, each named as `dev/ciharden-scan.sh` takes it:

- `ref`: this R, reference BLAS and LAPACK (the Windows runner's).
- `ob0.3.26-blas-t4`: cifix's emulator, OPENBLAS_NUM_THREADS=4.
- `ob0.3.26-lapack-t4`: ubuntu-24.04 as the runner has it, 4 threads
  (the runner has 4 cores and sets no thread count).
- `ob0.3.32-lapack-t4`: ubuntu-26.04.
- `ob0.3.32-lapack-t1`: ubuntu-26.04, one thread
  (OPENBLAS_NUM_THREADS = OMP_NUM_THREADS = 1).
- `ref` a second time (`b-ref2`), the noise floor: fresh processes on
  this hybrid CPU can differ by a few ulps (`dev/rtmb-pitfalls.md`
  item 21), and wall-clock tests move with load.

## Item 1: what depends on the configuration

Every run: all 354 test files of the eight suites, one R process each,
NOT_CRAN, FRMTMB_BRMS_FIT_TESTS, FRMTMB_DRMTMB_FIT_TESTS,
FRMTMB_FUZZ and FRMTMB_SCALE_TESTS set, sampler gates on. "Before"
runs the base commit's own test files (an export of 4f5ea39f,
`TREE=dev/ciharden-out/base-tree`) on rellib-r6; "after" runs this
worktree on the lane library. Every log's `lib:` line names the
library it claims (checked for all 8 packages; 0 lines name another
frmtmb). The comparison is `dev/ciharden-scan-diff.R`
(`b-diff.out`, `a-diff.out`); the table is `dev/ciharden-table.R`
(`table.md`), pasted verbatim.

Noise floor first: the two reference runs `b-ref` and `b-ref2`
differ in no expectation and no pass count (`b-noise`), so every row
below is the configuration.

Counts are pass/fail/error/skip/warn.

| File | arm | ref | ob0.3.26-blas-t4 | ob0.3.26-lapack-t4 | ob0.3.32-lapack-t4 | ob0.3.32-lapack-t1 |
|---|---| ---|---|---|---|--- |
| frmtmb test-bcm-latent-mixtures | before | 52/0/0/0/0 | 51/1/0/0/0 | 51/1/0/0/0 | 51/1/0/0/0 | 51/1/0/0/0 |
|  | after | 52/0/0/0/0 | 52/0/0/0/0 | 52/0/0/0/0 | 52/0/0/0/0 | 52/0/0/0/0 |
| frmtmb test-cumulative-cs | before | 52/0/0/0/0 | 47/0/1/0/0 | 47/0/1/0/0 | 47/0/1/0/0 | 47/0/1/0/0 |
|  | after | 52/0/0/0/0 | 52/0/0/0/0 | 52/0/0/0/0 | 52/0/0/0/0 | 52/0/0/0/0 |
| frmtmb test-ordinal-mixture | before | 99/0/0/0/0 | 99/0/0/0/1 | 99/0/0/0/1 | 99/0/0/0/1 | 99/0/0/0/1 |
|  | after | 99/0/0/0/0 | 99/0/0/0/0 | 99/0/0/0/0 | 99/0/0/0/0 | 99/0/0/0/0 |
| frmtmb test-se-check | before | 123/0/0/0/0 | 123/0/0/0/0 | 119/0/0/1/0 | 119/0/0/1/0 | 119/0/0/1/0 |
|  | after | 119/0/0/0/0 | 119/0/0/0/0 | 119/0/0/0/0 | 119/0/0/0/0 | 119/0/0/0/0 |
| frmtmb test-v11 | before | 22/0/0/0/0 | 22/0/0/0/0 | 23/0/0/0/0 | 23/0/0/0/0 | 23/0/0/0/0 |
|  | after | 22/0/0/0/0 | 22/0/0/0/0 | 23/0/0/0/0 | 23/0/0/0/0 | 23/0/0/0/0 |
| frmtmb.eam test-sampling | before | 97/0/0/0/0 | 97/0/0/0/2 | 97/0/0/0/2 | 97/0/0/0/2 | 97/0/0/0/2 |
|  | after | 97/0/0/0/0 | 97/0/0/0/0 | 97/0/0/0/0 | 97/0/0/0/0 | 97/0/0/0/0 |
| frmtmb.sample test-brms-shapes-draws | before | 76/0/0/0/0 | 76/0/0/0/1 | 76/0/0/0/1 | 76/0/0/0/1 | 76/0/0/0/1 |
|  | after | 78/0/0/0/0 | 78/0/0/0/0 | 78/0/0/0/0 | 78/0/0/0/0 | 78/0/0/0/0 |

Suite totals, every file of the eight suites:

| Run | files with a RESULT | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| b-ref | 354 of 354 | 24199 | 0 | 0 | 0 | 0 |
| b-ob0.3.26-blas-t4 | 354 of 354 | 24193 | 1 | 1 | 0 | 4 |
| b-ob0.3.26-lapack-t4 | 354 of 354 | 24190 | 1 | 1 | 1 | 4 |
| b-ob0.3.32-lapack-t4 | 354 of 354 | 24190 | 1 | 1 | 1 | 4 |
| b-ob0.3.32-lapack-t1 | 354 of 354 | 24190 | 1 | 1 | 1 | 4 |
| b-ref2 | 354 of 354 | 24199 | 0 | 0 | 0 | 0 |
| a-ref | 354 of 354 | 24225 | 0 | 0 | 0 | 0 |
| a-ob0.3.26-blas-t4 | 354 of 354 | 24225 | 0 | 0 | 0 | 0 |
| a-ob0.3.26-lapack-t4 | 354 of 354 | 24226 | 0 | 0 | 0 | 0 |
| a-ob0.3.32-lapack-t4 | 354 of 354 | 24226 | 0 | 0 | 0 | 0 |
| a-ob0.3.32-lapack-t1 | 354 of 354 | 24226 | 0 | 0 | 0 | 0 |

Per package, after, every configuration alike except the one count
(`a-*.sum`): frmtmb 215 files 17267 or 17268 passes; coupling 11
files 548; eam 29 files 1761; latent 10 files 366; learn 15 files 515;
ode 11 files 551; sample 48 files 2625; spline 15 files 592. All with
0 fail, 0 error, 0 skip, 0 warn.

What the configurations say about Ubuntu 26.04: 0.3.32 as BLAS and
LAPACK gives exactly what 0.3.26 as BLAS and LAPACK gives, with 4
threads and with 1. LAPACK matters (`test-se-check.R` and
`test-v11.R` move between 0.3.26 as BLAS only and as both); the
OpenBLAS version and the thread count did not, on this suite.

The rows, and what was done:

- `test-bcm-latent-mixtures.R:322`, "Malingering_2 matches its Stan
  program": FAIL with every OpenBLAS, not named by the task (the Stan
  identity tier is gated, and no CI job runs it). `phi2` runs to the
  binomial limit, and where nlminb stops on that flat ridge is
  rounding: 7.19e8 with the reference BLAS, 2.13e9 with OpenBLAS
  (`dev/ciharden-bcmmal.R`). There the log-gamma terms the two
  programs cancel sum to S = 1.23e12 and 3.83e12, and the two log
  densities disagree by 2.18e-5 and 1.30e-4: 0.080 and 0.152 of
  eps S, rounding of differently written formulas. The bound was
  1e-6 of the density, 5.9e-5. Fix: the bound is now the larger of
  that and 4 eps S, computed from the fitted point.
- `test-cumulative-cs.R:132`: ERROR "NA/NaN gradient evaluation" in
  the `mixture(cumulative, sratio)` fit with every OpenBLAS (5 of the
  test's 7 passes lost). The test asserts how often brms's cs()
  warning is given, and the frame gives it (`R/frame.R`), so the two
  mixture cases now build the frame only (`dry_run = "frame"`). A
  cs() mixture's fit is fragile in its own right (the backlog's
  "cs() on a cumulative component of an ordinal mixture often stops
  in the optimizer": 7 NaN-gradient errors in 20 seeds), which says
  nothing about the warning.
- `test-ordinal-mixture.R:751`: "singular convergence (7)" escapes
  with every OpenBLAS, and on the ubuntu devel runner of 0.68.1.
  Three components on two classes have a ridge; the convergence
  verdict there is the optimizer's, so it is allowed (not required)
  beside the required degenerate-boundary warning.
- `test-se-check.R:575`: SKIP with OpenBLAS as LAPACK ("this fit kept
  every standard error here"), 4 assertions lost: the smooth-plus-gp
  fit loses its smoothing sds with the reference LAPACK and keeps them
  with OpenBLAS's. On the runner's configuration those 4 assertions
  never ran. Fix: the platform-dependent half is gone from that test,
  and `frm_joint_cov()`'s NaN rows for a lost parameter are asserted in
  `test-se-lost-re-predict.R` on the OLRE fit, which loses sigma and
  the id sd by construction (one observation per id): "frm_joint_cov()
  shows the lost parameters as NaN", 6 assertions, passing on every
  configuration (it passes on rellib-r6 too; it moves coverage, it
  pins no defect).
- `test-v11.R`, "toep matches a hand-rolled reference (glmmTMB when it
  converges)": one more assertion with OpenBLAS as LAPACK, because
  glmmTMB's toep fit converges there and the comparison runs; never a
  failure. Where it runs the difference is 1.19e-5 against its bound
  of 1e-4 (`dev/ciharden-v11toep.R`). Left as is; it is the one
  difference the after runs keep.
- frmtmb.eam `test-sampling.R:29`: two ESS warnings escape with every
  OpenBLAS, and on the 0.68.1 runner. The chain is one chain of 200
  draws, and the test asserts that it runs and sits at the mode, so
  the ESS and R-hat warnings are allowed (not required).
- frmtmb.sample `test-brms-shapes-draws.R:196`: an R-hat warning
  escapes with every OpenBLAS (1.07) and on the 0.68.1 runner (1.11).
  The file's `short_chain` allowance had the two ESS warnings and not
  R-hat; it has all three now.

No fragile test asserted bitwise identity. The new `expect_identical()`
calls of this lane rest on construction: the same arithmetic per entry
(`gp_krig_cov()` scaled), a copied row (one position twice), and
string keys.

## How to rerun the scan on a merged tree

The scan is three scripts; other lanes change numerical code, so the
list of fragile tests is meant to be rerun on the merged tree.

1. Build the R copies once (about a minute each; output in
   `dev/ciharden-out/`, gitignored):

       bash dev/ciharden-openblas.sh 0.3.26 lapack
       bash dev/ciharden-openblas.sh 0.3.32 lapack

2. Run every test file of the eight suites, one R process per file,
   every gate on (NOT_CRAN, FRMTMB_BRMS_FIT_TESTS,
   FRMTMB_DRMTMB_FIT_TESTS, FRMTMB_FUZZ, FRMTMB_SCALE_TESTS; the
   sampler gates at their default, on), once per configuration:

       bash dev/ciharden-scan.sh m-ref ref <lib> 16
       bash dev/ciharden-scan.sh m-ref2 ref <lib> 16
       bash dev/ciharden-scan.sh m-26 ob0.3.26-lapack-t4 <lib> 16
       bash dev/ciharden-scan.sh m-32 ob0.3.32-lapack-t4 <lib> 16
       bash dev/ciharden-scan.sh m-32t1 ob0.3.32-lapack-t1 <lib> 16

   `<lib>` is the library holding the merged build, put ahead of the
   round's base library; set `ROUND_BASE` to that base library when it
   is not rellib-r6, and `TREE` to run another checkout's test files.
   The runner (`dev/ciharden-run1.R`) prints every non-passing
   expectation with its source line, and each test's pass count.

3. Compare:

       Rscript dev/ciharden-scan-diff.R dev/ciharden-log/m-diff \
         m-ref m-ref2 m-26 m-32 m-32t1

   It lists every expectation (kind, file:line, test, message with
   numbers masked) present in some runs and absent in others, and
   every test whose pass count differs, and writes the same as a TSV.
   A row that differs between the two `ref` runs is load or core-type
   noise, not BLAS; read it before calling anything a BLAS effect.

## What the 0.68.1 runners already show

The 0.68.1 runs on 4f5ea39f are green, but they let warnings escape
(`gh run view <id> --log`, saved as `dev/ciharden-log/ci-0681-*.log`
and `ci-sample-0681.log`):

| Job | testthat line | Escaped |
|---|---|---|
| R-CMD-check, ubuntu-latest (devel) | FAIL 0, WARN 1, SKIP 199, PASS 14689 | `test-ordinal-mixture.R:751` |
| check-frmtmb.eam | FAIL 0, WARN 2, SKIP 6, PASS 1732 | `test-sampling.R:29`, twice |
| check-frmtmb.sample | FAIL 0, WARN 3, SKIP 18, PASS 2274 | `test-brms-shapes-draws.R:106` and `:196`, `test-compat-preflight.R:95` |

Every other job has WARN 0. Each of these is fixed below.

## Item 2: the two frmtmb.sample warnings, both package defects

Neither is fixed by a wrapper alone, because each is the package doing
something wrong.

`test-brms-shapes-draws.R:106`. The runner's backtrace shows the call
reaching `brms::posterior_samples()`, brms's generic. brms's generic
warns "Method 'posterior_samples' is deprecated" BEFORE it dispatches,
and `posterior_samples.frmtmb_draws()` warned again, so one call gave
two warnings whenever brms was loaded (in R CMD check every file runs
in one process, and an earlier file loads brms). `expect_warning()`
absorbs one. Measured with `dev/ciharden-psamples.R` on rellib-r6:
frmtmb.sample's own generic, 1 warning; `brms::posterior_samples(ds)`,
2 warnings, and 1 escapes `expect_warning()`. Fix
(`R/methods-draws.R`): the method warns unless the generic directly
above it is brms's (`ps_via_brms_generic()`), as `parnames()` already
avoids its double warning. The test now uses `allow_warnings(...,
require = "is deprecated")`. New test "posterior_samples() warns once,
through brms's generic too": on rellib-r6 FAIL, `actual: 2, expected:
1` (`dev/ciharden-log/smp2-base-ref.sum`: pass=77 fail=1); on the lane
pass=78 fail=0.

`test-compat-preflight.R:95`. `dev/ciharden-tcltk.R` traces the load:
`sample_preflight()` -> `preflight_route()` -> `base_r_function()` ->
`asNamespace("tcltk")`. To decide whether a registry name is also a
function of a base-priority package, it LOADED all fourteen
namespaces, tcltk among them, and tcltk's `.onLoad()` warns "no
DISPLAY variable so Tk is not available" on a headless Linux machine.
A user calling `frm_sample()` on a server saw the same warning. Fix
(`R/sample.R`, `base_r_objects()`): a namespace that is loaded is
read; one that is not is read from its lazy-load database into an
environment of promises, which runs no `.onLoad()`. New test "the
base-function test reads packages without loading them": on rellib-r6
FAIL, `isNamespaceLoaded("tcltk")` TRUE (`smp-base-ref.sum`: pass=82
fail=1 skip=1); on the lane pass=85. Windows shows no warning either
way (it has a display), which is why the test pins the load and not
the warning.

## Item 3: test-perf.R counts instead of timing (bounds and reasons revised in punch round 1)

The scaling test failed under load in this lane's first scan too
(`b-ref1`, 20 processes: `t_large` 1.31 s against a bound of 1.00).
Both tests of the file were wall-clock bounds; the tape canary was an
absolute 20 s. Both now count, with `perf_census()`: the tape's
operations (`nrow(RTMB::GetTape(obj)$data.frame())`) and the bytes R
allocates (`Rprofmem(threshold = 0)`) while the frame is built, the
objective taped and evaluated once with its gradient, at n = 1000 and
n = 100000. Linear cost `a + b n`, `a >= 0`, bounds each ratio by 100.

`dev/ciharden-perfcount2.R`, the test's own `perf_census()`, seeds as
in the test; identical under the reference BLAS and OpenBLAS 0.3.32
with 4 threads, run after run:

| Model | nodes, n = 1e3 | nodes, n = 1e5 | ratio | bytes, n = 1e3 | bytes, n = 1e5 | ratio |
|---|---|---|---|---|---|---|
| Poisson GLM | 8341 | 829055 | 99.40 | 596200 | 35060272 | 58.81 |
| Poisson GLMM, 500 groups | 11829 | 912105 | 77.11 | 3132552 | 113108496 | 36.11 |

The node bound is 100: TMB shares operations between equal inputs, so
the count grows no faster than n. (Replicating the small data 100
times instead of drawing new data gives a node ratio of 12.9, because
TMB merges the duplicated rows; that construction was dropped.) The
allocation bound is 200, twice linear, because R rounds some sizes,
hash tables among them, up to a power of two. The instrument catches
the regression it guards: an objective with `dev/rtmb-pitfalls.md`
item 14's elementwise sub-assignment loop allocates 1.658e7 bytes at
n = 1e3 and 1.602e9 at n = 1e4, a ratio of 96.6 for a TENFOLD n
(`dev/ciharden-perfcount.R`). The whole fit is not counted: its
allocation follows the optimizer's evaluation count (17 at n = 1e3, 10
at n = 1e5 on these seeds), and its ratio was 96.21, too near the
bound to mean anything.

## Item 4: the gp() position key (the exact key was reverted in punch round 1; the record below is the first round)

brms decides first. `dev/ciharden-brmsgp.R` runs brms 2.23.0's
`standata(y ~ gp(x))` on `x = c(1/3, 1/3 * (1 + 2^-52), ...)`:
`Xgp_1` has 9 rows for 10 distinct doubles and `Jgp_1` is
`1 1 2 3 ...`. brms's `match_rows()` pastes the coordinates, so brms
itself makes rows equal to 15 significant digits one latent position.
So the key is split, not replaced:

- `gp_fit_poskey()` keeps brms's rule for the fitted positions: the
  frame's `gr = TRUE` grouping and `pred_design()`'s match of a new
  row to a fitted position. A row of the fitting data must find its
  own position again; with an exact key, `1/3 * (1 + 2^-52)` would
  have been kriged with the nugget's variance instead of read.
- `pos_rowkey()` is exact (`sprintf("%a")` per coordinate, -0 made 0)
  for everything that groups UNSEEN rows: `pred_design()`'s
  one-kriging-per-position, `gp_krig_cov()` and `gp_krig_factor()`.
  That is the line frmtmb.spline's `sp_row_key()` draws (equal as
  doubles under `match()`), so the two layers now agree on what one
  row is.

What changes, `dev/ciharden-poskey.R` (test-gp-multidim.R's fit, seed
17; new rows `u, u (1 + 2^-52), u`, `u = 10/3`): on rellib-r6 the
kriging covariance of the first two rows equals their variance,
7.870121455004321e-07, the gap `(S11 - S12) / (gp_nugget sd^2)` is 0
and the two design rows are identical; on the lane the covariance is
2.806895359325411e-07, the gap is 1.00000000003 (two positions, each
with its own nugget) and the design rows differ by at most 1.64e-10.
Tests: "the kriging key is exact and the fitted key is brms's"
(against `brms:::match_rows()` itself) and "two unseen positions one
bit apart are two positions" fail on rellib-r6 (`gpm-base-ref.sum`:
fail=2 error=1; the second is the behavioral one, a gap of 0 against
a bound of 1 +/- 64 eps / gp_nugget); "rows of the data one bit apart
find their fitted position" is the guard on the other side, and passes
on both.

`gp_krig_cov()` and memory. `dev/ciharden-krigmem.R`, a
`gp(x, by = w)` fit and a 2000-row grid, Rprofmem's allocations of at
least n^2 / 2 doubles inside the call: rellib-r6 makes 3 (the matrix,
a copy inside `diag<-`, and `outer(w, w)`; the product reuses
`outer`'s buffer), total 3.00 n^2; the lane makes 1, the matrix. The
copy in `diag<-` happened on every exact gp(), not only with `by`:
`diag<-` is a closure that duplicates its argument. Fix: the diagonal
is written by index, and a numeric `by` scales a column block at a
time with `outer(w, w[cc])`. The result is the same matrix to the bit:
each entry is still `S[i, j] * (w[i] * w[j])`, and a product of two
doubles does not depend on their order, so it stays symmetric. New
test "a numeric by scales the kriging covariance entry by entry" (850
rows, two column blocks, 50 repeated rows): identical to
`gp_krig_cov()` with `w = 1` times `outer(w, w)`, to its transpose,
and its diagonal to `extra_var`. It passes on rellib-r6 too; it guards
the rewrite, and the allocation count is the evidence for the change.
gc()'s "max used" moved only from 3.01 to 2.82 n^2, because it counts
garbage not yet collected; it is not the instrument.

## Item 5: the ported brms suite, seeded

`dev/brmsport-gen.R` now writes `withr::local_seed(<line>L)` as the
first line of every ported block, the line of the block's
`test_that()` in brms's file, so each block draws the same data at
every run and leaves the caller's RNG alone. Regenerated: the diff of
the 14 generated files is the seed lines and nothing else (the helper
copy is unchanged). `dev/rel068-msgdiff.R` is not changed; seeding
removes the cause.

`dev/ciharden-brmsstable.sh` records every generated file twice
(`dev/brmsport-run.R`, record mode, one process per file):

| Tree | rows | rows that differ between two runs |
|---|---|---|
| base (unseeded), rellib-r6 | 865 | 3: standata:310, :315, :316 |
| lane (seeded), lane build | 865 | 0 |

Seeding changes no outcome: the `held` column of all 865 rows is the
same on the base and the lane (`dev/ciharden-brmsheld.R`). No ledger
row flips.

## Item 6: the CI guard

Decision: pin, and watch the new image on a schedule.

- All 15 `ubuntu-latest` entries, in 13 workflows, are now
  `ubuntu-24.04`, each with a dated comment. The move to 26.04 then
  happens when someone chooses it, not on a date GitHub chooses.
- New `.github/workflows/ubuntu-next.yaml`: weekly (Monday 05:17 UTC)
  and on demand, the core R CMD check on `ubuntu-24.04` and
  `ubuntu-26.04` side by side, with a step that prints the BLAS,
  LAPACK and compiler. It does not run on push, so it cannot turn a
  push red; it reports. When it is green on 26.04, change the pins
  (search the workflows for "Pinned 2026-10-06").

Why not a matrix entry on every push: it doubles the Linux minutes of
the core check, and a red 26.04 column would mark every push failed
until it is fixed. The pin is the cheaper safe option; the scheduled
job is how the surprise is found before it is forced.

Proposed, then decided: the Stan-compiling job (`brms-likelihood.yaml`)
is where 26.04's newer GCC can break something the core check cannot
see (rstan compiles against StanHeaders with the runner's g++), and its
compiled-model cache key does not include the image, so a 26.04 run
could save programs that a 24.04 run would restore. User decision of
2026-10-07: try the Stan job on 26.04, in `ubuntu-next.yaml` only,
never on push. Done: a second job `stan-next` there runs
brms-likelihood.yaml's steps on ubuntu-24.04 and ubuntu-26.04, prints
g++, R's CXX17, BLAS, LAPACK, rstan, StanHeaders, brms and Stan
versions, and keys its cache `brms-stan-next-<image>-<rstan>-<hash>`
with the same restore prefix, so the images never share programs and
neither shares with the push job (prefix `brms-stan-<rstan>-`).
`brms-likelihood.yaml` is not changed. actionlint 1.7.12 passes on
all 14 workflows once `ubuntu-26.04` is declared as a label (its
built-in label list predates that image, which is its only finding
otherwise; shellcheck not run).

The 26.04 image itself could not be run from here: nothing is pushed
from a lane. What the emulator cannot see is the runner's CPU kernels
and the new compiler; `ubuntu-next.yaml`'s first run is that evidence.

## R CMD check --as-cran (first round; punch round 1 numbers are at the top)

Built and checked in `dev/ciharden-check/<pkg>/` (`dev/ciharden-check.sh`,
a copy of cifix's with this lane's library), NOT_CRAN=true,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, R_LIBS = lane library,
rellib-r6, user library:

| Package | Status | tests |
|---|---|---|
| frmtmb 0.68.1 + lane | 1 NOTE (V8 math rendering) | FAIL 0, WARN 0, SKIP 198, PASS 14712 |
| frmtmb.sample 0.16.0 + lane | OK | FAIL 0, WARN 0, SKIP 12, PASS 2392 |
| frmtmb.eam 0.11.2 + lane (tests only) | 1 NOTE (V8 math rendering) | FAIL 0, WARN 0, SKIP 6, PASS 1732 |

## Not done, and why

- `test-v11.R`'s toep comparison against glmmTMB runs only where
  glmmTMB converges, which is a LAPACK matter; its result is the same
  everywhere and its margin is 8x (1.19e-5 against 1e-4). Its absolute
  tolerances are older than this lane, like many in the suite
  (`expect_lt(abs(...), 1e-4)`); the scan shows none of them failing
  under any configuration, so none was changed here.
- frmtmb.spline is not changed. Its `sp_row_key()` already compares
  exactly, and core now draws the same line for unseen rows; nothing
  in the spline had to move for the two to agree.
- `dev/rel068-msgdiff.R` keeps its exact comparison: seeding removed
  the moving messages, and masking numbers would also hide a real
  change of a recorded message.
- The ledger (`dev/brmsport-*`) was not regenerated, per the round's
  rule; no verdict changes (all 865 `held` values equal).
- `dev/test-backlog.md` is not edited. This lane closes, from "Filed
  at 0.68.1": the two frmtmb.sample runner warnings, the two
  rounding-fragile core tests (`test-cumulative-cs.R:132`,
  `test-ordinal-mixture.R:751`), the gp() position key and the
  `gp_krig_cov()` `outer(w, w)`; from "Filed at the 0.68.0 release":
  the unseeded ported suite and the `test-perf.R` wall-clock bound.

## Defects found, not fixed

- `parnames.frmtmb_fit()` and `parnames.frmtmb_draws()` decide whether
  to warn from which generic the `parnames` BINDING resolves to, not
  which generic dispatched; a call through a generic the binding does
  not name (for example `brms::parnames()` while another package's
  binding wins) can still warn twice or not at all. Not measured to
  happen; `posterior_samples()` now asks the dispatching frame instead,
  and the two could share that test.
- The cs() mixture fit's NaN-gradient stop (backlog, "cs() on a
  cumulative component of an ordinal mixture often stops in the
  optimizer") is unchanged; the test no longer depends on it.
- The Malingering_2 fit stops wherever the flat phi2 ridge lets it
  (7.19e8 or 2.13e9). That is the documented infinite MLE, not a
  defect, but every quantity read at that point is rounding-sized
  noise in phi2.

## Version and floor

frmtmb: patch (bug fixes in prediction keys and memory, no API).
frmtmb.sample: patch (a double warning and a load side effect). Neither
needs a new floor: frmtmb.sample uses nothing new from core, and the
core change is internal to prediction. frmtmb.eam: tests only, no
version change needed.

## Files touched

Package code, after punch round 1: `R/utils.R` (`pos_rowkey()`'s
roxygen only; its body is the base's), `R/predict.R` (`gp_krig_cov()`
diagonal by index and blockwise `by` scaling; the doc line),
`extensions/frmtmb.sample/R/methods-draws.R` (`ps_via_brms_generic()`),
`extensions/frmtmb.sample/R/sample.R` (`base_r_objects()`, and
`exists()` alone for an unloaded package). `R/frame.R` is unchanged.
NEWS: `NEWS.md`, `extensions/frmtmb.sample/NEWS.md`.
Tests: core `test-gp-multidim.R`, `test-gp-by.R`, `test-perf.R`,
`test-cumulative-cs.R`, `test-ordinal-mixture.R`,
`test-bcm-latent-mixtures.R`, `test-se-check.R`,
`test-se-lost-re-predict.R`, the 10 regenerated
`test-brms-suite-*.R`; sample `test-brms-shapes-draws.R`,
`test-compat-preflight.R`, the 4 regenerated `test-brms-suite-*.R`;
eam `test-sampling.R`.
CI: 13 workflows pinned, `ubuntu-next.yaml` new.
dev: `brmsport-gen.R` (seed line), `test-backlog.md` (m5); new
`ciharden-*` scripts listed above, and this file.

Likely overlaps with other lanes: `R/predict.R` (gp kriging),
`R/utils.R`, `NEWS.md`, frmtmb.sample's `R/methods-draws.R` and `R/sample.R`, and the
workflows (any lane that edits a `runs-on:` line).
