# rv9 review findings

## ENVIRONMENT INCIDENT (report to user)
During the mid-session interruption the scratchpad was WIPED and the user's R library
at C:/Users/adf44/AppData/Local/R/win-library/4.6 was DAMAGED: 212 of 350 package
directories are now empty shells (dir present, no DESCRIPTION). Confirmed independently
via Bash and PowerShell. Gutted: RTMB, TMB, RTMBdist, reformulas, testthat, pkgload,
numDeriv, knitr, rmarkdown, brms, Rcpp, Matrix(absent), rstan, lme4, ggplot2, ... .
Intact: 138, incl. roxygen2. I did NOT touch the user library; I rebuilt the whole
closure into the private rv9-lib instead, which shadows the broken copies. The user's
library still needs repair (suggested: reinstall the 212 into their own lib).

## Setup
- main a233c3c, clean, unmoved. Worktrees ddmvar/gddm @ e684eb2, lba @ a233c3c.
- e684eb2..a233c3c touches ONLY docs/, tests/testthat/test-lkj.R, vignettes/case-studies.Rmd.
  `git diff --stat e684eb2 a233c3c -- extensions/frmtmb.ddm` is EMPTY => the two bases are
  identical for this extension. No rebase risk.
- All three worktrees: 0 files changed outside extensions/frmtmb.ddm (tracked and untracked).
- rv9-lib: frmtmb 0.49.1 built from MAIN; rtdists 0.11.6; RWiener 1.3.3; GLMMadaptive 0.9.7;
  RTMB 1.9; all 14 toolchain packages load clean.

## A. PLAIN WIENER PATH -- claim substantiated, wording wrong
Lane wording "a separate, literally unchanged route" is FALSE: R/wiener-density.R
ddm_lpdf_lower IS edited:
  -  u <- t / (a * a)
  +  ur <- t / (a * a); u <- ddm_floor(ur, ddm_u_floor)
The SUBSTANTIVE claim holds. Grid t{1e-8..30} x v{-6..9} x a{0.3..4} x w{0.05..0.95} x up{0,1}
= 11520 points, ddm_lpdf_both() under pkgload::load_all from each tree, dumped as raw doubles:
  positive segment BYTE-FOR-BYTE equal (cmp: first differing byte 92167, i.e. past the
  11520*8=92160-byte positive block); identical() TRUE; max abs diff 0; 11520/11520 equal.
Only difference: t<=0 (360 pts) main NaN -> ddmvar -Inf. Deliberate, documented, an
improvement for the mixture path.
=> gddm's accuracy table is NOT against a moved target.
Caveat: +ddm_u_floor is a rounding no-op only down to ~1e-300; at u<=1e-307 it perturbs.
u=t/a^2 is that small only where the density has already underflowed. Not a defect.

## B. FLOOR IDIOM -- all clean
Measured (lo=1e-300):
  x        ddm_floor   lba_atleast  BAD form  gd_relu
  0.3      0.3         0.3          0.3       0.3
  0        1e-300      1e-300       1e-300    0
  -1e-17   1e-300      1e-300       0  <--    0
  -1e-300  1e-300      1e-300       1e-300    0
  -1e-3    1e-300      1e-300       0  <--    0
The annihilating form is reproduced EXACTLY as the lane described. Both shipped helpers
are safe AND bit-exact passthroughs at 0.3/1/123.456/1e-8.
gd_relu (gddm) is max(x,0) with no floor, but its only call sites are the truncated-power
B-spline basis (gddm.R:52,58) where 0 is the correct value; never passed to log(). OK.
DUPLICATION: ddmvar ships ddm_floor(), lba ships lba_atleast(); same job, two spellings.
Punch-list item.

Pinned negative-span test (test-defects.R:215-244): first assertion HOLDS.
  span = -2.77555756156289135e-17  (< 0 TRUE); t0-st/2 = 0.251 > cap = 0.152 TRUE
  lpdf = -1.61184e9 finite; exp(lpdf) = 0; gradient (-0.7, 0, 0) all finite.
Negative-span FREQUENCY is wider than the lane's "9 to 14 percent":
  ndt=.40 st=.10 -> 5.0/9.4/11.0% over three seeds
  ndt=.35 st=.14 -> 14.6% ; ndt=.45 st=.06 -> 17.0%
  ndt=.30 st=.20 -> 40.1% ; ndt=.50 st=.02 -> 0.0%
Lane's figure is a point estimate on its own configuration, not a bound. Wording only.

+Inf overflow fix VERIFIED by reconstructing the pre-fix single-node reference:
  nodes=21: drift 118 PRE-FIX -1482.09 = FIXED ; drift 121 PRE-FIX +Inf, FIXED -1556
  nodes=41: drift 118 PRE-FIX -1481.63 = FIXED ; drift 121 PRE-FIX +Inf, FIXED -1555.48
  drift 150/200 also +Inf pre-fix, finite after. Threshold matches the lane's "about 120".
Fixed path finite at drift 60..200, nodes 7/21/41, and with sz+st both on.

## C. SCRATCH MERGE -- order ddmvar, gddm, lba
Clone at $SP/rv9-merge (rebuild script: rv9-mkmerge.sh). Lanes replayed as commits by
copying each worktree's extensions/frmtmb.ddm wholesale onto its own base. No worktree
index was touched.
MERGE 1 ddmvar -> main : CLEAN, 0 conflicts, 18 files.
MERGE 2 gddm : 3 conflicts
  DESCRIPTION            Description field only -> union of both clauses.
  man/frmtmb.ddm-package.Rd  same sentence -> union (roxygen regenerates).
  NEWS.md                ddmvar wrote "# frmtmb.ddm 0.2.0", gddm kept
                         "# frmtmb.ddm (development version)" -> folded into ONE 0.2.0.
  NAMESPACE auto-merged.
MERGE 3 lba : 3 conflicts
  DESCRIPTION            two hunks: Description text -> union; Suggests -> take lba's
                         `rtdists` but DROP lba's `RTMB (>= 1.9)`, because merge 2 already
                         moved RTMB to Imports. Keeping both would be a duplicate dep.
  NAMESPACE              disjoint export blocks -> union (roxygen regenerates).
  man/frmtmb.ddm-package.Rd -> union.
  NEWS.md AUTO-MERGED; vignettes/ddm.Rmd (ddmvar + lba) AUTO-MERGED.

MERGED DESCRIPTION:
  Depends: frmtmb (>= 0.49.0)      [main 0.47.0; ddmvar 0.49.0; gddm 0.47.0; lba 0.49.0]
  Imports: RTMB (>= 1.9), stats    [main+ddmvar+lba Suggests; gddm moved it to Imports]
  Suggests: brms, knitr, numDeriv, rmarkdown, rtdists, RWiener, testthat (>= 3.0.0)
  Version: 0.2.0 (all lanes agree, auto-merged)
DEFECT FOUND: gddm alone declares frmtmb (>= 0.47.0) and did NOT raise the floor.
The merge hides it (ddmvar/lba raise it). See punch list.

MERGED NEWS shape:
  # frmtmb.ddm 0.2.0
    ## Across-trial variability in `wiener()`
    ## `gddm()`, the generalized drift-diffusion model
    ## Defects fixed
    ## What frmtmb 0.49.0 let this package delete
  # frmtmb.ddm 0.1.0
  (lba's bullets auto-merged in; see below for where they landed)

## C. MERGE RESULTS (continued)
roxygenise (roxygen2 8.1.0, matching Config/roxygen2/version): produced ZERO changes on the
merged tree. NAMESPACE unchanged by roxygen => my hand union was exactly right.
  NOTE: my first roxygen run produced spurious `\code{}` -> backtick churn because
  `commonmark` (roxygen's markdown engine) was one of the gutted user-lib packages.
  After installing roxygen2's closure into rv9-lib the run is a no-op. Environment
  artifact, NOT a lane defect.

COMBINED SUITE (one process, load_all, NOT_CRAN=true), merged + bracket fix:
  files=134 tests=818 PASS=806 FAIL=0 ERROR=0 WARN=12 SKIP=0  (469s)
  All 12 warnings are test-gddm-recovery.R, exactly the deliberate ones gddm declared:
    10 "the estimator recovers the parameters it was given"
     1 "renormalizing the defective density is not optional"
     1 "a fitted model supports the surface a user reaches for next"
  (Before brms was reinstalled the run showed SKIP=4: 3 test-brms-parity + 1 test-sampling.
   Those were library-damage artifacts, not lane issues; both clear once brms is present.)

R CMD build (with vignettes): OK; both ddm.Rmd (the ddmvar+lba AUTO-MERGE) and gddm.Rmd knit.
R CMD check --no-manual on frmtmb.ddm_0.2.0.tar.gz against installed core 0.49.1: Status: OK
  (0 ERROR, 0 WARNING, 0 NOTE). "checking tests ... OK" - full suite passes inside check too.

## *** BLOCKER FOUND: core bracket lint fails on the merged tree ***
tests/testthat/test-bracket-access.R, second test ("no extension file reads a hazard
container with `$`"), run so test_path() resolves INTO the tree under test:
  FAIL: expected `found` identical to NULL, got
    gddm.R:964: dp$ndt
    gddm.R:969: dp$lapse
    lba.R:463 : fam$links
Attribution by running the scanner per tree:
  main            clean
  wt-ddmvar       clean
  wt-gddm         gddm.R:964 dp$ndt | gddm.R:969 dp$lapse
  wt-lba          lba.R:463 fam$links
So gddm and lba EACH fail core's suite on their own. Neither lane ran this test.
`dp` and `fam` are both on core's hazard_containers list (22 names).
Semantically these three are WRITES (`$<-` sets the exact name, so partial matching does
not actually bite), but the lint is a blanket reserved-name rule and core's suite is red.
FIX (verified in the scratch clone ONLY, 3 lines):
  gddm.R:964  dp$ndt        -> dp[["ndt"]]
  gddm.R:969  dp$lapse      -> dp[["lapse"]]
  lba.R:463   fam$links$ndt -> fam[["links"]][["ndt"]]
After the fix: bracket test 3 tests / 7 assertions PASS; core-boundary 2/4 PASS;
combined suite 818 tests 0 fail 0 error. Behavior-neutral.

METHOD NOTE: my first policing run was wrong - test_file() on MAIN's copy resolves
test_path() relative to the FILE's directory, so it scanned main's extensions and
reported a false PASS. Re-run against the merged clone's byte-identical copies.

## test-core-boundary.R
PASSES on the merged tree (2 tests, 4 assertions, no skips). Core R/ names none of
hmm, lca, mix_g, hmm_g, frm_ode, RTMBode. The three lanes add nothing to core.

## PER-LANE: ddmvar
DV1 sv CLOSED FORM vs my own adaptive quadrature (900 pts: t{.05..3} x vbar{-2..3} x
a{.6..2.5} x w{.3...7} x sv{.1..2.5}):
  naive integrate() of exp(logf)*dnorm gave max rel 1.16e-3 -- that was MY conditioning.
  Re-done with the peak factored out (integrate exp(logf - max logf), rel.tol 1e-13):
    max rel = 1.004e-13, MEDIAN = 0 exactly, 99th pct = 4.7e-15, 1/900 worse than 1e-13.
  Lane claimed 1.8e-15; my grid is more extreme (sv to 2.5) so my tail is looser, but the
  residual is the quadrature's, not the formula's. CLAIM SUBSTANTIATED.
DV2 DEGENERATE IDENTITY: sv=0 vs plain over 225 pts: max abs 3.55e-15, max rel 1.89e-15.
  Full var path with sv=sz=st=0 vs plain: max abs 2.84e-13. (Lane claimed 1.4e-14.) OK.
DV3 *** THE TIMING CLAIM IS REFUTED ***
  Density level, 1000 rows vectorized: plain 0.0002 s/call, sv 0.0002 s/call -- identical.
    var path st(21 nodes) 0.0170 s/call, sz(7 nodes) 0.0056 s/call. Variability costs MORE.
  FIT level, 1000 rows, 3 seeds:
    rep1 plain=1.29 sv=1.34 sz=1.69 st=4.31   (sv/plain 1.04x)
    rep2 plain=0.16 sv=0.16 sz=1.60 st=4.27   (sv/plain 0.97x)
    rep3 plain=0.18 sv=0.20 sz=1.62 st=5.45   (sv/plain 1.07x)
  Fresh processes: plain-first -> plain 2.62 then sv 2.86 ; sv-first -> sv 2.42 then plain 2.41.
  There is NO configuration in which sv is 4x faster than plain. sv/plain = 0.97-1.07x.
  The DEFENSIBLE restatement: the sv closed form is FREE relative to the plain density
  (it adds a sqrt and a few flops), which is the real and notable result.
  Likeliest source of "0.8 vs 3.5": 3.5 s matches the st path (4.3-5.5 s here), so the
  comparison probably mislabeled which family was timed; and there is an 8x in-process
  warm-up swing (plain 1.29 -> 0.16 s between rep1 and rep2) that will manufacture any
  ratio you like if the two arms are not both warm.
DV4 fitted()/simulate() REWEIGHTING -- CONFIRMED against simulation:
  40000 trials, sv=1.5, UNBIASED start point w=0.5:
    empirical mean RT | upper = 0.6528 (n=27819) ; | lower = 0.7282 (n=12181); SEP = 0.0754 s
    same model at sv=0: sep = 0.0030 s  => the separation is purely a variability effect,
    i.e. the plain closed form really would return the same mean for both boundaries.
    fitted() under variability: sep = 0.0726 s vs 0.0816 s empirical on the same 4000 rows.
    simulate() from the fit: mean 0.6795 vs data mean 0.6840.
  Lane claimed "0.06 s apart"; I get 0.075 s at sv=1.5. Direction and magnitude confirmed.
DV5 MOVED ASSERTION (test-surface.R:170) -- JUSTIFIED.
  Old test pinned core's REFUSAL of dec(). ddmvar registers `dec` through
  frmtmb_register_aterm (core 0.49.0), so that refusal no longer exists and the old
  assertion could not pass. The replacement pins a STRONGER property: dec() and vint()
  are one model (logLik equal to 1e-10, fixef to 1e-6). No coverage lost: the
  missing-indicator refusal is still pinned at test-defects.R:158 and test-family.R:113.

## PER-LANE: lba
LB1 single-accumulator density vs rtdists::dlba_norm, 810 pts (792 with a nonzero ref):
  max REL diff = 2.373e-13. (Lane claimed 2.5e-12.) CONFIRMED, mine is tighter.
LB2 race vs rtdists::n1PDF, 300 random points each:
  n=2 max abs diff in LOG = 2.309e-14 ; n=3 = 3.020e-14 ; n=4 = 2.309e-14
  (Lane claimed 1.6e-14.) CONFIRMED, same order.
LB3 UNTRUNCATED SURVIVAL: `s <- if (posdrift) sf/q else sf + (1 - q)` with q = Phi(v/s).
  Verified numerically: sF - sf == 1 - q == Phi(-v/s) at v = -1, 0.5, 2, 5, 8. So YES,
  it ADDS the never-arriving mass rather than dividing. CLAIM CONFIRMED.
  *** NEW DEFECT the lane did not report ***: it is spelled `1 - q`, not Phi(-v/s).
    v/s= 6 : 1-pnorm(v)=9.86588e-10  pnorm(-v)=9.86588e-10  rel err 5.6e-08
    v/s= 8 : 1-pnorm(v)=6.66134e-16  pnorm(-v)=6.22096e-16  rel err 7.1e-02
    v/s=10 : 1-pnorm(v)=0            pnorm(-v)=7.61985e-24  rel err 1.0e+00 (total loss)
  Does it matter? It co-occurs with the tiny-survival rows lsurv was written FOR:
    t=1.5 v/s=8  sf=9.55e-14  survival rel err 4.6e-04
    t=2.0 v/s=8  sf=2.62e-14  survival rel err 1.6e-03
    t=2.5 v/s=10 sf=3.28e-22  survival rel err 2.3e-02
    t=3.0 v/s=10 sf=1.68e-22  survival rel err 4.3e-02
  So up to ~4% error in the untruncated survival. posdrift=FALSE only (a deliberately
  defective model), so severity is low, but the fix is one token:
    lba.R:316   `sf + (1 - q)`  ->  `sf + RTMB::pnorm(-p$v / p$s)`
  This is the SAME defect class the lane documented for lba_phidiff, but in code the lane
  WROTE rather than inherited from rtdists.
LB4 IDENTIFICATION / scale invariance: multiplying (A, b, v, s) by c leaves logf exactly
  unchanged (max|diff| = 0.000e+00 at c=0.5 and c=2; 1.33e-15 at c=7). So the scale is
  unidentified and sd_v MUST be fixed. lba(sd_v=) is the right design. CONFIRMED.

LB5-LB8 THE pnorm-TAIL LIMIT -- the one open question. h = (b - A - v t)/(s t);
lba_phidiff(g, h) loses the v*(Phi(g)-Phi(h)) term once h passes ~8.3.
  VIGNETTE FIT (N=1500, seed 20), at the converged optimum:
    fitted A=0.4167 k=0.4382 b=0.8549 ndt=0.1920 (bound min(rt)=0.2764)
    max h = 4.668 ;  rows with any h>8.3 = 0 / 1500 (0.00%)
    (g>8.3 in 10/1500 rows, but g is the `hi` argument and its saturation is benign.)
  OCCUPANCY AS ndt APPROACHES THE BOUND (this is the regime that exists):
    ndt/bound 0.50 -> 0.00% (max h 3.2)      0.90 -> 1.00% (max h 15.3)
    ndt/bound 0.80 -> 0.00% (max h 7.4)      0.95 -> 2.33% (max h 31.2)
    0.99 -> 4.60% (max h 158)   0.999 -> 5.07% (max h 1585)  0.9999 -> 5.13% (max h 15854)
  DOES THE MISMATCH MOVE THE OPTIMUM?  NO.
    Nelder-Mead on the VALUE surface, started at the gradient optimum, run to convergence
    twice: improvement 1.199e-09, max|dpar| 4.624e-06. The gradient optimum IS the value
    surface's optimum to 6 decimal places in every coordinate.
  BUT the AD-vs-value gradient mismatch is real and grows toward the bound:
    at the optimum          0.36%
    ndt logit +1            1.7%
    ndt logit +2           16.5%   <-- worse than the lane's reported 6%
    +3 / +4 / +5            8.4% / 9.4% / 7.9%
  STRESS TEST (fast drifts, so min(rt) is close to ndt): three configurations reaching
  ndt/min(rt) = 0.866, 0.895, 0.904 -- max h stayed 4.2, 4.5, 4.4 and 0.00% of rows were
  in the regime; NM improvement 2.9e-10, 7.2e-09, 2.6e-09.
  WHY: the likelihood itself repels ndt from min(rt) (the fastest row's density collapses),
  so the optimum is structurally kept out of the defective regime.
  CONCLUSION: the limit is real and correctly described, but on every dataset I could
  construct it does NOT move the point estimate. The residual exposure is (a) predict()
  on new data with an rt below the training minimum, (b) a user pinning ndt manually,
  and (c) Hessian-based standard errors, which come from the derivative of the WRITTEN
  function and so inherit the 0.36% mismatch even at the optimum.

## *** LBA: THE LANE'S REASON FOR PINNING THE TAIL DEFECT IS WRONG ***
The lane argues (R/lba.R:341-355) that the defect cannot be fixed because the two
candidate forms fail in mirror-image regimes and "neither can be selected per
evaluation: RTMB refuses comparison on AD types, so there is no tape-safe conditional".
That is a FALSE DICHOTOMY. There is a third form that needs no selection:

  logsp(hi, lo):  la <- RTMB::pnorm(-lo, log.p = TRUE)
                  lb <- RTMB::pnorm(-hi, log.p = TRUE)
                  exp(la) * -expm1(lb - la)

It is algebraically Phi(-lo) - Phi(-hi) = Phi(hi) - Phi(lo), but every step keeps full
RELATIVE precision: pnorm(., log.p=TRUE) does not saturate, and expm1 removes the
cancellation. FEASIBILITY VERIFIED:
  - RTMB::pnorm has a log.p argument and it WORKS ON THE TAPE (MakeTape round-trip OK).
  - expm1 works on the tape.
  - AD derivatives through it are correct at hi/lo = (1,-1), (9,8.4), (12.5,11), (-8.4,-9):
    max rel err vs the analytic (dnorm(hi), -dnorm(lo)) is 7.6e-10, mostly ~1e-15.
  - NO comparison, NO TapeConfig, NO smooth blend needed.

ACCURACY, against a 200-bit Rmpfr reference (erfc form), on the lane's own density grid
(810 rows; 756 with a positive true value):
  BULK  (h < 8, n=696)   current max rel 7.246e-03, 15 rows worse than 1e-6
                         log-space max rel 1.019e-14, 0 rows worse than 1e-6
  TAIL  (h >= 8, n=60)   current max rel 1.000, 57 rows return EXACTLY 0
                         log-space max rel 3.571e-04, 0 rows return exactly 0
The log-space form is strictly better EVERYWHERE, by ~12 orders of magnitude in the bulk.
Note also that the 8.3 threshold UNDERSTATES the problem: the current form is already
0.7 percent wrong on some h < 8 rows.
In the density itself, at A=0.5 k=0.4 v=1.5 s=1 and t=0.040 (h=8.50), the current form
understates by 14.8 percent -- the lane's "13 percent at h = 11" reproduced in kind.

RECOMMENDATION: the log-space upper-tail reformulation is the fix. TapeConfig(comparison
= "tape") is NOT needed here and should not be adopted for this. The one real cost is
that the lane's "agrees with rtdists to 1e-11 everywhere" test becomes false by design,
because it would now be more accurate than rtdists in exactly the rows rtdists gets
wrong; that test should be restricted to h < 8 and a new Rmpfr-referenced test added for
the tail. Caveat I could not remove: 3 of 756 rows still exceed 1e-6 (max 3.6e-4) under
log-space, so it is a large improvement rather than a proof of correctness to 1e-13.
