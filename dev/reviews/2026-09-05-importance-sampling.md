# Review: `frm(importance =)` (lane wt-islap, base a233c3c)

Reviewer notes, written incrementally (so the numbered sections are in
the order they were finished, not in brief order: 2-6, then 1, then
7-11, then 12 and the verdict). Every number below was re-measured in
this review unless it says otherwise.

**VERDICT: GO WITH FIXES.** Jump to "VERDICT" and "Punch list" near the
end for the seven items; P1 (three wrong `frm_compat()` rows), P2
(`confint(method = "profile")` unguarded on a corrected fit) and P3
(every round re-optimizes from a cold start) are the substantive ones.
The reported suite discrepancy is a LIBRARY ARTIFACT, not a disabled
test: section 1.

## Review environment
- Private library `<scratchpad>/ri-lib`, worktree core installed with
  `R CMD INSTALL --library=ri-lib`. `R_LIBS` = ri-lib; user library;
  system library. `R_LIBS_USER` = ri-lib only.
- The SHARED user library `C:/Users/adf44/AppData/Local/R/win-library/4.6`
  is HEALTHY as of this review: all 40 Suggests present and loadable
  (brms 2.23.0, GLMMadaptive 0.9.7, lme4 2.0.6, numDeriv 2016.8.1.1,
  rstan 2.32.7, tmbstan 1.2.0, glmmTMB 1.1.14, ...). Nothing was
  missing, so nothing had to be installed beyond the core.
- Main checkout `C:/Users/adf44/source/r/frmtmb` verified at a233c3c,
  clean apart from one pre-existing untracked file
  `dev/brms-likelihood-tests.md` that is not this lane's and was not
  touched.

## 0. Surface
`git diff --name-only a233c3c` in the worktree:
NEWS.md, R/allfit.R, R/compat.R, R/confint.R, R/fit.R, R/influence.R,
R/methods-fit.R, R/objective.R, R/sandwich.R, R/sugar.R, man/frm.Rd,
man/frmtmb_control.Rd, vignettes/diagnostics.Rmd.
Untracked: R/importance.R, tests/testthat/test-importance.R,
dev/importance-findings.md, dev/laplace-bias-probe.R.
Nothing outside the declared surface. DESCRIPTION Version untouched at
0.49.1 as stated.

## 2. The gaussian identity: REPRODUCED
Design: 50 groups of 6, `y ~ x + (1|g)`, gaussian, seed 5.
Laplace objective at its optimum 470.102104039837.

| N | corrected objective at the anchor | abs diff | min ESS | MCSE |
|---|---|---|---|---|
| 2 | 470.102104039837 | 0.000e+00 | 1.0000000000 | 0.000e+00 |
| 10 | 470.102104039837 | 0.000e+00 | 1.0000000000 | 6.45e-09 |
| 1000 | 470.102104039837 | 0.000e+00 | 1.0000000000 | 0.000e+00 |

Away from the anchor on the SAME frozen tape (N = 200), the identity
goes as advertised: theta shift 0.05 gives IS - Laplace = -3.33e-02,
shift 0.20 gives -1.46e-01.

WHY it is exact. For a gaussian response the joint negative log-density
is exactly quadratic in `u`, so the Laplace Gaussian proposal IS the
conditional. Every log weight `-f_g(theta, u_i) + |z_i|^2/2` is then the
same constant for every draw i, the log-mean-exp collapses to that
constant, and the estimator returns the Laplace value algebraically,
for any N and any seed. That cancellation only holds where the proposal
matches the conditional, i.e. at the anchor; move theta and the
conditional moves while the frozen proposal does not, so the weights
spread and the estimate becomes an ordinary noisy Monte Carlo one.

WHAT THE REPORTED logLik MEANS, and the tape question. Read
`importance_fit()` (R/importance.R:668-712): each round optimizes the
current tape, sets `par <- opt$par`, then ALWAYS refreezes at that
`par` before testing the break conditions, and after the loop
`opt$objective <- fz$obj$fn(par)`. So the reported logLik is the
corrected objective on the FINAL tape, anchored at the reported theta,
evaluated at that theta. `logLik()` reads `-fit$opt$objective`
(R/methods-fit.R:235) and `sdreport()` runs on `fit$obj`, which is
`imp$obj`, the same final tape, at `fit$opt$par` (R/autoscale.R:143).
So logLik, vcov and confint(wald) all come from ONE tape - good.
The theta, however, is the argmin of the PREVIOUS round's tape, not of
the reporting tape. It is not a stationary point of the tape that
reports it: on the probe at N = 2000 the reporting tape's gradient
there is 7.94e-04. This is the deliberate design (the anchor is the
only place the Monte Carlo error is not optimistic) and the residual is
stored in `fit$importance$grad`, so I do not call it a defect. It is
worth one sentence in the docs that the reported covariance is a
Hessian at a point that is stationary to `grad` but not exactly.

## 3. The probe: REPRODUCED to the digit
`frm(bf(y ~ x + (x|g)) + bernoulli(), importance = 2000)`, 60 groups of
8, seed 11.

| fit | nll | intercept | slope |
|---|---|---|---|
| Laplace | 299.3017 | -0.6915 | 0.3592 |
| importance 2000 | 297.9147 | -0.7073 | 0.3648 |
| GLMMadaptive nAGQ 15 / 21 / 25 | 297.9185 | -0.7047 | 0.3617 |

draws 2000, rounds 3, capped FALSE, grad 7.94e-04, MCSE 0.0308,
min ESS 0.9312, median ESS 0.9701, moves 0.14938 / 0.01914 / 0.00127.
SE(intercept) 0.17918.
|297.9147 - 297.9185| = 0.0038 = 0.12 MCSE. Intercept off by 0.0026 =
0.015 SE. Laplace misses the nll by 1.387. nAGQ 21 gives the same value
as 15 and 25 to four decimals, so the reference is converged.
Wall clock 80.5 s against 0.4 s for the Laplace fit (200x).

## 4. The divergence refusal: REPRODUCED, guard and all
Design (from the test at tests/testthat/test-importance.R:426): 40
groups of 6 binary rows, seed 11, `y ~ x + (x|g)`, importance = 500.

WITHOUT the guard (scratch copy at `<scratchpad>/ri-noguard`, the
single `if (is.finite(worse) && ...)` at R/importance.R:726 replaced by
`if (FALSE)`; the worktree was NOT touched):
- Laplace: nll 153.0853, `singular convergence (7)`, SD(x) 0.0127,
  correlation -1.
- corrected: nll 156.0558, i.e. WORSE than Laplace by 2.97, after 5
  rounds. Rise over the corrected value at the Laplace estimates:
  +3.2267.
- both variance components collapsed: SD 5.48e-05 and 6.88e-06, theta
  -9.81 / -11.89 / -2.44.
- diagnostics blind exactly as claimed: min ESS 1.000000, median ESS
  1.000000, MCSE 0.000000. An anchor-ESS guard cannot see this.
- and the ESTIMATOR is fine there: the corrected objective AT the
  Laplace estimates is 152.8291 against GLMMadaptive nAGQ = 25's
  152.8506, better than Laplace's 153.0853. It is the fixed-point
  iteration that diverges, not the estimator. The lane's diagnosis is
  correct.

WITH the guard (the worktree as shipped), the same call refuses with
the named message:
"`importance` did not converge on this model: the corrected negative
log-likelihood ROSE from 152.82912 at the Laplace estimates to
156.05582 after 5 rounds, so the iteration moved away from the answer
instead of toward it. ..."

The nine-design separation, all nine re-measured through the guard-free
build so every one reports its actual rise:

| design | rise | 3*MCSE vs 0.1 floor | verdict |
|---|---|---|---|
| 40x4 seed11 (unidentified) | +0.3523 | 0.1000 | REFUSE |
| 40x6 seed11 (unidentified) | +3.2267 | 0.1000 | REFUSE |
| 50x6 seed3 | -0.0095 | 0.1000 | fine |
| 60x8 seed11 | -0.0709 | 0.1889 | fine |
| 40x5 seed7 wide | -0.0281 | 0.3254 | fine |
| 50x5 seed2 wide | -0.1279 | 0.2119 | fine |
| gaussian N=50 | +0.0103 | 0.1000 | fine |
| gaussian N=200 | +0.0069 | 0.1000 | fine |
| gaussian N=1000 | +0.0000 | 0.1000 | fine |

Worst legitimate rise +0.0103, smallest divergence +0.3523, floor 0.10:
9.7x above the one and 3.4x below the other. The lane's table
reproduces to the third decimal. The threshold is placed by
measurement, not by taste.

## 5. The anchor-versus-optimizer bias: REPRODUCED
Gaussian design, proposal frozen at the exact Laplace optimum
(470.10210404), then optimized ONCE on that frozen tape:

| N | optimizer's value | below the Laplace optimum | exact value there, above it |
|---|---|---|---|
| 2 | 469.904862 | -1.9724e-01 | +1.6836e-01 |
| 10 | 470.088938 | -1.3166e-02 | +1.2409e-02 |
| 200 | 470.095123 | -6.9807e-03 | +7.0622e-03 |

The lane's 1.97e-1 / 1.32e-2 / 6.98e-3 are exact. The mechanism is
confirmed too: the optimizer's value is optimistic by very nearly the
amount by which the exact objective at the point it stopped is
pessimistic, which is the signature of exploiting an error that is
zero at the anchor. "Report from the anchor" is the right fix.

## 6. The ESS threshold: partly reproduced, one number I could NOT pin
Displaced proposal (probe, both log-SDs moved by -0.5, weights read at
the true optimum):

| N | min ESS at the optimum | min ESS displaced | MCSE displaced |
|---|---|---|---|
| 1000 | 0.9467 | 0.0325 | 0.533 |
| 2000 | 0.9553 | 0.0066 | 0.552 |

The lane reports 0.0089 at N = 2000; I measure 0.0066. Same order,
same conclusion (30-40x below the 0.25 floor), but not the same digit.

The "hardest legitimate design, 40 groups of 3 binary rows, sd 2.0,
min ESS 0.429" I could NOT reproduce: the lane's measurement script was
lost in the scratchpad wipe and the findings do not record the seed.
Rebuilding the design as described and sweeping seven seeds at
N = 1000 gives min ESS 0.9703 / 0.7331 / 0.7732 / 0.8396 / 0.4802 /
0.9433 / 0.9434 (seeds 1, 2, 3, 7, 11, 101, 601). The hardest I found
is 0.4802.
This does not change the verdict on the threshold - the floor of 0.25
still clears every legitimate design I measured by at least 1.9x, and
the vignette design's 0.6426 reproduces exactly - but the specific
figure 0.429 and the "1.7x margin" sentence in the `imp_ess_floor`
comment (R/importance.R:775) rest on a measurement no surviving
artifact pins. PUNCH LIST: either re-measure and record the seed, or
soften the comment to the margin that is reproducible.

## 1. THE SUITE COUNT: RECONCILED. Not a blocker.
Full core suite, one file per process, `pkgload::load_all()`,
`NOT_CRAN=true`, ri-lib with every Suggests available.

**102 files, failed = 0, errors = 0, warnings = 0, skipped = 4,
passed = 5347.** No file failed to run.

The four skips are both opt-in tiers, not capability gaps:
- test-brms-agreement.R x3: "set FRMTMB_BRMS_FIT_TESTS=true to run
  brms fit tests"
- test-fuzz.R x1: "set FRMTMB_FUZZ=true to run the grammar fuzz tier"

This matches the ~5280 / 4-to-6-skips baseline the other runs report
(5347 - 89 new test-importance.R assertions = 5258).

THE GAP IS MISSING SUGGESTS IN THE LANE'S PRIVATE LIBRARY, NOT TESTS
THE CHANGE DISABLED. Per-file, every one of the lane's 39 extra skips
sits in a file that GAINED passes in my run, and no file lost a single
pass:

| file | lane | mine |
|---|---|---|
| test-brms-agreement.R | s=20 p=0 | s=3 p=168 |
| test-api-spellings.R | s=3 p=15 | s=0 p=22 |
| test-brms-port.R | s=3 p=7 | s=0 p=18 |
| test-mvn-mixture.R | s=1 p=165 | s=0 p=277 |
| test-lkj.R | s=1 p=152 | s=0 p=212 |
| test-mm.R | s=2 p=98 | s=0 p=116 |
| test-prior-compat.R | s=2 p=93 | s=0 p=105 |
| test-ordinal-fitted.R | s=2 p=81 | s=0 p=86 |
| test-review-v29.R | s=1 p=122 | s=0 p=136 |
| test-sandwich.R | s=1 p=50 | s=0 p=57 |
| test-pooled-anova.R | s=1 p=59 | s=0 p=64 |
| test-famgaps.R | s=1 p=94 | s=0 p=99 |
| test-portability.R | s=1 p=87 | s=0 p=91 |
| test-multiple-pooling.R | s=1 p=29 | s=0 p=33 |
| test-case-studies.R | s=1 p=25 | s=0 p=30 |
| test-effects.R | s=1 p=48 | s=0 p=50 |
| test-importance.R | s=0 p=84 | s=0 p=89 |

`test-importance.R` runs 89 assertions, 0 skipped, 0 failed. NO TEST IS
SILENTLY DISABLED BY THIS CHANGE. The lane's own suite figure is an
artifact of its `is-lib` and should not be quoted; it should be re-run
against a complete library before the lane's findings are archived.

## 7. trunc() and cens(): supported, and HARDER than the lane tested
Reproduced (40 groups of 6, seed 21):
- left-truncated Poisson, `trunc(lb = 1)`, n = 211: Laplace 448.8567,
  importance 1000 gives 448.7131, fixef 1.18187 / 0.39859, min ESS
  0.9353, median 0.9727, MCSE 0.0329, 2 rounds. `quadrature = TRUE`
  refuses the same model by name. Matches the lane exactly.
- right-censored lognormal, `cens(cen)`: min ESS 0.9753, median 1.0000,
  MCSE 0.0103, finite objective, estimates within 1e-3 of Laplace.
  NOTE: my run has 52 of 240 censored and a Laplace nll of 311.1358,
  where the lane's findings say 40 of 240 and 347.8633. The TEST's
  construction (test-importance.R:368) reproduces MY numbers, so the
  findings quote a different, lost scratch script. Cosmetic, but the
  findings' cens() row should be re-transcribed from the test.

NEW CASE the lane did not run - a truncation bound walked up to and
past the response median. Untruncated Poisson median 3.

| lb | n kept | median | pct at the bound | min ESS | median ESS | MCSE |
|---|---|---|---|---|---|---|
| 1 | 211 | 3 | 17.1 | 0.9353 | 0.9727 | 0.0329 |
| 2 | 175 | 4 | 21.1 | 0.9048 | 0.9636 | 0.0391 |
| 3 | 138 | 5.5 | 24.6 | 0.8356 | 0.9597 | 0.0461 |
| 4 | 104 | 6 | 16.3 | 0.5956 | 0.9573 | 0.0454 |
| 5 | 87 | 7 | 20.7 | 0.6398 | 0.9669 | 0.0404 |

And a gaussian RIGHT truncation at the median: `ub` = median + 0.25
keeps 61 percent of the rows and still gives min ESS 0.9723, MCSE
0.0185.

IT DOES NOT DEGENERATE. Even discarding 64 percent of the data at a
bound above the median, the worst group holds 0.60 of its draws, more
than twice the 0.25 floor, and every fit converged in 2-3 rounds
without being capped. The documentation does NOT need a warning here;
the lane's trunc()/cens() compat rows are honest and if anything
understate the margin.

## 8. imp_verify(): FIRES
Read at R/importance.R:545-573. Scratch copy `<scratchpad>/ri-break`
with one weight perturbed inside `amat` (the worktree untouched):

| perturbation on one weight | imp_verify | frm(importance =) |
|---|---|---|
| 0 | silent | fit succeeds |
| 1e-8 | silent | fit succeeds |
| 1e-3 | FIRES | fit REFUSED |

The 1e-3 message: "`importance` could not reproduce this model's joint
log-density from its per-group pieces (draw 1: -358.9903629 against
-358.9913629). ..." and the refusal propagates out of `frm()`. The
silence at 1e-8 is correct: the tolerance is `1e-6 * max(1, |joint|)`,
about 3.6e-4 here.

ONE REAL WEAKNESS, from reading it: `ours <- sum(a[, j] -
plan[["zsq"]][, j])` compares only the TOTAL over groups. An error that
adds eps to one group and subtracts it from another passes the pin
untouched - and a per-group error is precisely what would corrupt the
ESS diagnostics, which are read off the same matrix per row. Comparing
the per-group vector against a per-group joint density would close
this.

MEASURED, not conjectured. With the same scratch copy, adding +eps to
group 1 and -eps to group 2 (leaving the total untouched):

| injected error | imp_verify |
|---|---|
| +1e-3 on one group | FIRES |
| +1e-3 / -1e-3 on two groups | SILENT |
| +5 / -5 on two groups | SILENT |

An error of FIVE log-density units per group passes the pin without a
word. Low severity today (no such failure mode is known, and the ESS
would look odd), but the pin is weaker than its own comment claims, and
per-group correctness is exactly what the ESS diagnostics rest on.

## 9. imp_prior_terms(): EXACT on all 13 whitelisted structures
For each structure, a one-level block at a random theta, probing the
extracted closure at 0, e_j and e_j+e_k, then comparing the recovered
precision against `solve(covstruct_registry[[k]]$vcov(theta, blk))` and
the recovered normalizer against -0.5(q log 2pi + log|V|). Block
dimension 3 (dimension 1 also checked).

| structure | max abs precision error | normalizer error | in scope |
|---|---|---|---|
| us | 4.72e-16 | 0 | yes |
| diag | 8.88e-16 | 0 | yes |
| homdiag | 0 | 4.44e-16 | yes |
| cs | 4.44e-16 | 0 | yes |
| homcs | 8.95e-16 | 8.88e-16 | yes |
| ar1 | 4.44e-16 | 0 | yes |
| hetar1 | 4.44e-16 | 0 | yes |
| toep | 3.91e-14 | 8.88e-16 | yes |
| homtoep | 1.33e-15 | 4.44e-16 | yes |
| ou | 8.88e-16 | 4.44e-16 | yes |
| exp | 5.27e-16 | 0 | yes |
| gau | 5.00e-16 | 0 | yes |
| mat | 1.09e-11 | 2.44e-14 | yes |
| us (dim 1) | 3.33e-16 | 2.22e-16 | yes |
| us_t (student) | **3.55e-01** | 0 | NO - refused |
| gr(g, cov = A) | closure errors | - | NO - refused |

Every structure the whitelist admits is recovered to machine precision
(`mat`'s 1e-11 is Bessel-function roundoff). The two refusals are
LOAD-BEARING, not cosmetic: on a Student-t latent block the extraction
silently returns a gaussian precision that is wrong by 0.355, which is
exactly the silent-wrong-answer the whitelist exists to stop. The
lane's reasoning for refusing us_t/diag_t is confirmed by measurement.

## 10. The lp_eta_fixed() refactor: BIT-IDENTICAL
Pristine a233c3c extracted with `git archive` to
`<scratchpad>/ri-pristine`; the same objectives evaluated at the same
template under both, printed to 17 significant digits:

| model | pristine a233c3c | worktree | equal |
|---|---|---|---|
| autoscale bad-scale poisson GLMM | 1140.2489840683327 | 1140.2489840683327 | yes |
| autoscale FIT objective | 900.36652944182924 | 900.36652944182924 | yes |
| mo() | 322.97155822523359 | 322.97155822523359 | yes |
| offset() | 1070.2962195667228 | 1070.2962195667228 | yes |
| mi() | 426.01893826078947 | 426.01893826078947 | yes |
| gaussian GLMM | 565.49816018478828 | 565.49816018478828 | yes |

`diff` of the two outputs is empty. In the suite above, test-autoscale.R
46/46, test-mo.R, test-mi.R and every offset test pass with 0 failures.
The refactor preserves the ORDER of the additions (X beta, then Z b,
then offset, then mo(), then mi()), which is why it is bit-identical
and not merely close.

ADoverload sweep of R/: a mechanical scan for a literal-first `c(0, ...)`
or `c(1, ...)` in any top-level function that does NOT bind
`ADoverload` found NOTHING outside the functions that already bind it
(R/ad-env.R, R/autocor.R, R/covstruct.R, R/families.R, R/importance.R,
R/objective.R). No second instance of the bug.
Minor: `frmtmb_ad_overload()` (R/ad-env.R) documents THREE bindings
("c", "[<-", "diag<-") as the project convention; `lp_eta_fixed()` binds
two. `diag<-` is not used there so nothing is wrong, but the package
already owns a helper for exactly this and the hand-rolled spelling is
what went wrong once already.

## 11. The API judgement calls

### The four "refit paths" - one is misdescribed, and it matters
The lane's findings say the correction is threaded through "confint
profile, influence, frm_allfit, refit". Three of the four are right.
The fourth is not:

- `refit()` (R/sugar.R:240): VERIFIED. `refit(fi, newresp = dd$y)`
  returns `importance$draws = 100` and an objective 2.84e-14 from the
  original. Genuine.
- `influence()` (R/influence.R:89): VERIFIED, and provably genuine, not
  just non-NULL - `influence()` on the Laplace fit takes 1.9-3.0 s
  and on the importance fit 62-79 s, a 21-42x ratio that only 40
  real importance refits can produce. 40 refits, all converged.
- `frm_allfit()` (R/allfit.R:94): VERIFIED. Four optimizers all return
  logLik -100.3876, which is the IMPORTANCE value; the Laplace fit on
  the same data is -101.0184. A dropped correction would have shown the
  Laplace number. logLik spread 5.12e-06.
- `confint(method = "profile")`: **THE CHANGED LINE IS NOT THIS PATH.**
  R/confint.R:1190 sits inside `anova_refit_ml()`, which is called only
  from `anova(..., refit = TRUE)` on REML fits (R/confint.R:1305). And
  `importance` is REFUSED with `REML = TRUE`, so an importance fit can
  never reach `anova_refit_ml()`: that line is unreachable in practice.
  The real `confint(method = "profile")` (R/confint.R:436) does not
  refit at all - it calls `TMB::tmbprofile(object$obj, ...)` on the
  fit's own tape, which for an importance fit is the FINAL FROZEN
  proposal.

  This is not automatically wrong, but nothing warns and nobody
  measured it, so I did. Proposal frozen at the reported estimate, then
  the ESS read as a parameter is walked away from it (500 draws, 40
  groups of 4 binary rows):

  | parameter moved to | min ESS | median ESS | MCSE |
  |---|---|---|---|
  | x = 0.7014 (the estimate) | 0.9279 | 0.9881 | 0.0368 |
  | x = 0.2556 (profile lower) | 0.8507 | 0.9717 | 0.0620 |
  | x = 1.2758 (profile upper) | 0.7262 | 0.9517 | 0.0791 |
  | x = 2.7014 | 0.0587 | 0.7279 | 0.3483 |
  | theta = 0.2951 (the estimate) | 0.9279 | 0.9881 | 0.0368 |
  | theta = -0.2049 | 0.6024 | 0.8202 | 0.1457 |
  | theta = +0.7951 | 0.0921 | 0.5835 | 0.3461 |

  VERDICT: the FIXED-EFFECT profile is safe here - over its own
  interval the worst group still holds 0.73 of its draws. The
  COVARIANCE-parameter profile is not clearly safe: `confint(fi, parm =
  "theta_1", method = "profile")` returns [-0.4149, 0.7243], walking
  theta by +0.43 from the estimate, and at +0.5 the minimum ESS is
  0.0921, well below the 0.25 floor. The fit's reported ESS (0.9279)
  describes the anchor and says nothing about that path.
  So `confint(method = "profile")` on an importance fit can compute a
  covariance-parameter bound in the degenerate regime and report it
  without a word. This needs either a guard, or a documented warning,
  or at minimum a corrected sentence in the lane's findings. It is the
  most substantive gap I found.

### vcov_cluster(): refusal correct
`vcov_cluster(fi, dd$g)` refuses with the intended text: "vcov_cluster()
cannot use a fit made with an importance correction: that objective is
a sum over GROUPS of reweighted integrals, not the per-row sum the
cluster scores are read off, so a cluster carrying no rows of its own
still moves it. Refit with importance = 0 (Laplace)". Matches the
existing quadrature and prior guards in style and placement
(R/sandwich.R:172).

### print() and summary(): correct
Both emit exactly one line, in the right place (after logLik/AIC), and
a Laplace fit emits nothing:
"Marginal likelihood: importance-corrected, 100 draws per group in 3
rounds (MCSE 0.11, min ESS 0.68 of 1)"

### frm_compat("importance"): THREE WRONG ROWS
`frm_compat()` is the documented answer to "what does this feature do
with that one", and the new vignette section sends users to it
explicitly ("`frm_compat("importance")` lists the rules"). Three rows
contradict what `frm()` actually does. Verified by building a model for
each structure and calling `frm(..., importance = 10)`:

| structure | frm_compat says | frm() actually does |
|---|---|---|
| rr | **conditional** | SCOPE-REFUSED by name |
| us_t | **untested** | SCOPE-REFUSED by name |
| diag_t | **untested** | SCOPE-REFUSED by name |
| gr_cov | refused | refused (correct) |
| us, diag, homdiag, cs, homcs, ar1, hetar1, toep, homtoep, ou, exp, gau, mat | conditional | accepted (correct) |

Cause: R/compat.R:245-253 adds `importance_blocks` and
`crosslevel_blocks`, but `rr`, `us_t` and `diag_t` are in NEITHER, so
they fall through to the generic "estimation mode" default and come
back as works/conditional/untested. The lane's own DELIBERATE OMISSIONS
list says us_t and diag_t are "refused with the reason", so the table
contradicts the lane's stated intent, not just the code. And item 9
above shows the us_t refusal is load-bearing: the precision extraction
is wrong by 0.355 there, so a user who trusts the "untested" row and
tries it would get a refusal, which is the good outcome - but a user
who reads "conditional" for `rr` is told a supported thing is supported
when it is not.
FIX: one more group (every covstruct not in `importance_blocks` and not
in `crosslevel_blocks`) with status "refused" and the whitelist reason.
I did NOT make this edit: it changes shipped documentation and the
compat coverage assertions, which is the lane's call.

### frm.Rd and frmtmb_control.Rd: match the roxygen, and are accurate
`man/frm.Rd` carries `importance = 0L` in \usage in the right position
and a full \item; `man/frmtmb_control.Rd` carries all three control
arguments. Both are regenerated, not hand-edited. The prose is accurate
against the code with one exception noted in the punch list (the 0.43
figure).

### vignettes/diagnostics.Rmd: renders, and its numbers are EXACT
Re-ran the section's chunks fresh:
- `c(laplace = -119.5584, importance = -118.4392)` - exact match.
- print line: "importance-corrected, 500 draws per group in 4 rounds
  (MCSE 0.11, min ESS 0.64 of 1)" - exact match.
- `summary(imp$importance$ess)`: Min 0.6426, Median 0.8973, Max 0.9924
  - exact match. `mcse` 0.1084579 - exact match.
The section is placed after "Are the standard errors trustworthy?", and
that list correctly went from three remedies to four.
brms-migration.Rmd is untouched, as the coordinator asked.

## A DEFECT NOT IN THE BRIEF: every round re-optimizes from a COLD start
`importance_fit()` calls `optimize_obj(fz$obj, control, bounds,
par_units, verbose = vb)` (R/importance.R:672) and never passes
`start_par`. `optimize_obj()` defaults `start_par = obj$par`
(R/fit.R:1675), and the freeze tape is built as
`RTMB::MakeADFun(io$fn, otpl, ...)` where `otpl` is the START TEMPLATE.
So the optimizer begins every round at the template start, not at the
current estimate.

Measured on the probe design: the freeze tape's `obj$par` is
(-0.5286, 0, 0, 0, 0) while the Laplace optimum is
(-0.6915, 0.3592, -0.0801, -0.4879, -0.0853) - a maximum difference of
0.4879. The Laplace warm start is used for the ANCHOR and thrown away
for the OPTIMIZER.

Consequences, in order of importance:
1. Cost. Every round pays a full cold optimization on the most
   expensive tape in the package (the N = 2000 gradient is ~284 ms).
   The probe fit takes 80.5 s against 0.4 s for Laplace; a warm start
   should remove a large part of the per-round iteration count.
2. Robustness. The fixed-point iteration is "anchor at theta_k,
   minimize from the template start". On a ridge - exactly the
   unidentified-covariance case the divergence guard exists to catch -
   a cold start each round is free to land in a different basin than
   the anchor. I cannot prove it causes the divergence, but it is a
   plausible contributor and it is the opposite of what the design
   comment describes ("refreezing is what keeps the weights from
   degenerating").
FIX: pass `start_par = par`. One argument. It WILL change the reported
numbers, so it is the lane's to make and re-measure, not mine.

## Determinism, RNG hygiene and argument validation: all hold
- Two fits at the same seed: `identical()` on both objective and par.
- `set.seed(99); rnorm(3)` before and after a fit: identical. The
  private stream neither reads nor disturbs the session state.
- `importance_seed = 2` moves the answer (100.387632 -> 100.409020),
  so the seed is genuinely wired to the draws.
- An odd count rounds UP: `importance = 51` gives `draws = 52`.
- Every refusal is in the package's own words:
  `importance_ess = 1.5`, `importance_rounds = 0`, `importance = -1`,
  `importance = 2.5`, `+ quadrature = TRUE`, `+ REML = TRUE` all
  refuse with a named message.
- The capped warning fires correctly: at seed 2 the loop used all 5
  rounds still moving by 0.0536 and said so. Worth knowing that the
  default round cap of 5 is NOT always enough at 100 draws; the
  warning is the safety net and it works.

## The three "untested" compat rows are actually fine - measured
Same design, `importance = 200`, against the plain corrected fit
(logLik -100.390971, min ESS 0.7972):

| combination | logLik | importance record | min ESS |
|---|---|---|---|
| autoscale = TRUE | -100.390971 | present | 0.7972 |
| sparse_x = TRUE | -100.390971 | present | 0.7972 |
| autoscale engaged (x scaled by 1e6) | -100.390971 | draws 200 | 0.7972 |
| mo() + importance | -118.708630 | present | 0.8925 |

All four work and preserve the correction. The compat rows for
`autoscale`, `sparse_x` and `mo()` can be upgraded from "untested" to
"works" citing these. (Note the autoscale pre-fit at R/autoscale.R:116
deliberately does NOT thread `importance` - correct, that pre-fit only
derives a scaling and a Laplace pre-fit is the cheap right answer.)

## The several-blocks-over-one-factor follow-up: the description is RIGHT
Checked against the code. What the lane says generalizes unchanged does:
`build_importance_objective()`'s `zu` multiplies the WHOLE `lp[["Z"]]`
by the whole draw matrix (R/importance.R:471), so extra blocks need no
change there provided `plan$u` is assembled in `b` order; the
block-diagonality gate is written against a per-position level vector
(R/importance.R:308) and concatenates; the precision extraction is
per block and becomes a sum.

Its list of what must change is correct and complete in substance. One
refinement: the contiguity assumption `(k - 1) * q + seq_len(q)` is not
in ONE place but FOUR, and the write-up names only the first two
explicitly -
  R/importance.R:308 `lev <- rep(seq_len(ng), each = q)` (the gate)
  R/importance.R:344 `ii <- (k - 1L) * q + seq_len(q)` (Hessian slice)
  R/importance.R:364 `z[(seq_len(ng) - 1L) * q + j, ]` (zsq)
  R/importance.R:477 `plan[["u"]][(seq_len(ng) - 1L) * q + j, ]` (u_slices)
All four need the same per-group index vector. The lane's phrase "the
same index has to drive the Hessian slice and the draw placement"
covers 344 and 477, and 308 is covered by its own bullet, but 364 is
not named. Worth adding so the follow-up does not miss it.
Its conclusion that blocks over DIFFERENT factors stay refused (the
nested integral does not factorize) is correct and should not be
softened.

## Version
Nothing here argues against 0.50.0. This is a new user-facing argument
plus new control arguments, entirely additive: `frm()` gains
`importance = 0L` in a position every internal call site already passes
by name past, `frmtmb_control()` gains three arguments with defaults,
and `vcov_cluster()` gains one refusal that can only trigger on a fit
that could not previously exist. No existing behavior changes - proved
bit-for-bit in item 10. A minor bump is right, and sharing it with the
`get_prior()` lane is fine. DESCRIPTION is still at 0.49.1 and the
integrator must bump it; NEWS.md already opens the 0.50.0 section.

## The ESS warning is NOT dead code - it fires end to end
A worry worth ruling out: the reported ESS is the FINAL ANCHOR's, and
the lane's own analysis says the anchor is where the weights are as
equal as the proposal can make them, so one might expect the 0.25
warning never to fire on a real fit. It does. Pushing the design past
the point where a binary GLMM is identifiable at all (500 draws,
seed 11):

| design | min ESS | median ESS | MCSE | warning fired |
|---|---|---|---|---|
| 40 groups x 3 rows, sd 2.0 | 0.7308 | 0.9217 | 0.0964 | no |
| 40 x 2, sd 3.0 | 0.3175 | 0.8674 | 0.1368 | no |
| 60 x 2, sd 4.0 | 0.0833 | 0.8537 | 0.2233 | YES |
| 30 x 2, sd 5.0 | 0.1258 | 0.7067 | 0.2165 | YES |
| 50 x 3, sd 3.5 | 0.0267 | 0.7510 | 0.3419 | YES |

The gradation is smooth and 0.25 sits in a sensible place in it. This
also explains part of my failure to reproduce the lane's 0.429: the
lane's table reports the ESS at each design's own LAPLACE optimum,
whereas the fit reports it at the final CORRECTED anchor, which is a
better anchor. On the 40 x 3 sd 2.0 design that is 0.4802 at the
Laplace optimum against 0.7308 at the fit's anchor. The two numbers are
not the same quantity and the findings do not distinguish them.

## 12. R CMD check --as-cran: Status OK
Tarball built from the worktree, checked into ri-lib with
`_R_CHECK_CRAN_INCOMING_=false`, `--as-cran --no-manual`,
`NOT_CRAN=true`.

**Status: OK. Zero ERRORs, zero WARNINGs, zero NOTEs.**
  R code for possible problems [25s] OK
  examples [44s] OK
  examples with --run-donttest [36s] OK
  tests [485s] OK
  re-building of vignette outputs [102s] OK

ONE THING THE LANE SHOULD KNOW: my tests phase is 485s where the lane
reports 130s. The difference is `NOT_CRAN=true`. `test-importance.R`
opens with `skip_on_cran()`, so under a plain `R CMD check` the whole
file - including the 2000-draw probe fit with `se = TRUE` - is skipped.
The lane's "tests [130s] OK" therefore did NOT exercise the importance
tests inside the check. Mine did. Both come back OK, but CI budget
should assume roughly +350s wherever NOT_CRAN is set.

### Surface
`git status` in the worktree is exactly the declared surface plus my
one review file:
 M NEWS.md, R/allfit.R, R/compat.R, R/confint.R, R/fit.R,
   R/influence.R, R/methods-fit.R, R/objective.R, R/sandwich.R,
   R/sugar.R, man/frm.Rd, man/frmtmb_control.Rd,
   vignettes/diagnostics.Rmd
 ?? R/importance.R, tests/testthat/test-importance.R,
    dev/importance-findings.md, dev/laplace-bias-probe.R,
    dev/review-importance.md  <- the only file I wrote
Nothing outside R/, tests/testthat/, man/, NEWS.md,
vignettes/diagnostics.Rmd and dev/. DESCRIPTION untouched at 0.49.1.
No commits.

### The main checkout: clean, but it MOVED (not by me)
At the start of this review main was at a233c3c. It is now at
**3a856e6**, still clean apart from the pre-existing untracked
`dev/brms-likelihood-tests.md`. I made no write of any kind to the main
checkout; other lanes merged during the review (wt-ropensci,
wt-readme-words, and the ddm extension).
This does NOT affect the lane: `git diff --name-only a233c3c 3a856e6`
touches `R/srr-stats-standards.R`, `README.md`, `vignettes/brms-migration.Rmd`,
`dev/`, `docs/` and `extensions/frmtmb.ddm/**` - and NOT ONE file on
this lane's surface. The branch still applies cleanly.
Two consequences for the integrator: the lane should be re-based onto
3a856e6 and the suite re-run once before merge, and someone should
check whether `R/srr-stats-standards.R` now wants an srr tag for the
new `importance` code path.

---

# VERDICT: GO WITH FIXES

The estimator is correct, the design decisions are the right ones and
every one of them is backed by a measurement I reproduced. The gaussian
identity is exact, the probe lands on AGHQ at 0.12 MCSE, the precision
extraction is machine-exact on all thirteen whitelisted structures and
demonstrably wrong (and correctly refused) on the one structure outside
them, the divergence guard catches a real failure the diagnostics
cannot see, the refactor is bit-for-bit, the suite is clean and
R CMD check is OK. The suite-count scare is a library artifact, not a
disabled test.

Nothing here is a NO-GO. The fixes below are two substantive ones and a
handful of accuracy repairs, none of which touch the estimator.

## Punch list

**P1. `frm_compat("importance")` reports three structures wrong.**
R/compat.R:245-253. `rr` reports "conditional", `us_t` and `diag_t`
report "untested"; all three are SCOPE-REFUSED by name at
R/importance.R:186-196. The new vignette section points users at this
table ("`frm_compat("importance")` lists the rules"), and the lane's own
DELIBERATE OMISSIONS list already says us_t/diag_t are refused, so the
table contradicts both the code and the lane's intent.
FIX: add a group holding every covstruct in neither `importance_blocks`
nor `crosslevel_blocks` and give it status "refused" with the
whitelist reason.

**P2. `confint(method = "profile")` on an importance fit is unguarded
and can report a bound from the degenerate regime.**
R/confint.R:436 calls `TMB::tmbprofile(object$obj, ...)` on the frozen
final tape; it does not refit, so the line added at R/confint.R:1190 is
NOT this path (it is `anova_refit_ml()`, unreachable for an importance
fit because importance refuses REML). Measured: over the fixed-effect
profile's own range the worst group still holds 0.73 of its draws
(fine), but `confint(parm = "theta_1", method = "profile")` walks theta
by +0.43 from the estimate, and at +0.5 the minimum ESS is 0.0921,
below the 0.25 floor. Nothing warns; `fit$importance$ess_min` describes
the anchor only.
FIX: at minimum a documented caveat; better, a guard or an ESS check
along the profile path. And correct the findings' claim that
"confint profile" is one of the four threaded refit paths.

**P3. Every importance round re-optimizes from a COLD start.**
R/importance.R:672, `optimize_obj(fz$obj, control, bounds, par_units,
verbose = vb)` with no `start_par`; `optimize_obj()` defaults to
`obj$par`, which is the START TEMPLATE (measured 0.4879 away from the
Laplace optimum on the probe). The Laplace warm start is used for the
anchor and thrown away for the optimizer. Costs a full cold
optimization per round on the package's most expensive tape, and lets
each round land in a different basin on a ridge - the very failure the
divergence guard exists to catch.
FIX: pass `start_par = par`. One argument, but it changes the reported
numbers, so the lane must re-measure the probe and the nine-design
table after it.

**P4. The 0.43 ESS figure is not reproducible and is quoted in three
places.** R/importance.R:775 (`imp_ess_floor` comment, "1.7x margin"),
vignettes/diagnostics.Rmd (the ESS bullet), NEWS.md. Sweeping the
design as described over seven seeds gives 0.4802 as the worst; the
lane's measurement script was lost in the scratchpad wipe and the seed
is not recorded. Also, the findings' ESS table is measured at each
design's LAPLACE optimum while the fit reports at its CORRECTED anchor
- 0.4802 vs 0.7308 on the same design - and the two are not
distinguished.
FIX: re-measure with a recorded seed, say which anchor, and restate the
margin. The threshold itself is fine; only the citation is loose.

**P5. `imp_verify()` pins the total, not the per-group decomposition.**
R/importance.R:566, `ours <- sum(a[, j] - plan[["zsq"]][, j])`.
Measured: +5 / -5 on two groups passes silently. Per-group correctness
is what the ESS diagnostics rest on.
FIX: compare the per-group vector, or at least say in the comment that
the pin is on the sum.

**P6. Findings accuracy, for the archive.** The cens() row (Laplace nll
347.8633, "40 of 240 censored") does not match the test's own
construction, which gives 311.1358 and 52 of 240; re-transcribe from
tests/testthat/test-importance.R:368. The suite figure "4903 passed /
43 skipped" is an artifact of the lane's incomplete `is-lib` - the real
figure is 5347 / 4 - and should be corrected before the findings are
archived. And note that the lane's R CMD check ran without NOT_CRAN, so
it never exercised test-importance.R.

**P7. Small, optional.** Upgrade the `autoscale`, `sparse_x` and `mo()`
compat rows from "untested" to "works" - all three measured working and
preserving the correction (item on the three untested rows above).
Consider using the package's own `frmtmb_ad_overload()` in
`lp_eta_fixed()` instead of the hand-rolled two-binding spelling; the
hand-rolled version is what went wrong once already.

## Edits I made to the worktree
ONE, and it is not code:
- CREATED `C:/Users/adf44/source/r/frmtmb-wt-islap/dev/review-importance.md`
  (this file). No other file in the worktree was created, modified or
  deleted. No commits. The main checkout was not written to at all.
All experiments that needed modified sources were run against COPIES in
the scratchpad: `ri-noguard` (divergence guard disabled), `ri-break`
(imp_verify perturbation), `ri-pristine` (a233c3c via `git archive`).

## The several-blocks-over-one-factor follow-up
The lane's description of what is needed is RIGHT, and its judgement
that blocks over different factors stay refused is right. One
refinement: the contiguity assumption lives in FOUR places
(R/importance.R:308, 344, 364, 477), and the write-up does not name 364
(`zsq`). Add it so the follow-up does not miss it.

## Version
0.50.0 is correct and nothing here argues otherwise. The change is
purely additive and proved bit-for-bit non-disruptive to every existing
objective. Sharing the bump with the get_prior() lane is fine.
DESCRIPTION still needs the integrator's bump from 0.49.1.
