# Short check of lane wcheck2's fix round (round of 2026-11-13)

Scope: the fixes of F1 to F3 of `dev/reviews/2026-11-13-wcheck2.md`, the
"Fix round" section of `dev/wcheck2-findings.md`, and the user's
decision of 2026-11-13 (`check_ald_laplace()`). Not a new review of the
lane. Tree `frmtmb-wt-wcheck2`, library `wt-wcheck2-lib`, base
rellib-r21 (ec8178d3). Scripts: `dev/reviews/wcheck2-short-ald.R`,
`dev/reviews/wcheck2-short-student.R`, `dev/reviews/wcheck2-short-table.R`;
the review's own scripts were rerun unchanged. Logs were kept outside
the tree and deleted; the numbers below are copied from them. At most 2
R processes ran at once. No lane file was edited.

## Verdict: FIX

One FIX item (S1). Everything else that was asked passes.

### S1. The refusal fires on a gaussian block beside a Student-t block, where the proof does not hold

`ald_laplace_blocks()` (`R/ald-laplace.R:117-122`) skips a Student-t
block (`dist_nu`) and a block of another package, but it still refuses
every listed block of the same mu beside them. The proof needs `H =
Sigma^-1` for the whole inner problem. A Student-t block adds `H_t(b_t)
= -d2 log t(b_t)`, which depends on the inner mode. Its `-1/2 log H_t`
term does not cancel, and the mode moves as the gaussian block's sd
changes. Thus the approximation can fall as that sd grows. The lane's
own roxygen ("The proof needs the zero curvature, so the check stays
open for: ... a Student-t or fixed-covariance block, a structure
another package registered") and the HIGH backlog entry ("What
remains: ... Student-t blocks") say that this case stays open. The code
does not do that.

Case: `asym_laplace`, `y ~ x + (1 | g) + (1 | gr(h, dist = "student",
dist_nu = nu))`, 20 groups of 8, seed 20261116.
`check_ald_laplace()` refuses `(1 | g)` in each arm
(`wcheck2-short-ald.R lane student`). TMB's Laplace objective along
g's log sd, from the optimum, the rest re-optimized:

| arm | log sd + -2, -1, 0, 1, 2, 4, 8 | largest fall |
|---|---|---|
| nu 2.5, h a copy of g | non-finite, -100.129, -110.526, -105.524, -105.455, -105.444, -105.444 | 10.4 |
| nu 5, h a copy of g | -102.987, -104.520, -102.751, -102.710, -102.644, -102.644, -102.644 | 1.53 |
| nu 2.5, h crossed (8 levels) | non-finite at every point | |

TMB's inner search is not reliable for this family. Thus I computed the
exact-mode approximation by hand (`wcheck2-short-student.R`, the same
data, h a copy of g, exact inner maximum, b0 and b1 from `lm()`, sigma
0.3, quantile 0.5). With the t scale held at 1.304, the approximation
along log sd_g = -3, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8 is NA (H_t <= 0)
for the first five values, then -131.813, -139.273, -136.167, -135.523,
-135.511. It falls by 7.46 from log sd_g 0.5 to 1.0. The value at 0.5 is
3.7 above the value at 8. The proof says that the value never falls
for every value of the other parameters. That is false here. I did not
find a finite maximum. The t scale's own direction is flat or improper
too, because the profile over the t scale ran to the bound of its
search at most points. Thus the model can have no maximum for another
reason. But the refusal's reason, "the approximation never falls as
their standard deviation grows", is false for this model, and the
round rule asks for a proof.

Fix: in `ald_laplace_blocks()`, return `character()` when any own block
of the response has a Student-t density or comes from another package.
Gaussian blocks of any covariance are safe beside a refused block,
because each block's precision is constant in b, so its constant and
curvature cancel. Add the case to the "stays open" test. In the compat
note and the NEWS bullet, add "a Student-t block in mu" to the
exclusions.

## Items

### 1. check_ald_laplace()

(a) zero_inflated_asym_laplace. The premise of the brief does not hold
for frmtmb. A zero row's log density is `log(zi)` alone
(`R/families.R:4068-4077`: `i0 * g$l + (1 - i0) * (g$l1m + ald)`), not
`log(zi + (1 - zi) ALD(0 | mu))`. brms 2.23.0 does the same
(`chunks/fun_zero_inflated_asym_laplace.stan`: `if (y == 0) return
bernoulli_lpmf(1 | zi)`). Thus a zero row has no dependence on mu, and
the proof holds. Measured (`wcheck2-short-ald.R lane zial`, 30 groups of
4):

| seed, zeros | frm() | TMB inner Hessian, max abs(H - I / sd^2) | TMB profile never falls | exact-mode profile never falls |
|---|---|---|---|---|
| 20261113, 100 of 120 (10 all-zero groups, 3 of 4 in the rest) | refused | 8.47e-22 | TRUE | TRUE |
| 20261114, 68 of 120 (4 all-zero groups, 2 of 4 in the rest) | refused | 0 | TRUE | TRUE |
| 20261115, 90 of 120 (3 of 4 in every group) | refused | 0 | TRUE | TRUE |

Example (seed 20261114), log sd + -2 to +8 about the optimum: TMB
-169.725, -153.516, -142.814, -119.009, -113.851, -113.104, -113.048;
exact modes -157.048, -133.356, -109.747, -105.837, -105.243, -105.152,
-105.150. The ZIAL refusal is a proof (subject to S1).

(b) Weights. Pass. They scale each row's piecewise-linear term, so the
curvature stays 0 (`wcheck2-short-ald.R lane weights`, w uniform in
(0.2, 3), 20 groups of 4). In both seeds, the result is refused,
max abs(H - I / sd^2) = 0, and the TMB and exact-mode profiles never
fall. For seed 20261113, the exact-mode values are -203.603 to -130.829.

(c) Exclusions, spot-checked by calling `check_ald_laplace()` on the
frame and by calling `frm()` (`wcheck2-short-ald.R lane excl`):
`trunc(lb = -20)`, `me(x, sdx)`, a nonlinear body, `quadrature = TRUE`
and an `exponential(1)` prior on class sd all return NULL (the check
stays open). The fits run to conv 1 with the kink's warnings, as on
rellib-r21. Quadrature errors with "the Gauss-Kronrod objective broke
down". The control `(1 | g)` is refused. Note: `ar(time = tt, gr = g)`
on asym_laplace is a latent residual, not `frame$autocor`. Thus the
autocor exclusion does not reach it, and the check refuses `(1 | g)`
beside it. This is harmless: the latent is gaussian and linear in mu,
so the proof holds, and `frm()` refuses the model first with the
weights check's own `(1 | obs)`-type refusal.

(d) Refused structures: us, diag, homdiag, cs, ar1, gr_cov, smooth
(`s()`, `t2()`) and hsgp. Pass. Each is a gaussian block whose
covariance is a free scale squared times a shape from the other
parameters (ar1's and cs's correlation, the gr matrix, and the HSGP
spectral weights at fixed lscale). Each enters mu as `Z b`, which is
linear in the effects. A correlation parameter is nonlinear in the
covariance but does not change the argument, because the proof holds
at every value of the other parameters. The HSGP's sdgp multiplies the
basis, which is the same as scaling the covariance (b to c b). The lane
profiled each structure (findings table, all "never falls").

(e) Message: plain. It names the term, the family and the reason
(linear in mu on each side of its mode, no curvature, never falls). It
points to `frmtmb.sample::frm_sample()`, and it gives the measured
quadrature breakdown "on 40 test data sets its rule broke down on 24".
That is 8 + 5 + 4 + 7 of the review's four designs, which is correct.
After S1, the clause "never falls" is true for every refused model.

### 2. F1 cost: pass

The review's scripts were rerun on the lane library. Times are of
`check_weights_unbounded()` alone, in seconds:

| case | rellib-r21 (review) | first lane (review) | lane, this rerun |
|---|---|---|---|
| unweighted, n 2000, 100 levels | 0.04 | 3.00 | 0.02 |
| unweighted, n 4000, 200 levels | 0.16 | 13.75 | 0.19 |
| w 0.95, n 6000, 300 levels, no prior / Intercept prior | 0.80 / 0.61 | 48.07 / 47.80 | 0.33 / 0.26 |
| w 0.95, n 10000, 1000 levels | 10.24 / 8.87 | stopped | 4.84 / 5.39 |
| w (0.9, 1.1) + prior, n 6000, 300 (`dev/wcheck2-time-pair.R`) | 4.44 (lane's record) | | 0.64 |
| same, n 10000, 1000 | 190.93 (lane's record) | | 7.07 |

Whole `frm()` unweighted n 4000: 3.34 s (review: 5.14 base, 17.65 first
lane). Case table: I reran the three parts of `dev/wcheck-cases.R` on
the lane library and compared them with the lane's recorded rellib-r21
logs by verdict and value. The result is "58 of 58 rows identical".

### 3. F2 saturated designs: pass

`Rscript dev/reviews/wcheck2-review-pair.R lane grid` gives the
following. The lines are identical to the lane's record:

    three        agree_ref 42, agree_pass 108, false_refusal 0, miss 0, lp_disagree 0
    four_unbal   agree_ref 77, agree_pass 73, false_refusal 0, miss 0, lp_disagree 1
    inter        agree_ref 73, agree_pass 77, false_refusal 0, miss 0, lp_disagree 5
    fac_x        agree_ref 47, agree_pass 103, false_refusal 0, miss 0, lp_disagree 2
    sum3         agree_ref 49, agree_pass 101, false_refusal 0, miss 0, lp_disagree 0
    three_bnorm  agree_ref 0, agree_pass 143, false_refusal 0, miss 7, lp_disagree 0
    three_w0     agree_ref 48, agree_pass 102, false_refusal 0, miss 0, lp_disagree 1

The grid has 0 false refusals in 1050 configurations. f1 * f2 and
contr.sum have 0 misses. The 7 misses of the held column are filed.

Saturation (`unbounded_row_patterns()`, `R/weights-unbounded.R:466`):
the distinct rows are counted exactly (`match()` on the double values,
column by column). No tolerance is used. The rank is `qr(D)$rank` at
qr's default tolerance of 1e-7, on the U distinct rows. Rounding cannot
flip it toward a refusal. A truly dependent D leaves a residual near
1e-16, far below 1e-7, so it is never called saturated. An
ill-conditioned independent D can be called unsaturated, which falls
back to the cell units (fails open). Rows that differ only by rounding
count as distinct patterns. That raises U (to more than p it returns
NULL), and each distinct row of a full-rank D can move alone, so the
proof does not depend on the grouping.

### 4. F3 filings: pass

Under "Filed in the round of 2026-11-13 / From lane wcheck2" there is
the HIGH item "asym_laplace's Laplace approximation has no maximum with
a random effect in mu". It has the review's table and says what ships
and what remains. The LOW items are there: all-censored asym_laplace,
the shape pair rule, the mislabeled arm, the pair rule's misses (a held
column, unsaturated designs), huber `(1 | g4)` and identity-link sigma
to 0. Its "What remains" lists Student-t blocks as open, which the code
does not do (S1).

### 5. NEWS and compat: pass, after S1

There are three BREAKING bullets (point-mass families, the dispersion
intercept prior, asym_laplace). There are no lane names, `dev/` paths
or version numbers. The compat rows (`R/compat.R`, compared with
ec8178d3 on its own lines 2128-2135) are `<structure> x asym_laplace`,
"conditional", for us, diag, homdiag, cs, ar1, gr_cov, smooth, hsgp,
`s()` and `t2()`, plus `gp_pred() x asym_laplace`, which separates the
approximate gp (refused) from the exact gp (open). They match the code,
except that the Student-t exclusion is missing from the note (S1). The
p072 note (asym_laplace x gr_cov, from "equivalent" to "defect, frmtmb
refuses", record not edited, rerun at consolidation) is in the findings
under "F3 and the user's decision".

### 6. Tests on the lane library: pass

One file per process, 2 at a time, `dev/wcheck2-tf.sh lane`; every
`lib:` line names `wt-wcheck2-lib`:

    frmtmb test-ald-laplace.R             pass=21  fail=0 err=0 skip=0 warn=0
    frmtmb test-weights-unbounded.R       pass=210 fail=0 err=0 skip=0 warn=0
    frmtmb test-diagnostics-ux.R          pass=134 fail=0 err=0 skip=0 warn=0
    frmtmb test-v14.R                     pass=81  fail=0 err=0 skip=0 warn=0
    frmtmb test-se-exact-hessian.R        pass=35  fail=0 err=0 skip=0 warn=0
    frmtmb test-simulate-density.R        pass=620 fail=0 err=0 skip=0 warn=0
    frmtmb test-numerical-robustness.R    pass=723 fail=0 err=0 skip=0 warn=0
    frmtmb test-optstart.R                pass=69  fail=0 err=0 skip=0 warn=0
    frmtmb test-mi-families.R             pass=194 fail=0 err=0 skip=0 warn=0
    frmtmb test-latres.R                  pass=305 fail=0 err=0 skip=0 warn=0
    frmtmb test-families.R                pass=252 fail=0 err=0 skip=0 warn=0
    frmtmb test-cdfs.R                    pass=234 fail=0 err=0 skip=0 warn=0
    frmtmb test-brms-suite-brmsterms.R    pass=5   fail=0 err=0 skip=0 warn=0
    frmtmb test-brms-lik-cdfs.R           pass=70  fail=0 err=0 skip=0 warn=0
    frmtmb test-compat.R                  pass=868 fail=0 err=0 skip=0 warn=0
    frmtmb test-compat-extensions.R       pass=2   fail=0 err=0 skip=0 warn=0
    frmtmb test-compat-register.R         pass=103 fail=0 err=0 skip=0 warn=0
    frmtmb.sample test-sds-prior.R        pass=28  fail=0 err=0 skip=0 warn=0
    frmtmb.sample test-mi-laplace-route.R pass=3   fail=0 err=0 skip=0 warn=0
    frmtmb.sample test-brms-suite-brmsterms.R pass=4 fail=0 err=0 skip=0 warn=0
    frmtmb.sample test-default-priors-dpar.R  pass=57 fail=0 err=0 skip=0 warn=0
    frmtmb.sample test-weights-unbounded.R    pass=5  fail=0 err=0 skip=0 warn=0
    frmtmb.eam test-family.R              pass=61  fail=0 err=0 skip=0 warn=0

Each count is equal to the lane's record
(`dev/release/wcheck2-suite-par2.log`). No test covers S1. The fix
needs a "stays open" case, a model with `(1 | g)` beside
`(1 | gr(h, dist = "student"))`.

Free space on C: 146.55 GB at the start of the check, 146.02 GB at the
end. The temporary logs are deleted. The scripts
`dev/reviews/wcheck2-short-*.R` are kept.
