# Lane setier: standard errors frmtmb reports that it should not

Base: 4f5ea39f (frmtmb 0.68.1), branch wt-setier. Library
`C:/Users/adf44/source/r/wt-setier-lib` (core only; the seven
extensions load from rellib-r6 against this core, which every `lib:`
line of the suite logs shows). "Base" is rellib-r6. "OpenBLAS" is
R 4.6.1 with OpenBLAS 0.3.26 built by `dev/setier-openblas.sh` (the
cifix emulator with this worktree's paths; the 2000 x 2000 product
check ran in 0.03 s); "ref" is the reference BLAS. Scripts are
`dev/setier-*`; logs are in `dev/setier-log/` (local, gitignored).
Every number below was printed by the script named beside it; the
seeds are in the scripts.

## Punch round 2b

Final check after punch round 2: one blocker, RQ1. Logs in
`dev/setier-log/p2b/` (local).

### RQ1: no boundary verdict under quadrature

Under `quadrature = TRUE`, `se_sd_gain_up()` returns 0, because the
quadrature objective moved up from a stuck sd gave wrong signs (round
2). So the upward check never ran there, and a fit stuck short of the
maximum got the boundary message. I found no cheap way to run that
check on a correct quadrature objective, so this round takes the
conservative carve-out. Under quadrature, `se_edge_sd()` and
`se_boundary_names()` return nothing, and `se_boundary_stop()` does
not stop. Such an sd keeps base's "flat" SE warning, and any
convergence warning stays. NEWS states the exception.

Test "under quadrature no sd is called a boundary" (ri20, seed 36,
y x 1e-3, `quadrature = TRUE`). On the reviewer's trial library
`setier-rev-lib` (punch-2 setier) it fails with 2 expectations: one
boundary message, and se_lost "boundary" or "short"
(`p2b/stf-round2.txt`, 88 pass, 2 fail). On the lane it passes
(test-se-tier.R 90 pass, 0 fail).

- `dev/setier-rev3-trapq.R` (`p2b/trapq.txt`): no boundary message.
  Seeds 19 and 36 at 1e-3 say nothing, and seeds 21 and 23 at 1e3 give
  the SE warning "flat". This is the same output as rellib-r6, so the
  silent short quadrature fits were there before the lane (backlog,
  "stops short of the maximum").
- `dev/setier-rev3-quad.R` (`p2b/quad.txt`): the same as the
  reviewer's base log. Arm A: "nothing 20" in all four cells. Arm B:
  SEWARN 5 for each family (the script runs 5 seeds a family). On
  punch 2, arm A at sd 0 gave BOUNDARY 9 and 8. Those sds are at zero,
  but under quadrature they now keep base's verdict.
- The singular study has no quadrature arm.
- test-quadrature-defects.R: 59 pass, 0 warnings escaped. It is the
  only test-quadrature* file. One fit there now keeps the "flat" SE
  warning, so it is wrapped in `allow_warnings()` with that warning
  named.

### Separation warning: the list of coefficients is capped

`separation_check()` names the first `sep_name_max` (5) coefficients
and the count. `dev/setier-rev2-sepmiss.R`, 2050-level case
(`p2b/sepmiss.txt`): "mu: f10, f13, f14, f22, f25 and 334 more (339
coefficients)". All 339 still lose their SE ("separation 339"). The
other cells are as in punch 2. NEWS says the warning names the first
five and the count.

### Suite and check

Core suite, reference BLAS, one file per process, 10 at a time
(`p2b-ref`, 216 files, every lib: line wt-setier-lib): 17346 pass,
1 fail, 0 err, 0 skip, 0 escaped warnings. The failure is
test-perf.R's wall-clock bound under load (1.07 against 1.00). Run
alone, it passes 3 of 3 on the lane and on rellib-r6
(`p2b-perf-alone.txt`). Against `p3-ref` (punch 2), 5 files differ:

    test-se-tier.R            pass 88 -> 90 (the new test)
    test-quadrature-defects   boundary 2 -> 0, SE warning 0 -> 1
    test-fuzz.R               boundary 119 -> 112, SE warning 131 -> 136
    test-id-kron.R            boundary 0 -> 1, not converged 1 -> 0
    test-perf.R               the load failure above

test-fuzz.R: the changes are all next to its quadrature fits. Their
boundary messages became "flat" SE warnings or nothing, which is what
base says for those fits.

test-id-kron.R does not depend on the build. Its merged-dpar fit ends
at one of two optima from run to run. Ten copies of the file at once:
8 of 10 end at objective 128.07435688 on rellib-r6, and 7 of 10 on
the lane. The others end at 128.074356883, with gradient 5.6e-4
against 1.6e-4. The fit alone (`dev/setier-idk.R`), 8 at once, is the
first one 8 of 8 times. The first optimum is a boundary stop and the
second keeps code 7. p3-ref got the second and p2b-ref the first;
both pass. Filed in `dev/test-backlog.md` (medium) with
`dev/setier-idk-test.R`.

R CMD check --as-cran, frmtmb, after the last package edit
(`check-frmtmb-p2b.txt`): Status: 1 NOTE (V8 math rendering).
Tests: [ FAIL 0 | WARN 0 | SKIP 376 | PASS 11571 ].

Files changed in this round: `R/se-check.R`, `NEWS.md`,
`tests/testthat/test-se-tier.R`,
`tests/testthat/test-quadrature-defects.R`, `dev/test-backlog.md`,
`dev/setier-idk.R`, `dev/setier-idk-test.R`.

## Punch round 2

Re-check: "Re-check after punch round 1" in
`dev/reviews/2026-10-07-setier.md` (NOT MERGEABLE on RB1 and RB2). The
reviewer's `dev/setier-rev2-*` scripts were run unchanged against the
lane library; logs in `dev/setier-log/p2/` (local). New scripts:
`dev/setier-trapdbg.R`, `dev/setier-quadup.R`.

### RB1: symmetric second differences in the curvature probe

`se_curvature_real()` now takes `c(t) = f(p + t) + f(p - t) - 2 f(p)`
at the full step (the one the curvature says loses 2 * `grad_tol` each
way) and at half of it, and keeps the direction when `c(full)` is at
least 2 * `grad_tol` and `c(full) / c(half)` is in (3, 5.5). The
gradient's linear term cancels, so a fit that is not exactly stationary
along a weak direction keeps it.

`dev/setier-rev2-window.R lane` (`p2/final/window.txt`), 20 seeds each,
"max |SE/ref - 1|" against glm() among fits losing none:

    pois   pair at cor 1 - 5e-11   lost 0 of 20 (punch 1: 20)   0.0191
    bern   rare events, same pair  lost 0 of 20 (punch 1: 20)   0.154
    ppoly  poisson degree 5        lost 0 of 20 (punch 1: 4)    0.0179
    bpoly  bernoulli degree 5      lost 0 of 20 (punch 1: 7)    0.0457
    ubnear degree 5 near a ub      lost 0 of 20                 0.0556
    expb   a * exp(b * x), x on [1, 1 + 1e-5]: 11 of 20 lose a and b
           (base 10): conservative, recorded

The reference differences are base's own (0.019 and 0.154 on base).
c0k seed 77 still loses c0, a, k and theta_1 (`p2/final/c0k.txt`); the
raw polynomials of degree 3 to 6 keep every SE, at lm()/lme4 ratios
0.9821 to 1.004 (`p2/final/collin.txt`); degree 7 with (1 | g) now
loses 3 (punch 1: 5).

Test "a weak but identified direction keeps its SE off stationarity"
(the pois pair, seed 1). It fails on the punch-1 build (the reviewer's
trial merge library `setier-rev-lib`, punch-1 setier with optima),
3 expectations (`p2/stf-tier-round2.txt`).

m-c: `se_relabel()` says "separation" only in a predictor whose
certificate fired (`fit$cache$se_explained_sep`). Test "only a proved
separation relabels a flat coefficient" fails on the punch-1 build.

### RB2: an sd left near zero short of the maximum

`se_sd_gain_up()` moves a flat log sd that passed the edge test up to
0.01, 0.03, 0.1, 0.3 and 1 times the scale of its term (the residual
sd where the response has one, else 1 on the link scale), the rest
held. At a boundary every interior value is worse, so any gain above
`se_edge_tol` (1e-6) refuses the boundary verdict, in `se_edge_sd()`
and in `se_boundary_names()`. Such an sd gets the new reason "short":
"the fit stopped short of the maximum: this standard deviation is near
zero, where the likelihood is flat, but a larger value raises the
log-likelihood, so it is not at a boundary. Refit from a larger
starting value". `se_boundary_stop()` then keeps the convergence
warning, since not every loss is a boundary.

`grad_tol` as the threshold, as first suggested, is too loose: on seed
23 at x1e3 (0.0007 below lme4, whose sd is 0.053 sigma) the best gain
with the other parameters held is 4.3e-4 at 0.03 sigma, and 0.1 sigma
already loses (`dev/setier-trapdbg.R`); hence the grid starting at
0.01 and the 1e-6 threshold.

Under quadrature the question is not asked: that objective is not
accurate above the sd it was fitted at. On test-quadrature-defects.R's
nested beta fit (an sd at exp(-11.2)) a move up by 2 lowered it by
0.009 and a move to 0.01 by 17, where the Laplace objective rose by
0.0015 at 0.01 (`dev/setier-quadup.R`); the first version of the check
turned that boundary into "short" and the suite caught it.

`dev/setier-rev2-trap.R` (`p2/final/trap.txt`):

    seed 19 x1e-3  code 0  -0.666  SE warning, "stopped short"
    seed 36 x1e-3  code 8  -0.825  "Optimizer did not report
                                   convergence: false convergence (8)"
    seed 21 x1e3   code 0  -0.0045 SE warning, "stopped short"
    seed 23 x1e3   code 0  -0.0007 SE warning, "stopped short"

`dev/setier-rev2-short.R` (24 short fits), what the user is told:

    build     convergence warning  SE warning  boundary msg  nothing
    base      3                    18          0             3
    punch 1   2                    16          4             2
    punch 2   3                    19          0             2

The two told nothing (seeds 3 and 20 at x1e3) stop 8.6e-4 and 5.9e-4
below lme4, within `grad_tol`; lme4 calls both singular, and frmtmb
gives no boundary message because its sd stops near 7e-4 sigma, where
the edge step still moves the log-likelihood by more than 1e-6.

Test "an sd the fit left near zero short of the maximum is no
boundary" (seeds 36 and 19 at x1e-3): fails on the punch-1 build, 4
expectations.

Coverage is unchanged (`dev/setier-singular.R`, `dev/setier-rev-
smallsd.R`, reference BLAS, `p2/sing-ref.tsv`, `p2/smallsd-ref.txt`):
267 of 269 lme4-singular fits get the message and 0 of 331 others;
s1 8 of 8, s2 1 of 2, p1 7 of 7, c95 9 of 13, c99 10 of 20 plus the 2
correlation-near-1 fits. `dev/setier-rev2-scale2.R`: at x1e-3 17 of the
18 lme4-singular fits are flagged, at x1 18 of 18, at x1e3 0 of 18 and
0 others (punch 1: 2 others, seeds 21 and 23); at x1e3 the optimizer
stops at about 7e-4 sigma where the log-likelihood still moves by more
than 1e-6 over the edge step, which is the trap below and not this
check.

### Minors

- m-a: NEWS says separation is named once the fit has run far enough,
  with the reviewer's budget counts (`dev/setier-rev2-sepmiss.R`,
  `p2/sepmiss.txt`: complete separation at eval.max 5, 10, 20, 40 named
  on 5, 16, 19, 20 of 20; factor quasi-separation at 5, 10, 20 on 0, 0,
  20).
- m-b: above `sep_p_max` columns, `sep_column_cert()` tries each column
  alone (+e_j or -e_j, O(nnz)); a column all of whose nonzero rows move
  toward their own outcome is a certificate by the definition, so it
  cannot raise a false alarm. The 2050-level factor (39 levels all 0,
  300 all 1) is now named, its 339 level coefficients lost as
  "separation" (fit 19.9 s).
- The separation study (`dev/setier-sep.R`, reference BLAS,
  `p2/final/sep-ref.tsv`): 162 of 162 LP-separated fits named, 0 of 338
  others.

### Recorded for the user: the optimizer trap (pre-existing)

Gaussian random-intercept fits with the response multiplied by 1e3
stop more than 1e-4 below lme4's log-likelihood on 20 of 40 seeds, by
up to 0.825 over the scales, and on 4 of 40 at 1e-3 (none at 1), 23 of
the 24 at optimizer code 0. rellib-r6 gives the same fits. It is a
silent wrong answer on base (3 of 24 told nothing); this lane now says
"stopped short" on the code-0 ones that lose the sd's standard error,
but the estimate is still wrong. Filed in `dev/test-backlog.md`
("Filed by lane setier, 2026-10-07") with the reviewer's scripts, for
the optimizer's lane.

### Suites and check

One file per process, every gate on, on the final R code (install18,
after the last R edit); every lib: line shows wt-setier-lib.
pass/fail/error/skip/escaped warnings:

    p3-ref  frmtmb           17345/0/0/0/0 (216 files)
    p3-ref  frmtmb.coupling    542/0/0/5/0 (11 files)
    p3-ref  frmtmb.eam        1743/0/0/3/0 (29 files)
    p3-ref  frmtmb.latent      360/0/0/2/0 (10 files)
    p3-ref  frmtmb.learn       501/0/0/2/0 (15 files)
    p3-ref  frmtmb.ode         549/0/0/1/0 (11 files)
    p3-ref  frmtmb.sample     2616/0/0/1/0 (48 files)
    p3-ref  frmtmb.spline      590/0/0/1/0 (15 files)
    p3-ob   frmtmb           17338/2/1/0/1 (216 files)

Reference BLAS: clean on all eight suites. OpenBLAS (core only): the
emulator-only failures of rellib-r6 (bcm-latent-mixtures,
cumulative-cs, ordinal-mixture's code 7) and test-perf.R's wall-clock
bound under load (the same on rellib-r6, dev/setier-perfcheck.R).

Conditions raised (traced, muffled or not):

    p3-ref  bnd 271  SE warn 206  sep 39  not conv 142
    p3-ob   (core only) bnd 247  SE warn 194  sep 39  not conv 102

    lib lines with wt-setier-lib, p3-ref: 355
    lib lines with wt-setier-lib, p3-ob: 216

### Files touched in this round

- `R/se-check.R`: se_curvature_real() (symmetric differences),
  se_relabel() (separation gated on the certificate; "short"),
  se_boundary_names() and se_edge_sd() (the upward check), new
  se_sd_gain_up(), se_lost_clauses() ("short"), sep_certificate() and
  new sep_column_cert().
- `NEWS.md`: the probe bullet (RB1), the boundary bullet ("stopped
  short", RB2), the separation bullet (m-a, m-b).
- `tests/testthat/test-se-tier.R`: 3 new tests.
- `dev/test-backlog.md`: the optimizer trap.
- dev: `setier-trapdbg.R`, `setier-quadup.R`, this section.
- `R/fit.R` is not edited in this round.

R CMD check --as-cran on frmtmb, tarball built after the last edit of
the package (`dev/setier-check.sh`, `dev/setier-check/`): Status: 1
NOTE (V8 math rendering). Tests `[ FAIL 0 | WARN 0 | SKIP 376 |
PASS 11569 ]`.

`test-se-tier.R` now has 18 tests: pass 88 on the lane (p3-ref);
fail 8 on the punch-1 build in the 3 new tests.

## Punch round 1

Review: `dev/reviews/2026-10-07-setier.md` (NOT MERGEABLE on B1 and B2).
Scripts `dev/setier-*` (new this round: `setier-sepcost.R`,
`setier-rb4.R`, `setier-grby.R`, `setier-p1-attacks.sh`,
`setier-p1-mo.sh`, `setier-p1-mo-sum.R`); logs in `dev/setier-log/p1/`
(local). The reviewer's scripts were run unchanged against the lane
library; the two drivers that write into the reviewer's log directory
were copied (`setier-p1-*`) so that no reviewer log is overwritten.

### B1: the curvature probe and curved ridges

`se_curvature_real()` now also takes the half step each way and keeps a
direction only when both full steps lose at least `grad_tol` and each
loses 3 to 5.5 times its half step (a quadratic gives 4; the fourth
power of a curved ridge gives 16). This is the reviewer's fix, kept as
measured; the band is `se_quad_lo`/`se_quad_hi`.

- c0k seed 77 (`dev/setier-rev-c0k.R lane`, `p1/final/c0k.txt`):
  c0_(Intercept), a_(Intercept), k_(Intercept) and theta_1 lost, all
  NaN, as base. The objective along the exact ridge moves at most
  5.2e-11 (s from 0.5 to 2).
- Raw polynomials (`dev/setier-collin.R`, `p1/final/collin.txt`):
    degree 3 no RE: SE/lm() 0.9899 0.9899; lost 0; warnings 0
    degree 3 (1 | g): SE/lme4 1 1; lost 0; warnings 0
    degree 4 no RE: SE/lm() 0.9874 0.9874; lost 0; warnings 0
    degree 4 (1 | g): SE/lme4 1 1; lost 0; warnings 0
    degree 5 no RE: SE/lm() 0.9849 0.9849; lost 0; warnings 0
    degree 5 (1 | g): SE/lme4 1 1; lost 0; warnings 0
    degree 6 no RE: SE/lm() 0.9821 0.9821; lost 0; warnings 0
    degree 6 (1 | g): SE/lme4 1.004 1.004; lost 0; warnings 0
    degree 7 no RE: SE/lm() NaN NaN; lost 8; warnings 1
    degree 7 (1 | g): SE/lme4 NaN NaN; lost 5; warnings 1
- New test "a straight step off a curved ridge is not taken as
  curvature" (`test-se-tier.R`). Seen to fail on the round-1 build:
  the reviewer's trial merge library (`setier-rev-lib`, round-1 setier
  plus optima) gives 3 failures in it (`p1/stf-tier-round1.txt`); it
  passes on rellib-r6, which gave NaN there too.
- NEWS now says what the code tests.

### B2: separation without a dense design

`sep_certificate()` keeps the design as the frame holds it. The first
candidate is the estimate's own direction, checked with `X %*% beta`
(no decomposition); the null space of the holding rows comes from the
eigenvectors of their p x p cross product, scaled to unit columns, and
only when some observation has left a tail of 1e-3 (the condition
under which the old code could find a certificate at all) and p is at
most 2000. Every candidate is verified on every row, so the squared
conditioning and the loose null cut (1e-10 of the largest eigenvalue)
can only produce true certificates.

`dev/setier-sepcost.R` (`p1/sepcost-*.txt`), the check alone, minimum
of 3 calls, against the fit:

    design                     n     fit       check    check / fit
    dense, 19 covariates       1e5   19.7 s    0.140 s  0.7 %
    sparse_x, 500-level f + x  1e5   53.0 s    0.200 s  0.4 %
    dense, 19 covariates       1e6   153.9 s   1.050 s  0.7 %
    sparse_x, 500-level f + x  1e6   718.7 s   0.400 s  0.06 %

The session's largest memory use during the check stayed below the
fit's own (259 and 252 MB against 288 and 297 MB at 1e5; 526 MB against 674 MB
at 1e6 sparse; 1278 MB against 1481 MB at 1e6 dense). The reviewer
measured 54 to 68 s and 2297 MB for the old check on the 1e5 sparse
fit. The controls the script also times (one objective evaluation,
`X %*% beta`) are below the 10 ms clock except `X %*% beta` at 1e6
dense, 0.03 s. The 1e5 sparse fit is truly quasi-separated (a level
of f with no successes) and is named.

The 500-fit study (`dev/setier-sep.R`, `p1/sep-*.tsv`), both BLAS
builds, identical verdicts: 162 of 162 LP-separated fits named, 0 of
338 others. Logit link only, as the reviewer notes; probit and
cloglog complete separation are named in `dev/setier-rev-sep2.R`.

### m1: boundary coverage on larger designs

`se_edge_sd()` asks `se_at_edge()` of every log sd whose curvature is
at most `se_edge_screen` (0.1 in the optimizer's units; near zero a
move of 2 changes the log-likelihood by about a quarter of it), before
tier 1 is accepted, and tier 3 takes a TRUE as a flat row. The one
correlation of a block whose sd is at its edge is taken as undefined
with it. The edge test is now two-sided within `se_edge_tol` = 1e-6
(round 1: within `grad_tol`): with `grad_tol` it called 5 fits that
lme4 does not call singular boundary fits (sds 0.014 to 0.074 in the
singular study), whose log-likelihood falls by 1e-5 or more over the
move; at a real boundary it moves by about 1e-8. Two-sided, so that a
correlation whose log-likelihood still rises toward its end is not
called settled (the gr(g, by = f) fixture's correlation rose 9.3e-5,
`dev/setier-grby.R`; that fixture is back to nanse's verdict, theta_2
and theta_3 flat, with the SE warning).

    arm        design     n lme4s  lost       msg  warn  conv  se100
    ref        ri6      100    72    72   72 (72)     0     0      0
    ref        ri20     100    51    51   51 (51)     0     0      0
    ref        ri20s    100    18    18   18 (18)     0     0      0
    ref        rs20     100    69    67   67 (67)     0     2      0
    ref        bin15    100    57    57   57 (57)     0     0      1
    ref        bin15s   100     2     2    2 ( 2)     0     0      0
    ob         ri6      100    72    72   72 (72)     0     0      0
    ob         ri20     100    51    51   51 (51)     0     0      0
    ob         ri20s    100    18    18   18 (18)     0     0      0
    ob         rs20     100    69    67   67 (67)     0     2      0
    ob         bin15    100    57    57   57 (57)     0     0      1
    ob         bin15s   100     2     2    2 ( 2)     0     0      0
    
    dev/setier-rev-smallsd.R, 40 seeds each:
    arm des  lme4s  msg on-lme4s on-others  SE-warn  code!=0
    ref s1       8    8         8         0        0        0
    ref s2       2    1         1         0        0        0
    ref c95     13    9         9         0        0       10
    ref c99     20   12        10         2        0       17
    ref p1       7    7         7         0        0        0
    ob  s1       8    8         8         0        0        0
    ob  s2       2    1         1         0        0        0
    ob  c95     13    9         9         0        0       10
    ob  c99     20   13        11         2        0       16
    ob  p1       7    7         7         0        0        0

### m2: smooths

`s(x) + s(x2)` with `x2 = x` (`dev/setier-rev-smooth2.R`): the two sds
are flat only jointly and each costs 0.99 when moved; they are now
"flat" with the SE warning on seed 1 (round 1: "boundary", silent). The
all-flat shortcut is kept only for blocks with two or more
hyperparameters and asks a joint move along the block's flattest
direction (`se_block_ridge_edge()`). On the gamSim study 91 of 140
fits now lose a smoothing sd silently (round 1: 86; base warned on 7
ref, 5 OpenBLAS), identical on both BLAS builds (`p1/sx2-*.tsv`).

### m3, m5, m6

- m3: the coefficients the separation warning names lose their
  standard error on every path (`sdr_sep_lost()`): the separated
  logistic at nlminb's limit gave x 6.22e132 with the reference BLAS;
  now NaN with reason "separation" (`p1/attacks/...-a-rev-cases-*`).
  "weights 0" and "offset" in `dev/setier-rev-sep2.R`: 0 of 3 finite.
- m5: the probit complete separation stopped by underflow at code 0 is
  named (one warning; the non-finite-gradient warning is not given
  beside it, since the separation explains it), and its three
  coefficients are lost (`p1/final/probit.txt`, `sep2.txt`).
- m6: `se_line_probe()` guards its `min()`.

### m4: bounds

`se_tier3_could_act()` sends any fit with a bound-held parameter to
tier 3 (`se_bound_held()`, a gradient only when the box has a finite
bound). `dev/setier-rev-bound.R`: degrees 2, 3 and 5 now all report
"x_k: a bound holds it" with the others conditional on it (base: degree
2 and 3 silent with the unconstrained standard errors).

### m9, m10

- The boundary message names brms's parameters (`sd_g__Intercept`,
  `cor_g__Intercept__x`; others keep `theta_k (sd of ...)`), and says
  "at zero" for an sd, "at the limit of its range" for a correlation or
  mixing parameter, and "undefined there, since a standard deviation of
  its term is zero" for the correlation of a block whose sd is at zero.
- `check_se = "stop"` builds the check at `frm()` even where it would
  wait, so `frm()` stops (`dev/setier-rev-bnd.R` B now stops at its
  frm() call, which is the change; test "check_se = 'stop' stops at
  frm() when the check would wait").

### Item 7: nl_flat_message()

It now asks the standard-error analysis (built at the fit when it has a
candidate, `se_fit_analysis(force = TRUE)`, cached for se_check()) and
names only coefficients the check lost as flat. At cor 1 - 1e-9
(`dev/setier-b1-rho.R`, `p1/final/b1rho.txt`) it names c and dd only,
where it named a_x1 and b_x2 too; their SE stays 1268, the reference.
Edited in R/fit.R: nl_flat_message() and two lines of fit_assembled()
(the call passes the fit, and `fit$cache$bounds` is set before it).

### Item 8: code 7 at a boundary

check_convergence() gives no convergence warning when the optimizer's
code is not 0, the largest gradient is within `grad_tol` and every
standard error the analysis loses is a boundary one
(`se_boundary_stop()`, which builds the analysis even where the check
would wait); se_check() then gives the boundary message. An analysis
built on a non-converged fit that is not such a stop is not used to
repair its covariance (the fit keeps sdreport()'s own, as before; an
earlier draft of this change repaired 15 code-7 mo() fits'
coefficient SEs, now 0, `dev/setier-p1-mo-sum.R`).

Measured (`dev/setier-rev-rs20.R`, `dev/setier-rb4.R`,
`dev/setier-singular.R`, `dev/setier-rev-smallsd.R`; both BLAS builds):

- rs20 (y ~ x + (1 + x | g), no slope variance), 100 seeds: 69 code-7
  fits, all lme4-singular, logLik within 2.6e-6 of lme4, largest
  gradient at most 4.8e-4. 67 get the boundary message; 2 keep
  "Optimizer did not report convergence" because their analysis also
  loses theta_2 as "flat" (seeds 17 and 86 on ref, which two depends
  on the BLAS: their endpoints differ).
- RB4 (y ~ x + f + (1 + x | g2) and y ~ x + (1 + x | g2), 20 seeds
  each): 22 fits at code 1; 21 get the boundary message, 1 keeps the
  convergence warning (f40 seed 5: nothing lost); lme4 calls all 40
  singular, the boundary message is on 39 of them (base: 0).
- c95 and c99 (correlation 0.95 and 0.99 with real sds): 27 code != 0
  fits on ref (26 OpenBLAS); 18 (16) get the boundary message and 9
  (10) the convergence warning.
- Silenced code != 0 fits whose loss is not all "boundary": 0 in every
  study (by construction, and counted).

### m12: provenance and the mo() seeds

- Round 1's suites (19:49) and R CMD check tarball (19:37) predate the
  20:01 edit of R/se-check.R; the reviewer diffed the tarball against
  the tree, and the only difference is a roxygen comment. This round's
  suites and check ran after the last edit.
- The five mo() seeds whose verdict the lane changed: 3, 15, 103, 142
  and 171. Each loses zeta2_1 and zeta2_2 with the SE warning where
  base gave silent standard errors of 207 to 7994; seeds 3, 15, 142 are
  near-vertex simplexes with a unit-diagonal eigenvalue of -0.6 to -1,
  flat to 1e-5; 103 and 171 have a flat common shift (3e-10 to 4e-10 of
  the largest). This round's 200-seed study (`dev/setier-p1-mo.sh`)
  reproduces the reviewer's lane row exactly: 48 lost-SE fits, 31 SE
  warnings, 17 code != 0, 0 silent, 0 false alarms, the same 5 seeds,
  coefficient SEs unchanged on every row.

### Attack set (`dev/setier-p1-attacks.sh`)

Lane nanse's attack set on this round's build, reference BLAS and
OpenBLAS, compared with the reviewer's logs of round 1
(`dev/setier-rev-log/attacks/`):

- spread, crossed, pred, ridge-re: identical on both BLAS builds and to
  round 1.
- near-collinear (`-rev2-b1`): identical to round 1 except the 2-d flat
  subspace's slope ratio (2.9e-7 ref, -4.3e-9 OpenBLAS), as before.
- cor 1 - 1e-9 (`-rev2-b1-rho`): a_x1 and b_x2 keep 1268 on both paths;
  the identification warning now names only c and dd (item 7).
- RB4 sweep: 0 identified fixed effects lost on f40 and plain, both
  BLAS builds; the boundary messages now carry brms names.
- the separated logistic (`-rev-cases`): x is NaN with reason
  "separation" where round 1 printed 6.22e132 (m3).

### Suites and check

All eight suites, one file per process, every gate on (NOT_CRAN, the
brms, drmTMB and fuzz gates), on the final R code (install14, after
the last R edit); every lib: line shows wt-setier-lib (355 of 355 per
BLAS build). pass/fail/error/skip/escaped warnings:

    p1-ref  frmtmb           17334/1/0/0/1 (216 files)
    p1-ref  frmtmb.coupling    542/0/0/5/0 (11 files)
    p1-ref  frmtmb.eam        1743/0/0/3/0 (29 files)
    p1-ref  frmtmb.latent      360/0/0/2/0 (10 files)
    p1-ref  frmtmb.learn       501/0/0/2/0 (15 files)
    p1-ref  frmtmb.ode         549/0/0/1/0 (11 files)
    p1-ref  frmtmb.sample     2616/0/0/1/0 (48 files)
    p1-ref  frmtmb.spline      590/0/0/1/0 (15 files)
    p1-ob   frmtmb           17328/2/1/0/2 (216 files)
    p1-ob   frmtmb.coupling    542/0/0/5/0 (11 files)
    p1-ob   frmtmb.eam        1743/0/0/3/2 (29 files)
    p1-ob   frmtmb.latent      360/0/0/2/0 (10 files)
    p1-ob   frmtmb.learn       501/0/0/2/0 (15 files)
    p1-ob   frmtmb.ode         549/0/0/1/0 (11 files)
    p1-ob   frmtmb.sample     2616/0/0/1/1 (48 files)
    p1-ob   frmtmb.spline      590/0/0/1/0 (15 files)

- p1-ref: the failure is test-perf.R's wall-clock bound (1.39 s
  against 1.00); the same fits take 0.97 to 1.01 s on rellib-r6 and
  the lane alike (dev/setier-perfcheck.R), and the file passes 3 of 3
  alone on both BLAS builds. The escaped warning was test-autocor-cond.R
  (a box holding ar[1], now reported on every fit, m4); the test now
  allows and requires it: pass 57 on both BLAS builds, and it fails on
  rellib-r6, which says nothing there.
- p1-ob: besides those two, the emulator-only failures of rellib-r6
  (bcm-latent-mixtures, cumulative-cs, ordinal-mixture's code 7, eam
  sampling, sample brms-shapes-draws).

Conditions raised (traced, muffled or not), files they fire in:

    p1-ref  bnd 273 (61 files)  SE warn 204 (39)  sep 39  not conv 140
    p1-ob   bnd 271 (60 files)  SE warn 206 (41)  sep 39  not conv 134

Verdicts across BLAS (dev/setier-verdict-diff.R):

    p1-ref vs p1-ob: conditions 516 vs 516; files differing 7; unmatched 8
    (round 1: 8 unmatched in 6 files; rellib-r6: 19 in 10)

`test-se-tier.R` now has 15 tests: lane pass 79 on both BLAS builds;
rellib-r6 pass 42, fail 18, error 3 (`p1/stf-tier-base.txt`).

R CMD check --as-cran on frmtmb after the last edit (`dev/setier-check.sh`,
`dev/setier-check/`): Status: 1 NOTE (V8 math rendering). Tests
`[ FAIL 0 | WARN 0 | SKIP 376 | PASS 11560 ]`. frmtmb.learn did not
change this round.

### Files touched in this round

- `R/se-check.R`: cov_from_hessian(), se_tier3_could_act(), se_tier3(),
  se_line_probe(), se_curvature_real(), sdr_rescue() (split into it and
  sdr_rescue_hessian()), se_check(), se_fit_analysis() (cached, with
  `force`), se_boundary_names(), se_at_edge(), se_boundary_message(),
  se_boundary_labels(), separation_check(), sep_certificate(),
  sep_null_space(); new se_bound_held(), se_edge_sd(), se_edge_screen,
  se_edge_tol, se_quad_lo/hi, se_block_ridge_edge(), sdr_sep_lost(),
  se_boundary_stop(), sep_p_max.
- `R/fit.R`: check_convergence() (separation before the gradient
  checks, `se_sep_named`, the boundary stop of item 8),
  nl_flat_message() (item 7), fit_assembled() (two lines: the call
  passes the fit, and `fit$cache$bounds` is set before it), and the
  roxygen of frmtmb_control()'s `check_se`; `man/frmtmb_control.Rd`.
- `NEWS.md` (the development section rewritten: separation under
  Breaking changes, the probe and boundary bullets per B1, m1, item 8).
- Tests: `test-se-tier.R` (5 new tests), `test-se-check.R` (the gr-by
  test back to nanse's verdict), `test-autocor-cond.R` (m4).
- dev: `setier-sepcost.R`, `setier-rb4.R`, `setier-grby.R`,
  `setier-perfcheck.R`, `setier-p1-attacks.sh`, `setier-p1-mo.sh`,
  `setier-p1-mo-sum.R`, this section.

### Not done, and why

- m7 (degree-7 polynomial with `(1 | g)` keeps three SEs far below
  lmer's) and m8 (curved ridges with an eigenvalue above 1e-9 keep
  silent finite SEs on base and lane alike) are pre-existing, as the
  reviewer records; not touched.
- The 2 rs20 fits and the c95/c99 fits that keep a convergence warning
  at code != 0 lose a parameter as "flat" or nothing at all; item 8
  leaves them as they were by design.
- Fits that stop at different points along a flat direction on the two
  BLAS builds (lkj, frailty, brms-pins, two rs20 seeds) still differ
  in verdict; no rule on two different Hessians can make them agree.

## Round 1

### Summary

1. **Tier 3 is asked first** (backlog "Filed at 0.68.1", high, item 1).
   `cov_from_hessian()` keeps sdreport()'s inverse (tier 1) or the
   scaled inverse (tier 2) only when tier 3's own tests find nothing to
   remove (`se_tier3_could_act()`): no Hessian row within its noise, no
   unit-diagonal eigenvalue at or below 1e-9 of the largest. When they
   find something and tier 3 then removes nothing, the tier 1 inverse
   is returned, so a fit that loses nothing keeps sdreport()'s
   covariance bit for bit. Tier 3 gained a curvature probe
   (`se_curvature_real()`) so that a small eigenvalue that is real
   curvature (a raw polynomial) is kept. `autoscale_sdreport()` builds
   the outer Hessian itself when the fit-time check did not, at the
   cost sdreport() pays for it, so its noise is always known.
2. **A covariance parameter at the edge of its parameter space is a
   boundary fit, said with a message.** A flat parameter of a random
   effect block that the likelihood ignores when it moves further toward
   its end (`se_at_edge()`: a log sd moved 2 toward zero, a correlation
   or mixing parameter moved 2 toward its nearer infinity, a change
   within `grad_tol`) gets reason `"boundary"`. On a block over grouping
   levels it is reported as lme4 reports a singular fit: one message,
   "Boundary (singular) fit: ...", class `frmtmb_boundary_fit`, under
   `check_se = "warning"`; `"stop"` stops, `"ignore"` is silent. A
   smooth, gp() or hsgp() block at its own limit loses the standard
   error with no message. VarCorr(), confint_varcorr(), predict()'s
   draws and the delta-method standard errors of fitted() do not repeat
   the boundary report.
3. **The nonlinear ridge with random effects** (high, item 2) loses its
   standard errors through rule 1: its eigenvalue is below 1e-9 of the
   largest and the curvature probe finds the likelihood flat along it.
4. **Separation is proved and named** (high, item 3) at the end of every
   binomial-type fit, whatever the optimizer's code
   (`separation_check()` in `check_convergence()`), from a certificate
   checked against the definition.
5. **`ranef(condVar = TRUE)`** (medium, item 4) reads the repaired
   joint covariance when a standard error was lost.
6. **Backlog item 5** (fixef() and summary() NaN without a warning; a
   linear ridge next to a coefficient): both were closed by lane nanse
   in 0.68.0 and stay closed (measurements below).

### Item 1: the noise row, and the verdict across BLAS builds

#### The repro

`dev/setier-probe.R`, `Reaction ~ Days + (1 | Subject/a)` on
sleepstudy with `a = factor(Days %% 3)`; theta_1 (the Subject:a sd) at
log sd -29.03, its Hessian row 7.1e-12 (`probe-*-*.txt`):

    arm         diag of theta_1   SE of theta_1   what the user is told
    base ref    +7.10543e-12      375,150         nothing
    base ob     -7.10543e-12      NaN             SE warning (flat)
    lane ref    +7.10543e-12      NaN             boundary message
    lane ob     -7.10543e-12      NaN             boundary message

The other four standard errors are the same in all four arms to the
digits printed (9.50619, 0.801735, 0.0555555, 0.178981).

#### The tests that pin it

`tests/testthat/test-se-tier.R` (new). Seen to fail on rellib-r6 (logs
`stf3-tier-base-ref.txt`, `-base-ob.txt`):

    rellib-r6 ref        pass 32 fail 12 err 3 skip 0
    rellib-r6 OpenBLAS   pass 31 fail 12 err 4 skip 0
    lane ref             pass 63 fail 0  err 0 skip 0
    lane OpenBLAS        pass 63 fail 0  err 0 skip 0

The first test sets the noise row to each sign by hand and asserts the
same verdict and the same covariance of the others
(`solve(H[-j, -j])` to 1e-8) for both; on rellib-r6 the positive sign
keeps tier 1. One error on each base arm is the weak form (the probe
`se_at_edge()` does not exist there) in the guard test whose other
assertions pass on base. Two tests are guards that pass on base too, by
design: "a healthy fit with random effects is untouched and silent"
(`identical()` to `RTMB::sdreport()`, sleepstudy `(Days | Subject)` and
cbpp silent) and "an ill-conditioned but identified design keeps every
SE" (below).

#### The eigenvalue gate and its false alarm, found and fixed

The first version of the gate sent every fit with a unit-diagonal
eigenvalue at or below 1e-9 of the largest to tier 3, whose 1e-9 floor
then removed the direction. `dev/setier-collin.R` (raw polynomials of
degree d on x in [1, 2], 200 rows, seed d; against lm() and lme4):

    degree  eigenvalue ratio  base      gate alone   with the probe
    3       2.7e-07           kept      kept         kept
    4       2.0e-09           kept      kept         kept
    5       1.1e-11           kept      LOST 6       kept
    6       6.3e-14           kept      LOST 7       kept
    7       4.2e-16           lost 8    lost 8       lost 8
    with (1 | g): degree 5 kept / lost 3 / kept; 6 kept / lost 4 / kept;
    7 lost 5 on all three

Base's standard errors are lm's times sqrt((n - p) / n) (0.9849 at
degree 5) and lme4's to 0.4 percent at degree 6 with (1 | g). With
(1 | g) at degree 6 the finite-difference eigenvalue is 2.1e-13 of the
largest against an asymmetry noise of 1.9e-12 (`dev/setier-deg6.R`),
so the noise bound cannot decide it. The probe does: one step each way
along the direction, sized so that the quadratic model loses 2 *
grad_tol; real curvature loses at least grad_tol both ways, a ridge
loses nothing. It runs only on directions tier 3 would otherwise
remove, and only above 100 machine epsilons of the largest eigenvalue,
where the eigenvalue still has digits. At degree 7 it is below that,
and base lost those standard errors too.

The probe also keeps a near-collinear pair that the 1e-9 floor alone
removed (lane nanse's "at 1 - 1e-9 both paths lose it, as designed"):
`dev/setier-b1-rho.R` (nanse's `-rev2-b1-rho.R` with this lane's
arms; `b1rho-*.txt`), the c + dd ridge beside a ~ 0 + x1, b ~ 0 + x2 at
cor(x1, x2) = 1 - 1e-9: base loses a_x1 and b_x2 on the exact and the
finite-difference path; the lane keeps both at the reference standard
error 1268 (3.8e-6 and 1.2e-6 relative) and loses only c and dd. No
test asserted the old loss. Lane fixes' identification warning
(nl_flat_message(), threshold 1e-9 without a probe) still names a_x1
and b_x2 there, so the two checks now disagree on that one design.

#### The verdict across the BLAS builds, all eight suites

`dev/setier-verdict-diff.R` compares, per test file, the multiset of
standard-error conditions raised (boundary messages, SE warnings,
separation warnings), each cut to the parameters it names, from the
`.cond` logs that `dev/setier-run1.R` writes by tracing frm_warning()
and frm_message() (so muffled conditions count):

    base-ref vs base-ob: conditions 125 vs 124; files differing 10;
    conditions not matched 19
    lane3-ref vs lane3-ob: conditions 463 vs 465; files differing 6;
    conditions not matched 8

What is left differs because the FIT differs, not the verdict on one
Hessian. `dev/setier-blasdiff.R` (test-lkj.R's `homcs(t + 0 | g)`,
seed 79): the correlation parameter stops at 14.94 with the reference
BLAS and at 14.47 with OpenBLAS (logLik -256.51448012 and
-256.514480141), its row is 3.1e-5 and 4.3e-5, and only the first is
removed. The weibull softplus fit of test-numerical-robustness.R loses
shape on OpenBLAS in base and lane alike. A parameter running off along
a flat direction stops where the platform's rounding lets it, which is
`dev/rtmb-pitfalls.md` item 21 in another form; no rule on the Hessian
at two different points can make their verdicts agree.

### Item 1: the message level, measured

`dev/setier-singular.R`, 100 seeds each, against lme4's isSingular()
on the same data (lmer with REML = FALSE, glmer):

    arm        design     n lme4s  lost       msg  warn  conv  se100
    base-ref   ri6      100    72     0    0 ( 0)     0     0     72
    base-ref   ri20     100    51     0    0 ( 0)     0     0     51
    base-ref   ri20s    100    18     0    0 ( 0)     0     0     18
    base-ref   rs20     100    69     0    0 ( 0)     0    69      7
    base-ref   bin15    100    57     0    0 ( 0)     0     0     58
    base-ref   bin15s   100     2     0    0 ( 0)     0     0      2
    base-ob    ri6      100    72     0    0 ( 0)     0     0     72
    base-ob    ri20     100    51     0    0 ( 0)     0     0     51
    base-ob    ri20s    100    18     0    0 ( 0)     0     0     18
    base-ob    rs20     100    69     0    0 ( 0)     0    69      7
    base-ob    bin15    100    57     0    0 ( 0)     0     0     58
    base-ob    bin15s   100     2     0    0 ( 0)     0     0      2
    lane-ref   ri6      100    72    69   69 (69)     0     0      3
    lane-ref   ri20     100    51    44   44 (44)     0     0      7
    lane-ref   ri20s    100    18    14   14 (14)     0     0      4
    lane-ref   rs20     100    69     0    0 ( 0)     0    69      7
    lane-ref   bin15    100    57    56   56 (56)     0     0      2
    lane-ref   bin15s   100     2     2    2 ( 2)     0     0      0
    lane-ob    ri6      100    72    69   69 (69)     0     0      3
    lane-ob    ri20     100    51    44   44 (44)     0     0      7
    lane-ob    ri20s    100    18    14   14 (14)     0     0      4
    lane-ob    rs20     100    69     0    0 ( 0)     0    69      7
    lane-ob    bin15    100    57    56   56 (56)     0     0      2
    lane-ob    bin15s   100     2     2    2 ( 2)     0     0      0

So with a warning, 185 of these 600 ordinary fits (31 percent), 185 of
the 269 that lme4 calls singular, would warn; lme4 gives the 269 a
message. The boundary report is a message. It fires on 0 of the 331
fits lme4 does not call singular. Base, on both BLAS builds, said
nothing on all 600, and gave a finite theta standard error above 100
(on the log scale) on 208 of them.

The rs20 design (`y ~ x + (1 + x | g)`, slope sd 0) is not this lane's:
on 69 of 100 seeds lme4 is singular (seed 1: correlation 1.000), and
on the same 69 frmtmb stops with nlminb code 7, "singular
convergence", and warns
"Optimizer did not report convergence" on base and lane alike
(`dev/setier-rs20.R`, seed 1: logLik -187.0055 in both, lme4's message
"boundary (singular) fit"). That is a convergence verdict on a boundary
fit, in lane optima's territory; filed below.

Smooths. `dev/setier-sx.R` (lane fixes' gamSim design, 140 fits):

    arm        fits  lost  warned messaged    fixef SE NaN
    base-ref    140     7       7        0               0
    base-ob     140     5       5        0               0
    lane-ref    140    86       0        0               0
    lane-ob     140    86       0        0               0

On base, 7 fits (ref) and 5 (OpenBLAS) warned "Standard errors are not
available" for a smoothing sd, by the sign of noise. In the lane 86 of
140 lose a smoothing sd's standard error, and none is told, as mgcv
tells nobody; summary() lists them. A message on 86 of 140 GAM fits
would be noise. One fit differs between the BLAS builds (seed 13,
`y ~ s(x0) + s(x1) + s(x2) + s(x3)`): the s(x1) sd stops at log -8.08
with ref and -7.81 with OpenBLAS (`dev/setier-sx13.R`), the same kind
of difference as above.

A gp()/hsgp() block whose hyperparameters are all flat (an hsgp length
scale at exp(-11.9), where only sd^2 * length scale matters,
`dev/setier-gpprobe.R`) is a limit of the term too, which no single
coordinate's probe sees, and is silent. Two flat hyperparameters of
three (a non-isotropic gp's sd and z length scale trading off when y
does not depend on z, `dev/setier-gpbyprobe.R`) keep the warning.

The edge probe's step is 2 on the internal scale (a factor 0.14 on an
sd). 20 left the region where a quadrature objective is accurate: on
test-quadrature-defects.R's nested beta fit an sd at exp(-11.2) moved
by -5 changed the objective by 1.8, by -20 gave NaN, by -2 by 3e-9
(`dev/setier-quadprobe.R`).

The size of the sd is not the test: a cut at diagnose()'s sd < 1e-4
called a binomial group sd of 1.35e-4, which lme4 reports singular, an
unidentified parameter (`dev/setier-singular.R`, bin15 seed 57, first
build). diagnose() now also lists, as singular, every sd the boundary
message named (`singular_from_se()`).

### Item 2: the nonlinear ridge with random effects

`dev/setier-probe.R`, `y ~ a + b, a ~ 0 + f, b ~ 1 + (1 | g)`, seed 1,
k = 10 (the review's R2):

    arm        SE of a_f*      prediction of a: SE    warnings
    base ref   827,056         822,571                identification only
    base ob    428,568         427,801                identification only
    lane       NaN (11 lost)   NaN, one prediction warning
                               on both BLAS builds

sigma and the g sd keep 0.096225 and 0.307799 in all arms. The 13
outer parameters exceed se_check_np_free and the check waited, so on
base sdreport() built the Hessian and sdr_rescue() had none to look at;
autoscale_sdreport() now builds it with its noise. Lane fixes'
identification warning is the family-level explanation
(`se_explained = "nl_flat"`), so the SE check adds no warning.

### Item 3: separation

`separation_check()` asks the fit for a certificate. The observations
whose fitted probability is still away from both ends hold the estimate
in place, so the runaway direction lies in the null space of their
design rows; the projection of the estimate onto that space (and each
basis vector, both signs) is checked as the definition reads: no
observation moved against its outcome, a row with both outcomes among
its trials not moved, at least one moved. Three cut points for "away
from the ends" (1e-3, 1e-6, 1e-9) are tried; the certificate is exact
at each. It is checked whatever the optimizer's code; when the code is
not 0 the separation warning replaces "Optimizer did not report
convergence", which it explains, and on a converged fit the SE check
still runs for the other parameters, with the predictor's coefficients
explained.

`dev/setier-sep.R`, 50 seeds per design, against Konis's linear
program for separation (brglm2's detect_separation(), not installed
here) solved as the projection of sum(A) onto the cone {b: A b >= 0}
with quadprog (exact: the projection is non-zero exactly when some b in
the cone moves an observation):

    design       n LPsep  named  named  other code!=0
                           base   lane   lane    lane
    plain       50     0      0      0      0       0
    strong      50     0      0      0      0       0
    rare        50     0      0      0      0       0
    cells       50     3      0      3      0       0
    trials      50     0      0      0      0       0
    glmm        50     0      0      0      0       0
    complete    50    50      4     50      0      46
    quasi       50    50      0     50      0       0
    quasi_x     50     9      0      9      0       5
    glmm_quasi  50    50      0     50      0       0

Per BLAS build (the two agree on every verdict): 500 fits, 162
separated by the LP, all 162 named; 0 of the other 338 named. That
covers complete, quasi-complete, factor-level and random-intercept
separation. glm() itself warns "fitted probabilities numerically 0 or
1" on all 50 complete fits but on none of the 100 factor-level
quasi-complete ones. Base named 4 of 162 (code-0 complete fits,
through the SE warning).

"Early": the fit still runs to the optimizer's own end (median 1999 to
2000 evaluations on complete separation, about 0.05 s at n = 50); the
name now reaches the user at that end on every platform. Stopping the
optimizer early when the certificate holds would change the reported
estimates of every separated fit and was not done.

The default-budget repro (seed 514): with OpenBLAS nlminb stops at code
9 and the user now reads "The data separate the outcomes (complete
separation): along a combination of mu: (Intercept), z every
observation that moves (240) moves toward its own outcome ... The
reported values are where the optimizer stopped (function evaluation
limit reached without convergence (9))"; with the reference BLAS the
same sentence, without the parenthesis, and no SE warning.

### Item 4: ranef(condVar = TRUE)

On the gr(g, by = f) fixture (seed 11) condsd^2 now equals the b block
of get_joint_cov()'s diagonal to 1e-10 (sorted, since the long form is
ordered by term). On rellib-r6 the ratio ran 0 to 99 and the test fails
(both BLAS builds). A fit that lost nothing keeps sdreport()'s own
conditional variances, `identical()`.

### Item 5: closed before this round

- `dev/setier-sx.R` on rellib-r6, both BLAS builds: 140 fits, 0 with a
  non-finite fixef() standard error (`sx-base-*.txt`). Lane fixes saw
  3 and 5 on 0.67.0; lane nanse's check closed it.
- `dev/nanse-item3.R` on rellib-r6 and the lane (`item3-base.txt`,
  `item3-lane.txt`, identical warnings and SEs): `a + b + log(c0)` with
  c0 at 5.86e-5 names a_*, b_* and c0 (flat) and keeps sigma's 0.05;
  the curved ridge names a_(Intercept) and c0.

### Composition with the other warnings

- A boundary component that a structural check already named
  (check_re_structure(): a one-level factor) gets no message: it is in
  `se_explained_pars()`. With `check_nlev_1 = "ignore"` the message is
  given (test-diagnostics-ux.R now asserts it).
- A family's own verdict (`fit$cache$se_explained`, as lane fixes'
  nl_flat warning) explains the boundary parameters too.
- The separation warning explains every coefficient of its predictor
  (`fit$cache$se_explained_sep`); on complete separation the whole
  predictor's Hessian vanishes.
- VarCorr() warns only for entries that move along a lost direction
  that is not a boundary one (`hyp_prop_var()`'s `warn`);
  confint_varcorr()'s "No interval" warning skips blocks the boundary
  message named (`se_boundary_blocks()`); predict()'s draws hold the
  boundary parameters at their estimates instead of falling back to the
  plug-in with a warning when every lost parameter is a boundary one;
  fit_fd_se() (fitted() on ordinal and categorical fits) reads the
  propagation covariance and gives NaN with one warning only to rows
  that move along a non-boundary lost direction.

`dev/setier-predcheck.R` on the sleepstudy boundary fit, warnings per
call (`predcheck.txt`; predict() with ndraws = 200, VarCorr(),
confint_varcorr(), confint(), conditional_effects(), fitted(),
summary()):

    base ref   confint_varcorr() 1 ("Uninformative interval": the
               noise SE), the rest 0
    base ob    predict() 1 ("covariance ... could not be recovered"),
               VarCorr() 1, confint_varcorr() 1 ("No interval"), the
               rest 0
    lane ref   0 on all seven
    lane ob    0 on all seven

### Suites

All eight suites, one file per process, NOT_CRAN,
FRMTMB_BRMS_FIT_TESTS, FRMTMB_DRMTMB_FIT_TESTS and FRMTMB_FUZZ on
(`dev/setier-suite.sh`; the seven test-scale.R files skip without
FRMTMB_SCALE_TESTS, as in the release tiers). Every `lib:` line of the
lane runs (355 of 355 per BLAS) shows `wt-setier-lib` for frmtmb.
pass/fail/error/skip/escaped warnings per package
(`dev/setier-suite-sum.R`):

    run        package          p/f/e/s/w       files
    base-ref   frmtmb           17247/0/0/0/0   215
    base-ref   frmtmb.coupling    542/0/0/5/0    11
    base-ref   frmtmb.eam        1743/0/0/3/0    29
    base-ref   frmtmb.latent      360/0/0/2/0    10
    base-ref   frmtmb.learn       501/0/0/2/0    15
    base-ref   frmtmb.ode         549/0/0/1/0    11
    base-ref   frmtmb.sample     2616/0/0/1/0    48
    base-ref   frmtmb.spline      590/0/0/1/0    15
    base-ob    frmtmb           17241/1/1/0/1   215
    base-ob    frmtmb.eam        1743/0/0/3/2    29
    base-ob    frmtmb.sample     2616/0/0/1/1    48
    base-ob    (the other five as base-ref)
    lane3-ref  frmtmb           17320/0/0/0/0   216
    lane3-ref  (the seven extensions as base-ref)
    lane3-ob   frmtmb           17313/2/1/0/1   216
    lane3-ob   frmtmb.eam        1743/0/0/3/2    29
    lane3-ob   frmtmb.sample     2616/0/0/1/1    48
    lane3-ob   (the other five as base-ref)

The reference-BLAS lane run is clean. The OpenBLAS failures and
warnings are the emulator-only ones of rellib-r6 (test-cumulative-cs.R
error, test-bcm-latent-mixtures.R failure, test-ordinal-mixture.R's
"singular convergence (7)", eam test-sampling.R 2, sample
test-brms-shapes-draws.R 1), plus test-perf.R's wall-clock bound under
the load of the parallel lanes (1.78 s against 1.00; 3 of 3 pass
alone, `perf-ob-rerun.txt`; the backlog already lists that test).

Conditions raised in the suites (traced, muffled or not), and the
number of files they fire in:

    run        boundary msg   SE warning   separation   prediction
    base-ref       0  (0)      125 (31)       0 (0)        29
    base-ob        0  (0)      124 (35)       0 (0)        29
    lane3-ref    235 (56)      192 (31)      36 (5)        31
    lane3-ob     236 (55)      193 (32)      36 (5)        31

### Cost

`dev/setier-cost.sh` (`dev/setier-cost.R`): base and lane alternated in
fresh processes, five rounds, each model fitted in a block grown past
1.2 s; seconds per fit, the minimum over rounds (`cost.txt`):

    model                          base     lane     lane/base
    nested (boundary, tier 3)      0.0709   0.0831   1.17
    slope (healthy, fit-time)      0.0900   0.0813   0.90
    glmm_f30 (check waits)         0.3050   0.3150   1.03
    control (no random effects)    0.0213   0.0184   0.86

The machine ran 74 to 92 other R processes from the parallel lanes
throughout, and the control, whose code path the lane does not change,
reads 0.86: differences inside about 15 percent are not resolved here.
What the lane adds is counted work. A healthy fit: one eigen
decomposition of the outer Hessian (values only) and, on a fit with
random effects whose check waited, nothing more, since
autoscale_sdreport() builds the Hessian that sdreport() used to build.
A fit that tier 3 looks at: one gradient (the bound check), one full
eigen decomposition, one objective pair per direction it probes (line
or curvature), and two objective evaluations per flat covariance
parameter (the edge probe).

### Test changes, file by file

New: `tests/testthat/test-se-tier.R` (10 tests, above), and
`allow_boundary()` in `tests/testthat/helper-warnings.R`, which muffles
the boundary message alone (optionally requiring it).

Changed, because a verdict changed on purpose (each is a fit whose
noise SE is gone, or whose report changed level):

- `test-se-check.R`: "a small loading on a downward direction" (the
  gr(g, by = f) slope sd is now a boundary component: the message, no
  SE warning, x's SE as before); "separated data are named as
  separation" (one warning, the separation one, no SE warning; the SE
  reasons are all "separation"); "a lost group sd does not take the
  other variance components" (g2 is a boundary block: no VarCorr()
  warning; an OLRE fit under `check_olre = "ignore"`, whose sd is flat
  but not at an edge, keeps it).
- `test-diagnostics-ux.R`: `check_nlev_1 = "ignore"` now gives the
  boundary message (asserted, no warning); the OLRE test's "ordinary
  grouping factor" line fits a g with no variance (boundary message
  allowed).
- `test-dates.R`, `test-edgecases.R`: data with no group variance; the
  boundary message is allowed, every other message still fails.
- `test-car-spde.R`: bym2's mixing parameter sits at its limit; its
  confint_varcorr() row has NA bounds and every other row brackets its
  estimate.
- Warnings that are now raised on every platform, wrapped in
  allow_warnings() with the reason in a comment: `test-brms-likelihood.R`
  (a mo() interaction simplex weight at 0), `test-covstruct.R` (homdiag
  against sigma on rho = 0 data), `test-famgaps.R` (a cox() baseline
  hazard coordinate), `test-open-issues.R` (a mo() simplex on a rising
  slope, concave), `test-gp-by.R` (a non-isotropic gp's sd and length
  scale trading off), `test-ordinal-mixture.R` (fitted() rows along the
  collapsed thresholds), and frmtmb.learn's `test-families.R` (the
  Kalman filter's sigma0, flat).

### Not done, and why

- An optimizer stop at the separation certificate (see item 3).
- A rule that makes a fit that stops at different points on two
  platforms reach one verdict (item 1, last table). The verdict on a
  given Hessian is platform-independent now; the point is not.
- The boundary message does not carry the sd's value; diagnose() does.

### Defects found, not fixed

- **A random-slope fit whose correlation goes to +-1 warns
  "Optimizer did not report convergence: singular convergence (7)"**
  where lme4 gives its boundary message: 69 of 100 seeds of
  `y ~ x + (1 + x | g)` with slope sd 0 (`dev/setier-singular.R`
  rs20; `dev/setier-rs20.R` seed 1, the same logLik as lme4). Lane
  optima's area (convergence on the fit path); filed for it.
- **summary() prints "(Number of levels: NA)" for a nested term**
  (`(1 | Subject/a)`: the `Subject:a` block), on rellib-r6 and the lane
  alike.
- **test-bcm-latent-mixtures.R fails with the OpenBLAS emulator** on
  rellib-r6 and the lane alike (pass 51, fail 1), besides the four
  emulator-only files the cifix review lists.

### Files touched

- Core R: `R/se-check.R` (cov_from_hessian(), se_tier3(),
  se_fit_analysis(), sdr_rescue(), se_report(), se_relabel(),
  se_explained_pars(), se_lost_clauses();
  new se_tier3_could_act(), se_curvature_real(), se_boundary_act(),
  se_boundary_names(), se_at_edge(), se_edge_step,
  se_boundary_message(), se_boundary_labels(), se_boundary_blocks(),
  se_par_grouped(), separation_check(), sep_links, sep_certificate(),
  sep_null_space()); `R/autoscale.R` (autoscale_sdreport());
  `R/fit.R` (check_convergence(), and the roxygen of frmtmb_control()'s
  `check_se`; fit_assembled() is NOT edited); `R/confint.R`
  (hyp_par_cov(), hyp_prop_var(), new se_null_boundary(),
  confint_varcorr(), diagnose_singular(), new singular_from_se());
  `R/methods-fit.R` (ranef.frmtmb_fit(), VarCorr.frmtmb_fit());
  `R/brms-shapes.R` (fit_draw_space(), fit_fd_se());
  `R/predict-brms.R` (predict_par_drawer()). `man/frmtmb_control.Rd`
  regenerated.
- Core docs: `NEWS.md` (development-version section).
- Core tests: new `test-se-tier.R`; `helper-warnings.R`
  (allow_boundary()); `test-se-check.R`, `test-diagnostics-ux.R`,
  `test-dates.R`, `test-edgecases.R`, `test-car-spde.R`,
  `test-brms-likelihood.R`, `test-covstruct.R`, `test-famgaps.R`,
  `test-open-issues.R`, `test-gp-by.R`, `test-ordinal-mixture.R`.
- frmtmb.learn: `tests/testthat/test-families.R` (one allow_warnings()).
- dev: `dev/setier-*` scripts and this file.

Overlaps with the parallel lanes: `R/fit.R` (check_convergence(), and
frmtmb_control()'s documentation; lane optima edits the fit path),
`NEWS.md`, `R/confint.R`, `R/predict-brms.R`.

### R CMD check --as-cran

`dev/setier-check.sh` (built and checked in `dev/setier-check/`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
R_LIBS the lane library, rellib-r6, the user library):

- frmtmb: Status: 2 NOTEs. Examples timing (frm_hazard_reads 7.44 s
  elapsed against 1.70 s user, dharma_residuals 5.44 against 1.36: the
  load of the parallel lanes) and the V8 math-rendering note. Tests
  `[ FAIL 0 | WARN 0 | SKIP 376 | PASS 11544 ]`.
- frmtmb.learn: Status: 1 NOTE. Examples timing (bandit2arm_delta
  14.00 s elapsed against 1.51 s CPU). Tests
  `[ FAIL 0 | WARN 0 | SKIP 20 | PASS 404 ]`.

`man/frmtmb_control.Rd` rendered with Rd2txt: the new `check_se`
paragraph reads as written.

### Version

Core: a minor bump (behavior changes marked BREAKING: verdicts and
message levels of the standard-error check, separation naming).
frmtmb.learn: a test-only change, no bump and no new floor: the
allow_warnings() wrap does not require the warning, so the test passes
against 0.68.1's core and this one alike.
