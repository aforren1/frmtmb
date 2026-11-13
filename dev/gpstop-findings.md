# Lane gpstop: exact gp() fits that stop short

Round opened 2026-11-11, base e8abed0a (frmtmb 0.81.0). Worktree
`frmtmb-wt-gpstop`, branch `wt-r82-gpstop`, library `wt-gpstop-lib`;
"base" is `rellib-r20`. Scripts are `dev/gpstop-*.R` and
`dev/gpstop-*.sh`; logs are in `dev/gpstop-out/` (gitignored). Every
number below is pasted from the log named beside it.

Status: fix round after the review (dev/reviews/2026-11-11-gpstop.md,
FIX) done; installed in wt-gpstop-lib. The section "Fix round" holds where it disagrees with the sections
after it, which are the first submission's record. Nothing committed.

## Items

| backlog line | item | result |
|---|---|---|
| 4822 | `a * exp(b * z)` with a field in `a` or `b` stops short, silent | fixed: the body is linearized again at the fit's optimum, and its descent kept only at a regular maximum |
| 5033 | `student()` on a step field, 11 and 6.07 short, silent | open, re-described as a limit of the Laplace approximation (fix round, item 2); no code change |
| 4485 | `gp(x, k = )` too small, silent | not fixed: not this change (see "Not done") |
| 4491 | exact field positions not named | not fixed: not this change |
| 4388 | SE check's Hessian from the marginal's tape | not fixed: not this change |

## Fix round

The review (`dev/reviews/2026-11-11-gpstop.md`): 3 MAJOR, 5 MINOR. Each
finding was reproduced on the first lane build before the fix
(`dev/gpstop-fx-repro.R`, the review's helpers and data, arms lane and
norelin in one process; log `fx-repro-old.txt`), then rerun on the fixed
build (`fx-repro-new.txt`):

    row                     first build                           fixed build              norelin (0.81.0's route)
    logis_a student 40 2    27.822186 c1 w2 g 1.2e13 inner 1.97e-7  22.393662 c0 w0 inner 1  22.393662 c0 w0 inner 1
    logis_a bernoulli 40 2  -20.567787 c0 w0 log lscale -323.267    -24.511135 c0 w0          -24.511135 c0 w0
    pow_b student 40 1      9.270671 c1 w2                          9.532225 c1 w2            9.532225 c1 w2
    mm_K bernoulli 120 1    -79.035244 c1, 4 warnings               -79.035250 c1, 2 warnings -79.035250 c1, 2 warnings

(c = code, w = warnings, g = max outer gradient, inner = smallest
eigenvalue of the inner Hessian at the optimum.) The two MAJOR rows and
MINOR 8's row now return 0.81.0's point exactly, and MINOR 4's row has
0.81.0's two warnings.

### 1. MAJOR 1 and 2: a relinearized descent is kept only at a regular maximum

`gp_relin_regular()` (`R/gp-start.R`), asked in `gp_second_descent()`
after `better()`. A relinearized descent is kept only if it ended at code
0, its largest absolute gradient is within `grad_tol`, the smallest
eigenvalue of the inner Hessian exceeds `gp_relin_inner_ev_min = 1e-2`
(`gp_inner_ev_min()`), and every length scale lies between a tenth of the
smallest nonzero distance between positions and ten times their range.
Otherwise the point before the descent is kept and the rounds stop. The
second starts that 0.81.0 already had are unchanged.

The cut, read from the numbers (`dev/gpstop-fx-sweep.R`, the review's
sweep with every call of the check logged; logs `fx-sweep-gauss.txt`,
`fx-sweep-pb.txt`, `fx-sweep-student.txt`, lines `REG`): 20 descents rose
and were checked. The 14 kept read an inner eigenvalue of 0.997 to 1.000,
a gradient up to 2.04e-4, and length scales 22.5 to 740 times the
smallest gap and 0.0136 to 0.125 of the range. The 6 refused:

    REG logis_a_bernoulli_n40_s2 conv 0 gmax 8.45e-07 ls/mingap 2.89e-139 ls/range 4.18e-142 inner_ev 1.000e+00 verdict FALSE
    REG mm_K_bernoulli_n40_s2 conv 1 gmax 8.34e+00 ls/mingap 0.0268 ls/range 2.01e-05 inner_ev 1.000e+00 verdict FALSE
    REG mm_K_bernoulli_n120_s1 conv 1 gmax 1.79e+09 ls/mingap 2.43e+03 ls/range 0.25 inner_ev 2.991e-04 verdict FALSE
    REG pow_b_bernoulli_n120_s1 conv 1 gmax 1.33e+06 ls/mingap 1.05e+03 ls/range 0.109 inner_ev 5.195e-04 verdict FALSE
    REG pow_b_student_n40_s1 conv 1 gmax 2.94e+06 ls/mingap 231 ls/range 0.142 inner_ev 5.159e-04 verdict FALSE
    REG logis_a_student_n40_s2 conv 1 gmax 8.83e+09 ls/mingap 224 ls/range 0.023 inner_ev 1.906e-05 verdict FALSE

So the eigenvalues fall in two groups, 0.997 to 1 and 1.9e-5 to 5.2e-4
(the review's spikes: 4e-8 to 7e-4); 1e-2 sits a factor 20 above the
second and 100 below the first. Two things to know about the tests: on
these rows no descent was refused by the eigenvalue alone (each of the
four low ones also ended at code 1), and MAJOR 2 is caught only by its
length scale (code 0, gradient 8.5e-7, inner eigenvalue 1.0). The
eigenvalue test stays for a code-0 spike, which this sweep did not reach.
A refusal costs only the gain: the fit keeps the point 0.81.0 returns.

**The review's sweep, 122 rows, lane against norelin interleaved**
(summary by the review's `gpstop-review-sum.R`, output
`fx-sweep-sum.txt`):

    rows 122 | relin ran 115 | relin read a new scale 68 | left closure 0
    up >1e-6 10 up >1e-3 10 down < -1e-6 0 max up 173.4096
    warn or code changes 0
    suspect point (evmin<1e-4, gmax>1e-2 or |theta|>8): lane 26 norelin 26 lane-only 0
    wall lane 4172.3 norelin 4159.9 ratio 1.003 | rows relin hit: ratio 1.015 | rows no hit: ratio 0.983

Every row that moves, all up, each at code 0 with no warning: exp_a
gaussian n400 +173.4096, exp_a student n120 s2 +37.1216, exp_a gaussian
n120 s2 +37.0925, pow_a gaussian n40 s2 +21.5264, mm_K gaussian n400
+20.4866, mm_V student n120 s2 +7.4023, exp_b student n120 s2 +5.5676,
pow_b student n40 s2 +4.6511 (confirmed by the review's importance
sampling), mm_K gaussian n40 s2 +0.4792, mm_V poisson n40 s1 +0.3105.
Against the first build (the review's logs), 6 rows no longer move: the
two MAJOR artifacts (+3.94, +5.43), MINOR 8 (-0.26, now 0), and three
small gains on fits that end warned at code 1 (exp_a student n40 s2
+0.33, mm_K bernoulli n40 s2 +0.83 at log length scale -8.58, pow_b
bernoulli n120 s1 +0.32). Those three keep 0.81.0's warned point.

**fitedge's 192 rows** (`dev/gpstop-gpfresh.R <set> lane` on the fixed
build; `gpsum-fx.txt`): the same table as the first build, fresh short
2 to 0 and lane96 6 to 5, 3 up, 0 down; each of the 192 log-likelihoods
equals the first build's to 1e-6 (join of `gpf-*-lane1.tsv` and
`gpf-*-lane.tsv`: 0 rows differ). Interleaved against norelin
(`dev/gpstop-cost-body.R`, `bcost-fx-*.txt`):

    BCOSTSUM fresh rows 96 lane 389.3 s norelin 329.0 s ratio 1.18
    BCOSTSUM lane96 rows 96 lane 256.6 s norelin 205.9 s ratio 1.25

with the same 4 moves as before (E wig 120 2 +36.50, L wig 120 2 +5.24,
D step 60 1 +0.0016, D step 60 3 +0.00013), all up.

### 2. MAJOR 3: 5033 stays open

The review's st80 importance-sampling segment (71.0 to 71.3 at the
default, 73.5 to 73.8 at the spike) and my 73.94 against 71.33
(`stis2-base.txt`) agree: the model's marginal there is about 2.5 above
the default, which stops at code 0 with no warning and a healthy inner
Hessian (0.999). No check at the stop separates it from a correct fit.
The NEWS Limits bullet now says the default is the best regular maximum
of the Laplace objective and that the marginal can be higher, with no
version number; the CLOSED line under 5033 is removed and the item is
filed open, with these numbers, under "Filed in the round of
2026-11-11". My item 2 below said "misdiagnosed"; that holds for the
optimizer on the Laplace objective and not for the model on st80.

### 3. MINOR 4: raw nlminb warnings

The relinearized descent runs in a handler that muffles a warning
containing "NA/NaN", as `fit_error_context()` does for the first
descent of a nonlinear fit with no start (relin runs only there).
mm_K bernoulli n120 s1: 2 warnings, as norelin (`fx-repro-new.txt`). The
second starts of 0.81.0 are not changed.

### 4. MINOR 5: NEWS

The Fitting bullet is marked BREAKING and names the regularity tests,
where the step runs and where it does not; the Limits bullet has no
version number. Cost as measured: on the review's 122 rows 1.003 in total
(4172.3 s against 4159.9 s; rows with a new scale 1.015, without 0.983);
on the 192 smaller rows 1.18 and 1.25; on the n = 400 row with no new
scale (`pow_a`), the review measured 1.16 (270.2 s against 232.5 s) and
this run 0.99 (227.75 s against 229.66 s) under load. NEWS says "about
as long in total", "1.2 times" on the small ones and "up to 1.16" on
that row. The bump is left to the consolidation.

### 5. MINOR 6: where it runs

`fit_assembled()` (`R/fit.R`, one new line `gp_cold <- is.null(template)`
after `prefit_seeded <- FALSE`, and one argument at the `gp_start_scan()`
call) passes `relin = gp_body && announce_start && (gp_cold ||
prefit_seeded)`, nlstop's cold-start rule (`cold_tpl <- is.null(template)
|| prefit_seeded`, asked with `announce_start`). So it runs on the user's
own fit from its cold start, including the fit that follows an autoscale
pre-fit, and on `update()`, which evaluates a new `frm()` call; it does
not run in the pre-fit itself, on `refit()`, `influence()`, `anova()`'s
refits or `frm_allfit()` (each calls `fit_assembled()` with
`announce_start` at its default FALSE, and all but `frm_allfit()` with
`start_from_mle = FALSE` hand a template), nor on the probe and saddle
refits (`run(st, FALSE)`). Test:
"the body is linearized again only for its field, with no start" now
also asserts `update()` runs it and `refit(newresp =)` does not. My lines
in `R/fit.R` are separate from nlstop's hunk (which adds `cold_tpl` at
the `make_start()` call); the consolidation may fold `gp_cold` into
`cold_tpl`.

### 6. MINOR 7, MINOR 8

MINOR 7 is filed (approximate gp() in a nonlinear body, degenerate length
scale and one row 1.18 below a theta start); the gate is not widened.
MINOR 8 (pow_b student n40 s1, -0.26 with both arms at code 1): the
relinearized descent ends at code 1 and is refused, so the fit equals
norelin, 9.532225 (`fx-repro-new.txt`, and the sweep).

### 7. Tests

New: "a second reading of the body keeps only a regular maximum" (both
MAJOR rows; never below the fit with the step mocked off, and where above
it, code 0, a length scale above a tenth of the smallest gap, and an
inner eigenvalue above 0.01 of the reference fit's; `skip_on_cran()`).
Changed: the gate test checks `update()` and `refit()`.

    test-gp-start.R lane (t-gp-start-lane-fx.txt)   pass=247 fail=0 err=0 skip=0 warn=0
    test-gp-start.R base (t-gp-start-base-fx.txt)   pass=238 fail=2 err=2 skip=0 warn=0
    mutant regular  (gp_relin_regular() TRUE, the first build; fx-mut-regular.txt)   pass=250 fail=3
    mutant template (relin on a fitted template; fx-mut-template.txt)               pass=246 fail=1

On rellib-r20 the failures are item 1's (35.8 <= 72.3, 100.4 <= 105.6)
and the errors the missing `gp_body_relin_applies` and `gp_body_relin`.
The new test cannot fail on rellib-r20, which has no second reading; it
fails on the first lane build: code 1 against 0, inner eigenvalue
0.0000002 <= 0.01, and a length scale of 4e-141 <= 0.0014 (mutant
`regular`, `dev/gpstop-fx-mut.R`). The mutant `template` fails the
`refit()` check (`n_relin` 1, not 0).

### 8. Runs on the fixed build

- Core suite, one file per process, 3 at a time, gated tier
  (`bash dev/release/run-par.sh gpstop-fx-core dev/gpstop-core.jobs 3`,
  `FRMTMB_RELLIB=wt-gpstop-lib`; log `dev/release/gpstop-fx-core.log`,
  "GPSTOP-FX-CORE ran 317 of 317"; `gpstop-out/suite-fx-sum.txt`):

        frmtmb files 317 pass 26546 fail 0 err 0 skip 15 warn 0

  317 of 317 `lib:` lines name wt-gpstop-lib for frmtmb; the 15 skips
  are the drmtmb, scale and fuzz gates as before. The extensions were
  not rerun: the fix touches only paths that need an exact gp() in a
  nonlinear body, which no extension test fits (first submission,
  "Suites"; all eight passed there).
- Sessions of the gp files (the 5 `test-gp-*.R` and `test-nl-family.R`
  with core's helpers, `dev/release/run-session.R frmtmb <order> <dir>`;
  logs `gpstop-out/sess-gp-*.txt`):

        SESSION frmtmb forward files=6/6 pass=858 fail=0 err=0 skip=0 warn=0
        SESSION frmtmb preload files=6/6 pass=858 fail=0 err=0 skip=0 warn=0
        SESSION frmtmb reverse files=6/6 pass=858 fail=0 err=0 skip=0 warn=0

  LEAK lines 3, untagged 0.
- OpenBLAS 0.3.26 with LAPACK (`openblas-build-0326-fx.txt`: Rblas 74 and
  Rlapack 616 symbols to OpenBLAS, dgemm 2000 in 0.05 s;
  `ob0326-fx-sum.txt`, every `blas:` line the copy, every `lib:` line
  wt-gpstop-lib): test-gp-ar-start 49, test-gp-by 86, test-gp-multidim
  67, test-gp-noncentered 76, test-gp-start 247, test-nl-family 333
  passes, each fail 0 err 0 warn 0. `dev/ciharden-out/` deleted.

- `R CMD check --as-cran` of core (`bash dev/gpstop-check.sh .`, kept as
  `gpstop-out/check-frmtmb-fx.txt`; check tree deleted): tests
  `[ FAIL 0 | WARN 0 | SKIP 272 | PASS 23240 ]` in 29 min; `Status: 1
  NOTE`, the V8 math-rendering one. The tarball was built before a
  comment-only edit of the `relin =` comment in `R/fit.R`; the library
  was reinstalled after it (`install4.txt`).

### 9. After the short check (review, "Short check (fix round)": FIX, S1)

**S1, the low end of the span rule.** Reproduced first: the review's
`gpstop-short-white.R` on the build with a tenth of the SMALLEST gap
returned `regular TRUE` on plain k 2, k 5, k 20 and tie k 2, k 20
(`dev/reviews/gpstop-short-out/short-white.txt`).

The margin, measured before the factor was fixed
(`dev/gpstop-fx-margin.R`, no frmtmb code, `gpstop-out/fx-margin.txt`;
length scale over the median nearest-neighbour gap of the distinct x):

    MARGIN sweep-kept  fits  14 clean  14 | clean ls/mednn min 1.86 1% 1.92 median 21.4 | below 0.1: clean 0, all 0 | below 0.5: clean 0
    MARGINMOVED fe-fresh E wig 120 2 ls/mednn 20.9 ls/range 0.06699
    MARGINMOVED fe-fresh L wig 120 2 ls/mednn 25.8 ls/range 0.06807
    MARGINMOVED fe-lane96 D step 60 1 ls/mednn 1.84 ls/range 0.02034
    MARGINMOVED fe-lane96 D step 60 3 ls/mednn 3.75 ls/range 0.02835
    MARGINWHITE plain ls/mednn 0.141 | tie ls/mednn 0.0789

So every relinearized descent kept on the 122 and the 192 reads 1.84 of
the median gap or more, and the review's white-noise stationary points
read 0.141 and 0.079. The review's suggested factor, a tenth, does NOT
refuse the plain point: the median nearest-neighbour gap of MAJOR 2's 40
positions is 0.0665 (not about 0.2; here 0.27 of range / n, since a
nearest neighbour takes the shorter of two gaps), so a tenth is 0.00665
against the point's 0.0094.
Seen on a build with the tenth (`t-gp-start-lane-tenth.txt`): the new
test fails, `regular(...)` TRUE on the plain point. Shipped: half the
median gap (`gp_relin_ls_min = 0.5`), 3.5 times above the plain point
and 3.7 times below the lowest kept descent. A refusal costs only the
gain. On the shipped build `gpstop-fx2-white.R` (the review's script with
the eigenvalue for its log computed locally; `fx3-short-white.txt`)
reads `regular FALSE` on all six polished points.

New test "a white-noise field is refused at any length scale below the
gaps": both of the review's points (plain k 5, tie k 20) must be refused,
and the fit's own regular maximum must pass (the absent case).

Margin side note (filed): 0.81.0's own first descent ends at a length
scale under a tenth of the median gap, code 0 and no warning, on 5 rows
of the two sets (`MARGINLOW`, both arms equal).

**MINOR, cost of the eigenvalue test.** `gp_inner_ev_min()` is replaced
by `gp_inner_pd(obj, par, cut)`: a sparse Cholesky of `H - cut I`
(`Matrix::Cholesky(LDL = FALSE)`, which warns and does not stop on a
matrix that is not positive definite), so it answers "every eigenvalue
above the cut" without a dense matrix. `dev/gpstop-fx-pdtime.R` (poisson
`y ~ gp(x) + (1 | s)`, 60 positions, built at its start; one process,
arms interleaved, minimum of 3; `fx-pdtime.txt`):

    PDTIME q 500 random 500 | sparse 0.020 s (max MB 178) verdict TRUE | dense 0.010 s (max MB 180) ev_min 1 agree TRUE
    PDTIME q 2000 random 2000 | sparse 0.010 s (max MB 199) verdict TRUE | dense 1.490 s (max MB 272) ev_min 1 agree TRUE
    PDTIME q 4000 random 4000 | sparse 0.030 s (max MB 228) verdict TRUE | dense 13.220 s (max MB 555) ev_min 1 agree TRUE
    PDEDGE q 4000 cut 0.99 ev TRUE (want TRUE) | cut 1.01 ev FALSE (want FALSE)

(MB is the session's peak, not the call's.) Both sides of the answer
agree with the dense eigenvalue at 500, 2000 and 4000 (`PDEDGE`).

**Nit:** the Limits bullet names `frmtmb.sample::frm_sample()`.

**Reruns on the shipped build** (half the median gap, sparse test):

- The review's 122 rows, lane against norelin interleaved
  (`dev/gpstop-fx2-sweep.R`, logs `fx2-sweep-*.txt`, summary
  `fx2-sweep-sum.txt`):

        rows 122 | relin ran 115 | relin read a new scale 68 | left closure 0
        up >1e-6 10 up >1e-3 10 down < -1e-6 0 max up 173.4096
        warn or code changes 0
        suspect point (evmin<1e-4, gmax>1e-2 or |theta|>8): lane 26 norelin 26 lane-only 0
        wall lane 3394.1 norelin 3443.3 ratio 0.986 | rows relin hit: ratio 0.986 | rows no hit: ratio 0.986

  The same 10 rows rise by the same amounts as on the previous fix
  build, 0 fall; every lane log-likelihood equals that build's to 1e-6
  (123 of 123 lane rows, joined by tag). The checks: 14 kept, 6 refused,
  as before; the lowest kept length scale is 1.87 of the median gap.
- fitedge's 192, interleaved (`bcost-fx2-*.txt`): `fresh` 1.18 (296.9 s
  against 252.7 s), `lane96` 1.24 (223.9 s against 179.9 s); the same 4
  rows rise (+36.50, +5.24, +0.0016, +0.00013), 94 of 96 equal at 1e-6
  in each set, 0 fall.
- `test-gp-start.R`: lane `pass=251 fail=0 err=0 skip=0 warn=0`
  (`t-gp-start-lane-fx3.txt`); rellib-r20 `pass=238 fail=2 err=3` (the
  two item-1 failures 35.8 <= 72.3 and 100.4 <= 105.6, and the missing
  `gp_body_relin_applies`, `gp_relin_regular` and `gp_body_relin`;
  `t-gp-start-base-fx3.txt`). The white-noise test fails on the build
  with a tenth of the median gap (`t-gp-start-lane-tenth.txt`: `pass=250
  fail=1`).

- Core suite, one file per process, gated tier, 2 at a time (beside the
  sweep's last process; `bash dev/release/run-par.sh gpstop-fx2-core
  dev/gpstop-core.jobs 2`, log `dev/release/gpstop-fx2-core.log`,
  "GPSTOP-FX2-CORE ran 317 of 317"; `suite-fx2-sum.txt`):

        frmtmb files 317 pass 26550 fail 0 err 0 skip 15 warn 0

  317 of 317 `lib:` lines name wt-gpstop-lib; the skips are the same
  gates as before.
- gp sessions, one at a time (the 5 `test-gp-*.R` and `test-nl-family.R`
  with core's helpers; `sess2-gp-*.txt`): forward, reverse and preload
  each `files=6/6 pass=862 fail=0 err=0 skip=0 warn=0`; LEAK lines 3,
  untagged 0.
- OpenBLAS 0.3.26 with LAPACK, 2 at a time (`OB_N=2 bash
  dev/gpstop-openblas-run.sh 0.3.26`; `openblas-build-0326-fx2.txt`:
  Rblas 74 and Rlapack 616 symbols to OpenBLAS, dgemm 2000 in 0.05 s;
  `ob0326-fx2-sum.txt`, 6 of 6 `blas:` lines the copy and 6 of 6 `lib:`
  lines wt-gpstop-lib): test-gp-ar-start 49, test-gp-by 86,
  test-gp-multidim 67, test-gp-noncentered 76, test-gp-start 251,
  test-nl-family 333 passes, each fail 0 err 0 warn 0. Copy deleted.
- `R CMD check --as-cran` of core (`bash dev/gpstop-check.sh .`, kept as
  `gpstop-out/check-frmtmb-fx2.txt`; check tree deleted): tests
  `[ FAIL 0 | WARN 0 | SKIP 272 | PASS 23244 ]` in 23 min; `Status: 1
  NOTE`, the V8 math-rendering one. `R/` and the tests are those
  installed (`install6.txt`); one NEWS sentence on the span rule was
  edited while the build ran.

The installed build (wt-gpstop-lib, `install6.txt`) holds the final
`R/`. Nothing committed.

## Item 1 (4822): a field in a multiplicative body

### Seen on rellib-r20

`dev/gpstop-probe.R base E120w2,L120w2,E60w2,F120w1`, log
`probe-base-nl.txt`. Rows from `dev/review-fitedge-designs.R`
(`set.seed(52000 + 100 seed + n + 10 design + curve)`); reference: the
best of 15 theta starts (sd 0.3, 1.5, 6 by length scale 0.1, 0.3, 1, 3, 8).

    PROBE E120w2 default ll 35.802033 conv 0 nwarn 0 wall 4.26
    PROBE E120w2 best of 15 theta starts 72.301894 gap -36.4999
    PROBE L120w2 default ll 100.401247 conv 0 nwarn 0 wall 9.22
    PROBE L120w2 best of 15 theta starts 105.637283 gap -5.2360

E is `a * exp(b * z)` with `a ~ gp(x)`; L the same body with `b ~ gp(x)`;
both the wiggle curve `2 sin(3x) + x`, n 120, seed 2. The default and the
best start agree on the coefficients (E: `a` 2.4296 against 2.4199, `b`
0.39584 against 0.39912) and differ in the field: log length scale 1.7095
against -0.40361 (E), 0.79337 against -0.4002 (L). The default took a long
field.

### Cause

`gp_body_start()` linearizes the body at the template, where the body
scales the field by a derivative that is not the fit's. In E, `J_a =
exp(b z)` is 1 at the template's `b = 0`, so the working marginal reads
`y ~ a + gp(x)` with the factor `exp(0.4 z)` as noise. In L, `J_b = a z
exp(b z)` is 0 at the template's `a = 0`, so every row fails the `|J|`
test and there is no body start at all. The descent finds the
coefficients either way, and then sits at a long field.

### Fix

`gp_body_relin()` (`R/gp-start.R`): at the fit's optimum, the body is
linearized again at the fit's own coefficients (`gp_start_working(tpl =
<the fit's parList>, mu = FALSE)`, which is `gp_body_start()` at that
point), and where that reads another length scale (more than a factor 1.1
in some block) a descent runs from the fit's point with that theta. The
better optimum is kept. After a gain it reads the body again, at most
`gp_body_relin_rounds = 3` times. It runs inside `gp_second_descent()`
after the second starts, so `R/fit.R` does not change; `gp_start_scan()`
leaves a closure in `obj$env$.frmtmb_gp_relin`, as it leaves
`.frmtmb_gp_alt`, and `gp_second_descent()` removes both.

It runs only where `gp_body_relin_applies()`: an exact gp() in a parameter
of a nonlinear body, with `beta` on the tape in natural coordinates, and
the caller gave no `start` (`gp_start_scan(body = TRUE)`, fitedge's rule
for the body start). `gp_body_relin()` also checks that the fit's
`parList()` maps back to its own outer vector before it uses it.

On the lane build (`dev/gpstop-probe.R lane ...`, log `probe-lane1-nl.txt`):

    PROBE E120w2 default ll 72.301894 conv 0 nwarn 0 wall 6.06
    PROBE L120w2 default ll 105.637283 conv 0 nwarn 0 wall 11.69
    PROBE E60w2 default ll 19.371871 conv 0 nwarn 0 wall 1.13
    PROBE F120w1 default ll 74.624667 conv 0 nwarn 0 wall 8.44

Each equals the best of the 15 theta starts. (Walls from separate
processes; the cost is measured below.)

### fitedge's 192 rows, rellib-r20 against the lane

`dev/gpstop-gpfresh.R <set> <arm>` (lane fitedge's `fitedge-r2-gpfresh.R`
pointed at this round: arms base, lane, and ref = the linear spelling
where one exists plus 12 theta starts, which turn every warm start off and
so are build-free), sets `fresh` (the fitedge review's 8 designs x 3
curves x n 60, 120 x seeds 1, 2) and `lane96` (fitedge's own 4 bodies x 4
curves x n 60, 120 x seeds 1 to 3). Summary `dev/gpstop-gpsum.R`, output
`gpsum.txt`; the reference is the best of every arm, short is more than
1e-3 below it, silent is no warning:

    GPS fresh all  rows 96 | short base  2 (silent  2) lane  0 (silent  0) | worst base  -36.500 lane    0.000 | lane worse than base 0 (max 0.000), better 2 | time base 346.3 lane 499.7
    GPS lane96 all  rows 96 | short base  6 (silent  6) lane  5 (silent  5) | worst base   -0.635 lane   -0.635 | lane worse than base 0 (max 0.000), better 1 | time base 254.9 lane 391.6
    GPSB fresh E wig 120 2 base 35.8020 lane 72.3019 gain 36.4999
    GPSB fresh L wig 120 2 base 100.4012 lane 105.6373 gain 5.2360
    GPSB lane96 D step 60 1 base -69.6944 lane -69.6928 gain 0.0016
    GPSS lane96 B wig 60 2 lane -118.6294 ref -118.6111 short -0.0183 base -118.6294 warn ''
    GPSS lane96 B wig 60 3 lane -125.0282 ref -124.4175 short -0.6107 base -125.0282 warn ''
    GPSS lane96 B exp 60 1 lane -118.7819 ref -118.7755 short -0.0064 base -118.7819 warn ''
    GPSS lane96 B step 60 3 lane -122.7158 ref -122.0808 short -0.6350 base -122.7158 warn ''
    GPSS lane96 B step 120 1 lane -264.3318 ref -263.9660 short -0.3658 base -264.3318 warn ''

Every row that moves: 3, all up (E wig 120 2 +36.50, L wig 120 2 +5.24,
D step 60 1 +0.0016); 0 rows down. Short: 8 to 5, worst 36.50 to 0.635.
The 5 left are body B, `a + b * z` on poisson, where the body is additive
and the relinearization reads what the body start read; each equals
rellib-r20. Lane fitedge found them short on the linear spelling `y ~ z +
gp(x)` too (dev/fitedge-findings.md, item 1): the working-response limit
of low-count Poisson that 0.79.0's Limits name.
Per-design lines are in `gpsum.txt`. The walls are sums from separate
processes, not side by side; the cost is measured below.

### Cost

`dev/gpstop-cost-body.R <set>`: one process on the lane library, the lane
against the mutant `norelin` (`gp_body_relin_applies()` returns FALSE,
which is 0.81.0's route), interleaved row by row, order alternating. Logs
`bcost-fresh.txt`, `bcost-lane96.txt`:

    BCOSTSUM fresh rows 96 lane 632.8 s norelin 531.2 s ratio 1.19
    BCOSTSUM lane96 rows 96 lane 411.8 s norelin 337.5 s ratio 1.22

The log-likelihoods of the two arms agree to 1e-6 on 93 of 96 and 94 of
96 rows; the rows that differ are the 3 of the table above and lane96 D
step 60 3 (+0.000132, under the table's 1e-3). On most rows the working
marginal at the optimum reads the fit's own length scale and no descent
runs: over the 192 lane fits, `fit$opt$gp_second["runs"]` (second starts
plus relinearized descents) was 0 on 6, 1 on 137, 2 on 48, 3 on 1. So the
cost is mostly the working marginal, about 20 percent of the fit.

### Tests

`tests/testthat/test-gp-start.R`, two new tests:

- "a multiplicative body is linearized again at the fit's optimum" (E wig
  120 2 and L wig 120 2, each against one theta start; `skip_on_cran()`).
- "the body is linearized again only for its field, with no start":
  counts calls to `gp_body_relin()`: more than 0 on a body field with no
  start (the absent case), 0 with a caller's `start`, 0 on `y ~ gp(x)`
  poisson.

On rellib-r20 (`dev/gpstop-run-test.R base frmtmb ...`, log
`t-gp-start-base.txt`): `pass=238 fail=2 err=1`, the failures
`35.8 <= 72.3` and `100.4 <= 105.6`, the error `object 'gp_body_relin' not
found`. On the lane (`t-gp-start-lane.txt`): `pass=243 fail=0 err=0 skip=0
warn=0`.

Roundings (`dev/gpstop-ulps.R lane`, log `ulps-lane.txt`): both
assertions on `y * (1 + 2e-15 k)`, k = 1 to 20: `ULPSUM rows 40 fail 0`;
default minus reference -1.47e-11 to -1.93e-12 (the test allows 1e-5 of
the value).

## Item 2 (5033): student() on a step field

### Seen on rellib-r20

`dev/gpstop-probe.R base st80,st160,st40,gs40,pm80`, log
`probe-base-ng.txt`. Rows from `dev/arstep-ng.R` (`set.seed(74000 + 31
seed + n)`, seed 1), `y ~ gp(x)`, `student()`, `y = 2 (x > 5) - 1 + 0.02
t_5`:

    PROBE st80 default ll 68.550327 conv 0 nwarn 0 wall 1.80
    PROBE st80 best of 15 theta starts 79.545003 gap -10.9947
    PROBEP st80 best rho0 1 sd0 6 conv 1 nwarn 2 | ...
    PROBE st160 default ll 201.293775 conv 0 nwarn 0 wall 10.56
    PROBE st160 best of 15 theta starts 207.365455 gap -6.0717
    PROBEP st160 best rho0 0.1 sd0 1.5 conv 1 nwarn 2 | ...

So the gaps reproduce (10.99 and 6.07), the default is silent, and each
higher value comes from a fit that itself ends at code 1 with two
warnings.

### The higher values are not maxima

1. At each optimum the inner start does not matter (`dev/gpstop-stinner.R
   base ...`, log `stinner-base.txt`): at the default's outer optimum the
   objective reads the same from the fit's inner mode, from zero and from
   the reference fit's mode (68.5503 all three on st80, 201.2938 on
   st160), and at the reference's the same (79.5450, 207.3655).
2. The reference is a spike, not a stationary point
   (`dev/gpstop-stseg.R base ...`, log `stseg-base.txt`): on the segment
   from the default's optimum to the reference's, 41 points, the st80
   objective reads 68.550 at the default, 67.856 to 70.838 over t = 0.45
   to 0.98, and at t = 1 75.108 by continuation from the default's end
   but 79.545 from the reference's own state; st160 reads 201.294,
   201.061 at t = 0.98, and 204.682 or 207.366 at t = 1. That last pair
   is the path dependence lane arstep reported (204.77 against 207.37):
   it is at the spike only. The outer gradient at the reference is 9.1e8
   to 1.3e12 (st80) and 1.6e9 to 2.0e11 (st160); at the default 2.6e-6 to
   2.5e-5 and 4.9e-7 to 8.2e-6.
3. At the reference, one eigenvalue of the inner Hessian is near zero
   (`dev/gpstop-stis2.R base ... 40000`, log `stis2-base.txt`): the field's
   coefficients have a unit prior, so the inner Hessian holds eigenvalues
   of 1 or more where the likelihood is log-concave; the student density is
   not, past `sqrt(nu) sigma`, and there one direction is nearly flat.
   `-log(ev_min) / 2` is the part of the Laplace value that this one
   eigenvalue supplies:

        STIS2 student step 80 1 default laplace 68.5503 is 71.3256 (halves 71.3904 71.2563) ess 33 of 40000 | ev_min 0.999, -log(ev_min)/2 0.000
        STIS2 student step 80 1 ref     laplace 79.5450 is 73.9404 (halves 73.9513 73.9293) ess 74 of 40000 | ev_min 4.17e-08, -log(ev_min)/2 8.496
        STIS2 student step 160 1 default laplace 201.2938 is 203.6227 (halves 203.7120 203.5246) ess 424 of 40000 | ev_min 1, -log(ev_min)/2 0.000
        STIS2 student step 160 1 ref     laplace 207.3655 is 201.9325 (halves 202.0009 201.8590) ess 67 of 40000 | ev_min 1.82e-05, -log(ev_min)/2 5.457

   The near-zero eigenvalue supplies 8.50 of the 10.99 nats on st80 and
   5.46 of the 6.07 on st160. The marginal by importance sampling (t
   proposal, df 4, at the inner mode with the Hessian's eigenvalues floored
   at the prior's 1, 40000 draws, `set.seed(9001)`) reads the reference 2.6
   above the default on st80 and 1.7 below it on st160; with effective
   sample sizes of 33 to 424 these are rough, and each is biased low.

So the 11 and 6.07 nats are measured against non-stationary points of the
Laplace objective where its approximation fails, and every fit that
reaches them warns (code 1, a gradient warning). The default stops at a
regular maximum with a healthy inner Hessian. There is no warning to add
at the default's stop: nothing there differs from a correct fit. No code
change. (Fix round, item 2: on st80 the model's marginal is still about
2.5 above the default, so 5033 stays open as a limit of the Laplace
approximation; "no silent short stop" holds for the optimizer only.)

A first prototype that reweighted the working marginal at the default's
stop with the t density's EM weights (`dev/gpstop-stproto.R`, log
`stproto-base.txt`) read the same theta as the default on 6 student rows
(the descent from it gained -3.42 to 0.00), which fits this reading.

NEWS: a Limits bullet (rewritten in the fix round).

## The brief's other tables

### arstep's 150 rows, rellib-r20 against the lane

`dev/gpstop-ng.R` (lane arstep's `arstep-ng.R` with a second reference,
`bestc`: the best start that ended at code 0 with no warning), seed 1, the
10 cases: `lane 1 <cases> 1` in two processes (logs `ng-lane-a.txt`,
`ng-lane-b.txt`, with the 15 starts) and `base 1 <cases> 0` (log
`ng-base.txt`). Summary `dev/gpstop-ngsum.R`, output `ngsum.txt`:

    NGS all        ref all   rows 150 | short base  5 (silent  4) worst  -10.995 | lane  5 (silent  4) worst  -10.995
    NGS all        ref clean rows 150 | short base  2 (silent  2) worst   -1.591 | lane  2 (silent  2) worst   -1.591
    NGS student    ref all   rows  15 | short base  3 (silent  2) worst  -10.995 | lane  3 (silent  2) worst  -10.995
    NGS student    ref clean rows  15 | short base  0 (silent  0) worst    0.000 | lane  0 (silent  0) worst    0.000
    NGS moves: lane below base 0, above 0; equal at 1e-6 150 of 150
    NGSS pois_mid two 80 1 lane -202.6025 conv 0 w 0 | best -201.0112 (-1.591) | clean -201.0112 (-1.591)
    NGSS student step 40 1 lane 47.5939 conv 1 w 2 | best 47.7966 (-0.203) | clean 47.5939 (0.000)
    NGSS student step 80 1 lane 68.5503 conv 0 w 0 | best 79.5450 (-10.995) | clean 68.5503 (0.000)
    NGSS student step 160 1 lane 201.2938 conv 0 w 0 | best 207.3655 (-6.072) | clean 201.2938 (0.000)
    NGSS gsigma step 40 1 lane 34.4124 conv 0 w 0 | best 35.8612 (-1.449) | clean 35.8612 (-1.449)

No row moves (by construction: none of these fits holds a nonlinear body,
so `gp_body_relin_applies()` is FALSE and the route is 0.81.0's). Against
the starts that end cleanly, no student row is short. The two short rows
left, `pois_mid two 80` (1.59) and `gsigma step 40` (1.45), are the two
0.80.0's Limits names beside the student rows; neither is in this
lane's items, and each reads as on rellib-r20 (filed, below).

### arstep's cost designs

`dev/gpstop-cost-ar.R` (lane arstep's review cost script, arms lane and
`norelin` interleaved in one process), logs `acost.txt`, `acost-a500.txt`.
Every log-likelihood equals the other arm's and the value lane arstep
recorded (p600 32.589404, p1000 173.765935, a300 227.607924, pois240
-459.715793, bern240 -134.893263); walls per design, lane against norelin,
2 reps: p600 2.50 (the process's first fit), 0.77 against 0.92, 0.96;
p1000 1.06, 0.75 against 0.78, 0.86; a300 30.00, 24.06 against 24.44,
25.89; pois120 1.71, 1.58 against 1.81, 1.42; pois240 17.77, 16.98
against 17.89, 17.83; bern240 18.78, 17.20 against 17.68, 17.23; a500
(1 rep, ll 269.968071 both) 133.22 against 117.81. The routes are the
same code but for one test that returns FALSE, so these differences are
load.

## Suites

Every test file of core and the eight extensions, one file per process,
3 at a time, `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`,
`FRMTMB_RELLIB=wt-gpstop-lib` ahead of rellib-r20
(`bash dev/release/run-par.sh gpstop-gated dev/gpstop-all.jobs 3`; log
`dev/release/gpstop-gated.log`, "GPSTOP-GATED ran 510 of 510"; tally from
each file's RESULT and `lib:` lines into `gpstop-out/suite-sum.txt`):

    frmtmb files 317 pass 26542 fail 0 err 0 skip 15 warn 0 badlib 0
    frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0 badlib 0
    frmtmb.eam files 33 pass 1854 fail 0 err 0 skip 3 warn 0 badlib 0
    frmtmb.latent files 12 pass 687 fail 0 err 0 skip 2 warn 0 badlib 0
    frmtmb.learn files 16 pass 545 fail 0 err 0 skip 2 warn 0 badlib 0
    frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0 badlib 0
    frmtmb.sample files 87 pass 4414 fail 0 err 0 skip 1 warn 0 badlib 0
    frmtmb.spatial files 8 pass 561 fail 0 err 0 skip 0 warn 0 badlib 0
    frmtmb.survival files 15 pass 556 fail 0 err 0 skip 0 warn 0 badlib 0

`badlib` counts files whose `lib:` line does not name wt-gpstop-lib for
frmtmb: 0. The skips are other gates: core's `test-drmtmb-agreement.R`
13, `test-curve-scale.R` 1, `test-fuzz.R` 1, and each extension's
`test-scale.R`. No extension was reinstalled: none changed, and each
loads this lane's core from the library path.

Sessions, three ways, for core and frmtmb.sample (`bash
dev/release/run-sessions.sh gpstop-sess 3 frmtmb frmtmb.sample`, same
environment; log `dev/release/gpstop-sess.log`):

    SESSION frmtmb forward files=317/317 pass=26542 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb reverse files=317/317 pass=26707 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb preload files=317/317 pass=26810 fail=0 err=0 skip=17 warn=0
    SESSION frmtmb.sample forward files=87/87 pass=4414 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample reverse files=87/87 pass=4414 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample preload files=87/87 pass=4414 fail=0 err=0 skip=1 warn=0
    GPSTOP-SESS ran 6 of 6

LEAK lines 173, untagged 0; every `lib:` line names wt-gpstop-lib for
frmtmb. Only core and frmtmb.sample ran as sessions: by construction
nothing else reaches the change. `gp_body_relin_applies()` is TRUE only
for an exact gp() in a nonlinear body's parameter, `gp_start_scan()`
returns at its first line on a frame with no exact gp(), and the
extensions' tests fit exact gp() terms only in frmtmb.sample (lane
arstep's grep, section 13); every extension's files ran one per process
above.

## R CMD check --as-cran (core)

`bash dev/gpstop-check.sh .` (built and checked in `dev/gpstop-check/`,
`R_LIBS` = wt-gpstop-lib, rellib-r20, user library; pandoc and TinyTeX
on PATH; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; `NOT_CRAN=true`). Kept as
`gpstop-out/check-frmtmb.txt`; the check tree is deleted.

    Running 'testthat.R' [48m]
    [ FAIL 0 | WARN 0 | SKIP 272 | PASS 23236 ]
    * checking HTML version of manual ... [32s] NOTE
    Skipping checking math rendering: package 'V8' unavailable
    Status: 1 NOTE

The NOTE is the expected V8 one. No extension changed, so no extension
check ran.

## OpenBLAS 0.3.26 with LAPACK

`bash dev/ciharden-openblas.sh 0.3.26 lapack` (log
`openblas-build-0326.txt`: Rblas 74 symbols and Rlapack 616 to OpenBLAS;
dgemm 2000 in 0.12 s), then `bash dev/gpstop-openblas-run.sh 0.3.26`
(4 threads; logs `gpstop-out/ob0.3.26/`), every `blas:` line naming the
copy and every `lib:` line wt-gpstop-lib:

    RESULT test-gp-ar-start.R pass=49 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-noncentered.R pass=76 fail=0 err=0 skip=0 warn=0
    RESULT test-gp-start.R pass=243 fail=0 err=0 skip=0 warn=0
    RESULT test-nl-family.R pass=333 fail=0 err=0 skip=0 warn=0

0.3.32 was not run: the new assertions sit on no rounding-decided branch
(each compares the default fit with a fit from a start and accepts any
point at or above it; on 20 roundings the two read within 1.5e-11 of
each other). `dev/ciharden-out/` is deleted.

## Not done

- 4485 (`gp(x, k = )` too small): a check of the fitted length scale
  against the basis's top frequency is a new diagnostic for the
  approximate gp(), not part of the exact field's start. Not this change.
- 4491 (position names): touches the block's `cnms`, which every label
  reader shares. Not this change.
- 4388 (the SE check's Hessian from the marginal's tape): the SE check's
  tiers, not the start. Not this change.
- The two short rows of the 150 set (`pois_mid two 80`, `gsigma step 40`)
  and the 5 low-count Poisson rows of the 192 set: not in this lane's
  items, unchanged, filed with a candidate fix (a relinearization of the
  mean's working marginal at the fit's optimum).
- No warning for the student rows: at the default's stop nothing differs
  from a correct fit (healthy inner Hessian, gradient 1e-5), so any rule
  there would fire on correct fits. The spikes themselves are warned.
- Mutants: two in the fix round (`regular`, `template`, fix round item 7).

## Files touched (final, after the fix round)

- `R/gp-start.R`: `gp_start_working()` (argument `mu`), `gp_start_scan()`
  (argument `relin`; sets `.frmtmb_gp_relin`, a list of `start` and
  `regular`), `gp_second_descent()` (runs it, muffles its "NA/NaN"
  warnings, keeps only a regular point); new `gp_body_relin_applies()`,
  `gp_body_relin()`, `gp_relin_regular()`, `gp_inner_pd()` (after the
  short check; it replaced `gp_inner_ev_min()`), `gp_relin_ls_min`,
  `gp_relin_inner_ev_min`, `gp_body_relin_rounds`.
- `R/fit.R` (shared): `gp_cold <- is.null(template)` after
  `prefit_seeded <- FALSE`, and the `relin =` argument at the
  `gp_start_scan()` call in `fit_assembled()`. Two hunks, apart from
  nlstop's.
- `tests/testthat/test-gp-start.R`: design L in `gs_rv_row()`,
  `gs_logis_row()`, `gs_inner_ev()`, four tests.
- `NEWS.md`: a development section, one Fitting bullet (BREAKING), one
  Limits bullet.
- `dev/test-backlog.md`: a CLOSED line under 4822; section "Filed in the
  round of 2026-11-11", lane gpstop, five items (0.81.0's white-noise
  stops, 5033 re-described, approximate gp() in a body, the two short
  rows, the student spikes).
- `dev/gpstop-*` (scripts, this file).

No man page changes (`@noRd` only; roxygenise moved nothing).

## Bump and floor

The bump is left to the consolidation (the Fitting bullet is BREAKING).
No extension calls a changed function, so no floor moves.
