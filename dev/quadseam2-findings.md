# Lane quadseam2: one-block quadrature seams

Round opened 2026-11-13, base ec8178d3 (frmtmb 0.82.0). Worktree
`frmtmb-wt-quadseam2`, branch `wt-r83-quadseam2`, library
`wt-quadseam2-lib`, base library `rellib-r21`. Brief:
`round-briefs/r083-quadseam2.md`.

Status: fix round done (review dev/reviews/2026-11-13-quadseam2.md, FIX).
The section "Fix round" at the end wins over every section before it.
Items 4836, 4861, 4867, 5222, 5620 fixed; 5408 closed as not a defect
(its change reverted); 561 closed as not reproduced; 4846 measured, route (b) closed,
route (a) filed; 4856, 5615 (two pairs) and the rest of cluster 5 left
open. Logs are in `dev/quadseam2-log/` (local,
gitignored); every number below is pasted from them, with the script
that made it. "base" is rellib-r21, "lane" is wt-quadseam2-lib.

## Item 5222: the stationary escape made one pass

Reproducer `dev/quadseam2-hmm.R` (a port of `dev/silent-hmm-probe.R`,
plus `mixture(skew_normal(), skew_normal())` of 400 rows, seeds 101 to
103). Every alpha started at 0. Logs `hmm-base.txt`, `hmm-lib.txt`:

    base HMM seed 1 free -998.2551 alpha 4.18 -5.451 | stall -1019.4522 alpha -8.365e-14 -5.607 escape 4.0 37.5
    base HMM seed 2 free -960.6004 | stall -985.6858 alpha 4.468 -1.449e-14
    base HMM seed 3 free -953.5887 | stall -976.1491 alpha 5.126 -5.035e-14
    base MIX seed 1 free -635.4970 | stall -644.7523 alpha 12.35 3.451e-14
    base MIX seed 2 free -579.4192 | stall -588.1827 alpha 7.194 -1.917e-13
    base MIX seed 3 free -614.0458 | stall -621.4442 alpha -1.255e-14 -6.635
    lane HMM seeds 1 to 3, MIX seeds 1 to 3: stall = free to the 4
         printed decimals, both alphas off 0, escape starts 6

Change (`R/fit.R`): `stationary_escapes()` records how many predictors
sit on their point (`n_on`), and `escape_stationary()` asks again from
the kept optimum while a pass gains, at most `n_on` passes. The record's
`starts` is the total over the passes.

## Items 4836 and 4861: an sd stuck near zero under quadrature

4836 was filed as "no Newton polish on the one-block route". It is not:
the stall fit's group sd is at zero. `dev/quadseam2-snref.R` (seed 6 of
the item, the exact marginal log-likelihood by `stats::integrate()` per
group, rel.tol 1e-12; log `snref-base.txt`):

    FIT default  alpha 5.26735 ll -244.100465993 exact -244.100465992 err -1.360e-09
    FIT stall    alpha 5.00776 ll -244.155229735 exact -244.155229731 err -4.116e-09
    EST default theta -2.54807   EST stall theta -20.389

The rule is right at both points; the stall stop is the integral's own
value with the sd at exp(-20.4), a point the Laplace fit leaves through
`sd_up_escape()` and the quadrature fit could not (item 4861). So both
items are fixed by one change.

Change (`R/fit.R`): `sd_up_start()` (the probe half of
`sd_up_escape()`) also runs on the one-block route, where
`se_sd_gain_up()` already reads rules taped at each larger sd
(`se_probe_fn()`). `quad_fit()` returns `restart(par)`: a rule taped at
`par`, optimized from there and settled by `quad_settle()`, whose
unsettled warning is left to the caller (`warn = FALSE`, the message in
`info$msg`). `quad_sd_up_escape()` keeps the restart when it gains, and
`fit_assembled()` then replaces the fit's tape, optimum, estimates
(conditional modes from the Laplace objective) and cache.

`dev/quadseam2-seed6.R`, logs `seed6-base.txt`, `seed6-lib.txt`:

    base quad TRUE start stall   ll -244.155229735 alpha 5.0078
    lane quad TRUE start stall   ll -244.100466369 alpha 5.2681
    both quad TRUE start default ll -244.100465993 alpha 5.2673

## Item 4867: a coefficient run off under one-block quadrature

`dev/quadseam2-limit.R` (the data of `dev/reviews/quadseam-review-limit.R`
plus a bernoulli separation and the two absent cases), logs
`limit-base.txt`, `limit-lib.txt`:

    base runoff QUAD TRUE  se fc 7.471e+03, no warning
    lane runoff QUAD TRUE  se fc NaN, "fc (at -22.53, toward -Inf)", as Laplace
    both absent QUAD TRUE  every SE finite, no warning, bitwise the same SEs
    base bernsep QUAD TRUE frmtmb_separation + "Some standard errors are not finite"
    lane bernsep QUAD TRUE frmtmb_separation only, as Laplace

Change (`R/se-check.R`, three hunks): the limit rule's gate in
`cov_from_hessian()` is `se_edge_probed()` (off only on the nested
route) in place of `!isTRUE(fit$quadrature)`; `se_edge_limit()` reads
`se_probe_fn()`; `se_limit_curv()` takes its central differences on a
rule taped at the point it reads (`attr(quad_reanchored_fn(fit),
"obj_at")`, new in `R/fit.R`), not on the AD Hessian, which is NaN on
Gauss-Kronrod tapes of poisson fits (dev/quadseam-findings.md).

## Item 5408: the escape ran twice under quadrature

Change (`R/fit.R`): the rule's escape in `fit_assembled()` runs only when
`quad_fit()`'s Laplace-side escape did not. `stationary_escape_keep()`
is removed, since the second record it merged no longer exists.

`dev/quadseam2-escdup.R` wraps `escape_stationary()` and logs each call
that has starts (logs `escdup-base.txt`, `escdup-lib.txt`). FORCED is a
gaussian custom family whose declared window on sigma holds every fit
(tol 1e6), with `(1 | g)`, 120 rows, seeds 1 to 20; SYM and STALL are
`dev/rel082-quad-escape.R`'s two skew-normal sets:

    base SUM FORCED fits 20 laplace 20 rule 20 rule secs 0.79 of 5.33 max rule gain 1.02e-12
    base SUM SYM    fits 40 laplace 4  rule 0
    base SUM STALL  fits 10 laplace 2  rule 0
    lane SUM FORCED fits 20 laplace 20 rule 0
    lane SUM SYM    fits 40 laplace 4  rule 0
    lane SUM STALL  fits 10 laplace 2  rule 0

The second escape gained at most 1.02e-12 (rounding) where it ran.
logLik base against lane: FORCED and SYM bitwise equal on all 60;
STALL equal on 9 of 10, seed 6 moving from -244.155229735 to
-244.100466369 (the sd escape of 4861).

## Tests, seen to fail on base

On rellib-r21 (logs `dev/quadseam2-log/tests/*-base.txt`):

    RESULT test-quad-escape-record.R pass=4 fail=1
      FAIL `ran` actual 2, expected 1
    RESULT test-stationary-passes.R pass=1 fail=3 warn=1
      FAIL all(abs(a) > 1); logLik gap 9.255287 >= 6e-6; starts 4, not 6
    RESULT test-quadrature-defects.R pass=80 fail=6 err=1
      FAIL sd escape gain; limit warnings length 0; se_lost NULL; SE not NaN

On the lane: test-quad-escape-record.R pass=5, test-stationary-passes.R
pass=4, test-quadrature-defects.R pass=89, each fail 0 err 0 warn 0.

## Items 4846 and 561: the rule's error at its anchor

Route (b), a smaller `dx` through RTMB's GK spec, is closed without an
upstream change. TMB's `MakeADFun()` passes the spec to
`TransformADFunObject()` as `config`, and `tmb_core.hpp` (TMB 1.9.25,
lines 1588 to 1596) reads only `adaptive` and `debug` from it into
`TMBad::gk_config`; `dx` (1), `ytol` (1e-2) and `nan2zero` keep their
defaults (`TMBad.cpp:3863`). RTMB 2.0's DLL holds no `ytol` string
(`dev/quadseam2-rtmbgk.R`, `dev/quadseam2-tmbgk.R`). Filed upstream as
TMB-7 in `dev/upstream-bugs.md` (not posted).

What the calibration costs, measured instead of route (a). The tape
can be given an integrand whose mode and scale are exact: integrate u
with b = bhat + s u, bhat and s the Laplace conditional mode and
1/sqrt(inner Hessian diagonal) at the anchor, constants on the tape, and
subtract sum(log s). TMB then calibrates a unit-width integrand at mode
0, where its differences of step 1 and its 1e-2 mode search are
accurate. `dev/quadseam2-white.R` tapes both rules at each fit's own
estimate on the 120 fits of the 0.82.0 table (bern500, pois200,
bern2sd3, poissd2, seeds 6 to 35, the base build's fits) and reads them
against `ref_ll()` of `dev/reviews/quadseam-review-acc.R`. Logs
`white-base-a.txt`, `white-base-b.txt`, summarized by an awk over the
ROW lines:

    DES bern2sd3 n 30 max|raw| 1.064e-03 max|white| 1.188e-03 better 4 worse 5
    DES bern500  n 30 max|raw| 9.408e-06 max|white| 9.423e-06 better 0 worse 0
    DES pois200  n 30 max|raw| 1.173e-05 max|white| 8.803e-06 better 17 worse 9
    DES poissd2  n 30 max|raw| 6.837e-05 max|white| 6.275e-05 better 8 worse 7
    ALL n 120 max|raw| 1.064e-03 max|white| 1.188e-03 better 29 worse 21
        raw > 1e-5 on 6, white > 1e-5 on 7   ("better/worse" by > 1e-6)

The exact mode and scale do not remove the error, so `dx` and `ytol`
are not what sets it. TMBad's `integrate()` runs with its default
control, `reltol = abstol = 1e-4` per integral (`integrate.hpp:993`,
`TMBad.cpp:4680`; `logIntegrate_t` calls `integrate(*this)`), and the
subdivision it settles on depends on the integrand's shape at the
anchor. That is the inferred source; it was not isolated by a run with
other tolerances, which RTMB cannot set. Whitening was therefore not
built.

What remains for 4846 is route (a), a Gauss-Hermite rule written on the
tape with enough nodes, as the nested route writes its rule (accuracy
6.5e-10 at k = 31 on the nested Bernoulli reference of 0.80.0). It
changes every one-block value and needs the 0.80.0 accuracy tables for
one block; size L, not started. Filed under "From lane quadseam2".

561's breakdown does not reproduce. `dev/quadseam2-break.R`: 15 groups
of 4, true group sd 0 and 0.05, seeds 1 to 20, poisson, Gamma(log),
Beta, negbinomial, bernoulli, lognormal, `quadrature = TRUE`, on base
(log `break-base.txt`): 240 fits, 0 refusals, 0 errors; 145 of them
end with the sd below 1e-3. CLOSED as not reproduced.

## Items 5620 and 5615: frmtmb_latres_laplace on a lognormal

`dev/quadseam2-latres.R` (the lognormal data of
`dev/fitedge-mi-latres.R`, seeds 40401 to 40405, each as `y | mi()` and
as its weight-0 spelling; plus poisson and bernoulli latents, seeds 501
to 503, where the warning is wanted). Logs `latres-base.txt`,
`latres-lib.txt`:

    base LN 40401 mi: latres 0 | w0: latres 1   (other 9 lognormal fits 0)
    lane LN 40401..40405, both spellings: latres 0 on all 10
    both BERN 501..503: latres 1 on each; POIS 501..503: 0 on each

On 0.82.0 the item's other two seeds (40403, 40405) no longer fire;
40401's weight-0 spelling does, the rounding-decided case of 5615.

Change (`R/autocor-latent.R`, arspat2's file; two hunks):
`latres_laplace_reason()` returns NULL when the new `latres_quadratic()`
finds every row's log density quadratic in eta: the third and fourth
differences at a step of the latent's marginal sd vanish to 1e-8 of the
row's log density. A first version reused the share path's cubic test
(`latres_cubic_floor`), and it dropped the warning on bernoulli seed 502,
whose rows near p = 0.5 have a third derivative near 0 and a fourth that
is not; the fourth difference is why the gate reads both.

5615's other two pairs (Gamma `ma(cov = TRUE)` seed 3, weibull
`ar(cov = TRUE)` seed 3) are not quadratic and are untouched; they were
not rerun (`dev/review-fitedge-mi.R` sources a helper in another lane's
tree). 5615 stays open.

## Cluster 5 of part B, read

- 5620: fixed (above). 5615: the lognormal pair fixed, the Gamma and
  weibull pairs left open.
- 5729 (nested quadrature, an unsettled fit carries an error of a few
  1e-2): not in reach. It is the nested rule's own error on hard
  designs (more nodes do not always help, `?frm`); a fix is a rule
  that converges there, not a seam. Left open.
- 5832 (fit the first-order Laplace correction's shift instead of
  warning): a feature, size M to L; not started.
- 5873 (importance-sampling bundle): four separate items (crossed
  factors, a heavy-tail guard, a poisson false convergence, MCSE by
  pairs), none a quadrature seam; not started.
- 5738 (`fitted()` under a draw on a periodic body; the lognormal
  mixture question Q1 for the user): a mean-integration rule, not
  quadrature of the likelihood; not started. Q1 still needs the user.
- 5593 (`fitted()` cost with a latent residual beside many levels):
  a sparse-pattern change in `Q_rr`, cost only; not started.
- 4513 (the SE check's outer Hessian on an exact gp field from the
  marginal's tape): touches the SE tiers (selimit2's file); not started.

## Item 4856: the nested route's boundary verdict

Not done. The probes need a nested rule taped at each probed sd, which
needs the Laplace objective for `quad_nest_anchor()`; the fit does not
keep it. `quad_reanchored_fn()` rebuilds a Laplace-free GK tape for the
one-block route; the nested route would rebuild the Laplace objective
from `build_objective()` in the same way, then anchor and tape the
nested rule at each point. Left open as filed.

## Acceptance: the 0.82.0 accuracy table and tight sweeps

The lane build against rellib-r21's own runs of the release record
(`frmtmb-wt-rel082/dev/rel082-log/kt/quadseam/acc3-lane.txt` and
`fr-sweep-lane.txt`, both on rellib-r21). Lane runs: `dev/quadseam2-acc3.R`
(a copy of `dev/reviews/quadseam-review-acc.R` pointed at the round's
libraries) and `dev/quadseam2-sweep.R` (a copy of
`dev/quadseam-fr-sweep.R`), driver `dev/quadseam2-acc-run.sh`, logs
`acc3-lane.txt`, `sweep-lib.txt`:

    ACC fits 120 logLik identical to base 120 lane worse >1e-6 0 better 0
        max|err| lane 1.064e-03 base 1.064e-03 warned 0
    grad_tol 1e-05  fits 50 warned 0 rounds>0 29 max|err| 5.331e-06 (base the same)
    grad_tol 0.0001 fits 50 warned 0 rounds>0 17 max|err| 5.354e-06 (base the same)
    SWEEP logLik identical 100 of 100

So the glmer(nAGQ = 25) comparison of the release record holds
unchanged. None of these 220 fits stops with an sd near zero short of a
maximum or with a coefficient run off, so neither new path moves them.

## ulp draws

`dev/quadseam2-ulp.R <lib> <case>`, k = 1..20, logs
`ulp-lib-<case>.txt`:

    ULP limit pass 20 of 20     (4867, the covariate moved)
    ULP sdesc pass 20 of 20     (4836/4861, the response moved)
    ULP mix pass 20 of 20       (5222, the response moved)
    ULP latln pass 20 of 20     (5620, both spellings, the response moved)
    ULP latbern pass 20 of 20   (5620's guard-present bernoulli)

## Runs

All on the lane build; every per-file log's `lib:` line names
`wt-quadseam2-lib` for frmtmb, and rellib-r21 for each extension (no
extension's R code changed; frmtmb.latent's test file did and ran from
this tree).

- Every test file of core and all eight extensions, one per process, 3
  at a time, `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`
  (`dev/release/run-par.sh quadseam2-all`, job list
  `dev/quadseam2-log/all.jobs`, log `dev/release/quadseam2-all.log`):
  `QUADSEAM2-ALL ran 531 of 531`, 0 without a RESULT line.

      frmtmb           files 332 pass 27471 fail 0 err 0 skip 15 warn 0
      frmtmb.coupling  files  11 pass   551 fail 0 err 0 skip 5 warn 0
      frmtmb.eam       files  33 pass  1857 fail 0 err 0 skip 3 warn 0
      frmtmb.latent    files  13 pass   699 fail 0 err 0 skip 2 warn 0
      frmtmb.learn     files  16 pass   545 fail 0 err 0 skip 2 warn 0
      frmtmb.ode       files  11 pass   549 fail 0 err 0 skip 1 warn 0
      frmtmb.sample    files  91 pass  4603 fail 0 err 0 skip 1 warn 0
      frmtmb.spatial   files   8 pass   612 fail 0 err 0 skip 0 warn 0
      frmtmb.survival  files  16 pass   575 fail 0 err 0 skip 0 warn 0

  Every skip is an opt-in tier: test-scale.R in six packages,
  test-curve-scale.R, test-fuzz.R, test-drmtmb-agreement.R. test-fuzz.R
  with `FRMTMB_FUZZ=true` (its specs include a quadrature mode): pass=2
  fail=0 err=0 skip=0 warn=0 (`dev/quadseam2-log/tests/fuzz-lib.txt`).
- test-se-tier.R gained two assertions in its quadrature test (the sd
  started at exp(-20) now reaches the Laplace fit's log-likelihood, and
  `sd_escape` records the gain). On rellib-r21 (`tests/setier-base.txt`):
  pass=112 fail=2, "0.00606179 >= 1e-08" and `sd_escape` gain FALSE; on
  the lane pass=114 fail=0.
- One session per suite, three orders, every package
  (`NOT_CRAN=true bash dev/release/run-sessions.sh quadseam2-sess 3`, log
  `dev/release/quadseam2-sess.log`): `QUADSEAM2-SESS ran 27 of 27`, fail
  0 err 0 warn 0 in all 27:

      frmtmb           forward/reverse/preload files=319/332 pass 24417/24582/24685 skip 240/240/242
      frmtmb.coupling  11/11 pass 551 skip 5 (each order)
      frmtmb.eam       33/33 pass 1851 skip 5
      frmtmb.latent    13/13 pass 699 skip 2
      frmtmb.learn     16/16 pass 474 skip 13
      frmtmb.ode       11/11 pass 549 skip 1
      frmtmb.sample    86/91 pass 4238 skip 18
      frmtmb.survival  16/16 pass 575 skip 0
      frmtmb.spatial    8/8  pass 529 skip 6

  46 LEAK lines, every one tagged `[with load of ...]`; 0 untagged. The
  ungated file counts (319/332, 86/91) are file-level gated skips, as in
  earlier records.
- frmtmb.latent's new test on rellib-r21 (`tests/hmm-base.txt`):
  test-hmm-stationary.R pass=10 fail=2 warn=1, "all(abs(alpha)) > 1"
  FALSE and a logLik gap of 21.19709; on the lane pass=12 fail=0.
- test-latres-quadratic.R on rellib-r21 (`tests/lq-base.txt`): pass=1
  fail=1 err=2, `c(r1$n, r2$n)` was `0 1` (the other two are the missing
  `latres_quadratic()`); on the lane pass=4 fail=0.
- OpenBLAS with LAPACK (`dev/quadseam2-openblas.sh`, which runs
  `dev/ciharden-openblas.sh <ver> lapack`; logs
  `dev/quadseam2-log/ob<ver>-*.txt`), the files this lane added or
  changed plus the neighbors that test the changed paths; each the same
  at 0.3.26 and 0.3.32, fail 0 err 0 skip 0 warn 0:

      test-quadrature-defects.R pass=89    test-stationary-passes.R pass=4
      test-latres-quadratic.R pass=4       test-quad-escape-record.R pass=5
      test-se-tier.R pass=114              test-skew-normal-start.R pass=70
      test-quadrature-nested.R pass=69     test-se-limit.R pass=240
      frmtmb.latent test-hmm-stationary.R pass=12

  `dev/ciharden-out/` removed by the script.

## Files touched

- `R/fit.R`: `escape_stationary()`, `stationary_escapes()` (5222);
  `fit_assembled()` (the rule's escape gate, 5408; the one-block sd
  escape, 4861); `quad_fit()` (`restart`), `quad_settle()` (`warn`,
  `info$msg`), new `quad_unsettled_warn()`; `quad_reanchored_fn()`
  (`obj_at` attribute); `sd_up_escape()` split into `sd_up_start()`,
  new `quad_sd_up_escape()`; `stationary_escape_keep()` removed; the
  `frm()` `@param quadrature` text.
- `R/se-check.R` (selimit2's file; three hunks, read against its tree,
  which had no change when this lane started): the limit-rule gate in
  `cov_from_hessian()`, `fn` in `se_edge_limit()`, `se_limit_curv()`.
- `R/autocor-latent.R` (arspat2's file; two hunks):
  `latres_laplace_reason()` gate, new `latres_quadratic()`.
- `man/frm.Rd` regenerated.
- Tests: `test-quadrature-defects.R` (two tests), `test-quad-escape-record.R`
  (rewritten: the function it tested is gone), new
  `test-stationary-passes.R`, new `test-latres-quadratic.R`,
  `test-se-tier.R` (two assertions), and
  `extensions/frmtmb.latent/tests/testthat/test-hmm-stationary.R` (one
  test).
- `NEWS.md` (a `# frmtmb (development version)` heading with three
  sections), `dev/test-backlog.md` (7 CLOSED lines, one filing, by
  `dev/quadseam2-backlog-edit.R` and `-wrap.R`), `dev/upstream-bugs.md`
  (TMB-7).
- New dev scripts `dev/quadseam2-*`.

## Bump and floor

Behavior changes marked BREAKING (new restarts move estimates; a new
NaN standard error and warning under quadrature; a second escape pass):
a minor bump of core. No new dependency. frmtmb.latent's R code is
unchanged, but its new test needs core's two-pass escape, so its floor
on frmtmb must rise to the core release that carries this (0.83.0) for
that test to pass against an installed core.
- `R CMD check --as-cran`, built and checked in `dev/quadseam2-check/`
  (`dev/quadseam2-check.sh`; R_LIBS lane, rellib-r21, user;
  `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`): core `Status: 1 NOTE`, and
  frmtmb.latent `Status: 1 NOTE`, each the V8 math-rendering NOTE on the
  HTML manual (logs `check-frmtmb.txt`, `check-frmtmb.latent.txt`; the
  check trees were deleted after reading).

## Merge notes for the consolidation

Read at the end of this lane (selimit2's and arspat2's working trees,
read-only):

- `R/se-check.R`: selimit2 has restructured the limit-rule block of
  `cov_from_hessian()` (its hunk at old line 388 keeps the old gate line
  as context and changes the two lines after it) and the head of
  `se_edge_limit()` (its hunk at old line 791 ends one line above this
  lane's `fn` line). Both of this lane's hunks there touch adjacent
  lines, so a textual conflict is likely. The intent to keep: the gate is
  `se_edge_probed(fit) && is.null(fit$importance)` in place of
  `!isTRUE(fit$quadrature) && ...`, and `fn` reads `se_probe_fn(fit)`
  (`pf`) in place of `obj$fn`. selimit2's new joint directions reuse
  `fn` and `se_limit_curv()`, so with both kept they read rules taped at
  each point on the one-block route; that combination was not run here.
  This lane's `se_limit_curv()` hunk does not overlap selimit2's.
- `R/autocor-latent.R`: arspat2 changes old lines 213 and 673 only;
  no overlap with this lane's two hunks (old 1045 and 1054).
- `R/fit.R`: nlstop2's hunks (old 213, 1106) do not overlap this lane's.

## Fix round (review `dev/reviews/2026-11-13-quadseam2.md`, verdict FIX)

This section wins over every section above. Where they differ (item
5408, the latres gate's step, the NEWS bullets), the numbers here are
the ones measured on the fix-round build.

### F1 (MAJOR): item 5408 was a silent wrong answer; reverted

The rule's escape was skipped whenever the Laplace-side escape left a
record. The rule's own optimization can walk back onto alpha = 0 after
that escape, and then nothing moved it. The reviewer's
`dev/reviews/quadseam2-review-esc2.R` (45 skew normals with `(1 | g)`
under quadrature, every alpha started at 0) found 3 such fits on the
first-round build, 1.45 to 6.67 below 0.82.0. My `escdup` sets had never
put the rule's optimum on the point, so "gained nothing" rested on the
FORCED set alone.

Change (`R/fit.R`): 0.82.0's behavior restored. The rule's escape runs
whenever its optimum is on a declared point, `stationary_escape_keep()`
is back, and the quadrature sd restart keeps the fit's
`stationary_escape` and `stationary_escape_laplace` records.

    esc2 on the fix-round build (dev/quadseam2-log/fr-esc2-lib.txt) against
    the reviewer's esc2-base.txt: fits 45, lane below base by > 1e-8: 0,
    differ beyond 1e-10: 0
    SYMSMALL 1 rule 2.291 ll -84.9253951576   (first-round build -87.2167801711)
    WEAK 4     rule 1.453 ll -74.1377898690   (first-round build -75.5906679346)
    ALPHAX 11  rule 6.667 ll -84.5783086103   (first-round build -91.2452361629)

5222's gains are kept (`hmm-lib-fr.txt`: HMM seeds 1 to 3 and MIX
seeds 1 to 3, stall = free to 4 decimals, escape starts 6).

`test-quad-escape-record.R` is back to 0.82.0's record test plus a
regression test on SYMSMALL seed 1 (the rule's escape gains more than
1e-3 of |logLik|, alpha is off 0, a Laplace record is kept, no
warning). On the first-round build (`tests/fr-qer-pre.txt`): pass=1
fail=2 err=1, alpha 7e-15, `stationary_escape_laplace` NULL, and the fit
warned (`frmtmb_warning`); on rellib-r21 pass=7 fail=0; on the fix-round
build pass=8 fail=0.

NEWS: the "escaped once" bullet is removed. Backlog: 5408's CLOSED line
now says "not a defect", with the reviewer's three fits as the reason.
The FORCED cost (the second escape ran on 20 of 20 fits whose window
holds every fit, 0.79 of 5.33 s, `escdup-base.txt`) is the price of a
window that always fires; on a fit off its point the escape returns at
once.

### F2 (MEDIUM): `latres_quadratic()` step

The differences were taken at the latent's marginal sd, so at a small
enough sd any family read as quadratic: FEWBERN seed 703 (a bernoulli,
4 groups of 5, sderr at exp(-10.12), marginal sd 7.6e-5) lost the
warning 0.82.0 gave. Change (`R/autocor-latent.R`): the step is
`pmax(h, 1)` on the link scale; a row finite at the fit but not finite
at a moved point makes the answer FALSE, and a row not finite at the fit
(an unread `mi()` row) is left out.

    dev/reviews/quadseam2-review-fewbern.R lib (fr-fewbern-lib.txt):
      gate FALSE on all 4 seeds at every step (0.001, 0.01, 0.1, 0.3, 1)
    dev/reviews/quadseam2-review-latres.R lib (fr-latres-lib.txt) against
    the reviewer's latres-lib.txt, 33 rows: one row differs,
      FEWBERN 3: latres 0 quad TRUE  ->  latres 1 quad FALSE (as 0.82.0)
    lognormal at y x 1e-8 and 1e8, seeds 40401 to 40403, both spellings:
      gate TRUE, latres 0 on all 12
    ulp draws (fr-ulp-*.txt): ULP latbern pass 20 of 20, ULP latln pass 20 of 20

New test in `test-latres-quadratic.R` on FEWBERN 703: on the first-round
build (`tests/fr-lq-pre.txt`) pass=4 fail=2 (`r$n` 0, gate TRUE); on the
fix-round build pass=6 fail=0. The NEWS bullet is marked BREAKING.

### Minor

`test-stationary-passes.R` accepts `starts` of 4 or 6: two passes (four
starts, then two from the kept optimum), or one pass of four when a
start that moves one alpha carries the other off its point on the way.
Both end at the same fit, which the logLik assertion checks.

### Merge note for the consolidation

After priors2 merges: priors2 adds `check_map_unbounded(opt, bounds)`
before `fit` is built, and this lane's quadrature sd restart replaces
`opt` after that. Call `check_map_unbounded(esc$opt, bounds)` in the
restart branch of `fit_assembled()` before it is kept. frmtmb.latent's
floor on frmtmb rises at consolidation (its new test needs this core).
The reviewer's three-way merges of `R/se-check.R` with selimit2 were
clean, so the conflict predicted above does not occur.

### Fix-round runs

- 120-fit accuracy table on the fix-round build (`fr-acc3.txt`) against
  rellib-r21's record: `ACC fits 120 logLik identical to base 120 max|err|
  lane 1.064e-03 warned 0`.
- Changed and neighboring test files on the fix-round build, fail 0 err 0
  warn 0 each: test-quad-escape-record.R 8, test-stationary-passes.R 4,
  test-latres-quadratic.R 6, test-quadrature-defects.R 89,
  test-se-tier.R 114, frmtmb.latent test-hmm-stationary.R 12.

### The F1 regression test, and a defect 0.82.0 has under OpenBLAS

The first version of the new test used SYMSMALL seed 1. It failed under
OpenBLAS 0.3.26 with LAPACK, `OPENBLAS_CORETYPE=Haswell`
(`fr-ob0.3.26-test-quad-escape-record.txt`): alpha 8e-15 and a warning.
`dev/quadseam2-symsmall.R` shows why: under OpenBLAS the rule's escape
runs (2 starts) and both restarts come back to alpha = 0, on rellib-r21
and on the lane alike (`fr-symsmall-ob-base.txt`, `-ob-lib.txt`:
logLik -87.2167801711 on both); under the reference BLAS the same escape
gains 2.29 (`fr-symsmall-ref-lib.txt`). The reviewer's esc2 under
OpenBLAS (`fr-esc2-ob-base.txt`, `fr-esc2-ob-lib.txt`): lane = base on
all 45; SYMSMALL 1 and WEAK 4 stop at alpha 1e-15 on both, 2.29 and 1.45
below their reference-BLAS values; ALPHAX 11 escapes on both (rule gain
6.685, logLik -84.5602809547). So this is 0.82.0's own behavior,
BLAS-dependent, and not this lane's change; filed under "From lane
quadseam2" (MEDIUM).

The test now uses ALPHAX seed 11 (`alpha ~ x`, both coefficients started
at 0), which escapes under both BLASes. It accepts the warnings a run-off
alpha can carry (`frmtmb_conv_code`; and on OpenBLAS also the
standard-error `frmtmb_warning`) and no others. Seen failing: the
first-round gate was put back into `R/fit.R` from a copy, installed, the
test run, then the fix-round file restored and installed again
(`tests/fr2-qer-firstround.txt`: pass=7 fail=2, the Laplace record NULL
and `all(abs(alpha) > 1)` FALSE). On the fix-round build pass=9 fail=0
(`tests/fr2-qer-lib2.txt`); on rellib-r21, which escapes as the fix
round does, pass=8 fail=0 on the file before its warning assertion was
written (`tests/fr2-qer-base.txt`).

OpenBLAS with LAPACK, `OPENBLAS_CORETYPE=Haswell` (`dev/quadseam2-fr-openblas.sh`,
then the rewritten test alone, `fr2-ob<ver>-qer.txt`), each fail 0 err 0
skip 0 warn 0 at both 0.3.26 and 0.3.32: test-quad-escape-record.R 9,
test-stationary-passes.R 4, test-latres-quadratic.R 6,
test-quadrature-defects.R 89, test-se-tier.R 114,
test-skew-normal-start.R 70, test-quadrature-nested.R 69,
test-se-limit.R 240, frmtmb.latent test-hmm-stationary.R 12.

Disk: drive C: had 2.6 GB free when the second OpenBLAS build ran out
of space ("No space left on device"); the builds were then made and
deleted one at a time. `dev/ciharden-out/` is gone.

### The step of 1 and an escaped warning

The full core suite on the fix-round build (`dev/release/run-par.sh
quadseam2-fr-core`, `dev/quadseam2-log/core.jobs`, one file per process,
3 at a time, `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`; log
`dev/release/quadseam2-fr-core.log`): `QUADSEAM2-FR-CORE ran 332 of 332`,
pass 27477 fail 0 err 0 skip 15 warn 1. The warning: test-latres.R, "a
link that does not hold the mean is read by its margin", "NaNs
produced". A step of 1 on the link scale takes some rows outside the
support (a Gamma inverse link at eta = 2), and the density then warned
from inside `latres_quadratic()`. Such a row reads non-finite, which is
the gate's answer, so the evaluation is now quiet (`suppressWarnings()`
around the density only). After it, every test file that reaches a
latent residual (grep for `cov = TRUE` or `latres`: 35 core, 12
frmtmb.sample, 1 frmtmb.learn; `dev/release/run-par.sh quadseam2-fr3`):
`QUADSEAM2-FR3 ran 48 of 48`, pass 6019 fail 0 err 0 skip 0 warn 0; and
the reviewer's 33 latres rows are identical to the fix-round rows above
(`fr3-latres-lib.txt`). The other 284 core files ran on the build before
this change, which differs from it only in that warning; the 15 skips
are the opt-in tiers.

Disk: by the time of the core suite drive C: had 115 to 134 GB free
again.
- OpenBLAS with LAPACK, `OPENBLAS_CORETYPE=Haswell`, after that change
  (`fr3-ob<ver>-*.txt`), each at 0.3.26 and 0.3.32 fail 0 err 0 skip 0
  warn 0: test-latres-quadratic.R 6, test-latres.R 305,
  test-quad-escape-record.R 9. `dev/ciharden-out/` removed.
- `R CMD check --as-cran` of core on the final fix-round build
  (`dev/quadseam2-check.sh`, `dev/quadseam2-check/check-frmtmb.txt`):
  `Status: 1 NOTE`, the V8 math-rendering NOTE. The check tree was
  deleted after reading. frmtmb.latent was not checked again: its only
  change, the test file, is the same as in the first round's check.
- The coordinator's 2-process limit arrived during the last job; by then
  only the check was running.

### Short check (review `dev/reviews/2026-11-13-quadseam2-short.md`, SHIP)

The ALPHAX regression test now uses `allow_warnings()` with the two
texts that fit gives ("Optimizer did not report convergence: false
convergence", seen under the reference BLAS by
`dev/quadseam2-alphax-warn.R`; "Standard errors are not available for",
the class seen under OpenBLAS, whose text was not rerun there), so any
other warning reaches the test. Lane, reference BLAS: pass=8 fail=0
warn=0 (`tests/sc-qer-lib.txt`, and `sc-qer-lib2.txt` after the
reinstall). First-round gate put back from a copy and installed:
pass=6 fail=2, Laplace record NULL and alpha at 0
(`sc-qer-firstround.txt`); the final `R/fit.R` was then restored and
installed.

The OpenBLAS text, observed (coordinator's follow-up): under OpenBLAS
0.3.26 and 0.3.32 with LAPACK, `OPENBLAS_CORETYPE=Haswell`,
`dev/quadseam2-alphax-warn.R` gives, besides the false-convergence code,
"Standard errors are not available for 4 of 6 parameters. (Intercept),
x, alpha_(Intercept), alpha_x: the likelihood is flat along them ..."
(`sc2-ob<ver>-warn.txt`); the reference BLAS gives only the code. The
allowance now pins that text up to "the likelihood is flat". It is true
for alpha_(Intercept) and alpha_x, which run toward the half-normal
limit (about 1e6). For (Intercept) and x it is not established: the
reference BLAS keeps their standard errors on the same fit, so that part
of the verdict depends on the BLAS's rounding of a near-flat Hessian;
filed (LOW) as a remark under the OpenBLAS item in "From lane
quadseam2". The file: pass=8 fail=0 warn=0 under the reference BLAS
(`tests/sc3-qer-lib.txt`) and under OpenBLAS 0.3.26 and 0.3.32
(`sc3-ob<ver>-qer.txt`). The count is 8, not 9, because the class
check `expect_true(all(w %in% ...))` was an expectation and
`allow_warnings()` without `require` adds none.
