# Lane quadseam: quadrature seams

Round opened 2026-11-11, base e8abed0a (frmtmb 0.81.0). Worktree
`frmtmb-wt-quadseam`, branch `wt-r82-quadseam`, library
`wt-quadseam-lib`, base library `rellib-r20`. Brief:
`round-briefs/r082-quadseam.md`.

Status: fix round done (review dev/reviews/2026-11-11-quadseam.md, FIX).
The section "Fix round" at the end wins over every section before it.
Items 828, 3466 (= 3614) and 3477 fixed; 561 recorded, not
started. All suites, sessions, OpenBLAS and the check: "Runs" at the end.

Logs are in `dev/quadseam-log/` (local, gitignored). Every number
below is pasted from them; the script that made it is named beside it.
"base" is rellib-r20, "lane" is wt-quadseam-lib.

## What changed

- `R/fit.R`
  - `quad_fit()`: escapes a family's declared stationary point on the
    Laplace warm start (`escape_stationary()` on `lap_obj`) before the
    Gauss-Kronrod rule is anchored, and carries the escape's record to
    the returned optimum. After the calibration loop it calls the new
    `quad_settle()`, and returns a `settle` closure for later moves.
  - New `quad_settle()`: tapes the rule again at the estimate; while the
    rule taped there has a gradient above `grad_tol`, optimizes it from
    there and tapes again. A round is kept only when the rule anchored
    at its end is no worse than the one anchored at its start (the
    nested route's acceptance rule), at most `quad_settle_rounds` (10)
    rounds; unsettled warns with class `frmtmb_quad_unsettled` (shared
    with the nested route). Returns the tape anchored at the estimate
    and records `opt$quad_settle = c(rounds, gain, grad)`.
  - `fit_assembled()`: after `escape_stationary()` and `mo_search()` on
    the one-block route, a moved estimate is settled again
    (`qf$settle()`). The GK `integrate` spec is now `quad_gk_integrate()`.
  - New `quad_reanchored_fn(fit)`: the objective at an outer vector by a
    rule taped there (template: the fit's estimates and conditional
    modes; nll rebuilt from `build_objective()` plus the fit's prior),
    memoized per point in `fit$cache`. NULL off quadrature and on the
    nested route.
  - New `se_edge_probed(fit)` (FALSE only on the nested route) and
    `se_probe_fn(fit)` (the re-anchored function, else `fit$obj$fn`).
  - `headroom_edge()`: gated on `se_edge_probed()` instead of
    `fit$quadrature`.
  - `frm()` docs (`@param quadrature`): the re-anchoring and its warning.
- `R/se-check.R`: `se_edge_sd()`, `se_boundary_names()`,
  `se_boundary_stop()` and `se_sd_gain_up0()` gate on
  `se_edge_probed()` instead of `fit$quadrature`; `se_sd_gain_up0()`
  and `se_at_edge()` evaluate through `se_probe_fn()`.
- `tests/testthat/test-quadrature-defects.R`: three new tests (one per
  item), each with the case where the guarded thing is absent; the
  non-Gaussian loop takes `allow_boundary()` (a one-block Beta fit with
  its sd at zero now gives the boundary message).
- `tests/testthat/test-se-tier.R`: the RQ1 quadrature test is retitled
  ("an sd stuck short of the maximum is no boundary") and its comment
  says why it still holds; assertions unchanged.
- `NEWS.md`: a `## Quadrature` section with three bullets under
  `# frmtmb (development version)`.
- `man/frm.Rd` regenerated.
- `dev/test-backlog.md`: CLOSED lines under 828, 3466, 3477 and 3614;
  four items filed under "Filed in the round of 2026-11-11" / "From
  lane quadseam" (`dev/quadseam-backlog-edit.R` made the edit).

How my se-check.R hunks relate to lane boundfam's (its working tree
read at the end of this lane, `git diff R/se-check.R`): my hunks are
the gates at the top of `se_edge_sd()`, `se_boundary_names()`,
`se_boundary_stop()` and `se_sd_gain_up0()` and the two `obj$fn` reads
in `se_sd_gain_up0()` and `se_at_edge()`. boundfam changes
`se_boundary_stop()` too, but 8 lines below my gate (the `ok <-` line
takes "limit"), so the hunks do not overlap; its other se-check.R hunks
(cov_from_hessian() at line 383, new code after se_sd_at_zero(),
se_report(), se_relabel(), se_lost_clauses()) touch none of my lines.
Its new family-parameter verdict in cov_from_hessian() keeps its own
`!isTRUE(fit$quadrature)` gate; after the merge, whether that verdict
should run on the one-block quadrature route through `se_probe_fn()` is
for the consolidation to decide. boundfam does not touch R/fit.R.

## Item 828: the escape on a frozen rule

Reproducer ported to the round's libraries: `dev/quadseam-punch1c.R`
(the first half of `dev/skewinit-punch1c.R`; wt-skewinit's acceptance
seed 1, n 200, plus a 20-level `(1 | g)`, started at alpha's stationary
point). Logs `punch1c-base.txt`, `punch1c-lib.txt`:

    base  ARM quadrature=FALSE default -259.338935424 rescued -259.338935434 d -1.019e-08 esc 24.272
    base  ARM quadrature=TRUE  default -259.338935428 rescued -259.348966632 d -1.003e-02 esc 24.262
          rescued: grad on its own tape 2.162e-04
    lane  ARM quadrature=TRUE  default -259.338935428 rescued -259.338935438 d -1.018e-08

The base fit is stationary on its own tape (2.2e-4 < grad_tol) and
1.0e-2 short: the defect as filed, seen on 0.81.0.

The order, measured (`dev/quadseam-escorder.R`, one process, arms
interleaved per seed, seeds 1 to 10 of the acceptance stream; log
`escorder-lib.txt`). Arm "settle" removes the Laplace-side escape from
`quad_fit()` (the escape then runs on the rule frozen at the stall, and
`quad_settle()` repairs after it); arm "lane" is the build:

    SEED  1 lane   gap -1.018e-08 esc_gain 24.272
    SEED  1 settle gap -3.118e-05 esc_gain 24.262
    SEED  5 lane   gap  2.842e-14 esc_gain 23.405
    SEED  5 settle gap -7.560e-12 esc_gain 23.392

Seeds 1 and 5 are the two where the escape fires. Re-taping alone
leaves seed 1 3.1e-5 short; escaping before the rule is taped is
needed. The other eight seeds stop off the stationary point and both
arms agree to the digit; seed 6 is 5.5e-2 short on both arms and on
base alike (`dev/quadseam-seed6.R`: Laplace reaches -244.098622627 from
both starts, quadrature -244.100465993 from the default and
-244.155229735 from the stall, the same on base). That one is not the
escape and is filed.

Decision: escape first (on the Laplace objective), and re-tape after
any later move. Both are needed: the escape on the rule still runs
(`fit_assembled()`, unchanged) and a move by it or by `mo_search()` is
settled again.

## Item 3466 (= 3614): a runaway stopped at a false maximum

The filed construction (`dev/rel073-wcens-probe.R`, censored gaussian
with an uncensored excess) is refused by the weights check on 0.81.0
for every excess from 2 to 32 (`dev/quadseam-runaway.R`, log
`runaway-base.txt`), so it no longer reaches the rule. A runaway the
check does not reach: two gaussian responses, `mvbf()`, each with
`(1 | obs)`, weights 1.5, 80 rows, seed 20261111. Same script:

    base  mvbf quadrature  logLik -24.1244  log sigma  -4.008  conv 0  warn 2 (the obs confound notes only)
          exact at est -18.4577  log sigmas -1: 61.5423  -3: 221.5423
    base  mvbf Laplace     logLik 2284.4888 log sigma -37.811  conv 1
    lane  mvbf quadrature  logLik 1072.5789 log sigma -25.296  conv 1  warn 5
          frmtmb_quad_unsettled, frmtmb_conv_code, frmtmb_conv_gradient

"exact" is the closed-form marginal log-likelihood written in the
script from the model alone. At base's estimate the frozen rule is
5.7 off the integral, and the integral rises by 80 and 240 as both log
sigmas fall by 1 and 3: no maximum, and base said nothing. The lane
build warns. The censored design without excess (E 0) settles in
silence on both builds (logLik -107.1113).

Decision between the two routes the brief names. One-block quadrature
does not move to the Gauss-Hermite route. What protects the nested
route from this is its re-anchor loop (`quad_nest_fit()`), not the
Gauss-Hermite rule: a frozen Gauss-Hermite rule cannot follow a spike
narrower than its nodes either. So the one-block route got the same
loop (`quad_settle()`), and keeps TMB's rule, whose accuracy the
existing tests pin against `glmer(nAGQ = 25)`. A separate runaway guard
was not needed: the loop is the guard, and it is not tied to weights.

Cost of the loop on fits that have a maximum: one tape. Rounds counted
from `fit$opt$quad_settle` on 130 fits that have a maximum:

- `dev/quadseam-acc.R`, log `acc-lib.txt`: 30 fits (Bernoulli groups
  of 4 sd 1.2, binomial and poisson 60 x 6), rounds 0 on 30 of 30; the
  largest gradient of the rule taped at the estimate 8.5e-4.
- `dev/quadseam-setier-quad.R`, log `setier-quad-lib.txt`: 100 fits
  (arm A 80, arm B 20, the latter started at a group sd of exp(-12)),
  "rounds>0 0" in all six cells.

Tapes per fit (`dev/quadseam-timing.R`, `RTMB::MakeADFun` counted with
`trace()`; logs `timing-base.txt`, `timing-lib.txt`): base 2 on every
design; lane 3 on the two designs with interior sds, 11 on the two with
an sd at zero (8 of them the boundary probes of item 3477). Time, min
of 3, base then lane, parallel processes:

    skewn     0.33  1.25   (sd at zero)
    bern_sd0  0.11  0.37   (sd at zero; the item-3477 fit)
    bern_sd1  0.25  0.28
    pois      0.26  0.34

The runaway costs about 10 s per fit with 10 rounds (`dev/quadseam-ulp.R
runaway`: 20 fits in 197.3 s, under the load of two other processes).

Accuracy. The reported log-likelihood is now the rule taped at the
estimate, not the one taped at the Laplace optimum, so it moves by the
rule's own error. Against `glmer(nAGQ = 25)` (`acc-base.txt`,
`acc-lib.txt`, 10 seeds each):

    base  SUM bern n 10 max|d| 3.692e-07   bin 1.104e-08   pois 2.510e-08
    lane  SUM bern n 10 max|d| 8.968e-07   bin 9.794e-09   pois 2.457e-08

The Bernoulli maximum moved, so I read the rule against an independent
reference at each build's own estimate (`dev/quadseam-bernref.R`,
`stats::integrate()` per group, rel.tol 1e-12; logs `bernref-*.txt`):

    seed      1        2        3        4        5        6        7        8        9       10
    base -1.3e-08 -1.4e-08  3.9e-07 -1.2e-08 -4.1e-08 -4.0e-08 -7.6e-09  1.4e-08 -2.8e-08 -7.5e-09
    lane  1.4e-08  1.3e-08 -1.7e-09  1.3e-08 -9.0e-07  2.6e-07  1.3e-08  9.3e-09  4.6e-09  1.3e-08

The estimates agree (seed 5: -0.546529, 0.753033, 1.364176 on both,
`dev/quadseam-bern5.R`); the rule's error at its anchor is 1e-8 typical
and up to 9e-7 at some anchors on both builds, at different seeds. This
is TMB's calibration (curvature by differences with a step of 1, mode
search stopped at a gain of 1e-2), not the re-anchoring. Filed.

Nested reference, `dev/nestquad-ref.R` through `dev/quadseam-nestnodes.R`
(a copy of `dev/nestquad-nodes.R` pointed at the round libraries; logs
`nestnodes-{base,lib}-{bern,zero}.txt`): the nested rule's values are
bitwise equal on both builds at every k (k 11 err 1.301e-04, k 31
6.468e-10, as recorded at 0.80.0). The "zero" case's one-block GK
reference moved by the re-anchoring: gk-ref -1.015e-09 (base) and
-2.070e-09 (lane); nested-gk at sd_inner 1e-4 -1.631e-08 and
-1.525e-08.

## Item 3477: an sd at zero under quadrature

Construction: seed 3 of `dev/fuzzfix-quad-boundary.R`, Bernoulli, 20
groups of 6. `dev/quadseam-boundary.R`, logs `boundary-base.txt`,
`boundary-lib.txt`, and `dev/quadseam-scratch1.R` for base's quadrature
arm (the boundary script stopped there on base at a `vcov()` row name
before it was fixed):

    base quadrature: sd(Intercept) Estimate 0, Est.Error 0.14,
                     interval 0 to Inf; cov.fixed[theta] 8.052334e+06;
                     no message
    lane quadrature: MSG frmtmb_boundary_fit "Boundary (singular) fit:
                     sd_g__Intercept is at zero..."; sd 4.773e-05,
                     Est.Error NaN, Q2.5 NaN; cov.fixed diag 0.03386
                     0.04893 NaN; se_lost theta_1=boundary
    lane Laplace:    the same message, the same NaN, the same
                     vcov(beta) diag 0.0338581 0.0489252

Why the rule was off: the fit's own tape is frozen at the estimate,
and read away from it it is garbage. `dev/quadseam-probe1.R` (base):
at log sd -11.95 (the move toward zero) frozen and re-taped agree to
2.7e-9, but at sds of 0.01, 0.03, 0.1, 0.3 and 1 the frozen tape is off
by +2.5, -6.8, +0.52, +20.6 and +40.7, while a rule taped at each
point costs 0.01 to 0.04 s and tracks the Laplace objective (82.2501
against 82.2501 at 0.01; 86.43 against 86.69 at 1). So the boundary
probes now evaluate a rule taped at each point (`se_probe_fn()`).

The absent case and false alarms, `dev/quadseam-setier-quad.R` (a copy
of the setier reviewer's `dev/setier-rev3-quad.R`, which pinned the
reason for the old skip; logs `setier-quad-base.txt`,
`setier-quad-lib.txt`). 20 seeds per cell of arm A, 5 of arm B:

    base A bin  sd 0.0 | glmer singular 8, no msg 8 | nothing 20
    base A pois sd 0.0 | glmer singular 8, no msg 8 | nothing 20
    lane A bin  sd 0.0 | boundary msg where glmer is not singular: 1 | glmer singular 8, no msg 0 | BOUNDARY 9, nothing 11
    lane A pois sd 0.0 | glmer singular 8, no msg 0 | BOUNDARY 8, nothing 12
    lane A bin/pois sd 0.3 | nothing 20 (both builds)
    lane B bin/pois sd 0.8 | 5 fits short of glmer | SEWARN 5 (both builds)

16 of 16 glmer-singular fits now get the message (0 of 16 on base). Arm
B, the fits stuck near sd = 0 short of the maximum that made the
setier reviewer block the boundary verdict under quadrature, still get
the "flat" warning and no boundary message: the upward probe on rules
taped at each sd finds the gain. The one message where glmer is not
singular (`dev/quadseam-fa.R`, log `fa-lib.txt`, seed 4) is the
Laplace fit's verdict too:

    quad sd 0.002339 ll -639.47737585 [frmtmb_boundary_fit]
    laplace sd 0.002339 ll -639.47737584 [frmtmb_boundary_fit]
    glmer25 sd 0.001346 ll -639.47737580 singular FALSE

The quadrature profile there is flat to 1e-9 from sd 1e-6 to 1e-3
(-639.4773758179 to -639.4773758101). The brief asks for the Laplace
rule under quadrature, and that is what it gives; lme4's isSingular()
threshold (1e-4 on theta) is a different rule.

The setier trap fits (`dev/quadseam-setier-trapq.R`, copy of
`dev/setier-rev3-trapq.R`): identical lines on both builds, all four at
lme4's maximum with no message.

## Tests, seen to fail on base

`tests/testthat/test-quadrature-defects.R` on rellib-r20 (log
`dev/release/quadseam-qdbase.log` and its `-files/`):

    RESULT test-quadrature-defects.R pass=66 fail=6 err=0 skip=0 warn=0
    FAIL [the stationary-point escape runs before the rule is anchored] Expected `gap/gain` < 1e-06. Actual comparison: 0.000413 >= 0.000001
    FAIL [quadrature does not stop a runaway at a false maximum] Expected `unsettled` to be identical to 1L. actual: 0
    FAIL [a quadrature sd at zero is a boundary fit, as under Laplace] Expected `seen` > 0L. Actual comparison: 0 <= 0
    FAIL [...] Expected `unname(sdr$se_lost)` to be identical to "boundary". `actual` is NULL
    FAIL [...] Expected `is.nan(v[[3L]])` to be TRUE.
    FAIL [...] Expected `is.nan(VarCorr(f)[["g"]]$sd[1L, "Q2.5"])` to be TRUE.

On the lane build: pass=72 fail=0 err=0 skip=0 warn=0; test-se-tier.R
pass=112 fail=0.

Each new assertion over 20 data sets moved by ulps (`dev/quadseam-ulp.R`,
v * (1 + 2e-15 k), k 1..20; the response for the skew normal and the
two-response runaway, the covariate for the Bernoulli; logs
`ulp-lib-*.txt`):

    ULP boundary pass 20 of 20
    ULP escape pass 20 of 20
    ULP runaway pass 20 of 20

The runaway ends at a different log sigma on different ulps (-23.2,
-20.7, -10.0 on k 1 to 3), which is why the test asserts only the
warning.

## Not done, and why

- **561** (single-block quadrature refused on a variance near zero,
  size L). Not started. What it would take: TMB's `marginal_gk`
  calibrates each integrand once, with finite differences of step
  `dx = 1` (`gk_config`, TMBad `graph_transform.hpp`,
  `logIntegrate_t::rescale_integrand()`), so an integrand narrower than
  1 is measured badly. Two routes: (a) the one-block rule written on
  the tape as the nested route writes its rule (Gauss-Hermite at the
  Laplace Gaussian of each level, `quad_nest_objective()` reduced to
  one block), whose scale comes from the inner Hessian and not from
  differences; this also removes the anchor-dependent error filed
  above, but changes every one-block value and needs the accuracy
  tables of 0.80.0 for one block, M to L; (b) pass a smaller `dx`
  through RTMB's GK spec, if RTMB forwards it (not checked; TMB's
  `gk_config` has the field). The re-anchoring of this lane makes (a)
  easier: the loop exists.
- **`sd_up_escape()` under quadrature.** The probe now works, but the
  restart would need `quad_settle()` after it; filed.
- **The nested route's boundary verdict.** Needs the Laplace objective
  kept on the fit for the anchor; filed.
- No Newton polish on the one-block route (seed 6 above); filed.

## Defects found, not fixed

All four are filed in `dev/test-backlog.md` under "Filed in the round of
2026-11-11", "From lane quadseam".

## Bump and floor

Behavior changes marked BREAKING in NEWS (reported log-likelihood moves
by the rule's error; a new warning; a new message and NaN SE): a minor
bump. No new dependency; no floor change.

## Runs

All on the lane build (`lib:` line of every per-file log names
`wt-quadseam-lib` for frmtmb; the extensions are rellib-r20's, which call
none of the changed functions: `grep` of `extensions/*/R` for them finds
only compat-table text).

- Every test file of core and all eight extensions, one per process,
  3 at a time, `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`
  (`dev/release/run-par.sh quadseam-all`, log
  `dev/release/quadseam-all.log`): `QUADSEAM-ALL ran 510 of 510`.

      frmtmb           317 files pass 26554 fail 0 err 0 skip 15 warn 0
      frmtmb.coupling   11 files pass   551 fail 0 err 0 skip  5 warn 0
      frmtmb.eam        33 files pass  1854 fail 0 err 0 skip  3 warn 0
      frmtmb.latent     12 files pass   687 fail 0 err 0 skip  2 warn 0
      frmtmb.learn      16 files pass   545 fail 0 err 0 skip  2 warn 0
      frmtmb.ode        11 files pass   549 fail 0 err 0 skip  1 warn 0
      frmtmb.sample     87 files pass  4414 fail 0 err 0 skip  1 warn 0
      frmtmb.spatial     8 files pass   561 fail 0 err 0 skip  0 warn 0
      frmtmb.survival   15 files pass   556 fail 0 err 0 skip  0 warn 0

  The 29 skips are opt-in tiers: test-scale.R in six packages
  (FRMTMB_SCALE_TESTS), test-curve-scale.R, test-fuzz.R (FRMTMB_FUZZ)
  and test-drmtmb-agreement.R (FRMTMB_DRMTMB_FIT_TESTS). The last two
  were run with their gates:
  - test-fuzz.R, `FRMTMB_FUZZ=true` (its specs include a quadrature
    mode): lane pass=2 fail=0 err=0 skip=0 warn=0, base the same
    (`quadseam-fuzzlib.log`, `quadseam-fuzzbase.log`).
  - test-drmtmb-agreement.R, `FRMTMB_DRMTMB_FIT_TESTS=true`: pass=131
    fail=0 err=0 skip=0 warn=0 (`quadseam-drm.log`).
- OpenBLAS with LAPACK (`dev/ciharden-openblas.sh <ver> lapack`; logs
  `dev/quadseam-log/openblas-*.txt`, `openblas32-*.txt`), on the three
  files that test the changed paths:

      0.3.26  test-quadrature-defects.R pass=72, test-quadrature-nested.R pass=69, test-se-tier.R pass=112
      0.3.32  test-quadrature-defects.R pass=72, test-quadrature-nested.R pass=69, test-se-tier.R pass=112
      (each fail 0 err 0 skip 0 warn 0)

- One session per suite, three orders, every package
  (`NOT_CRAN=true bash dev/release/run-sessions.sh quadseam-sess 2`, log
  `dev/release/quadseam-sess.log`): `QUADSEAM-SESS ran 27 of 27`, fail 0
  err 0 warn 0 in all 27:

      frmtmb           forward/reverse/preload files=304/317 pass 23504/23669/23772 skip 239/239/241
      frmtmb.coupling  11/11 pass 551 skip 5 (each order)
      frmtmb.eam       33/33 pass 1848 skip 5
      frmtmb.latent    12/12 pass 687 skip 2
      frmtmb.learn     16/16 pass 474 skip 13
      frmtmb.ode       11/11 pass 549 skip 1
      frmtmb.sample    82/87 pass 4050 skip 18
      frmtmb.survival  15/15 pass 556 skip 0
      frmtmb.spatial    8/8  pass 478 skip 6

  46 LEAK lines, every one tagged `[with load of ...]` (0 untagged): the
  extension registrations of test-families-extreme.R, test-cox-home.R,
  test-sampling.R and the like, and GDAL_DATA / rxTempDir. The ungated
  file counts (304/317, 82/87) are file-level gated skips, as in earlier
  records (0.80.0: 243/253 ungated).
- `R CMD check --as-cran` of core, built and checked in
  `dev/quadseam-check/` (`dev/quadseam-check.sh`, log
  `dev/quadseam-check/check.txt`; R_LIBS lane, rellib-r20, user;
  `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`): `Status: 2 NOTEs`. The V8
  math-rendering NOTE (expected), and the examples timing NOTE for
  `frm_hazard_reads` (user 3.26 s, elapsed 6.92 s) and
  `dharma_residuals` (1.25, 6.16), neither of which this lane touches;
  elapsed far above user time is load (lane-rules.md). The check tree
  and tarball were deleted after reading; check.txt stays. No extension
  changed, so no extension was checked.
- `dev/ciharden-out/` (925 MB) deleted after the OpenBLAS runs.

Process count: `dev/ciharden-openblas.sh` ends by running its own
`ciharden-blascheck.R`; the first build ran it while three of my test
processes were running, so for about 20 s I had four R processes.

## Fix round (review `dev/reviews/2026-11-11-quadseam.md`, verdict FIX)

This section wins over the sections above. Where it and they differ
(the 3466 case, "130 of 130 settle at once", the 9e-7 accuracy claim,
the cost claim), the numbers here are the ones measured on the final
build. Every finding was reproduced on the first-round lane build before
it was fixed. Logs `dev/quadseam-log/fr-*`.

### MAJOR 1: false `frmtmb_quad_unsettled` on fits at the maximum

Reproduced on the first-round build:

    dev/reviews/quadseam-review-gradtol.R lane 1:3 (fr-gradtol-before.txt)
    FIT seed 1 settle 10/0/2.56e-04 secs 32.3 nwarn 1  frmtmb_quad_unsettled "... rose by 0 over them (0 in the last step)"
    FIT seed 2 settle 0/0/6.36e-05  secs  9.4 nwarn 0
    FIT seed 3 settle 10/0/4.14e-04 secs 26.1 nwarn 1
    dev/reviews/quadseam-review-acc.R lane (fr-acc-before-a.txt, -b.txt), default grad_tol
    bern2sd3 35  sgain 0 sgrad 2.85e-3  warn frmtmb_quad_unsettled
    poissd2 29   sgain 0 sgrad 1.06e-3  warn frmtmb_quad_unsettled,frmtmb_conv_gradient

What changed in `quad_settle()` (R/fit.R):

- A round whose optimizer converges (code 0) without a step of 1e-8 in
  the optimizer's units ends the loop as settled. The gradient of the
  rule taped there is then the rule's calibration noise.
- A round that converges and whose two rules (taped at the round's two
  ends) differ by at most `quad_settle_tol` = 1e-5 (absolute, MINOR 3)
  ends the loop as settled, and the estimate is the stop of the rule the
  round optimized, which is stationary there. Reporting the old point
  instead left poissd2 seed 29 with a gradient of 1.06e-3 on its tape
  and a `frmtmb_conv_gradient` warning.
- A round that loses more than that is rejected and retried in a box of
  a quarter of its step (the nested route's trust region); three
  rejections in a row warn.
- The loop warns only when it is still gaining after 10 rounds, after
  three rejections, or when the optimizer stops without converging and
  without a step the rule can resolve. The last case is needed: on
  hurdle_lognormal (`dev/quadseam-runaway2.R`) the runaway's rule
  stopped at "false convergence" with no step and a gradient of 3e5,
  and the no-move rule alone called that settled; on one ulp draw in
  five of the `mi()` runaway the stop moved 1.39e-8 with the two rules
  5.6e-7 apart and a gradient of 3.7e5 (`dev/quadseam-fr-ulp3.R`, log
  `fr-ulp3.txt`). Both now warn.
- `grad_headroom()` and `conv_headroom()` (R/fit.R): a Gauss-Kronrod
  tape's AD Hessian is NaN in every entry on poisson fits, on base and
  lane alike (`dev/quadseam-fr-pois29.R`, `fr-pois29-*.txt`), and the
  gradient verdict then warned "the curvature there is unusable". A
  non-finite `he()` now falls to `fit_outer_hessian()`'s differences,
  as the SE check already does. This took the last 4 false
  `frmtmb_conv_gradient` warnings of the tight sweep (pois200 seeds 1
  and 2 at both tolerances) to 0.

After, on the final build:

    fr-gradtol-lane-final.txt  seed 1 secs 24.4 nwarn 0, seed 2 10.8 0, seed 3 14.9 0
    fr-gradtol-base.txt        seed 1 secs 20.6 nwarn 0, seed 2  8.7 0, seed 3 11.9 0
    fr-acc3-lane-final.txt     bern2sd3 35 and poissd2 29: warn -

The reviewer's tight sweep (10 designs x seeds 1 to 5, `gen()` and
`ref_ll()` of `dev/reviews/quadseam-review-acc.R`, run uninstrumented by
`dev/quadseam-fr-sweep.R`; logs `fr-sweep-lib-final.txt`,
`fr-sweep-base-a.txt`):

    lane grad_tol 1e-4: 50 fits, any warning 0, rounds > 0 on 17, max |err| 5.354e-6, 216 s
    base grad_tol 1e-4: 50 fits, any warning 0, max |err| 1.821, |err| > 1e-3 on 5, 158 s
    lane grad_tol 1e-5: 50 fits, any warning 0, rounds > 0 on 29, max |err| 5.331e-6, 310 s
    base grad_tol 1e-5: 50 fits, any warning 1, max |err| 3.779, |err| > 1e-3 on 6, 252 s

(The reviewer counted 49 and 48 fits; this run fitted all 50 at each
tolerance with no error.) The first-round build warned on 13 and 25.

The default-tolerance table, remeasured on the final build
(`dev/reviews/quadseam-review-acc.R lane bern500,pois200,bern2sd3,poissd2
6:35`, log `fr-acc3-lane-final.txt`, against the reviewer's
`acc3-base.txt`):

    120 fits; max |err| lane 1.064e-3, base 4.483e-2
    lane worse than base by > 1e-6: 19; better by > 1e-6: 17
    warnings: lane 0, base 0; rounds > 0: 9
    seconds summed: lane 283, base 212 (separate runs, not interleaved)

New test "a quadrature fit at its maximum settles in silence"
(test-quadrature-defects.R): bin500 seed 1 and binz seed 4 at
`grad_tol = 1e-5`, and the Bernoulli of the first round at the default,
no warning of any class. The first-round build warned
`frmtmb_quad_unsettled` on both tight fits (reviewer's
`settle2-lane.txt`: rounds 10 each). It passes on rellib-r20, which has
no loop. Over 20 ulp draws: `ULP atmax pass 20 of 20`
(`ulp-lib-fr-atmax.txt`).

### MAJOR 2: the runaway case belongs to wcheck's 3447

The two-response gaussian case is refused on wcheck's build. The test
and the NEWS example now use the case wcheck's review gave
(coordinator's addendum): a gaussian `mi()` response with 3 rows
missing, weights 1.5, `(1 | obs)`, the data of
`frmtmb-wt-wcheck/dev/reviews/wcheck-review-misc.R` part "quad".
`dev/quadseam-runaway3.R`, logs `runaway3-{base,wcheck,lib}.txt`:

    base    quadrature logLik -11.5195 log sigma -4.057 conv 0, no warning but the obs note
            exact (observed rows) at est -3.1786, log sigma -1: 35.3214, -3: 112.3214
    wcheck  the same, bitwise (the check skips mi() by design)
    lane    quadrature logLik 763.1509 log sigma -23.948 conv 1, settle 5 rounds,
            frmtmb_quad_unsettled + frmtmb_conv_code + frmtmb_conv_gradient

The exact log-likelihood of the observed rows is written in the script
from the model alone; at base's estimate it rises by 38.5 per unit of
log sigma. Test "quadrature does not stop a runaway at a false maximum"
rebuilt on it; over 20 ulp draws `ULP runaway pass 20 of 20`
(`ulp-lib-fr-runaway.txt`; 16 of 20 before the code-0 condition on the
noise band, above). hurdle_lognormal (`dev/quadseam-runaway2.R`) also
runs away silently on base and wcheck; not used, since wcheck may add it
to its table.

The runaway test depends on lane wcheck leaving `mi()` out of its
weights check; if a later round adds `mi()` to that check, the test's
case is refused before the rule runs and needs another runaway.

CLOSED lines under 3466 and 3614 reworded: the loop fixes the rule's
false stop on runaways the check does not reach, and the weighted
`(1 | obs)` cases are refused by the weights check.

### MINOR 3: acceptance tolerance

Absolute `quad_settle_tol` = 1e-5 in place of `1e-8 * |f|`; a round
that loses more is rejected and retried in a box (above). The roxygen
no longer says "the kept values fall and the loop cannot cycle". Why
1e-5 and not 1e-6: poissd2 seed 29's first round lost 7.83e-6 between
two rules taped 1.56e-4 apart (`fr-pois29-lib.txt`, verbose stages),
which is the rule's own calibration noise; at 1e-6 it was rejected, the
retry ended on the box's edge, and the fit warned. Rounds whose gains of
2.6e-7 and 2.8e-7 the reviewer saw lower the exact value now end the
loop, at the frozen rule's stop.

### MINOR 4, 5, 6, 7

- NEWS rewritten with the measured accuracy (worst 4.48e-2 to 1.06e-3;
  19 worse and 17 better by 1e-6 of 120), the unclaimed fix (1.82 and
  3.78 at tight tolerances, now 5.4e-6), the SE move (up to 6.8e-4
  relative, the reviewer's `misc-*.txt`), and the cost at the default
  and at tight tolerances. The backlog item on the anchor error carries
  the 120-fit range.
- The 3477 NEWS bullet says "standard deviation at zero" and that a
  coefficient run-off still gets no verdict under quadrature.
- MINOR 7 filed under "Filed in the round of 2026-11-11 / From lane
  quadseam" as a cross-lane item for boundfam's `se_edge_limit()`; not
  implemented.
- Not done: the Laplace-side escape in `quad_fit()` still passes no
  `tally` (the review's cosmetic note).

### Fix-round runs

All on the final build (installed after the last edit of R/; every
`lib:` line names wt-quadseam-lib).

- Changed test files. On rellib-r20 (`dev/release/quadseam-fr-qdbase.log`):
  test-quadrature-defects.R pass=65 fail=6, the six behavioral failures
  of the first round (the runaway test's now on the `mi()` case:
  "`unsettled` actual 0 expected 1"). The new at-maximum test passes on
  base, which has no loop; it failed on the first-round build (above).
  On the lane (`quadseam-fr-lib1.log`): test-quadrature-defects.R
  pass=71 fail=0, test-se-tier.R pass=112 fail=0.
- Full core suite, one file per process, 3 at a time,
  `NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true`
  (`dev/release/quadseam-fr-core.log`): `QUADSEAM-FR-CORE ran 317 of
  317`, pass 26553 fail 0 err 0 skip 15 warn 0; the skips are the
  opt-in tiers (test-scale.R, test-curve-scale.R, test-fuzz.R,
  test-drmtmb-agreement.R). test-fuzz.R with `FRMTMB_FUZZ=true`:
  pass=2 fail=0 (`quadseam-fr-fuzz.log`). The extensions were not
  rerun in this round (the brief asked for core); R/fit.R's changes are
  inside `quad_settle()` and the two headroom functions, which no
  extension calls.
- Core sessions on the 25 files that fit `quadrature = TRUE` plus
  test-se-tier.R (a copy of tests/testthat with the other test files
  removed, deleted after; logs `fr-sess-*.txt`):

      SESSION frmtmb forward files=25/25 pass=3089 fail=0 err=0 skip=2 warn=0
      SESSION frmtmb reverse files=25/25 pass=3199 fail=0 err=0 skip=2 warn=0
      SESSION frmtmb preload files=25/25 pass=3364 fail=0 err=0 skip=3 warn=0

  8 LEAK lines, all `[with load of frmtmb.eam]` from
  test-families-extreme.R; 0 untagged.
- OpenBLAS with LAPACK (`dev/quadseam-fr-openblas.sh`, logs
  `fr-ob<ver>-*.txt`): 0.3.26 and 0.3.32 each test-quadrature-defects.R
  pass=71, test-quadrature-nested.R pass=69, test-se-tier.R pass=112,
  fail 0 err 0 skip 0 warn 0. `dev/ciharden-out/` removed.
- `R CMD check --as-cran` of core (`dev/quadseam-check/check.txt`):
  `Status: 1 NOTE`, the V8 math-rendering NOTE.
- ulp draws on the final build (`ulp-lib-fr-*.txt`): escape 20 of 20,
  boundary 20 of 20, runaway 20 of 20, atmax 20 of 20.
- Item reproducers on the final build: `punch1c-lib.txt` quadrature
  rescued gap -1.018e-08 (esc gain 24.272); `boundary-lib.txt` the
  boundary message, sd 4.773e-05, Est.Error NaN; `runaway3-lib.txt`
  `frmtmb_quad_unsettled`.

Files touched in the fix round: R/fit.R (`quad_settle()`, its
constants, the `frm()` docs, `grad_headroom()`, `conv_headroom()`),
man/frm.Rd, NEWS.md, tests/testthat/test-quadrature-defects.R,
dev/test-backlog.md, dev/quadseam-findings.md, and new dev scripts
`dev/quadseam-fr-*`, `dev/quadseam-runaway2.R`,
`dev/quadseam-runaway3.R`, an updated `dev/quadseam-ulp.R`.

Bump and floor unchanged: minor, no floor change.
