# Lane wt-phase3b: items 3.3, 3.4 and 3.5

Worktree `frmtmb-wt-phase3b`, off main 51bfaa4e (frmtmb 0.62.0,
frmtmb.eam 0.10.0, frmtmb.learn 0.6.0). Nothing under core's `R/` was
changed and none was needed. Every number below was produced by a
script named beside it; the generated blocks are pasted verbatim from
`dev/phase3b-log/`.

## Punch round 2 (2026-09-24): what changed, and what it supersedes

Where this section and the text below disagree, this section holds.
The statements it overturns are marked WITHDRAWN where they stand.

### The drift-intercept under-coverage: a seed artifact, WITHDRAWN

Round 1 filed the drift intercept's 84.8 percent (196 of 231, pooled
over six arms) as a finding about the hierarchical Wiener's Wald
interval. That attribution is WITHDRAWN. The reviewer's numbers:

- An oracle that knows every parameter but the population mean covers
  on 197 of 231 on the same seeds, one more than the fit
  (`dev/phase3b-review2/arm-ideal.log`). The 231 fits rest on 80
  distinct seeds, because every arm draws from seeds 1 upward.
- Seeds 4, 11, 14 and 17 miss in every arm and for the oracle
  (`paired.log`): the realized subject effects there put the sample
  mean far from 0.4, and no interval from those data covers.
- sd of the realized mean of the subject effects is 1.21 times its
  theoretical value over seeds 1 to 40, and 1.004 over seeds 1 to
  10000 (`seed-u.log`).
- On fresh seeds 2001 to 2070, on the base build, the drift intercept
  covers on 65 of 70, 92.9 percent, Wilson [84.3, 96.9], mean SE
  0.0647 against an sd of 0.0657 (`cov-summary-fresh.log`).

The sentence "Read the drift intercept's interval from a hierarchical
fit as narrow" and the pooled figure are removed from `?wiener`.

Rule for this lane's records from now on: a coverage claim uses seeds
that are independent across arms, or states that they are shared.
`?wiener` now says the censoring arms draw from the same seeds.

The same fresh-seed run found one quantity low that is not an
artifact of shared seeds: log sd(log bs | s) covers on 60 of 70,
Wilson [75.7, 92.1], on the base build. It is not caused by items 3.4
or 3.5 and is filed in `dev/test-backlog.md`.

### The survival table holds on its grid

`?wiener` now says the by-band table holds on the 700-bit grid (|v| a
to 120, w 0.001 to 0.999) and gives the worst measured anywhere: the
review's 798 points at 1200 bits (`ref-analyse.log`) reach 1.7e-11
where |v| a is at most 1, 1.46e-9 at |v| a = 120 and 1.6e-9 at 250, all
at w = 0.9999, u = 0.01. The table was not widened: no new reference
was computed for this round.

### The interval mass: a difference of logs in every regime

The round-1 `ddm_rt_linterval_b()` blended an early and a late route
on the lower edge's normalized time. On the review's 511 Rmpfr
intervals (`interval.log`) it was 25 units out on the log where a
weight of 1e-6 fell on a route that had cancelled to the 1e-300 floor.

It now chooses among three spellings by which is accurate:
F(t2) (1 - F(t1) / F(t2)); T(t1) (1 - T(t2) / T(t1)) with T the mass
still to come; and 16-node Gauss-Legendre quadrature of the density
where both of those hold under a thousandth of the mass on their side.
T is taken from the eigenfunction tail wherever its 30 terms have
converged (u past 0.008) and from P - F only before that. The first
rewrite kept the time blend for T and was still 6.2e-3 out at u = 0.17
(the P - F share leaked through); the gate on convergence fixed it.
Measured (`dev/phase3b-interval-check.R`):

| set | points | worst abs error of the log mass | log |
|---|---|---|---|
| the review's, u1 0.05 to 3, \|v\| a to 24 | 511 | 5.7e-14 | `interval-check.txt` |
| this lane's, u1 1e-3 to 0.05, \|v\| a to 80 | 240 | 4.9e-11 | `interval-check-small.txt` |

The 240 are from `dev/phase3b-interval-ref-small.R`, with the
review's 1200-bit tail. Masses go down to 1e-99 of F(t1). The tape
gradient is finite at all 751. "Checked against quadrature to 1e-9" is
replaced by these numbers, and the test block now reads all 751 from
`fixtures/wiener-interval-ref.csv` and asserts 1e-10 on the log.

### Smaller items

- A response outside `contaminant_range` is refused by name, class
  `frmtmb_eam_contaminant_range_error`, the interval upper edge
  included (`ddm_cont_check_range()`; test "a response outside
  contaminant_range is refused, classed").
- Stale messages: `wiener()`'s range validation no longer offers NULL
  for the observed range; `lba(contaminant = NA)` and `rdm()` now say
  "`contaminant` must be FALSE" instead of reporting TRUE.
- The refusal test covers `trunc()` on `lba()`, and `cens()` and
  `trunc()` on `gddm()`.
- Spelling: "centre" and "behaviour" corrected in the R sources, this
  file, the dev scripts, and the two validation logs that print them
  (`cdf-validate.txt`, `cdf-validate3.txt`, header text only); the
  file `seenfail-cdf-single-centre.txt` is renamed to `-center`. Test
  lines over 80 columns are rewrapped.
- `?bandit2arm_delta`'s Sessions section, which every family points
  to, says why `a b a b` is refused and files interleaved contexts as
  a possible later extension.
- Two upstream reports are drafted, not filed, with the minimal
  reproduction `dev/phase3b-rtmb-repro.R` (`rtmb-repro.txt`).
  `dev/phase3b-rtmb-report-pnorm.md` is addressed to RTMB: its
  `pnorm` method calls `distr_log_pnorm()` for `log.p = TRUE`, which
  only loops over `log_pnorm_both()`, RTMB's own atomic (in neither
  TMB's nor RTMB's installed headers). The defect named is that
  atomic's derivative; that it is `exp(log dnorm - log pnorm)` is
  inferred from the x^2 eps / 2 error growth, not read. Fix proposed:
  the Mills ratio for the lower tail.
  `dev/phase3b-rtmb-report-tanh.md` is addressed to TMB: `TanhOp::reverse()` in `TMBad/global.hpp` (line 3534 in 1.9.25)
  divides by cosh^2, which overflows past 710; the fix proposed is
  `dy * (1 - y * y)`. The truth for pnorm is the Mills-ratio
  expansion: the review's own "truth" column,
  `exp(dnorm(log) - pnorm(log.p))`, has the same cancellation as the
  tape. Dropped: the NaN claimed at x = -8.3e9 (did not reproduce) and
  `log(cosh(711))` (plain overflow, the same in base R).
- `dev/test-backlog.md`: the core seam for a log-difference slot, the
  core seam for refusing NA `dec()` on censored rows, and the base
  under-coverage of log sd(log bs | s).

### Timeline of round 1, corrected

Round 1 said `phase3b-lib` was installed at "14:3x" and that every log
postdated the last edit at 14:33. Both are wrong: `phase3b-lib` was
installed at 14:29:50, and the library match and bitwise battery ran
before the last edit, `wiener-family.R` at 14:33:38 (roxygen text
only). The checks below were rerun after this round's final edits.

### Checks, punch round 2

Final build: `phase3b-lib`, both packages installed at 19:59:24 and
19:59:29 (`install-r2-*.log`). Every source file under
`extensions/frmtmb.eam` and `extensions/frmtmb.learn` predates it
except `frmtmb.learn/NEWS.md` (20:03:12, one sentence), which no test
reads. The R CMD check of frmtmb.learn built its tarball after that
edit.

- Library match, 19:59:42 (`libmatch-r2.txt`): 220 of 220 frmtmb.eam
  and 50 of 50 frmtmb.learn source objects are identical to the
  installed build.
- Bitwise battery, 19:59:51 to 20:00:01 (`bitwise-compare.txt`): 68 of
  68 quantities over 19 models are identical to the base build.
- Suites, 20:00:29 to 20:37:46 (`dev/phase3b-final-suites3.sh`, one
  file per process, labels `plain-r2` and `gated-r2`):

| run | files | expectations | pass | fail | error | skip |
|---|---|---|---|---|---|---|
| eam plain, new | 28 | 1730 | 1726 | 0 | 0 | 3 |
| eam gated, new | 28 | 1730 | 1726 | 0 | 0 | 3 |
| learn plain, new | 15 | 443 | 430 | 0 | 0 | 13 |
| learn gated, new | 15 | 507 | 501 | 0 | 0 | 2 |
| eam plain and gated, base | 28 | 1730 | 1659 | 67 | 0 | 3 |
| learn plain / gated, base | 15 | 443 / 507 | 409 / 480 | 21 | 0 | 13 / 2 |

  The base build fails only in `test-wiener-cdf.R` (36),
  `test-contaminant.R` (28), `test-surface.R` (2), `test-family.R` (1)
  and `test-session.R` (21).
- Seen failing (`seenfail-punch2-base.txt`,
  `seenfail-punch2-round1.txt`): on the round-1 build, only the new
  interval block, the `lba(contaminant = NA)` assertion and both
  assertions of the range refusal fail in the two files. On the
  released build, every block in both files fails.
- `R CMD check --as-cran` with vignettes and the manual: frmtmb.eam
  `Status: 1 NOTE` (V8 math rendering), 20:17:51, in-check tests
  `FAIL 0 | WARN 1 | SKIP 6 | PASS 1715`. The warning is the
  large-gradient warning of the `wiener_gng` recovery test, which this
  lane does not touch. frmtmb.learn `Status: OK`, 20:23:49,
  `FAIL 0 | WARN 0 | SKIP 15 | PASS 422`. `testthat.Rout` is copied to
  `check-*-testthat.Rout`; round 1's copies are kept as `r1-check-*`.

One test changed for a reason other than a new assertion: "lambda
recovers on contaminated data" had `contaminant_range = c(0.35, 4)`
over data with 11 diffusion rows below 0.35, and the new refusal
caught that. The window now starts at the non-decision time, 0.3.

## Punch round 1 (2026-09-24): what changed, and what it supersedes

Where this section and the round-0 text below disagree, this section
holds; the round-0 statements it overturns are marked WITHDRAWN where
they stand. Final build: `phase3b-lib`, checked against the worktree
sources by `dev/phase3b-libmatch.R` (217 and 50 source objects, 0
differ; `dev/phase3b-log/libmatch-final.txt`).

### B1, the staircase in lambda: fixed

`exp(lc[["l"]] + d - m)` became `exp(lc[["l"]] + (d - m))` in
`ddm_cont_lpdf()` and `ddm_cont_lccdf()`. The reviewer's diagnosis
reproduced: on the old spelling, a central difference at h = 1e-8
disagrees with one at h = 1e-4 by 4.76 on a 400-row objective with
rows below ndt (`test-contaminant.R`, the "smooth in lambda" block;
`dev/phase3b-log/seenfail-b1-old-spelling.txt`), and passes on the fix.
The round-0 "False convergence, and what it is not" section is
WITHDRAWN: it was this defect. At 30 x 400, arm `contfix` (window
given): nlminb code 0 on 40 of 40 fits against 14 of 60 on the round-0
build; `lambda` 0.0499, covering on 39 of 40. Arm `cont` (the old
observed-range window, now only reachable by passing it): code 0 on
34 of 34, `lambda` 0.0451 covering on 19 of 33, which is the window's
bias and not B1 (`dev/phase3b-log/summary2-all.txt`).

### M1, the survival at large |v| a: rewritten

The review was right: the small-time survival was `1 - F` floored at
1e-300. `R/wiener-rtcdf.R` replaces the whole response-time
distribution: every quantity is a log, the survival's small-time route
is the method of images applied to S directly (the killed driftless
density integrated against the Girsanov weight, each image a normal
mass taken from the tail that keeps its digits), the large-time route
is the eigenfunction series with the k = 1 term and the larger boundary
exponent factored out, and the defective function of each boundary is
built the same way. My own reference over the reviewer's ranges
(`dev/phase3b-cdf-reference-v3.R`, 700 bits, 3150 rows: v -20..20,
a 0.3..6, w 0.001..0.999, u 1e-3..10; images against eigenfunctions
agree to 2.2e-181 on S, and F_l + F_u + S = 1 to 3e-210):

```
  band        S        F       Fl       Fu rows
   <=1 1.65e-12 2.84e-14 7.69e-13 7.71e-13  630
   <=5 2.92e-12 2.84e-14 7.10e-13 7.09e-13  700
  <=12 3.76e-12 2.84e-14 9.38e-13 9.38e-13  560
  <=24 7.41e-12 2.84e-14 6.25e-13 6.25e-13  420
  <=48 1.84e-11 2.84e-14 4.19e-13 4.19e-13  420
  <=72 4.57e-11 4.26e-14 1.14e-13 1.14e-13  140
 <=120 7.98e-11 2.13e-14 1.99e-13 2.27e-13  280
smallest reference log S reproduced: -72062.7 (got -72062.7)
```

Errors are on the log (the relative error of the quantity). The
survival's error grows with |v| a and is 8.0e-11 at 120; the
unqualified "2.4e-11" is WITHDRAWN. `?wiener` carries the table by band.
Blend centers, from the sweep in `dev/phase3b-log/cdf-validate3.txt`:
F 0.1, S 0.05, defective 0.2 (no single center serves all three).
`RWiener::pwiener()` per boundary is 3.5e-9 from the reference and this
package 2.4e-14; `WienR::pWDM()` 2.0e-13. Tape gradients match central
differences to 1.2e-9 at 60 random points.

Three defects the recovery arm found in the rewrite, all fixed and
pinned in `test-wiener-cdf.R` (seen failing on the lib5 build,
`seenfail-punch1-cdf-lib5.txt`): RTMB's `pnorm(log.p = TRUE)` has a
wrong and then NaN tape derivative past |x| of about 1e5 (4.85e8 at
x = -2e8 where the truth is 2e8), so decision times are held at
u = 1e-10; `ddm_lsinhc()`'s regularizer has a second derivative of
1e10 at zero drift, every fit's start, which made the inner Hessian
-1e15, so the new code uses a Taylor branch (`ddm_lsinhc_s()`;
`wiener_gng()` keeps the old one and does not move); and RTMB's second
derivative of `tanh` is NaN past about 700, so the blends clamp it at
40 (`ddm_tanh_s()`). An interval late in the distribution had F(y2)
and F(y) equal to their error and the difference went negative, so an
interval's mass is formed directly (`ddm_rt_linterval_b()`), checked
against quadrature to 1e-9 (SUPERSEDED in round 2: see the interval
mass above).

Memory: frmtmb calls `lcdf` and `lccdf` on every row. Formed on all
12,000 rows the tape did not fit: every left/interval fit and many
right-censored ones ended in `std::bad_alloc`. Both slots now work on
the rows frmtmb reads (`ddm_on_rows()`), which changed no value
(logLik to 1e-13 across the three builds,
`dev/phase3b-log/cleft-lib7-vs-lib8.txt`, `cleft-lib8-vs-lib9.txt`;
pinned by a test seen failing on lib7). One fit now peaks at 3.0 to
3.7 GB (`dev/phase3b-peakmem.ps1`).

### Decision (b), the known boundary on left and interval censoring

A left- or interval-censored row is scored with the defective
distribution function of the boundary `dec()` names; a right-censored
row with the survival over both, its `dec()` unread. Left or interval
censoring with `trunc()` is refused by name, because frmtmb forms both
from the one `lcdf` slot. The two defective functions sum to the
marginal F to 1.4e-14. A test flips `dec()` on left-censored rows and
the log-likelihood moves; a test flips it on right-censored rows and it
is `identical()`. All four codes match the hand-written RWiener
likelihood within 1e-10.

Not implementable from the family: "refused by name when `dec()` is
missing". frmtmb's `na.action` drops a row with an NA `dec()` before
the family sees it, with the message "1 row removed because of missing
values"; the family cannot refuse it. That is a core seam.

Recovery at 30 x 400, arm `cleft`: trials faster than 0.45 s
left-censored (29 percent), every fifth other trial interval-censored
into its 100 ms bin (14 percent), 60 replicates, code 0 on 60 of 60,
coverage 90.9 to 100 percent on every quantity but ONE: the drift
intercept, 47 of 55, 85.5 percent, whose Wilson interval [73.8, 92.4]
EXCLUDES 95. Arm `cens` (right censoring at 1 s, rerun on the new
survival), 40 replicates: code 0 on 40 of 40, and again the drift
intercept alone is low, 33 of 38, [72.7, 94.2]. WITHDRAWN in round 2
from here to the end of the paragraph (a seed artifact; see above).
It is not the
censoring: the same quantity is low in every arm this round, including
the two with no censored row at all, `contfix` 34 of 40 [70.9, 92.9]
and `collapse` 23 of 29 [61.6, 90.2], and the round-0 `cens` arm had it
at 88 of 96 [84.4, 95.7]. It is the population intercept of a drift
with a subject deviation of 0.35 over 30 subjects, the quantity
`?wiener`'s hierarchical recovery table already reports at 54 of 60,
90.0 percent. Each arm's count is small (29 to 55); pooled over the
six arms of this round the drift intercept covers on 196 of 231, 84.8
percent. That is a finding about the hierarchical Wiener's Wald
interval on the drift intercept, not about items 3.4 or 3.5, and it is
filed rather than fixed. Full blocks in
`dev/phase3b-log/summary2-all.txt`.

`rdm()`: the analogue, P(T <= t, winner j), is an integral of one
accumulator's density times the others' survivals with no closed form,
so it is NOT built; `?rdm` now says left and interval rows discard the
known winner.

### Decision (a), the contaminant's window

`contaminant_range =` given: used. Not given: a model with
`trunc(ub = )` uses the fastest response to the deadline; any other
model is refused with class `frmtmb_eam_contaminant_range_error`; a
`trunc(ub)` that differs between rows is refused. Top of the window
under a deadline, measured on my own design (`dev/phase3b-cont-window.R`,
one subject, 4000 trials, 25 seeds, contaminants uniform on [0, 5],
5 s deadline): the deadline gave 0.0508 and the slowest response
0.0518, truth 0.0496, both 25 of 25; the deadline is the one that does
not move with the sample, so it is used. Without a deadline, the
observed range gave 0.0212 and 2 of 25, the true window 0.0498 and
25 of 25. Rerun on the final build (`contwin2.txt`): same numbers. At
30 x 400 with a 3 s deadline and the default window (arm `contdl`,
44 replicates): `lambda` 0.0306 against a recorded share of 0.0307,
covering on 35 of 36. Precedents: I did not check HDDM, DMAT or
Ratcliff and Tuerlinckx against their sources; `?wiener` says the
attributions are the review's. Filed, not changed: `gddm(lapse =
"uniform")` spans [0, t_max].

Collapse at 30 x 400 with the window `c(0, 5)` and no contaminant,
30 replicates: `lambda` at the link edge on 25, `diagnose()` named 27
(all 25 at the edge among them), log-likelihood the plain family's to
8.9e-07 at the edge.

### Minor 1: a session label reused in non-adjacent runs is refused

Where the trial column orders a subject's sessions (numbers unique
across them), a label in two runs that are not adjacent is refused by
name. With numbering restarted per session there is no order between
sessions and nothing to check. Seen failing on the round-0 build
(`seenfail-minor1-session-runs.txt`).

### Minor 2: refusals that name the family

`lba(contaminant = TRUE)`, `rdm(contaminant = TRUE)`, and `cens()` or
`trunc()` on `lba()` and `gddm()` are refused with "<family>(): ... is
not built for this family". Both packages' NEWS carry the scope gap.

### Minor 3: the round-0 first-pass errors, corrected

Of the 30, 11 were `std::bad_alloc` and 19 NA/NaN gradients
(`dev/phase3b-log/errtimes.txt`); the round-0 text lumped them and
tied both to memory without evidence. The NaN gradients were not all
B1: 3 were in the `cens` arm, which has no contaminant. None
reproduced: every one refitted cleanly on a sequential rerun of the
same seed on the same build, which says the failure was not a property
of the data or the model; the cause is not established. They fall in
the same two windows as the bad_alloc errors (00:01 to 00:47 and 02:09
to 02:52), and in this round the machine went down at 09:24 with
eleven of this lane's R processes and other lanes' running, physical
memory at 0 GB free; that is the evidence for memory pressure, and it
is circumstantial.

### Checks, punch round 1

All on the final build (`phase3b-lib`, installed 14:29:50 on
2026-09-24; CORRECTED in round 2 from "14:3x";
0 of 217 and 0 of 50 source objects differing from the worktree). The
last source edit is `wiener-family.R` at 14:33:38 (roxygen text). The
suite and check logs postdate it; the library match and bitwise
battery do NOT (CORRECTED in round 2, where both were rerun).

| run | files | expectations | pass | fail | error | skip | log |
|---|---|---|---|---|---|---|---|
| eam plain, new | 28 | 1725 | 1720 | 0 | 0 | 3 | 14:51 |
| eam gated, new | 28 | 1725 | 1720 | 0 | 0 | 3 | 15:14 |
| learn plain, new | 15 | 443 | 430 | 0 | 0 | 13 | 14:53 |
| learn gated, new | 15 | 507 | 501 | 0 | 0 | 2 | 15:19 |
| eam plain/gated, base | 28 | 1724 | 1659 | 61 | 0 | 3 | |
| learn plain/gated, base | 15 | 507 / 443 | 480 / 409 | 21 | 0 | 2 / 13 | |

Base against new differ only in `test-wiener-cdf.R`, `test-contaminant.R`,
`test-surface.R`, `test-family.R` and `test-session.R`. The gated runs
skip only the scale tier (`FRMTMB_SCALE_TESTS`), which
`dev/release/run-gated.ps1` does not set; the Stan identity tier ran,
StanHeaders 2.32.10. Bitwise: 68 of 68 quantities over 19 models
identical to the base (`bitwise-compare.txt`), now including `lba()`,
`rdm()` with right censoring and `wiener_gng()`.

`R CMD check --as-cran`, vignettes and manual: frmtmb.eam `Status: 1
NOTE` (V8 math rendering), 15:46, in-check tests `FAIL 0 | WARN 2 |
SKIP 6 | PASS 1709` (both warnings are large-gradient warnings from
recovery tests, one of them the untouched `wiener_gng` one);
frmtmb.learn `Status: OK`, 15:50, `FAIL 0 | WARN 0 | SKIP 15 |
PASS 422`. `testthat.Rout` copied to `dev/phase3b-log/check-*`.

### The crash at 09:24, and the recovery arms

Every recovery RDS read back after the crash (0 unreadable,
`rds-check.txt`). Missing seeds were refitted by
`dev/phase3b-recovery-driver4.sh` with three workers, each waiting for
5 GB free before a seed; no fit on the final build failed. Each record
names the library it was fitted with; the arms mix builds whose values
agree (above), and `summary2-all.txt` prints the count per build.
Scaled down from the first plan: cleft 60, contfix 40, contdl 44,
cont 34, collapse 30, cens 40 replicates.

## Answer

- **3.3, `session =` on every frmtmb.learn family: done.** The value
  store restarts at each session's first trial; the subject stays the
  unit random effects and `frm(importance =)` group on. A two-session
  objective equals the sum of the two single-session objectives at the
  same parameters with relative difference 0. Every single-session fit
  is bit-identical to 0.6.0. Over 202 replicates at 40 subjects by
  2 x 100 trials, no Wald interval's Wilson bound excludes 95 percent;
  ignoring the boundary instead covers `tau` on 74.8 percent.
- **3.4, censoring on `wiener()`: done for `wiener()` only.** Against a
  400-bit reference, 4.0e-14 relative on the distribution function
  (down to 2.9e-68) and 2.4e-11 on the survival (down to 2.5e-133)
  [WITHDRAWN: on that grid only; see punch round 1, M1];
  against `RWiener::pwiener()`, 4.5e-11 absolute, all of which is
  RWiener's own error. At 30 x 400 with 16 percent right-censored, 100
  replicates, every coverage's Wilson interval contains 95. `gddm()`
  and `lba()` are NOT done.
- **3.5, `wiener(contaminant = TRUE)`: done for `wiener()` only, and
  the plan's range is the wrong default.** A uniform over the OBSERVED
  range lets the slowest diffusion trial set the contaminant's
  density: at 30 x 400 with 5 percent contaminants `lambda` came back
  at 0.0451 against 0.05 and its Wald interval covered on 57 of 99
  fits. With the new `contaminant_range =` set to the range the
  contaminants came from, 0.0499 and 58 of 60. The collapse case is
  asserted: on clean data at 30 x 400, `lambda` ran to its link edge
  on 49 of 60 fits, `diagnose()` named every one, and each
  log-likelihood was the plain family's to within 1e-06 on the 47 where
  the plain fit to the same data also returned.
- WITHDRAWN in punch round 1 (see Minor 3): One environmental hazard, not a defect: 30 of 340 first-pass eam
  recovery fits stopped with `std::bad_alloc` or an NA/NaN gradient
  while the machine had 2.5 GB of 31.7 GB free. All 30 refitted
  without error when rerun sequentially on the same build, and six
  successful seeds rerun as controls reproduced their log-likelihoods
  to the last printed digit and their coefficients `identical()`.

## 3.3 `session =`

### Construction

`ln_pack_sessions()` (`frmtmb.learn/R/family.R`) lays the recursion out
over one SEQUENCE per (subject, session) pair, subject level first and
session within it, so `init()` runs at each sequence's first trial.
The engine is unchanged: it always walked sequences, and a sequence was
a subject. The block keeps `group` and `subject` as the subject and
adds `seq_subject`. `ln_group_sum()` sums per-sequence log-likelihoods
into per-subject ones, replicate-major, for `loglik_group`; without
`session =` it returns its input untouched, which is why a
single-session fit is bit-identical. `loglik_row` needs nothing.

Trial numbers need be unique within a (subject, session) only.
`frm_value_trace()` gains `session` between `subject` and `trial`;
`frm_task_simulate()` reads `session =`. All eight constructors take
`session = NULL` after `trial`, by name.

### Validation (`dev/phase3b-session-check.R`, seed 331; `test-session.R`)

```
session= logLik        -604.949529869049
subject=id:session     -604.949529869049
identical logLik: FALSE
max |coef diff|: 4.440892e-16
objective, 2-session fit at its optimum: 604.949529869049
sum of the two single-session objectives at the same pars: 604.949529869049
relative difference: 0
no session (store carries), logLik: -611.362436851688
first-trial rows: 20; all q1 == 0 and q2 == 0: TRUE
```

The relabelled fit (`subject = id:session`) is the same model and not
`identical()`: it sums the same terms in another order.

`test-session.R` also asserts that `loglik_group` returns one value
per SUBJECT in level order, each the sum of that subject's rows over
both sessions; that the stacked (importance) layout stays
replicate-major; that `frm(importance =)` runs with `(1 | id)` and
returns one effective sample size per subject; that
`frm_task_simulate()` draws the same data as the relabelled design;
and that `rlddm()` restarts too.

### Recovery and coverage (`dev/phase3b-session-recovery.R`)

Seeds 1 to 202; 40 subjects x 2 sessions x 100 trials;
`bandit2arm_delta`; `alpha = plogis(qlogis(0.35) + u_id)`,
`u_id ~ N(0, 0.5)` shared by a subject's sessions; `tau = 3`; each
session drawn from a fresh store. The same data fitted with and without
`session =`. `dev/phase3b-log/summary-session.txt`:

```
session: 202 files, seeds 1 to 202, gaps: 

== with_session: 202 fits, 0 errors, code 0 on 202
quantity                   truth      mean     mcse   covered    rate  Wilson 95
logit alpha              -0.6190   -0.6047   0.0072  191/202    94.6%  [90.5, 96.9]
log tau                   1.0986    1.0987   0.0016  192/202    95.0%  [91.1, 97.3]
log sd(alpha | id)       -0.6931   -0.7634   0.0144  197/202    97.5%  [94.3, 98.9]

== without: 202 fits, 0 errors, code 0 on 202
quantity                   truth      mean     mcse   covered    rate  Wilson 95
logit alpha              -0.6190   -0.5893   0.0075  184/202    91.1%  [86.4, 94.3]
log tau                   1.0986    1.0709   0.0018  151/202    74.8%  [68.3, 80.2]
log sd(alpha | id)       -0.6931   -0.7330   0.0151  197/202    97.5%  [94.3, 98.9]
```

The log standard deviation is 0.070 low with `session =` (4.9 Monte
Carlo standard errors) and still covers at 97.5 percent. ML shrinkage
at 40 groups, `0.5 log(1 - 1/40) = -0.013`, is a fifth of that; the
rest is consistent with the Laplace caveat `?frmtmb.learn` documents
for binary choices, and this study does not separate the two.
Without `session =`, `tau` is 2.7 percent low and covers 74.8 percent.

## 3.4 censoring on `wiener()`

### Construction (`frmtmb.eam/R/wiener-cdf.R`)

The distribution function is of the response time over BOTH
boundaries, `F = F_lower + F_upper`, the upper one the lower reflected.
That is the only reading under which core's `cens()` arithmetic is
right: core scores a right-censored row as `log(1 - F)` or from
`lccdf`, and a boundary-specific defective distribution function tends
to the boundary probability rather than to one. A censored trial
reached no boundary, so its `dec()` is not read; `test-wiener-cdf.R`
flips it on every censored row and gets an `identical()`
log-likelihood.

- Small time, `ddm_cdf_small()`: the existing image series twice.
  Gives F directly.
- Large time, `ddm_surv_large()`, new: the large-time density
  integrated term by term from t to infinity, the two boundaries'
  remaining masses collapsed into one sum. Gives S directly.
- `ddm_rt_lcdf2()` blends the logs in `log(u)` as the density does,
  with a SEPARATE center per output: F at `u0 = 0.1`, S at `0.02`.
- `lcdf` returns `exp(lF)`, because core's slot is a probability;
  `lccdf` returns `lS`. Both pass through `ddm_ndt_install()`'s wrapper,
  so `ndt_group()` bounds reach them.
- Refused under `variability =`, by name, from `valid_y`. The slots are
  filled with the refusal because core checks that an `lcdf` exists
  before it runs the family's validator, and its generic message names
  neither the family nor the reason.
- `wiener` accepts `cens` and `trunc`; the two compat rows read
  `works`. Two existing assertions pinned the old refusal and were
  changed: `test-family.R` (`trunc()` status) and `test-surface.R` (the
  refusal block now fits `trunc()` and asserts the variability
  refusal).

### The reference

`dev/phase3b-cdf-reference-v2.R`, merged by `dev/phase3b-cdf-merge.R`:
the large-time series in Rmpfr at 400 bits, summed until a BOUND on
the next term is below 2^-400 of the total, over the 1125-row grid
`test-density.R` pins the density on and a 192-row tail extension
(w = 0.1 and 0.9, v = +-5, t to 15). `dev/phase3b-log/cdf-reference-400.txt`:

```
400-bit reference: 1317 rows (1125 density grid, 192 tail)
rows where the small-time image series was also summed: 969
worst relative disagreement of the two 400-bit series: 3.31e-16
smallest F: 2.87e-68; smallest S: 2.49e-133
```

The two series are separate derivations. The 3.3e-16 is the check's,
not the reference's: the image series received `1 - w` rounded to
double, and the ten worst rows are all `w = 0.2`.

**A defect in my first reference, caught before use.** The first
version stopped when a TERM fell below tolerance, and at `w = 0.5`
`sin(k pi w)` is exactly zero at every even k, so it stopped at k = 6
and reported survivals above one (F down to -0.143). Its merge printed
that. Its outputs are kept as
`dev/phase3b-log/cdf-ref-400-chunk*-DEFECTIVE-earlystop.csv`. The
1000-bit script `dev/phase3b-cdf-reference.R` has the same defect and
had not finished when this was written; nothing reads its output.

### Validation (`dev/phase3b-cdf-validate.R`, `dev/phase3b-log/cdf-validate.txt`)

```
== blend center sweep, max relative error over all rows ==
(one center for both columns here; the shipped function takes F at
 ddm_rtcdf_u0F and S at ddm_rtcdf_u0S, each from its own column)
    u0    F max rel    S max rel  F>1e-10  S>1e-10
 0.005    1.000e+00    3.155e-08       92       11
 0.010    2.164e-04    9.556e-11       50        0
 0.015    2.515e-07    9.330e-11       17        0
 0.020    2.080e-09    2.431e-11        7        0
 0.030    4.291e-11    3.530e-11        0        0
 0.050    3.327e-13    1.759e-07        0        4
 0.070    1.347e-13    4.793e-05        0       10
 0.100    3.997e-14    1.813e-02        0       12
 0.150    3.798e-14    1.000e+00        0       27
 0.200    3.798e-14    1.000e+00        0       45
 0.350    3.798e-14    1.000e+00        0       72

== the shipped function, ddm_rt_lcdf2() ==
density grid: F max rel 4e-14, S max rel 2.43e-12, F max abs 2.29e-14
tail grid: F max rel 2.82e-14, S max rel 2.43e-11, F max abs 2.81e-14

== established implementations, against the same 400-bit truth ==
frmtmb.eam vs RWiener::pwiener, density grid, max abs: 4.46e-11
frmtmb.eam vs RWiener::pwiener, all rows, max abs: 4.46e-11
RWiener vs the 400-bit truth, max abs: 4.46e-11; max rel: 2.01e+55
frmtmb.eam vs WienR::pWDM(precision = 1e-12), all rows, max abs: 7.53e-14
WienR vs the 400-bit truth, max abs: 7.53e-14; max rel: 7.53e-14
frmtmb.eam vs the 400-bit truth, max abs: 2.81e-14

== F integrates the density ==
40 random (v, a, w, t1, t2): max |F(t2) - F(t1) - integral| / integral = 4.39e-14
```

The plan asked for `RWiener::pwiener()` to 1e-10 on the density's
grid: met at 4.46e-11, and the reference shows that difference is
RWiener's. The sweep is why there are two centers: no single center
from 0.005 to 0.35 is under 1e-10 on both F and S. The first build had
one center at 0.1, which is 1.8 percent wrong on S = 8.8e-19 (`t = 3,
v = -5, a = 4, w = 0.1`); `test-wiener-cdf.R`'s tail block fails on
that build (`dev/phase3b-log/seenfail-cdf-single-center.txt`).

Against the hand-written likelihood (RWiener as density and
distribution function, at the fitted parameters), within 1e-10
relative: right censoring on 400 rows, all four codes on 300 rows, and
an upper truncation bound. A per-group `ndt_group()` bound reaches
censored rows: 1e-12 against the family's own functions at
`ndt_time()`.

### Recovery and coverage (`dev/phase3b-eam-recovery.R`, arm `cens`)

30 subjects x 400 trials, `mu ~ cond + (1 | s)`, `bs ~ 1 + (1 | s)`,
`ndt ~ 1`, bias 0.5; truth as `?wiener`'s recovery section. Deadline
1.0 s: a slower trial is recorded at 1.0, right-censored, and its
`dec()` set to 0, which is wrong half the time and unread.
`dev/phase3b-log/summary-cens.txt`:

```
arm cens: 100 files, seeds 1 to 100, gaps: 
first-pass errors: 13; refitted without error on rerun: 13
control seed 1: logLik first pass -5321.0389785807, rerun -5321.0389785807, identical: TRUE
control seed 2: logLik first pass -4256.6162520131, rerun -4256.6162520131, identical: TRUE
errors after the rerun: 0
convergence code 0: 100 of 100; diagnose() reported a finding: 0; diagnose() failed: 2
   64  No convergence problems detected
max |gradient| at the optimum: median 0.00082, max 0.0039
fits with a confint(): 96 of 100

quantity                   truth      mean     mcse   covered    rate  Wilson 95
mu intercept              0.4000    0.3971   0.0074   88/96     91.7%  [84.4, 95.7]
mu condition              0.9000    0.9035   0.0030   93/96     96.9%  [91.2, 98.9]
log bs                    0.3365    0.3381   0.0039   89/96     92.7%  [85.7, 96.4]
ndt: the true 0.25 is below the fitted bound on 100 of 100
ndt (link scale)          2.7261    2.7267   0.0251   92/96     95.8%  [89.8, 98.4]
ndt natural scale: mean 0.2500 (truth 0.25)
log sd 1 (mu | s)        -1.0498   -1.0797   0.0151   89/96     92.7%  [85.7, 96.4]
log sd 2 (log bs | s)    -1.6094   -1.6413   0.0123   93/96     96.9%  [91.2, 98.9]
censored share: mean 0.1636, range 0.1166 to 0.2203
seconds per replicate: median 303
```

The 36 fits with neither a verdict nor a finding from `diagnose()` had
a maximum gradient between its two thresholds, where it prints the
three standard lines and nothing more. The four fits without a
`confint()` are seeds 53 and 63, where `diagnose()` also hit
`bad_alloc`, and seeds 12 and 39, where `confint()` failed and the
harness kept NULL rather than the message, so the cause is not
recorded.

This arm ran on the build with ONE blend center. The build check,
`dev/phase3b-eam-recovery4.R` on the final build, refitted seeds 3, 5
and 7: the log-likelihoods move by at most 5.0e-11 and the
coefficients by at most 1.6e-09 (`dev/phase3b-log/buildcheck.txt`).
The `contfix` seeds 2 and 3, refitted the same way, are `identical()`.
The `cont`, `contdef` and `collapse` arms ran on the build before
`contaminant_range` existed, whose contaminant code path with the
argument unset is the final one.

## 3.5 `wiener(contaminant = TRUE)`

### Construction (`frmtmb.eam/R/wiener-contaminant.R`)

Row density `(1 - lambda) f(t, b) + lambda g`, `g = 1 / (2 (hi - lo))`
on `[lo, hi]`, 0 outside. The 2 is the coin-flip boundary: the response
is the pair (time, boundary), and without it the mixture integrates to
`1 + lambda`. `lambda` is a dpar on a logit link, read through
`dpar_log_complement()`. On the log scale,
`lf + M + log((1 - lambda) e^{-M} + lambda e^{d - M})` with
`d = lg - lf`, `M = (d + |d|) / 2`. `lcdf` is
`(1 - lambda) F + lambda G`; `lccdf` is the same anchored log-sum over
`(1 - lambda) S + lambda (1 - G)`. `fitted()` and `simulate()`
condition on the row's boundary: the contaminant's share of a row at
boundary b is `0.5 lambda / (0.5 lambda + (1 - lambda) P(b))`.

`[lo, hi]` is `contaminant_range =`, or `range(y)` of the fitted data,
kept on the family like the ndt bound. With `contaminant = TRUE`,
`max_ndt` may exceed the fastest response, and the Wiener part then
uses the floored density `allow_unreachable` selects; at the default
bound the plain density is kept. `contaminant` with
`allow_unreachable` is refused, and so is `contaminant_range` without
`contaminant`.

### Validation (`test-contaminant.R`)

Integrates to one over time and both boundaries, split at the range's
edges, within 1e-9. At `lambda = 0` the density is `identical()` to the
plain family's on every row whose diffusion density is at least the
contaminant's, and within 64 ulp elsewhere. A censored row under the
contaminant matches the hand-written mixture within 1e-12. A row
outside a stated range scores `log(1 - lambda) + lf`. The conditional
mean agrees with 5000 draws within four Monte Carlo standard errors.

### Recovery: the default range is wrong

Arm `cont`: 5 percent of rows replaced by contaminants uniform on
[0.1, 5] s with a coin-flip boundary, `max_ndt = 0.5`, default range.
Arm `contfix`: the same draws (same seeds), `contaminant_range = c(0.1,
5)`. `dev/phase3b-log/summary-cont.txt` and `summary-contfix.txt`:

```
arm cont: 100 files, seeds 1 to 100, gaps: 
first-pass errors: 11; refitted without error on rerun: 11
control seed 1: logLik first pass -9129.1343256111, rerun -9129.1343256111, identical: TRUE
errors after the rerun: 0
convergence code 0: 29 of 100; diagnose() reported a finding: 2; diagnose() failed: 1
   17  No convergence problems detected
    3  Non-finite objective at 1 trial point; the optim
max |gradient| at the optimum: median 0.0092, max 0.32
fits with a confint(): 99 of 100

quantity                   truth      mean     mcse   covered    rate  Wilson 95
mu intercept              0.4000    0.3938   0.0071   90/99     90.9%  [83.6, 95.1]
mu condition              0.9000    0.8942   0.0030   91/99     91.9%  [84.9, 95.8]
log bs                    0.3365    0.3409   0.0039   90/99     90.9%  [83.6, 95.1]
ndt: the true 0.25 is below the fitted bound on 100 of 100
ndt (link scale)          0.0000   -0.0047   0.0010   88/99     88.9%  [81.2, 93.7]
ndt natural scale: mean 0.2494 (truth 0.25)
log sd 1 (mu | s)        -1.0498   -1.0925   0.0151   90/99     90.9%  [83.6, 95.1]
log sd 2 (log bs | s)    -1.6094   -1.6439   0.0121   96/99     97.0%  [91.5, 99.0]
lambda (logit)           -2.9444   -3.0607   0.0135   57/99     57.6%  [47.7, 66.8]
lambda natural: mean 0.0451, sd 0.0055 (truth 0.05)
realized contaminant share: mean 0.0501
seconds per replicate: median 424

arm contfix: 60 files, seeds 1 to 60, gaps: 
first-pass errors: 1; refitted without error on rerun: 1
control seed 1: logLik first pass -9131.9991556241, rerun -9131.9991556241, identical: TRUE
errors after the rerun: 0
convergence code 0: 14 of 60; diagnose() reported a finding: 1; diagnose() failed: 0
    9  No convergence problems detected
    1  Non-finite objective at 1 trial point; the optim
max |gradient| at the optimum: median 0.0092, max 0.16
fits with a confint(): 60 of 60

quantity                   truth      mean     mcse   covered    rate  Wilson 95
mu intercept              0.4000    0.4058   0.0098   54/60     90.0%  [79.9, 95.3]
mu condition              0.9000    0.9010   0.0039   55/60     91.7%  [81.9, 96.4]
log bs                    0.3365    0.3369   0.0044   56/60     93.3%  [84.1, 97.4]
ndt: the true 0.25 is below the fitted bound on 60 of 60
ndt (link scale)          0.0000   -0.0018   0.0011   57/60     95.0%  [86.3, 98.3]
ndt natural scale: mean 0.2498 (truth 0.25)
log sd 1 (mu | s)        -1.0498   -1.0813   0.0205   54/60     90.0%  [79.9, 95.3]
log sd 2 (log bs | s)    -1.6094   -1.6476   0.0151   60/60    100.0%  [94.0, 100.0]
lambda (logit)           -2.9444   -2.9481   0.0073   58/60     96.7%  [88.6, 99.1]
lambda natural: mean 0.0499, sd 0.0027 (truth 0.05)
realized contaminant share: mean 0.0502
seconds per replicate: median 261
```

The mechanism, `dev/phase3b-cont-bystatus.R`
(`dev/phase3b-log/cont-bystatus.txt`), on `cont`:

```
logit lambda: sd over 99 fits 0.1348, mean Wald se 0.0584, ratio 2.31
estimate minus realized share (logit): mean -0.1168, sd 0.1321
observed range top:  4.981 4.999 5.413 5.928 6.350 9.928 
cor(log width of range, lambda - realized) = -0.949
```

The top of the observed range is a diffusion trial on most replicates,
up to 9.93 s against contaminants that stop at 5. On the same draws
with the range given, `lambda` recovers and covers. Recorded in
`?wiener` and NEWS with the advice to give the task's response window.
In `cont`, `ndt` covers at 88.9 percent, Wilson [81.2, 93.7]; with the
range given, 95.0. The plan's acceptance criterion names the observed
range; I built it as the default and measured that it fails the
criterion's own recovery test, and did not change the default without
the user.

**The bound.** Arm `contdef`, the same draws at the default `ndt`
bound, 20 seeds: the fastest contaminant pulls the bound under the
true 0.25 on all 20, `ndt` comes back at 0.111 s, the drift effect
covers on 0 of 20, and `diagnose()` names `ndt` at the end of its link
on all 20 (`summary-contdef.txt`). Hence `max_ndt` above the fastest
response is allowed under the contaminant.

### The collapse case

On clean data. Small design, `dev/phase3b-collapse-scan.R`, one subject
of 600 trials, seeds 1 to 40: `lambda` stopped inside its link on 36,
at 0.0011 to 0.039, with log-likelihood gains over the plain family of
0.007 to 4.8; the default range always covers the slowest trial.
`diagnose()` named the 4 that collapsed and none of the 36.

At 30 x 400, arm `collapse`, seeds 1 to 60,
`dev/phase3b-log/summary-collapse.txt` and `dev/phase3b-collapse-detail.R`:

```
convergence code 0: 60 of 60; diagnose() reported a finding: 49; diagnose() failed: 0
   49  Distributional parameter at the end of its link:
   11  No convergence problems detected
lambda at the link edge (logit below -15): 49 of 60
diagnose() names lambda at the end of its link: 49; of those at the edge: 49
interior lambda: median 0.00019, min 5.7e-05, max 0.00044
logLik(contaminant) - logLik(plain): min -9.97e-07, median -6.04e-07, max 1.55
  at the edge: max |difference| 9.97e-07
```

At the edge the logit is -19 to -23.5 with a standard error of 1001 to
1639; the other 11 stop at logits of -7.7 to -9.8 with standard errors
of 0.86 to 5.2 and gains up to 1.55, which is the small design's
behavior again. `test-contaminant.R` asserts the collapse on seed 1
of the small design: log-likelihood within 1e-6 relative of the plain
fit, logit below -8, `diagnose()` naming `lambda`, and the drift within
1e-3 of the plain fit's.

### False convergence, and what it is not (WITHDRAWN: it was defect B1; see punch round 1)

In `cont` 71 of 100 fits and in `contfix` 46 of 60 return nlminb code 1
(false convergence), where `cens`, `collapse` and `contdef` return
code 0 on every fit. It does not cost coverage: in `cont`, `lambda`
covers on 57.1 percent of code-0 fits and 57.8 percent of code-1 fits.
On one subject of 4000 trials (`dev/phase3b-cont-simple.R`, 30 seeds
each; `cont-simple-*.txt`), code 0 on 30 of 30 with contaminants on
[0.3, 4] and `max_ndt = 0.5`, on 28 of 30 on [0.1, 5] with the bound
below the truth, and on 11 of 30 on [0.1, 5] with the truth inside the
bound; `lambda` covered on 30, 27 and 30 of 30. So it needs
contaminants at and below an interior non-decision time. Two causes
were measured and refuted (`dev/phase3b-gradnoise.R`): the tape
gradient for a row just above `ndt` matches a central difference, and
a row below `ndt` gives the same value to the bit as the separation
moves by 1e-9. The cause is not identified.

## Nothing else moved

`dev/phase3b-bitwise.R`: 16 models fitted on the base build
(`rellib-r3`) and on the final build, compared with `identical()` on
the log-likelihood, the optimizer's parameter vector, the fixed
effects, `fitted()` and a seeded `simulate()`. Plain `wiener()` with a
covariate, a `bias` formula, `variability = "sv"`, `ndt_group()`,
`max_ndt`, `vint()`, a lognormal `mixture()`, `wiener_lpdf()`, and the
learn families `bandit2arm_delta` (with and without `(1 | id)`),
`bandit2arm_dual`, `prl_fictitious`, `igt_pvl_delta`, `rlddm`, a value
trace and a `ts_par7` task draw. `dev/phase3b-log/bitwise-compare.txt`:
`IDENTICAL 59 of 59`.

frmtmb.learn imports seven functions from frmtmb.eam; none changed.

## Checks

Every run below is one test file per R process
(`dev/phase3b-suite.sh`), on the base build (`rellib-r3`) and on this
lane's final build (`phase3b-lib`, installed 04:24 on 2026-09-24),
with the log in `dev/phase3b-log/suite-<pkg>-<plain|gated>-<base|new>.txt`.
`gated` sets every variable `dev/release/run-gated.ps1` sets
(`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS`, `FRMTMB_FUZZ`,
`FRMTMB_STAN_CACHE`), plus `R_MAKEVARS_USER`; each file's log prints
the values it saw and `StanHeaders 2.32.10`. The scale tier
(`FRMTMB_SCALE_TESTS`) is not in that script and was not run.

| run | files | expectations | pass | fail | error | skip | log finished |
|---|---|---|---|---|---|---|---|
| eam plain, new | 28 | 1697 | 1693 | 0 | 0 | 3 | 05:00 |
| eam gated, new | 28 | 1697 | 1693 | 0 | 0 | 3 | 05:00 |
| eam plain, base | 28 | 1697 | 1653 | 40 | 0 | 3 | 05:00 |
| eam gated, base | 28 | 1697 | 1653 | 40 | 0 | 3 | 05:00 |
| learn plain, new | 15 | 441 | 428 | 0 | 0 | 13 | 05:03 |
| learn gated, new | 15 | 505 | 499 | 0 | 0 | 2 | 05:04 |
| learn plain, base | 15 | 441 | 409 | 19 | 0 | 13 | 05:03 |
| learn gated, base | 15 | 505 | 480 | 19 | 0 | 2 | 05:04 |

Base against new, file by file, differs ONLY in the files this lane
added or changed: on the base build `test-wiener-cdf.R` fails 17 of 19,
`test-contaminant.R` 20 of 22, `test-surface.R` 2 of 50, `test-family.R`
1 of 51, and `test-session.R` 19 of 20. Those are the tests seen
failing, kept as `dev/phase3b-log/seenfail-*.txt`. What passes on base
in them is preconditions (a fixture's row count, a censored count, the
data selection) and refusals that already existed. Every other file
has identical counts on both builds. Every new value in the new test
files is computed inside `tryCatch()`, so on base each assertion fails
rather than the block stopping at an error: base shows 0 errors.

`R CMD check --as-cran`, with vignettes and the manual
(`dev/phase3b-check.ps1`; logs `dev/phase3b-log/check-*.log`, the
in-check `testthat.Rout` copied beside them):

- frmtmb.learn: `Status: OK`, log 05:08. In-check tests
  `FAIL 0 | WARN 0 | SKIP 15 | PASS 420`.
- frmtmb.eam: `Status: 1 NOTE`, log 05:25, the environmental V8 NOTE
  ("Skipping checking math rendering: package 'V8' unavailable"). In-check tests
  `FAIL 0 | WARN 1 | SKIP 6 | PASS 1682`; the warning is the existing
  `wiener_gng` recovery test's large-gradient warning, in a family and
  test this lane did not touch.

Both check logs postdate every source file (the last, `wiener-family.R`,
04:43). The eight suite runs started before that last edit, which
changed roxygen text only; the in-check test runs above are the ones
that postdate it.

## What I did not do

- **3.4 for `gddm()` and `lba()`**, which the row names. The brief
  scoped the item to `wiener()`, and every new likelihood here is held
  to a reference and a 30 x 400 recovery study; I spent that on
  `wiener()`. `lba()` is the race identity, the three lines `rdm()`
  already has. `gddm()` needs the absorbed mass off its own grid.
- **3.5 for `lba()` and `rdm()`.** Same reason.
- **Censoring under `variability =`**: refused by name.
- **The default contaminant range.** Measured wrong on the plan's own
  design; left as the plan specified, documented, and the fix offered
  as an argument. Changing the default is the user's call.
- **The false-convergence cause** under fast contaminants.
- The v2 reference's cross-check could form `1 - w` at 400 bits; it
  validates the reference to double precision as it stands.

## Needs the user

- Whether `contaminant_range` should default to something other than
  the observed range, for example `c(0, max(y))` or a required
  argument.
- The machine: free space on C: fell to 738 MB during this lane (16 GB
  at the end) and Storage Sense is the recorded mechanism of the
  seventh library loss. Nothing was deleted.
