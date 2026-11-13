# Lane ordstart: a multi-start for ordinal mixtures (backlog 2349)

Base e8abed0a (frmtmb 0.81.0), branch `wt-r82-ordstart`, uncommitted.
Before: `rellib-r20`. Lane library: `C:/Users/adf44/source/r/wt-ordstart-lib`.
Scripts `dev/ordstart-*.R`; logs in `dev/ordstart-out/` (local, gitignored).
Every number below is pasted from the log named beside it.

Status: FIX ROUND DONE (2026-11-11), after review dev/reviews/2026-11-11-ordstart.md.
Sections 4 to 10 describe the first build; where they differ, the "Fix
round" section at the end is the record. Everything is uncommitted.

## 1. The 40-fit table on rellib-r20

`dev/optima-ordmix-starts.R` read rellib-r6 and the optima lane library,
so it was ported as `dev/ordstart-starts.R` (same data, same jitter: N(0,
0.5) on every outer coordinate about the default start, `set.seed(1000 +
s)`; it also records each run's `mixture_ord_degeneracy()` row).
Summary `dev/ordstart-table.R`, log `dev/ordstart-out/base-table.txt`:

```text
cum_cum : seeds 20 ; best of 10 above the default by > 1e-3: 9 ; of those degenerate (warning): 5 ; non-degenerate: 4 ; default degenerate (warning): 6
cum_sratio : seeds 20 ; best of 10 above the default by > 1e-3: 11 ; of those degenerate (warning): 7 ; non-degenerate: 4 ; default degenerate (warning): 6
ALL: seeds 40 ; best above default: 20 ; of those degenerate: 12 ; kinds of the degenerate better points:
          step+big step+collapsed+big    step+runoff+big
                 6                  4                  2
```

The optima record (9/5/4/6 and 11/7/4/6) reproduces exactly on 0.81.0.
Every one of the 12 degenerate better points is a STEP FUNCTION (a
component with `sharpen > -0.1`); none is degenerate by collapsed
thresholds or run-off thresholds alone. The rule I reimplemented agrees
with the fit-end warning on 440 of 440 runs (identity check of my
reimplementation, not a measurement of the rule).

## 2. What "degenerate" is, measured on the 40 fits

The model's rule is `mixture_ord_degenerate_check()` (0.68.0 lane
ordmix): a component is a step function (`sharpen > -0.1`: doubling its
latent distances changes the log-likelihood by less than 0.1), has a
threshold no row is near (`near > 50`), or has two thresholds that are
the same double. `R/se-check.R` has no mixture-specific boundary test;
its "theta" is the covariance vector. Weight at 0 or 1: no run of 1640
had a mean mixing weight below 1e-3, so that kind does not occur here;
`disc` is held at 1 in these models, so a scale at a boundary shows up
as a step function (all latent distances to infinity).

40 starts per seed (`dev/ordstart-starts.R ... 40`, summary
`dev/ordstart-classes.R`, log `dev/ordstart-out/base40-classes.txt`),
kinds over 1640 runs:

```text
            cum_cum cum_sratio
  collapsed      42         15
  proper        574        588
  runoff         24         14
  step          180        203
```

Two findings decide the definition used for ranking:

- The flat kinds are path-dependent. cum_cum seed 5: the default
  (collapsed, warns) and run 21 (`cover` 2.5e-16, two thresholds 1e-15
  apart but not the same double, does not warn) have the same logLik to
  1e-10 (`dev/ordstart-out/q3.R`). cum_sratio seed 16: a threshold of
  the sratio component at -39.8 with the log-likelihood flat to 1e-10
  over a shift of 10 (`dev/ordstart-probe16.R`); some runs drift past
  `near = 50` (runoff, warns) and some do not, at one logLik. Ranking
  "proper first" with the full rule would swap a warned point for an
  unwarned one of the same likelihood.
- The flat kinds buy almost no likelihood: over 40 starts the best
  candidate without a step function and the best candidate the full
  rule passes agree within 1e-3 on 36 of 38 seeds; the two others
  differ by 0.17 (cum_cum 14, collapsed) and 0.044 (cum_sratio 5,
  runoff). A step function buys up to 10.4 (cum_cum 14).

So the ranking uses three tiers: proper (the rule passes), flat boundary
(collapsed or runoff, not a step), step function.

## 3. brms 2.23.0 on the same data

`dev/ordstart-brms.R <family> <seed>`: brms's default priors are
`dirichlet(1)` on theta, flat `b`, and `student_t(3, 0, 2.5)` on each
component's Intercept (thresholds). Log posterior (no Jacobian) and
log-likelihood at frmtmb's candidates, then `rstan::optimizing()` from
10 random inits on the default-prior program (posterior mode) and on
the flat one (maximum likelihood).

cum_cum seed 4 (`dev/ordstart-out/brms-cum_cum-4.txt`):

```text
CAND default          ll -388.0955 brms_ll -388.0955 lpost -405.1954 maxb 6.52
CAND best_proper      ll -388.0955 brms_ll -388.0955 lpost -405.1954 maxb 6.52
CAND best_degenerate  ll -382.0119 brms_ll -382.0119 lpost -458.0703 maxb 2252.23
OPT default 1..10: ll -389.27 (7 of 10), -389.33 (2), -397.11 (1); lpost -401.60 at best; maxb <= 7.3
OPT flat: ll -388.0955 (4 of 10, maxb 6.6), -383.2 / -383.8 / -384.4 (maxb 1728, 705, 112), -393.28 (3)
```

cum_sratio seed 14 (`dev/ordstart-out/brms-cum_sratio-14.txt`):

```text
CAND default          ll -380.2392 brms_ll -380.2392 lpost -464.1137 maxb 2132.57
CAND best_proper      ll -389.4253 brms_ll -389.4253 lpost -408.3046 maxb 14.58
CAND best_degenerate  ll -379.4302 brms_ll -379.4302 lpost -464.7739 maxb 2499.49
OPT default 1..10: ll -388.16 (2), -388.44 (8); lpost -401.63 at best; maxb 3.5
OPT flat: ll -380.24 / -380.37 / -380.24 (maxb 1627, 1417, 6619), -385.49 (4, maxb 17-18), -384.53, -387.38, -389.31
```

The translation is exact (`brms_ll` equals frmtmb's logLik at each
candidate). Under brms's default priors the step function's posterior
density is 52.9 and 55.8 units below the proper maximum's, and Stan's
posterior-mode optimizer reaches a proper point from all 20 random
inits. Stan's L-BFGS on the flat likelihood reaches step functions from
3 of 10 and 4 of 10 inits: the likelihood itself has the same suprema.
So brms's answer is near the proper maximum, which decides the choice:
return the best proper maximum and warn about the higher step function,
not the reverse.

Sampling (4 chains, 1000 + 1000, brms's random inits, seed 1; per chain
the median log-likelihood of 100 thinned draws and the median largest
|b| or |Intercept|):

```text
cum_cum 4 (3291 s): CHAIN 1 ll median -399.24 max|b| median 63816 | CHAIN 2 -396.34, 158955 | CHAIN 3 -394.69, 10542983 | CHAIN 4 -399.43, 389023 ; Rhat max 8.04 ; divergent 17
cum_sratio 14 (2935 s): CHAIN 1 -399.77, 102669 | CHAIN 2 -398.96, 69583 | CHAIN 3 -399.65, 2727 | CHAIN 4 -399.40, 677269 ; Rhat max 4.59 ; divergent 0
```

brms's sampler does not find a maximum at all: with brms's flat prior
on `b`, the slopes of one component run to 1e3 to 1e7 in every chain
(its thresholds stay finite under the student_t prior, so the component
becomes a step function of x at a finite cut), the chains do not mix
(Rhat 4.6 and 8.0), and the log-likelihood of the draws (-394.7 to
-399.8) is below both frmtmb maxima. The posterior is improper along the
step-function direction. So the brms sampler is no reference for which
maximum to report; brms's posterior MODE is, and it is proper. A user
who wants the step function held needs priors on `b`. This is a
consequence of brms's documented flat default prior on `b`, not a brms
defect, so nothing was filed in `dev/upstream-bugs.md`.

## 4. The build

`R/mixture-starts.R` (new): `mixture_starts()` runs after the fit, beside
`sd_up_escape()` in `fit_assembled()`, for a family that declares
`post$start_degenerate` (only `mixture()` of ordinal families, through
`mixture_ord_start_class()` in `R/families.R`). Each start is
`optimize_obj()` from `obj$par + 0.5 * u * z`, `z` standard normal from a
private stream with a fixed seed (as `imp_draws()`), `u` one over the
root mean square of the coefficient's design column (1 elsewhere). The
pick (`mixture_starts_pick()`):

- proper default: replaced only by a proper point better by more than
  `grad_tol`;
- flat-boundary default: replaced by a better point (by `grad_tol`) that
  is not a step function;
- step default: replaced by the best proper point, else the best
  flat-boundary point, whatever its likelihood;
- a step function above the point kept by more than `grad_tol`: warning
  of class `frmtmb_mixture_step_start`.

The winner is Newton-polished as the default was. `frmtmb_control(starts
= NULL)`: 10 on a cold fit (no `start`, no refit template) without random
effects; `0` off; a number runs that many on every fit.

## 5. The 40-fit table before and after

`dev/ordstart-compare.R`, log `dev/ordstart-out/compare-v2.txt` (lane
build without the random-effect rule, which these models do not reach):

```text
fits 40 ; logLik moved by > 1e-6: 13 ; up: 6 ; down: 7
kind before -> after:
before      collapsed proper step
  collapsed         1      0    0
  proper            0     28    0
  runoff            1      0    0
  step              0      7    3
degenerate warning before 12 after 5 ; step-start warning after 11
proper before and down after: 0
elapsed seconds, sum before 3.16 after 26.48
```

Moved up: cum_cum 3 (+0.122), 10 (+2.784, runoff to collapsed), 11
(+0.174), 16 (+2.409); cum_sratio 3 (+0.118), 11 (+0.156). Moved down,
step function to proper: cum_cum 14 (-10.371), 15 (-7.365); cum_sratio
4 (-4.450), 9 (-3.731), 14 (-9.186), 15 (-5.824), 18 (-2.691). Left at a
step function (no start found a point without one): cum_cum 9, 12;
cum_sratio 12. No fit that 0.81.0 passed lost likelihood or gained the
degenerate warning. The first build (step functions disqualified, flat
kinds free) moved cum_sratio 10 from a proper -394.4668 to a runoff
point 0.23 higher that warns (`compare-v1.txt`); the tier rule exists
for that.

## 6. Count and cost

Count (`dev/ordstart-count.R` on the 40-start runs, the pick rule
applied to the default plus the first n runs; log
`dev/ordstart-out/base40-count.txt`):

```text
n  5: at the n = 40 choice 31 of 40; moved from the default 13; step answer  3
n 10: at the n = 40 choice 35 of 40; moved from the default 17; step answer  2
n 20: at the n = 40 choice 37 of 40; moved from the default 19; step answer  2
n 40: at the n = 40 choice 40 of 40; moved from the default 21; step answer  2
```

Cost (`dev/ordstart-cost.R`, one process, arms interleaved, min of 3;
log `dev/ordstart-out/cost.txt`):

```text
COST seed 1 n300    off 0.06s    82 evals | on 0.46s   980 evals | ratio time 7.7 evals 12.0
COST seed 1 n3000   off 0.59s   144 evals | on 4.86s  1475 evals | ratio time 8.2 evals 10.2
COST seed 1 re1000  off 4.50s   248 evals | on 103.41s  5154 evals | ratio time 23.0 evals 20.8
COST seed 2 n300    off 0.06s    65 evals | on 0.44s   796 evals | ratio time 7.3 evals 12.2
COST seed 2 n3000   off 0.50s   124 evals | on 5.19s  1449 evals | ratio time 10.4 evals 11.7
COST seed 2 re1000  off 2.47s   121 evals | on 37.09s  1421 evals | ratio time 15.0 evals 11.7
COST seed 3 n300    off 0.11s   160 evals | on 0.59s  1389 evals | ratio time 5.4 evals 8.7
COST seed 3 n3000   off 0.83s   159 evals | on 5.24s  1599 evals | ratio time 6.3 evals 10.1
COST seed 3 re1000  off 2.79s   142 evals | on 41.26s  1521 evals | ratio time 14.8 evals 10.7
```

(`re1000` is `y ~ x + z + (1 | g)`, 50 groups, measured under the
default before the random-effect rule existed; under the final build
the same 10 starts need `starts = 10`.)

How often the default starts move the answer at larger n (lane default
fit against rellib-r20's, 20 seeds per family, `dev/ordstart-starts.R
<arm> 1:20 <out> <family> 0 <n>`, `dev/ordstart-compare.R`; logs
`dev/ordstart-out/compare-1000.txt`, `compare-3000.txt`):

```text
n = 1000: fits 40 ; logLik moved by > 1e-6: 9 ; up: 9 ; down: 0 ; every fit proper before and after ; elapsed seconds, sum before 14.88 after 97.07
n = 3000: fits 40 ; logLik moved by > 1e-6: 1 ; up: 1 ; down: 0 ; every fit proper before and after ; elapsed seconds, sum before 42.97 after 302.1
```

Moves at n = 1000: +0.033 to +0.809 (cum_cum 2, 5, 10, 13; cum_sratio 2,
7, 12, 13, 17); at n = 3000: cum_sratio 8, +0.715.

Decision: default ON (10 starts) for a cold fit without random effects,
because the starts change the answer on 13 of 40 fits at n = 300 and 9
of 40 at n = 1000, every change to a higher proper maximum or off a step
function, for 5 to 10 times the optimization. With random effects the
cost is 15 to 23 times (each start is a Laplace fit), so the default is
off there and `starts = n` opts in. The user may prefer a size cutoff
(at n = 3000 only 1 of 40 moved for 7 times the time); I did not add one.

## 7. Tests

New file `tests/testthat/test-mixture-starts.R`:

- the pick rule on constructed candidates (13 cases, every tier pair);
- cum_sratio seed 4, x perturbed by 0, 2 and 4 ulps (k = 0..2): the fit
  has no step function, or the fit-end warning says no start escaped
  one (both valid endings accepted);
- cum_cum seed 4 (proper default, higher step functions): with the
  starts the fit never loses more than 1e-3 against `starts = 0` and
  never gains the degenerate warning; the step-start warning has its
  class;
- the session's random state is untouched, and two fits after different
  seeds are identical;
- the starts run only for ordinal mixtures, not with a `start`, not with
  `starts = 0`, and do with an explicit count;
- random effects: no starts by default;
- `frmtmb_control(starts = )` refuses -1, 1.5, "a", c(1, 2), NA.

Seen to fail on rellib-r20 (`dev/ordstart-out/t-mixture-starts-base.txt`):

```text
RESULT test-mixture-starts.R pass=6 fail=6 err=4 skip=0 warn=0
DETAIL FAIL [a step-function default gives way to a proper maximum] Expected `any(grepl("None of 10 jittered starts", r$msgs, fixed = TRUE))` to be TRUE. ... k 0 (and k 1, k 2)
```

The behavioral failure is that one: 0.81.0 returns the step function
(-383.6344, component 2) without an escape. The other errors are the
weak form (`starts` and `mixture_starts_pick` do not exist there).

Every rounding (`dev/ordstart-ulps.R`, k = 0..23, logs
`dev/ordstart-out/ulps-lane.txt`, `ulps-base.txt`):

```text
lane: SUMMARY step-default pass 24 of 24; proper-default pass 24 of 24
      STEP: 24 of 24 at -388.0844, no step function, step-start warning
      PROPER: 24 of 24 off -388.0955 on -388.0955, neither degenerate
base: SUMMARY step-default pass 0 of 24; proper-default pass 24 of 24
      STEP: 24 of 24 at -383.6344, component 2 a step function
```

## 8. Defects met in the build and fixed before the record

- The first spelling caught errors of the family's classifier as "no
  degenerate component", which would have called an unread point
  proper. A classifier error now drops that candidate, and on the
  default start's optimum stops the starts.
- Runs from a jittered start leaked nlminb's "NA/NaN function
  evaluation" 10 times on the cs() mixture of `test-linpred-layers.R`
  (warn=10 there). A start where the objective is not finite is now
  skipped, and a start's run is quiet: the point kept is judged by the
  convergence and fit-end checks afterwards.

## 9. Existing tests the change reaches

- `test-ordinal-mixture.R`: the fit-end degenerate-check tests (seed 19
  and 3, the margin designs, the weights test, three components) now
  fit with `frmtmb_control(starts = 0)` (`omx_one_start`): their subject
  is the check at the default start's point, and whether a start leaves
  a step function is a rounding path. "cs() and thres(gr = ) reach each
  component" allows the step-start warning (it fired there: 6 of 10
  starts reached -532.8668 against -534.2310).
- `test-ordinal-thres-names.R`: four mixture fits through `otn_frm()`,
  which allows the step-start warning (fired twice, at gains of 0.30 and
  0.98).

## 10. Suites

Core, one file per process (`dev/release/run-par.sh ordstart-core`, 3 at
a time, NOT_CRAN=true, lane library first; log
`dev/release/ordstart-core.log`; every file's `lib:` line names
wt-ordstart-lib):

```text
ORDSTART-CORE ran 318 of 318 ; pass=23547 fail=0 err=0 skip=239 warn=0
```

(run on the build before the "log-likelihood plus log prior" wording of
the step-start warning; the files that wording reaches were rerun after,
below.)

Extensions, one file per process, final build (log
`dev/release/ordstart-ext.log`; every `lib:` line loads core from
wt-ordstart-lib and the extension from rellib-r20, no extension code
changed):

```text
ORDSTART-EXT ran 193 of 193
frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=33 pass=1848 fail=0 err=0 skip=5 warn=0
frmtmb.latent files=12 pass=687 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=87 pass=4050 fail=0 err=0 skip=18 warn=0
frmtmb.spatial files=8 pass=478 fail=0 err=0 skip=6 warn=0
frmtmb.survival files=15 pass=556 fail=0 err=0 skip=0 warn=0
```

Reached files on the final build, gated tier on (NOT_CRAN=true,
FRMTMB_BRMS_FIT_TESTS=true; log `dev/release/ordstart-t3.log`):

```text
RESULT test-linpred-layers.R pass=83 fail=0 err=0 skip=0 warn=0
RESULT test-mixture-starts.R pass=60 fail=0 err=0 skip=0 warn=0
RESULT test-optima.R pass=107 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture-brms.R pass=24 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture.R pass=99 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-thres-names.R pass=46 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture-draws.R (frmtmb.sample) pass=34 fail=0 err=0 skip=0 warn=0
```

The 40-fit table on the final build (`dev/ordstart-out/compare-final.txt`)
is the table of section 5 exactly (same logLik, kinds and warnings per
seed); only the elapsed sum differs (48.16 s, two processes at once).

OpenBLAS with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack` and
`0.3.32 lapack`, OPENBLAS_NUM_THREADS=4; logs `dev/ordstart-out/ob326/`,
`ob332/`; `dev/ciharden-out/` deleted after):

```text
0.3.26: test-mixture-starts 60/0/0/0/0, test-ordinal-mixture 99/0/0/0/0, test-ordinal-thres-names 46/0/0/0/0, test-optima 107/0/0/0/0, test-linpred-layers 83/0/0/0/0 (pass/fail/err/skip/warn)
0.3.26 dev/ordstart-ulps.R: SUMMARY step-default pass 24 of 24; proper-default pass 24 of 24 (all at -388.0844 and -388.0955)
0.3.32: test-mixture-starts 60/0/0/0/0, test-ordinal-mixture 99/0/0/0/0
0.3.32 dev/ordstart-ulps.R: SUMMARY step-default pass 24 of 24; proper-default pass 24 of 24
```

Equivalence survey: the one pair with an ordinal mixture, p131
(`cs_pred() x mixture`, `mixture(sratio(), sratio())`), remeasured with
`dev/equiv-measure.R p131 dev/ordstart-equiv-out`
(FRMTMB_EQUIV_LIB = lane library; rellib-r20): verdict equivalent, L1,
"every level agrees", as recorded; `l1_const` 0 as recorded; `l1_grad`
23.64 to 25.87 (the fit point moved). No verdict changes. No brms-port
ledger row reaches a fitted ordinal mixture (the ledger's only ordinal
mixture rows are constructor refusals), so no row flips; no compat row
changed.

`R CMD check --as-cran` of core (`dev/ordstart-check.sh`, built with
vignettes in `dev/ordstart-check/`, R_LIBS lane, rellib-r20, user;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; log kept as
`dev/ordstart-out/check.log`, the check tree deleted after):

```text
* checking examples ... [119s] NOTE
Examples with CPU (user + system) or elapsed time > 5s
frm_hazard_reads 3.25   0.13    5.55
* checking tests ... Running 'testthat.R' [45m] OK
* checking HTML version of manual ... [33s] NOTE
Skipping checking math rendering: package 'V8' unavailable
Status: 2 NOTEs
```

Both NOTEs are the known ones: the V8 note, and an examples-timing note
on an example this lane did not touch, under the load of two session
runs (lane-rules.md: such a note measures load). No extension changed,
so no extension was checked.

Sessions (`dev/release/run-sessions.sh ordstart-sess 2`, NOT_CRAN=true,
lane library first; per-session logs
`dev/release/ordstart-sess-sessions/`, summary
`dev/release/ordstart-sess.log`), ORDSTART-SESS ran 27 of 27, untagged
LEAK lines 0:

```text
SESSION frmtmb forward files=305/318 pass=23547 fail=0 err=0 skip=239 warn=0
SESSION frmtmb reverse files=305/318 pass=23712 fail=0 err=0 skip=239 warn=0
SESSION frmtmb preload files=305/318 pass=23815 fail=0 err=0 skip=241 warn=0
SESSION frmtmb.coupling forward|reverse|preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward|reverse|preload files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward|reverse|preload files=12/12 pass=687 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward|reverse|preload files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward|reverse|preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward|reverse|preload files=82/87 pass=4050 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.survival forward|reverse|preload files=15/15 pass=556 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial forward|reverse|preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
```

(305/318 and 82/87: the files whose every test skips without the gated
tier leave no row, as at earlier releases, e.g. 228/238.)

## Files touched

- `R/mixture-starts.R` (new): `mixture_starts()`, `mixture_starts_pick()`,
  `mixture_starts_count()`, `mixture_starts_draws()`,
  `mixture_starts_units()`, `mixture_starts_warn()`.
- `R/fit.R`: `cold_start` flag at the top of `fit_assembled()`; one hunk
  after `sd_up_escape()` that calls `mixture_starts()`; `frmtmb_control(
  starts = NULL)` (argument, check, list entry, `@param`); one paragraph
  in `frm()`'s `@param start`.
- `R/families.R`: `post$start_degenerate` on the ordinal mixture family;
  `mixture_ord_start_class()`; the degenerate warning's last sentences
  (names `frmtmb_control(starts = )`, or says no start escaped a step
  function); the `mixture()` doc paragraph.
- `man/frm.Rd`, `man/frmtmb_control.Rd`, `man/mixture.Rd` (roxygenised).
- `NEWS.md`: a top `# frmtmb (development version)` section.
- `vignettes/brms-migration.Rmd`: one sentence on mixture starts.
- `tests/testthat/test-mixture-starts.R` (new);
  `tests/testthat/test-ordinal-mixture.R`,
  `tests/testthat/test-ordinal-thres-names.R`,
  `tests/testthat/test-ordinal-mixture-brms.R` (allow the new warning,
  or fit the default start alone where the subject is the fit-end check).
- `dev/test-backlog.md`: CLOSED line under item 2349; three items under
  "Filed in the round of 2026-11-11 / From lane ordstart".
- `dev/ordstart-*.R`, `dev/ordstart-check.sh`, this file.

Release: a MINOR bump (a new `frmtmb_control()` argument and a BREAKING
change of default fits). No new dependency; no floor change.

## Open questions and not done

Questions for the user:

1. Default on, or behind the switch? I chose ON without random effects
   (section 6). The cost is 5 to 10 times the optimization; the answer
   moves on 13 of 40 fits at n = 300, 9 of 40 at n = 1000, 1 of 40 at
   n = 3000. A size cutoff (for example off above some n) is an option
   I did not build.
2. The new warning fires on proper fits that were silent at 0.81.0 when
   a start reaches a higher step function (4 of the 40 fits: cum_cum 4,
   18; cum_sratio 3, 6). It says the reported logLik is a local maximum
   below a boundary supremum. If that is too loud, it could be a
   message; I kept a warning because `logLik()`, AIC and LR tests read a
   point that is not the supremum.

Not done:

- No size or time cutoff; no adaptive count (stop early when starts
  agree).
- Random effects default off; not measured how often starts would move
  such fits.
- The weight-at-0 kind is not in the rule; it never occurred (0 of
  1640 runs had a mean weight below 1e-3).
- The degenerate-check rule itself is unchanged (backlog items filed).

## Fix round (review dev/reviews/2026-11-11-ordstart.md: FIX, 3 MAJOR, 4 MINOR)

Logs of this round: `dev/ordstart-fix-out/` (local). Each finding was
reproduced on the first build before it was fixed.

### F1 (MAJOR): the fit depended on RNGkind()

Reproduced (`dev/reviews/ordstart-review-rngkind.R`, log
`dev/ordstart-fix-out/repro-rngkind.txt`): "fits whose logLik differs by
RNG kind: 11 of 40"; the draws under L'Ecuyer-CMRG differ from those
under the default kind.

Fix: `mixture_starts_draws()` draws under `set.seed(20261111, kind =
"Mersenne-Twister", normal.kind = "Inversion", sample.kind =
"Rejection")` and puts back the caller's three kinds and `.Random.seed`;
a `.Random.seed` that did not exist is removed. The default kind is
Mersenne-Twister with Inversion, so the draws, and every fit of the
tables, are those of the first build under the default kind.

After (`fix-rngkind.txt`): "draws, default kind: 0.9977 -0.108 0.325 |
L'Ecuyer-CMRG: 0.9977 -0.108 0.325"; "fits whose logLik differs by RNG
kind: 0 of 40". Tests: "the starts draw under a named generator,
whatever the session's" (three kinds; kinds and seed put back; an absent
seed stays absent) and "the same call gives the same fit under any RNG
kind" (cum_cum 14, the review's -389.68 against -379.31; seen failing in
the review's `rngkind.log` on the first build).

What a fixed seed hides, on the record: the answer is the best of the
points 10 starts reach. Under another set of draws, 11 of the 40 fits of
the first build differed, by up to 10.4.

### F2 (MAJOR): under autoscale the starts ran twice and the warning was lost

Reproduced (`ordstart-review-paths.R autoscale`, `repro-autoscale.txt`):
15 of 15 rescaled fits ran the starts twice (picks 2); no rescaled fit
gave the step-start warning; cum_cum 10 lost its degenerate warning when
rescaled.

Fix, in three parts:
- `autoscale_prefit()` sets `starts = 0`: the pre-fit only seeds.
- A fit warm-started by a template jitters about the start the plain
  fit would use (`outer_from_template(make_start(...))`), not about the
  warm start. About the pre-fit's optimum, no start left its basin.
- The autoscale "choose" path (pre-fit not converged) passes the call's
  cold flag to both of its fits (`fit_assembled(starts_cold = )`).
  `autoscale_choose()` compares the tier of the points they keep
  (recorded as `mixture_starts["tier"]`) before the objective; comparing
  objectives alone handed the step function back.

After (`fix-autoscale.txt`, x * 1e-4, 15 fits):
- Every rescaled fit that runs the starts runs them once in the fit it
  reports: picks 1 on the warm-start path. On the 5 "choose" fits, picks
  is 2, one per fit compared, as for any choose fit, which fits twice.
- Of the 7 plain fits that give the step-start warning, the 7 rescaled
  fits give it too. 1 rescaled fit (cum_cum 18) gives it where the plain
  fit does not, because the optimizer steps in other coordinates.
- cum_cum 10 warns degenerate in both, at 0.81.0's -394.2865 (see F3b).

Test: "autoscale runs the starts once, in the fit that is reported".
On cum_cum 4 rescaled, the starts run once and warn as unscaled. On
cum_cum 10, the degenerate warning comes with and without the rescaling.
The review's `autoscale.log` shows the first build failing these.

### F3 (MAJOR): unmeasured and unconverged points were ranked

Reproduced (`dev/reviews/ordstart-review-cum3.R`, `repro-cum3.txt`):

```text
SEED 1 lane ll -460.5832 code 1 (function evaluation limit reached without convergence (9)) max|par| 1177 | starts=0 ll -455.1767 code 0
SEED 4 lane ll -456.0897 code 0 max|par| 1.054e+04 | starts=0 ll -451.6403 code 1
SEED 6 lane ll -454.2496 code 0 max|par| 2.008e+04 | starts=0 ll -451.4528 code 0
SEED 7 lane ll -448.8986 code 1 (singular convergence (7)) max|par| 2229 | starts=0 ll -453.6446 code 0
```

Fix:
- `mixture_ord_start_class()` counts an NA sharpening (doubled
  likelihood not finite) as a step function.
- Only runs with convergence code 0 are candidates.
- A default with a non-zero code runs no starts and stays as 0.81.0 had
  it.

After (`fix-cum3.txt`, the final build's ranking code):
- Seeds 1, 4, 6 and 7 return exactly the `starts = 0` fit: the same
  logLik, code and max|par|.
- No fit changed from code 0 to code 1.
- No point with coefficients of 740 or more was chosen.

Tests: "only converged runs are ranked, and an unconverged default
stays" (seeds 1, 4 and 7, every valid ending accepted) and "an
unmeasured sharpening counts as a step function".

F3b, found while checking F2: cum_cum 10's default is at a flat boundary
(a threshold run off; it warns). It moved to a point 2.78 higher. Under
the rescaling, that point had a gap of 1e-10 where the plain fit had a
gap of 0, so the degenerate warning went (`dev/ordstart-fix-out/q10.R`).
A flat-boundary default now stays, and runs no starts.

### What the three fixes cost in moves (40-fit table, `compare-fix.txt`)

```text
fits 40 ; logLik moved by > 1e-6: 8 ; up: 5 ; down: 3
before      collapsed proper runoff step
  collapsed         1      0      0    0
  proper            0     28      0    0
  runoff            0      0      1    0
  step              0      3      0    7
degenerate warning before 12 after 9 ; step-start warning after 6
proper before and down after: 0
```

- Moved up: cum_cum 3 (+0.122), 11 (+0.174), 16 (+2.409); cum_sratio 3
  (+0.118), 11 (+0.156).
- Moved from a step function to a proper maximum: cum_cum 14 (-10.371),
  15 (-7.365); cum_sratio 9 (-3.731).
- No longer moved:
  - cum_cum 10 (a flat-boundary default).
  - cum_sratio 4, 14, 15 and 18. Their step-function defaults stopped at
    code 1 (singular convergence), so they stay as 0.81.0 had them, with
    the degenerate warning.

The first build moved 13 fits. 4 of the 5 lost moves are the price of
F3's rule that an unconverged default stays; cum_cum 10 is the price of
F3b.

The review's 145-fit sweep, both arms rerun (`ordstart-review-sweep.R`
sets sweep, wide and steep; summary `dev/ordstart-fix-sweep-sum.R`; log
`sweep-all-sum.txt`):

```text
fits 145 ; passed by 0.81.0 119 ; logLik moved 18 (up 15 , down 3 )
passed fits: logLik down 0 ; gained the degenerate warning 0
degenerate warning: base 26 lane 21 ; lost on 5 ; step-start warning 10 (on passed fits 7 )
convergence code 0 -> not 0: 0 ; not 0 -> 0: 0
moved down to a point with max|par| >= 740: 0
time ratio lane/base: median 7.03 range 0.83 to 15.28
```

- The review counted 116 passed fits from the old warning text. This
  count reads "ended at a degenerate boundary".
- The 3 moves down are acat_cumcll_K7 seeds 4 (-4.32) and 6 (-4.89) and
  theta_formula seed 4 (-2.09). Each left a step function for a point
  with max|par| from 10.9 to 30.4.

Every rounding: `dev/ordstart-ulps.R` now uses cum_sratio seed 9, a
code-0 step default (seed 4's default is code 1 and stays).
- Lane: "SUMMARY step-default pass 24 of 24; proper-default pass 24 of
  24". All 24 end at -378.6102 with the step-start warning
  (`ulps-lane.txt`).
- Base: "step-default pass 0 of 24". All 24 end at -374.8792, with
  component 2 a step function (`ulps-base.txt`).

### F4 (MINOR): brms agrees on direction, not on the point

The NEWS and help text now say what was measured: the posterior density
at the step functions is far below that at the proper maximum, and
brms's posterior mode is near the proper maximum, not at it.

20 starts against 10 on the 40 (`dev/ordstart-starts.R lane 1:20 ... 0
300 20`; summary `starts20-sum.txt`):
- 20 starts are above 10 by more than 1e-3 on 6 of 40: cum_cum 15
  (+0.380), and cum_sratio 2 (+0.209), 6 (+2.155), 10 (+3.066), 13
  (+0.415) and 17 (+0.765). They are below on 0.
- Elapsed: 39.22 s against 21.46 s, 1.8 times.
- cum_cum 15 reaches -384.8006 with 20 starts, still below the brms
  mode's -383.58.

I kept 10. The gain is real on 6 of 40, but the cost nearly doubles,
and the choice of count is the user's (below).

### F5 (MINOR): the warning's text

- Grammar: "components 1 and 3 are step functions of their predictors";
  a single component "is a step function of its predictors".
- "the highest -455.18 against ..." replaces "at most".
- The remedy names the prior classes as the degenerate warning does
  ("on class Intercept and b with dpar = "mu<k>""), and says that
  brms's default prior on b is flat.
- Only `mixture_starts_signal()` raises the condition, so changing the
  warning to a message is a one-line change.

### F6 (MINOR): cost

- `frm_allfit()` sets `starts = 0` for every refit. Test: "frm_allfit()
  refits run no starts" (a counter on `mixture_starts_pick()` sees 0
  calls).
- Documented under `starts`: bootstrap and `refit()` replicates
  warm-start from the chosen point and run no starts unless `starts` is
  a number.
- Moves and cost by n, two-class process, 40 fits each, lane default
  against rellib-r20 (`compare-fix.txt`, `compare-1000.txt`,
  `compare-3000.txt`, `compare-20k.txt`; single runs, two processes at
  once):

```text
n       moved (up/down)   degenerate warning before -> after   elapsed sum before -> after   ratio
300     8 (5/3)           12 -> 9                               3.16 s -> 24.31 s             7.7
1000    9 (9/0)            0 -> 0                               14.88 s -> 67.77 s            4.6
3000    1 (1/0)            0 -> 0                               42.97 s -> 235.64 s           5.5
20000   0                  0 -> 0                              228.45 s -> 1836.07 s          8.0
```

  The moves at n = 1000 are +0.033 to +0.809; at 3000, cum_sratio 8
  moved +0.715. In the review's steep set, one fit with n = 20000 moved
  +8.82 (slope 8, K = 6).
- Proposal for the user: a size cutoff would save 8 times the time at
  n = 20000, where 0 of 40 two-class fits moved, but it would lose the
  steep set's +8.82. If you want one: off by default above 10000 rows
  (the starts are then one `starts = 10` away). I did not build it.

### F7 (MINOR, record only): merge note for consolidation

- nlstop's `prior_conflict_escape()` sets `opt$prior_escape` (and the sd
  escape sets `opt$sd_escape`). On the user's call, nlstop's message
  says that "frm() returns that fit".
- When a start wins, `mixture_starts()` keeps only `stationary_escape`
  and `mo_search` from the old optimum. `prior_escape` and `sd_escape`
  are then dropped after their message was given.

At consolidation, do one of these, the same way for both fields:
- (a) Carry `prior_escape` and `sd_escape` into the winning optimum's
  record in `mixture_starts()` (the loop over `c("stationary_escape",
  "mo_search")`), and reword the message to "the fit starts from that
  point".
- (b) Clear them when a start wins, and emit nlstop's message after
  `mixture_starts()`, only when its point is the one returned. (b) keeps
  the message true, but it needs the message moved below the starts
  hunk in `fit_assembled()`.

### Files touched in the fix round

- `R/mixture-starts.R`.
- `R/families.R`: `mixture_ord_start_class()`. The degenerate warning's
  "Compare several starts (frm(start = ))" is back to 0.81.0's words.
- `R/fit.R`: `fit_assembled(starts_cold = )`, the `center` closure and
  the starts doc.
- `R/autoscale.R`: the pre-fit's `starts = 0`, and tiers in
  `autoscale_choose()`.
- `R/allfit.R`: `starts = 0`, and its doc.
- `man/frmtmb_control.Rd`, `man/frm_allfit.Rd`, `NEWS.md`.
- `tests/testthat/test-mixture-starts.R`.
- `dev/ordstart-starts.R` (a starts argument), `dev/ordstart-ulps.R`
  (seed 9), `dev/ordstart-fix-sweep-sum.R` (new).

`R/autoscale.R` and `R/allfit.R` are outside the first round's files.

### Fix-round test runs

- `test-mixture-starts.R` on the lane: pass=84 fail=0 err=0 skip=0
  warn=0 (`t-ms-lane.txt`). On rellib-r20: pass=9 fail=3 err=10
  (`t-ms-base.txt`). The behavioral failure is the step-default test:
  0.81.0 returns cum_sratio 9's step function at code 0 with no start
  run, on k = 0, 1, 2. The rest are the weak form (the new functions and
  argument do not exist there). The RNG-kind and autoscale tests pass on
  0.81.0, which has no starts; they were seen failing on the first lane
  build, in the review's `rngkind.log` and `autoscale.log`.
- Full core suite, one file per process (`dev/release/ordstart-fix-core.log`,
  every `lib:` line wt-ordstart-lib): ORDSTART-FIX-CORE ran 318 of 318;
  pass=23571 fail=0 err=0 skip=239 warn=0.
- Gated (NOT_CRAN and FRMTMB_BRMS_FIT_TESTS true;
  `dev/release/ordstart-fix-t4.log`): test-ordinal-mixture-brms 24/0/0/0/0,
  frmtmb.sample test-ordinal-mixture-draws 34/0/0/0/0, test-linpred-layers
  83/0/0/0/0, test-optima 107/0/0/0/0 (pass/fail/err/skip/warn).
- Core sessions on the changed test files (`dev/release/run-session.R
  frmtmb <order> <dir>`, the dir holding the helpers and
  test-mixture-starts, test-ordinal-mixture, test-ordinal-thres-names and
  test-ordinal-mixture-brms; `sess-*.txt`): forward, reverse and preload
  each files=4/4 pass=229 fail=0 err=0 skip=3 warn=0; no LEAK line. The 3
  skips are the gated brms file.
- OpenBLAS 0.3.26 with LAPACK (`ob326-*.txt`): test-mixture-starts
  pass=82 fail=0 err=0 skip=0 warn=0 (2 fewer expectations than on the
  reference BLAS: a valid ending with fewer assertions was taken),
  test-ordinal-mixture 99/0/0/0/0. `dev/ciharden-out/` deleted.
- Final build rechecks: `final-cum3.txt` (seeds 1, 4, 6, 7 as in F3) and
  `final-rngkind.txt` ("fits whose logLik differs by RNG kind: 0 of 40").
- No extension code changed; the extension suites were not rerun in this
  round (they were clean on the first build, section 10).
- `R CMD check --as-cran` of core on the fix build (`dev/ordstart-check.sh`;
  log kept as `dev/ordstart-fix-out/check.log`, the check tree deleted):
  "Status: 1 NOTE", the known V8 math-rendering note; tests (27 min) OK;
  examples OK (59 s and 67 s with --run-donttest).

### Questions for the user after the fix round

1. Default on, a size cutoff, or behind the switch: the moves-and-cost
   table above (F6). My recommendation is on, without a cutoff, because
   of the steep +8.82 at n = 20000; a cutoff above 10000 rows is the
   cheap alternative.
2. 10 or 20 starts (F4): 20 moved 6 of 40 higher, for 1.8 times the time.
3. Warning or message for `frmtmb_mixture_step_start`: the reviewer
   recommends a message. It fired on 7 of the 119 fits of the sweep that
   0.81.0 passed, and on 6 of the 40. Changing it is
   `mixture_starts_signal()` alone.
4. "An unconverged default stays" keeps 4 of the 40 fits (cum_sratio 4,
   14, 15, 18) at 0.81.0's step function with code 1, where the first
   build moved them to proper maxima. Under autoscale those same fits do
   move (the warm-started fit converges). If you prefer, a converged
   proper start could replace an unconverged step-function default; the
   review asked for the stricter rule, so I built that.

### Q4, decided by the coordinator: an unconverged step-function default

Rule now (`mixture_starts()`, `mixture_starts_pick(converged = )`):
- A default the optimizer did not call converged runs the starts only
  if it is a step function (or its sharpening cannot be measured). Then
  the best converged proper start replaces it, when one exists; flat and
  step candidates do not.
- If no converged proper start exists, it stays as 0.81.0 had it.
- An unconverged default that is proper or at a flat boundary stays, as
  in the fix round, and runs no starts.
- Candidates are still code-0 runs only.

Logs: `dev/ordstart-q4-out/` (local).

40-fit table (`compare-40.txt`, `changes-40.txt`):

```text
fits 40 ; logLik moved by > 1e-6: 12 ; up: 5 ; down: 7
  step              0      7      0    3
degenerate warning before 12 after 5 ; step-start warning after 10
proper before and down after: 0
code base->q4 worse (0 to non-0): 0 ; better: 4
```

Changes against the fix round (each was code 1 at a step function, and
is now code 0 at a proper maximum, the first build's values):

```text
cum_sratio  4  -383.6344 code 1  ->  -388.0844 code 0  max|par| 6.36
cum_sratio 14  -380.2392 code 1  ->  -389.4253 code 0  max|par| 14.83
cum_sratio 15  -377.8161 code 1  ->  -383.6402 code 0  max|par| 5.19
cum_sratio 18  -378.8494 code 1  ->  -381.5406 code 0  max|par| 22.79
```

The 145-fit sweep, both arms rerun (`sweep-all-sum.txt`,
`changes-sweep.txt`; the base arm is identical to the fix round's, in
logLik and code on every case):

```text
fits 145 ; passed by 0.81.0 119 ; logLik moved 23 (up 15 , down 8 )
passed fits: logLik down 0 ; gained the degenerate warning 0
degenerate warning: base 26 lane 16 ; lost on 10 ; step-start warning 15 (on passed fits 7 )
convergence code 0 -> not 0: 0 ; not 0 -> 0: 5
moved down to a point with max|par| >= 740: 0
time ratio lane/base: median 7.45 range 1.04 to 18.18
```

The changes against the fix round are 5 cases. Each was code 1 at a step
function, and is now code 0 at a proper maximum with no degenerate
warning:

```text
cs_sr_acat_seed1    -399.3975 code 1 -> -404.9916 code 0  max|par| 17.95
cs_sr_acat_seed3    -401.1663 code 1 -> -401.2247 code 0  max|par| 38.08
cs_sr_acat_seed4    -391.3687 code 1 -> -397.2175 code 0  max|par| 325.5
dpar_formula_seed4  -398.3574 code 1 -> -401.3254 code 0  max|par| 32.45
dpar_formula_seed7  -394.0352 code 1 -> -400.1462 code 0  max|par| 2.15
```

Confirmed:
- 0 passed fits lose likelihood, in the 40 and in the 145.
- 0 convergence codes change for the worse; 4 (of 40) and 5 (of 145)
  change from 1 to 0.
- No point with coefficients of 740 or more is chosen.
- cs_sr_acat_seed4's new point has a cs() coefficient of 325.5. The
  sharpening rule calls it proper, and the fit-end check passes it.

Autoscale (`autoscale.txt`):
- All 15 rescaled fits now have the plain fit's logLik, including
  cum_sratio 4, 14, 15 and 18, which disagreed in the fix round.
- All 10 plain fits that give the step-start warning also give it
  rescaled; cum_cum 10 warns degenerate in both, at -394.2865.

cum_cum 18's rescaled-only warning remains (`q18.R`, `q18.txt`). Plain,
`mixture_starts_pick()` sees 10 candidates, all at 381.39661 and proper:
one of the 10 starts did not end at a converged run, so it is not
ranked. Rescaled, it sees 11 candidates, and that start converges, at a
step function with objective 378.84939 (logLik -378.8494, 2.55 higher).
Both fits keep the same proper maximum (keep 1). Only the rescaled fit
saw the higher step function, so only it says so. The cause is the
optimizer's coordinates under autoscale, not the rule.

Tests (`t-*.txt`):
- `test-mixture-starts.R`: pass=95 fail=0 err=0 skip=0 warn=0. New:
  pick cases with `converged = FALSE`, and "an unconverged step-function
  default gives way to a proper one" (cum_sratio 4, k = 0 and 1, every
  valid ending accepted). The three-component test now accepts an
  unconverged default replaced by a converged proper point, and still
  fails a fit that goes from code 0 to non-0.
- `test-ordinal-mixture.R`: pass=99 fail=0 err=0 skip=0 warn=0.

The NEWS bullet and the `starts` help say: "when it did not converge,
only by a proper maximum".

`R CMD check --as-cran` of core on the Q4 code (`dev/ordstart-check.sh`;
log `dev/ordstart-q4-out/check.log`, the check tree deleted): "Status:
1 NOTE", the V8 math-rendering note; tests (25 min) OK; examples OK
(56 s and 70 s with --run-donttest). After the check, only the numbers
in the `starts` help text and NEWS changed (12 of 40, 12 to 5, 18 and 23
of 145). They were roxygenised and installed, and no code changed.

Unchanged and left with the user: 10 starts, the default on with no
size cutoff, and the warning class `frmtmb_mixture_step_start`.

### After the short check (review "Short check", SHIP with minors)

Logs: `dev/ordstart-short-out/` (local).

**1. A stricter "proper" test, measured, not built.**
`dev/ordstart-strict.R` refits the 40 with `starts = 0` and with the
default. At each point it measures the smallest eigenvalue of the exact
outer Hessian (`obj$he()`, as in the review's `ordstart-short-moves.R`)
and the smallest category cover. `dev/ordstart-strict-sum.R` scores two
floors against the six flat points of the review: cs_sr_acat seeds 4
and 1, dpar_formula seed 4, cum_sratio 14 and 18, and acat_cumcll_K7
seed 6. Their numbers come from the review's `q4pts.log`, `profile.log`
and `moves.log`. Log `strict40-sum.txt`:

```text
six flat points: rejected by the eigenvalue floor 5.1e-08 : 6 of 6
six flat points with a cover measured: rejected by cover < 1e-08 : 1 of 4
passed defaults the eigenvalue floor would call flat (they would stay): 11
of those, fits whose upward move would be lost: 2   (cum_cum 16, cum_sratio 3)
passed defaults the cover floor would call flat: 8
of those, fits whose upward move would be lost: 2   (cum_cum 16, cum_sratio 3)
step-to-proper moves whose new point the eigenvalue floor rejects: 4   (cum_cum 15, cum_sratio 9, 14, 18)
upward moves to a point the eigenvalue floor would not call proper: 5   (cum_cum 3, 11, 16; cum_sratio 3, 11)
```

- The eigenvalue floor rejects all six. But it changes fits 0.81.0
  passed: 11 of the 28 passed defaults would count as flat, and the 5
  upward moves land on points it rejects. It also rejects cum_cum 15
  and cum_sratio 9, whose largest standard errors are 1.33 and 7.55.
  On these data a raw-scale eigenvalue does not separate a flat
  direction from a nearly collapsed category: the six have 9.8e-9 to
  5.07e-8, and sound fits have 3.0e-8 to 5.4e-8.
- The cover floor misses cum_sratio 14 (cover 6.6e-7). It also calls 8
  passed defaults flat.
- So the build condition fails on the 40-fit table alone: the floor
  must reject all six and change no fit that 0.81.0 passed or that I
  classed as clean proper. I did not also run the sweep.
- Filed with these numbers, as an extension of the backlog item "The
  degenerate check's 'collapsed' test is exact-double equality". The
  next candidate filed there is a standard-error test in natural units:
  above 5e3 on the flat points, at most 43.6 on the other moves. It is
  not measured.

**What the help, NEWS and CLOSED line now say about "12 to 5".**
- The 7 step-to-proper points on the 40 are now described as "a point
  the degenerate check passes", not as proper maxima.
- Of the 7, 2 are clean maxima (cum_sratio 4 and 15) and 3 have a
  category near 1e-8 to 1e-7 (cum_cum 14 and 15, cum_sratio 9).
- 2 lie on a flat direction with a standard error above 5000 (cum_sratio
  14 and 18).
- Of the 8 such points of the 145, 4 are flat in the same way, and one
  of the 4 is not a maximum (cs_sr_acat seed 4: a refit gains 0.46).
- The help adds: "A point the degenerate check passes is not always a
  clean maximum ... Check the standard errors of a fit that moved."
- The backlog CLOSED line adds that 2 of the 7 new points are flat
  points the check does not see.

**2. Warning text.** It now reads "so the fit reports the best converged
point (the default start's) without one. logLik() and the information
criteria are those of this point." The class is unchanged.

**3. Backlog.** The filed item now lists the six cases with the review's
numbers and the measurement in item 1.

**4. NEWS.** The bullet is rewritten to 80 columns (no line of the
section is longer).

**5. Runs on the final build (NOT_CRAN and FRMTMB_BRMS_FIT_TESTS true;
`dev/release/ordstart-short-t5.log`):** see below. The `R CMD check` of
core was run again, because the docs changed after the last tarball.

```text
RESULT test-mixture-starts.R pass=95 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture.R pass=99 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-thres-names.R pass=46 fail=0 err=0 skip=0 warn=0
RESULT test-linpred-layers.R pass=83 fail=0 err=0 skip=0 warn=0
RESULT test-optima.R pass=107 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture-brms.R pass=24 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture-draws.R (frmtmb.sample) pass=34 fail=0 err=0 skip=0 warn=0
ORDSTART-SHORT-T5 ran 7 of 7
```

Every `lib:` line loads core from wt-ordstart-lib (frmtmb.sample from
rellib-r20).

`R CMD check --as-cran` of core on the final build (log
`dev/ordstart-short-out/check.log`, the check tree deleted): "Status: 1
NOTE", the V8 math-rendering note; tests (25 min) OK.

Unchanged and left with the user: 10 starts, the default on, and
warning or message for `frmtmb_mixture_step_start`.
