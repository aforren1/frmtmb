# Lane nlstart: nonlinear fits that stop short with no warning

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-nlstart`, library `wt-nlstart-lib`, base library `rellib-r18`.
Status: DONE, with a fix round after review (last section). Numbers are pasted
from the logs named beside them (`dev/nlstart-out/`, gitignored).

## The three defects, reproduced on rellib-r18

`dev/nlstart-repro.R`, log `dev/nlstart-out/repro-before.txt`.

1. `b1 * (1 - exp(-b2 * x))`, gaussian, the 450-row scan's data (seed
   `s * 100 + 1`, slope 0.6): default -601.7851644 at code 0, no warning,
   b1 = 4.8845e+04, b2 = 1.9427e-05 (the linear limit); from a start of
   (-2, -0.5): -530.4714139 at b1 = -0.73239, b2 = -1.11954. Seeds 201 to
   501: -632.06 (code 1, warned), -628.95, -617.71, -632.50 against
   -557.23, -555.31, -547.05, -563.44. With slope -0.6 the default reaches
   the maximum on all 5.
2. `a * x` with x = 1e-20 * rnorm (seed 6201): -299.3282639 at code 0, no
   warning, a = 9.07e-19; `y ~ 0 + x` -272.503545.
3. `a + b * x`, y = 1e6 + 5e4 x + noise (seed 9901): -3045.621732, code 0,
   a = 0.166, b = -0.180, with the standard-error warning; `y ~ x`
   -1658.76266.

## Causes, measured

1. **Start (sign pattern).** The body is constant at its zero start, so
   nl_fit_off_saddle() restarts from moved coefficients with the sign
   patterns (+,+), (-,+), (+,-) at magnitude about 0.16. All three run to
   the linear limit; the maximum needs (-,-). `dev/nlstart-probe-starts.R`
   (log `dev/nlstart-out/probe/gaussian.txt`): at magnitude 0.15 only
   (-,-) reaches it on 5 of 5 seeds; at magnitude 1, (-,+) and (-,-) do.
   The optimizer's stopping rule is not the cause: the fit stops at a
   point where it has converged along the ridge toward the linear model.
   The other short rows of the scan are the same class: bernoulli
   `b1 * exp(b2 * x)` and `b1 * exp(b2 * x) + b3 * z` stop at a true
   second maximum (the other sign of b2), which no restart tries because
   that fit is not on a saddle at all (no message, no warning).
2. **Scale of the nonlinear parameter.** The gradient in `a` at the start
   is about 1e-18, and nlminb steps in absolute units, so it stops at
   once. The linear spelling is stepped in units of u / sd(x) by
   predictor_par_units(). The body rewritten as `a * (xt * 1e20)` reaches
   -272.503545 (`dev/nlstart-scale-probe.R`, log
   `dev/nlstart-out/scale-probe-before.txt`). Not the start: no start of
   order 1 helps.
3. **Start and scale together.** sigma starts at log(sd(y)) = 10.79 and mu
   at 0, 20 sigma below every row; a step of order 1 in `a` changes
   nothing, so sigma takes the residual. A start of `a` at mean(y) alone
   reaches -1658.76266 (`dev/nlstart-scale-probe.R`). So does the unit
   below, without a change of start.

## What changed

All in core, `R/par-template.R` (new functions after nl_local_message(),
away from lane fitedge's hunks in nl_fit_off_saddle()) and two hunks of
`R/fit.R`. nl_fit_off_saddle() itself is unchanged (it holds `<<-`, which
lane superassign rewrites).

- **nl_par_units()**: each nonlinear coefficient gets a natural unit,
  the response's unit (the sd of `y` for an identity-link location-scale
  `mu`, as predictor_unit_input() decides; 1 otherwise) over the root mean
  square of the body's derivative in the coefficient at the start, to a
  power of 2. Multiplied into `par_units` in fit_assembled() before
  predictor_par_units() (a separate hunk from lane bigx's, which inserts
  after it). Fixes defects 2 and 3. Two limits, each from a measured
  failure of a broader first version:
  - A unit within a factor of 16 of 1 is not applied. Without it, a user
    start of (-2, -0.5) on the slope -0.6 sat data went from -539.108463
    to -606.198 (`dev/nlstart-out/repro-after1.txt`, first build): a
    factor-2 unit changed the path to the other maximum.
  - Only a coefficient the body reads as an added term gets a unit
    (nl_body_additive()). The second build gave every coefficient one, and
    the full suite (`dev/nlstart-out/all1.log`) failed frmtmb.ode's
    test-ode-tv.R 4 times: `lk ~ 1 + phase` inside an ODE rate got a unit
    of 16 or more on `phaselate`, from the start's derivative, and ran to
    9.1 against the truth log(2), where the late phase is 0 and flat
    (logLik -193.8 against the flat model's -179.8). test-ode-nlf.R failed
    once (the nlf() and direct spellings 1.6e-8 apart). Both pass on the
    restricted build (`dev/nlstart-out/re1.log`). Defects 2 and 3 are both
    added terms (`a * x`, `a + b * x`), so the restriction keeps them.
- **nl_sign_probe()**, called by **nl_fit_probed()**, which wraps
  nl_fit_off_saddle() at frm()'s one call site. From the default start
  only (`start = NULL`): for every coefficient the body does not read as
  an added term (nl_body_additive(): the derivative is the same at two
  generic points), the sign patterns (all 2^k for k <= 2; else all
  positive, all negative, the first two flipped: 4 at most), each at one
  natural unit over the nlpar column's half range, the rest of the model
  at the fit's own start (`obj$par`). Each runs optimize_obj() on the
  fit's own tape (no new tape, no standard errors); the objective's state
  is put back (obj_settle()). If the best is higher than the fit by more
  than max(1e-4, 1e-9 |ll|), frm() refits from that point with `run()` and
  returns that fit, with a message of class `frmtmb_nl_multistart`. The
  first fit's warnings and its `frmtmb_nl_local` message are dropped then
  (they describe a fit not returned); other messages are kept. Skipped
  under quadrature, importance, profile, a whitened beta, a map on a
  body coefficient, a custom optimizer, an early-stopping optCtrl, and
  `options(frmtmb.nl_sign_probe = FALSE)`. Fixes defect 1 and the
  bernoulli class.

Not changed: nl_domain_start() (defect 3 is fixed by the unit alone, so
the start rule, which the record says is deliberately narrow, stays).

## The 450 rows

`dev/nlstart-scan.R` (the rows of `dev/review-nlfix-d7.R`: 9 families x 5
bodies x slopes +-0.6 x seeds 101..501, n = 300; reference the best of 8
random N(0, 1) starts, seed `seed * 7 + 1`), driver `dev/nlstart-scan.sh`,
summary `dev/nlstart-scansum.R`. The reference of a row is the best any
arm reached. Short = more than 1e-4 below it; silent = no warning of any
class. Logs `dev/nlstart-out/scan-before/` and, on the final build,
`scan-after3/`, summary `dev/nlstart-out/scansum-after3.txt`:

```
scan-before: rows 450, default NA 0, short > 1e-4: 16 (silent 11, warned 5), short > 1: 15 (silent 10)
scan-after3: rows 450, default NA 0, short > 1e-4: 0 (silent 0, warned 0), short > 1: 0 (silent 0)
scan-after3 minus scan-before default logLik: higher by > 1e-4: 16  lower by > 1e-4: 0  within 1e-4: 434
```

The "16 of 450" of the record reproduces exactly. The 16 before:

```
            family seed slope body   default conv nwarn         wclass nmsg      best   short
          gaussian    2   0.6  sat -632.0585    1     1 frmtmb_warning    1 -557.2304 74.8281
          gaussian    3   0.6  sat -628.9511    0     0              -    1 -555.3073 73.6438
          gaussian    1   0.6  sat -601.7852    0     0              -    1 -530.4714 71.3138
          gaussian    4   0.6  sat -617.7062    0     0              -    1 -547.0531 70.6531
          gaussian    5   0.6  sat -632.5010    0     0              -    1 -563.4429 69.0581
             gamma    5   0.6  sat -489.4618    1     1 frmtmb_warning    1 -471.8400 17.6218
         bernoulli    3  -0.6  exp -200.3857    0     0              -    0 -190.7695  9.6162
         bernoulli    3  -0.6  two -197.0672    0     0              -    0 -187.8296  9.2376
              acat    2   0.6  two -354.9559    1     2 frmtmb_warning    0 -345.8564  9.0995
 hurdle_cumulative    5   0.6  two -436.6934    1     2 frmtmb_warning    0 -427.8905  8.8029
              acat    4   0.6  two -383.3198    1     2 frmtmb_warning    0 -377.5364  5.7834
         bernoulli    2  -0.6  exp -195.8425    0     0              -    0 -192.0260  3.8165
         bernoulli    2  -0.6  two -191.6263    0     0              -    0 -188.9074  2.7189
         bernoulli    5   0.6  exp -195.3480    0     0              -    0 -193.5546  1.7934
         bernoulli    5  -0.6  exp -198.4821    0     0              -    0 -196.8912  1.5909
         bernoulli    5   0.6  two -193.5242    0     0              -    0 -193.0972  0.4270
```

Caveat on the reference: it is 8 random starts per arm. A row can be
short of a maximum no arm found; "0 of 450" is against this reference.

### Which patterns, measured

`dev/nlstart-patterns.R` runs every pattern on each row's tape with the
probe off; `dev/nlstart-patsum.R` replays strategies offline (logs
`dev/nlstart-out/pat1/`, base = the fit's start; `pat-opt/`, base = the
fit's optimum for the other parameters). Time ratio counts the fit plus
the probes, not the refit of a winner.

```
start of the rest at the fit's start (shipped):
all            short   0  time ratio 1.749
neg_first      short   4  time ratio 1.290
pos_neg_first  short   2  time ratio 1.604
other_orthant  short   0  time ratio 1.666
start of the rest at the fit's optimum:
all            short   2  time ratio 1.417  (acat 2 0.6 two; acat 4 0.6 two)
```

Shipped: all patterns from the fit's start. A cascade (stop when the
first one or two patterns agree with the fit) loses 2 to 4 bernoulli
rows. Skipping the fit's own sign orthant lost none here, but on
`dev/nlstart-probe-starts.R` bernoulli seed 201 slope -0.6 `exp` is
reached from (-,-) at magnitude 1, the orthant of its short fit, so it
was not taken. An earlier build varied every coefficient, 8 patterns for
`two` and `prod`: 0 short as well, wall ratio 2.35 (scan-after1, not
interleaved); nl_body_additive() halves those to 4.

### Probe-off equals 0.79.1

`dev/nlstart-offcmp.R`, final build, cost run `cost4`'s probe-off arm
against scan-before (log `dev/nlstart-out/offcmp4.txt`):

```
rows 450  probe-off vs 0.79.1 default: |diff| > 1e-8: 0  max |diff| 0 
```

On the second build (units on every coefficient) 7 rows differed by up
to 6e-6, all `logi` under `cumulative()` and `hurdle_cumulative()`
(`dev/nlstart-offchk.R`): from nl_fit_off_saddle()'s small restart starts
the derivative in b2 is about 0.04, so b2 got a unit of 32, and both
builds stopped on the ridge (|b1| 600 to 1200) a few 1e-6 apart. With
units on added terms only, none differs.

## Cost

Interleaved in one process, probe off / on / off again in rotating
order, per row; off2 / off is the control (`dev/nlstart-cost.R`,
`dev/nlstart-costsum.R`). Two runs, at two loads:

```
cost3 (second build; the probe is the same as the final one):
seconds: off 73.79  on 140.07  off2 74.54
ratio on / off: 1.898   control off2 / off: 1.01
by body, on / mean(off, off2):  exp 2.428  logi 2.346  prod 1.160  sat 1.822  two 2.031
cost4 (final build, machine about 2.5 times as loaded):
seconds: off 182.48  on 280.53  off2 177.26
ratio on / off: 1.537   control off2 / off: 0.971
by body, on / mean(off, off2):  exp 2.450  logi 1.694  prod 1.014  sat 1.522  two 1.684
```

So 1.5 to 1.9 times on these rows, depending on load; the `exp` body,
whose default fit is cheapest, 2.4 to 2.5. (The scans' own totals,
scan-before 49.5 s against scan-after3 240.2 s, are separate runs at
different loads; their ratio is not a measurement.)

Larger fits, min of 3 interleaved rounds (`dev/nlstart-cost-big.R`, seed
7301, n = 20000, log `dev/nlstart-out/cost-big1.txt`; this was the build
with up to 8 patterns, but every body there has k = 2, so 4 patterns
either way):

```
gauss_exp_2e4 min s: off 0.16 on 0.49 off2 0.19 | on/off 3.062 control 1.187
gauss_sat_2e4 min s: off 1.10 on 2.91 off2 1.11 | on/off 2.645 control 1.009
pois_re_2e4   min s: off 0.58 on 1.88 off2 0.64 | on/off 3.241 control 1.103
```

Bodies with 3 and 4 varied coefficients, blocks of 30 fits, min of 3
rounds (`dev/nlstart-repro2.R`, seed 7304, n = 300, log
`dev/nlstart-out/repro2-after.txt`): prod3 on/off 1.629 (control 0.933),
two_exp 2.103 (control 1.140).

The reason the cost is paid on every curved body: nothing in a short fit
marks it. The bernoulli rows stop at a true second maximum with code 0, a
clean gradient and finite standard errors, and they are not on a saddle,
so no existing restart runs. The only test is to look. Skip: `start`, or
`options(frmtmb.nl_sign_probe = FALSE)`.

## Tests changed because the old fixture's trap is gone

The full suite on the second build (`dev/nlstart-out/all1.log`, 500 files,
every `lib:` line reads the lane's core) failed 5 files. Besides the two
ODE files above:

- test-conditions-census.R: nl_fit_probed() re-signals held conditions
  with `warning(cnd)`, `message(cnd)` and `stop(e)`, as nl_fit_off_saddle()
  does. Three exemptions added in `dev/conditions-census-gen.R`, the test
  regenerated (the extension copies come out byte-equal; `git diff` shows
  none).
- test-diagnostics-ux.R (9 fails) and test-se-flags.R (5 fails): their
  fixture is a default fit of a Gaussian bump with `pk` at 0 while w runs
  1 to 45, so the bump is 0 on every row and three directions are flat.
  The probe's starts of pk at +-1 put the bump on the data and the fit
  leaves that point. The tests are about diagnose() and the SE flags on a
  flat fit, so they hold the probe off with
  `withr::local_options(frmtmb.nl_sign_probe = FALSE)`. One of them, "a
  theta in the flat set ...", had skipped on the second build ("this
  fixture did not drive the group sd flat"); with the probe off it runs
  again (skip 1 to 0).
- test-nl-family.R: "a point no restart leaves is a warning" mocks the
  saddle restarts away to keep the fit on the point; it now also holds
  the probe off. "a constant body at a local maximum is said in a
  message" reads that message with the probe off, and pins the new
  default: higher than the thresholds-only model by more than 1e-3 of
  |logLik|, no `frmtmb_nl_local` message, a `frmtmb_nl_multistart` one.

An early run of the nl files (`dev/nlstart-out/nl1.log`) read rellib-r18
by mistake (the runner needs `FRMTMB_RELLIB`); it counts only as a base
control. Every run cited below checks its `lib:` lines.

## Other measurements

- The flat-peak fixture of those tests (`dev/nlstart-flatpeak.R`, seed 5,
  log `dev/nlstart-out/flatpeak.txt`): probe off 3455.731789, code 0, the
  SE warning, pk -0.90; probe on 3510.85196, code 0, no warning, pk 10.18
  (the truth draws centres from N(10, 1.1)). 55.1 higher. The old default
  fit was one of the stops this lane is about.
- Equivalence survey, the 14 pairs with an nl feature (`dev/nlstart-equiv.sh`,
  `dev/nlstart-equiv-cmp.R`, log `dev/nlstart-out/equiv-cmp.txt`): no
  verdict changes. 12 pairs keep every verdict column; p014 (`cs_pred() x
  nl`) and p166 (`nbinom1 x nl`) differ from the record in `level` (NA
  against an empty string) and p166 in `agree` (FALSE against NA). The same
  two read the same on rellib-r18 (`dev/nlstart-out/equiv-base/`), so it is
  the record's adjudicated text against a raw remeasure, not this change.
- The nl x family sweep (`dev/equiv-nl-sweep.R` into
  `dev/nlstart-out/nlsweep/`, compared by `dev/nlstart-nlsweep-cmp.R`):
  48 rows, every logLik column within 1e-6 of `dev/equiv-data/nl-sweep.tsv`,
  no outcome changed. Its body `b0 + b1 * x` is all added terms, so the
  probe does not run there, as designed.
- Gated (FRMTMB_BRMS_FIT_TESTS=true, `dev/nlstart-out/gated1.log`):
  test-nl-brms.R 91 pass, test-nl-stan-funs.R 162 pass, frmtmb.sample
  test-nl-body-draws.R 15, test-prior-nlpar-intercept.R 50,
  test-brms-suite-brmsformula.R 16; 0 fail, 0 skip.

- The help's bernoulli example (`dev/nlfix-bern.R`'s data, seeds 9202 and
  9205; `dev/nlstart-repro2.R`): rellib-r18 -193.4966858 at b1 -0.36,
  b2 -0.87; this build -189.9989833 at b1 0.0017, b2 5.01, with the
  message. The help of `par_template()` and `frm(start =)` is rewritten.
- `b^2 * x * (1 - b^2)` under cumulative (test-nl-family.R, seed 11005):
  the default stopped at the thresholds-only model with the
  `frmtmb_nl_local` message; the probe now leaves it, higher by more than
  1e-3 of |logLik| (the record's 28.8 below brms's mode). The test is
  rewritten to read the message with the probe off and pin the new fit.
- Defects after (`dev/nlstart-out/repro-after1.txt`, first build; the
  final build's tests pin the same): sat seeds 101..501 slope 0.6 reach
  -530.4714139, -557.2304184, -555.3073014, -547.0530689, -563.4428877 at
  code 0 with no warning; slope -0.6 unchanged (-539.108463 and so on).
  Tiny x -272.503545 (a = 5.2021e+19, the linear coefficient). Large y
  -1658.76266.

## Tests

New file `tests/testthat/test-nl-start.R` (8 tests). Seen to fail on
rellib-r18 (`dev/nlstart-out/newbefore-files/`): `pass=11 fail=13 err=2`,
the failures behavioral: the sat body at -601.8 against -530.5 and -632.1
against -557.2 with b2 > 0 and no multistart message, a warning on seed
201; the option test (no gain); the probe-count test (0 optimizations
against 4); tiny x -304.1 against -281.4 (this file's own data); large y
-3045.7 against -1680.2 with a warning; the two internal-function tests
error (no such function). On the lane build `pass=30 fail=0`
(`dev/nlstart-out/newafter-files/`, `re1.log`). Guard-absent cases in it:
an all-added body (`a + b * x`) is not probed and equals `y ~ x`; a body
the probe cannot improve returns the probe-off fit bit for bit (opt$par
and logLik identical, vcov equal to 1e-8 relative, after 4 counted
optimizations), which also checks that the objective's state is put
back; a caller's `start` is not probed; units of 1 on an order-1 body.

Full suite (500 files) on the second build (`dev/nlstart-out/all1.log`,
495 clean); the 5 files that failed there and the 13 nl files
rerun on the final build all pass (`re1.log`). The session runs below
are the whole suite of every package on the final build. Gated as above.

One session per suite, three orders, final build
(`dev/nlstart-run-sessions.sh`, log `dev/nlstart-out/sess1.log`), every
`lib:` line frmtmb from `wt-nlstart-lib`:

```
SESSION frmtmb forward files=296/309 pass=22544 fail=0 err=0 skip=239 warn=0
SESSION frmtmb reverse files=296/309 pass=22694 fail=0 err=0 skip=239 warn=0
SESSION frmtmb preload files=296/309 pass=22787 fail=0 err=0 skip=241 warn=0
SESSION frmtmb.coupling forward/reverse/preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam      forward/reverse/preload files=33/33 pass=1848 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent   forward/reverse/preload files=12/12 pass=683 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn    forward/reverse/preload files=16/16 pass=474 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode      forward/reverse/preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample   forward/reverse/preload files=80/85 pass=3962 fail=0 err=0 skip=18 warn=0
SESSION frmtmb.survival forward/reverse/preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial  forward/reverse/preload files=8/8 pass=478 fail=0 err=0 skip=6 warn=0
SESS1 ran 27 of 27
```

(extension lines folded: the three orders read the same numbers.) Every
LEAK line carries `[with load of ...]`; none is this lane's.

OpenBLAS 0.3.26 with LAPACK routed (`dev/ciharden-openblas.sh 0.3.26
lapack`, OPENBLAS_NUM_THREADS=4, a 2000 x 2000 crossprod 0.25 s), logs
`dev/nlstart-out/ob/`: test-nl-start.R 30, test-nl-family.R 323,
test-nl.R 92, test-diagnostics-ux.R 130, test-se-flags.R 68,
frmtmb.ode test-ode-tv.R 27; 0 fail. 0.3.32 not run: no test here sits
near a rounding-decided branch (each margin is 1e-3 of |logLik| or more).
`dev/ciharden-out/` deleted.

## Bump and floor

Minor: the default fit of a curved nonlinear body changes (BREAKING in
NEWS). No extension code changed, no export changed, and every extension
suite passes against this core from rellib-r18's extension builds, so no
floor moves.

## Not done, and why

- No change to nl_domain_start(): the unit alone fixes defect 3, and the
  record keeps that start rule narrow on purpose.
- No detector for a short fit that the probe also misses. None exists
  that is cheap: the bernoulli rows are true second maxima with code 0, a
  small gradient and finite standard errors. The probe reaches them
  instead. "0 of 450" is against 8 random starts per arm.
- No `frmtmb_control()` argument for the probe: lane edits to
  frmtmb_control() in this round would conflict; an option is the switch.
- The brmsport ledger is not regenerated (consolidation does). The
  probe can change a ported row whose default nl fit stopped short; I did
  not identify such rows.
- `dev/test-backlog.md` not edited: the consolidation closes the three
  filed items ("Filed at the 0.79.0 release": MAJOR `b1 (1 - exp(-b2 x))`
  item, MINOR `a * x` near 1e-20; "From the lanes' records", nlfix: `a + b
  * x` at y about 1e6).

## Defects found, not fixed

- First build, both fixed in the fix round below: the probe start of
  magnitude 1 overflowed a body on a wide covariate (review major 2), and
  the cap of 4 patterns by position made the result depend on term order
  (review major 1).
- Left: 14 of the review's 335 fresh rows are still short, 4 silently
  (Michaelis-Menten under binomial and bernoulli, a power law under
  cumulative); no sign pattern reaches them (fix-round section).
- With five or more curved coefficients only the four the body moves most
  by get every pattern.

## R CMD check

`dev/nlstart-check.sh` (built and checked in `dev/nlstart-check/`, R_LIBS
lane lib, rellib-r18, user lib; pandoc and TinyTeX on PATH;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`), log `dev/nlstart-check/check.log`,
`dev/nlstart-out/check-00check.log`: `Status: 1 NOTE`, the expected one
(HTML manual: "Skipping checking math rendering: package 'V8'
unavailable"). Tests 41 min, OK. Only core changed, so only core is
checked. The `.Rcheck` tree is deleted; its log is kept.

## Files touched

- `R/par-template.R`: new nl_par_units(), nl_body_additive(),
  nl_sign_patterns(), nl_sign_probe(), nl_fit_probed() after
  nl_local_message(); the "default starting values" section of
  par_template()'s help.
- `R/fit.R`: frm()'s call site (nl_fit_off_saddle() to nl_fit_probed()),
  the `nu` hunk before predictor_par_units() in fit_assembled(), the
  `@param start` text.
- `man/frm.Rd`, `man/par_template.Rd` (roxygen), `DESCRIPTION`
  (roxygen's `Config/roxygen2/version` 8.1.0 to 8.1.1, from the updated
  user library; every lane that roxygenises gets it).
- `NEWS.md`: two bullets under `# frmtmb (development version)`.
- Tests: new `test-nl-start.R`; `test-nl-family.R`,
  `test-diagnostics-ux.R`, `test-se-flags.R` (probe held off on fixtures
  that rely on the old stop; one test pins the new default);
  `test-conditions-census.R` regenerated from `dev/conditions-census-gen.R`
  (three exemptions).
- `dev/`: `nlstart-*` scripts and this file.

## Fix round (review dev/reviews/2026-11-06-nlstart.md, verdict FIX)

### Merge with lane loglik (MAJOR 3): what consolidation must change

Not changed here: `logLik_penalized()` does not exist in this tree. After
the merge with loglik, in `R/par-template.R`, `nl_sign_probe()`:

- the two lines (1342 and 1343 in this tree)
  `ll0 <- as.numeric(stats::logLik(fit))` and
  `ll1 <- suppressWarnings(as.numeric(stats::logLik(r)))` become
  `ll0 <- as.numeric(logLik_penalized(fit))` and
  `ll1 <- suppressWarnings(as.numeric(logLik_penalized(r)))`, so that the
  probe compares what the optimizer, nl_restart_round() and
  nl_fit_off_saddle() compare. Wrap neither in a call that can stop: under
  REML with a prior on fixed effects loglik's `logLik()` stops, and the
  probe must not.

In `nl_fit_probed()`, the `frmtmb_nl_multistart` message (lines 1397 to
1403 in this tree) says "log-likelihood" twice. Take loglik's word:
`llw <- if (fit_penalized(pr$fit)) "penalized log-likelihood" else
"log-likelihood"`, as loglik's nl_fit_off_saddle() does at its line 843,
and write `"reached ", llw, " ", ll[1L]` and `"the best reached ", llw, " ",
ll[2L]`. On the reviewer's trial merge (`rev-nlstart-lib/merge-loglik`,
sat seed 101 with `normal(0, 2)` on b1) the kept restart said "penalized
log-likelihood -616.7576" and the probe "log-likelihood -611.7488" of the
same fit.

Lane bigx: the reviewer's trial merge is clean and does not double the
scale (`opt_coords()` replaces `par_units` in a basis block). Keep a test of
`a * xt`, `a ~ 1 + wb` at consolidation: `fit$par_units` on the merge still
reports nl units the optimizer did not use there.

### What changed in the fix round

- **Blocker 1, prior in the unit.** nl_par_units() takes the prior
  entries (fit_assembled() passes `prior_entries`) and adds the curvature
  of a normal, Student t or logistic prior on the internal scale
  (nl_prior_precision()), per row: `u = 1 / sqrt(r^2 / resp_u^2 + p / n)`.
  nl_resp_unit() factors out the response unit. The review's 63 prior rows
  (`dev/review-nlstart-prior.R`, seeds 6201 to 6203, log
  `dev/nlstart-out/fix/prior-after2.log`, summary
  `dev/nlstart-out/fix/priorsum.R` against the review's `prior-before.log`):
  `rows 63; lane body lower than base by > 1e-6: 0, higher: 0`. Codes and
  warning counts equal base's (1 row at code != 0, 7 with a warning, on
  both). The blocker's row (seed 6201, s = 1e-20, sd 1e-3) reads
  -298.07016213 at code 0, the linear spelling's value.
- **Major 1, order.** nl_sign_patterns(k, lead): every pattern of the
  coefficients at `lead`, at most four, the rest positive. nl_sign_probe()
  ranks the probed coefficients by the RMS of the body's derivative
  (nl_probe_scale()'s `key`), so the choice is the body's, not the
  formula's; k <= 4 runs all 2^k patterns. `dev/review-nlstart-order.R`
  (log `dev/nlstart-out/fix/order2.log`): all three spellings equal on all
  6 rows: bernoulli 33 -162.0191, bernoulli 22 -156.5785, binomial 33
  -531.9234, poisson 11 -407.5344, gamma 33 -294.9961 (was -299.1039 and
  -299.1040), cumulative 22 -331.7084, each with 0 warnings.
- **Major 2, start size.** nl_probe_scale(): the body's derivative in the
  coefficient where it is 0 and the other probed coefficients are 1; start
  size `2 * resp_u / rms`, to a power of 2, used where it is 16 or more
  times from 1, else the old 1 over the column's half range. On the
  review's bigx body (`dev/nlstart-bigx-dbg.R`, seeds 22, 44, 55, output on
  screen, numbers here): one unit (`resp_u / rms`, 2^-11 = 4.9e-4) stayed
  at -404.4046, -412.2400, -405.8144; two units (2^-10 = 9.8e-4), 1e-3 and
  2e-3 reached -385.4386, -397.55, -392.98. So 2 of 5 seeds at one unit,
  5 of 5 at two. `dev/review-nlstart-bigx.R` (log
  `dev/nlstart-out/fix/bigx2.log`): defaults -375.6325, -385.4386,
  -381.3314, -397.5533, -392.9816 against the review's best -375.6322,
  -385.4386, -381.3345, -397.5503, -392.9815, each at code 1 with 1 or 2
  warnings (as the started fits are), each with the multistart message.
  0.79.1 was 12.8 to 33.0 short.
- **Minors.**
  - Threshold: nl_probe_gain() = max(0.01, 1e-9 |value|), for both the
    tape's objective and the refit's log-likelihood. The smallest gain on
    the 450 rows was 0.43; on the fresh rows after this round 0.106
    (`freshsum-fix.txt`, gain quantile 0%); the review's three ridge slides
    were 9.6e-4 to 2.0e-3. 0 of the 89 messages are for a gain below 0.01.
  - The kept restart's message has the class `frmtmb_nl_restart` (one
    argument added to nl_fit_off_saddle()'s frm_message() call, away from
    lanes fitedge's and superassign's hunks) and is dropped when the probe
    wins; pinned in "a body whose maximum has both signs negative is
    reached".
  - `?frm` (`@param start`) names `options(frmtmb.nl_sign_probe = FALSE)`
    and says a located prior does not skip the probe; so does
    `?par_template`. nl_start_message()'s error text ends "Only (1) stops
    frm() from also trying the sign patterns of the body's coefficients".
  - NEWS: both bullets BREAKING, every line 80 columns or fewer; "now 0 do"
    is scoped to the 450 rows, and the fresh rows are given with what is
    left.
  - `dev/test-backlog.md`: the review's two "Filed, not this lane's" items,
    under "Filed by lane nlstart, round of 2026-11-06".

### Fix-round measurements

450 rows, final build (`scan-after4`, `dev/nlstart-out/scansum-after4.txt`):

```
scan-after4: rows 450, default NA 0, short > 1e-4: 0 (silent 0, warned 0), short > 1: 0 (silent 0)
scan-after4 minus scan-before default logLik: higher by > 1e-4: 16  lower by > 1e-4: 0  within 1e-4: 434
```

Cost there, interleaved (`cost5`): `ratio on / off: 1.845   control off2 /
off: 0.996`; by body exp 2.445, logi 2.333, prod 1.129, sat 1.758, two
2.095. Probe off equals 0.79.1 on all 450 (max |diff| 0).

The review's fresh rows (425, of which 90 fail from the default start on
every build, leaving 335): `dev/review-nlstart-fresh.R` and `fresh2.R`
rerun on this build without the reference (`dev/nlstart-out/fresh-after-fix/`,
and `fresh2-after-fix3/`, the review's three fresh2 bodies), summarized by
`dev/nlstart-freshsum.R` against the review's before logs, with the
reference from the review's after logs (best of 16 to 17 starts and every
arm). Log `dev/nlstart-out/freshsum-fix.txt`:

```
rows 425, default error (after) 90, (before) 90
short > 1e-4: after 14 (silent 4), before 97 (silent 38)
short > 1: after 5 (silent 2), before 86 (silent 32)
after - before: higher > 1e-4 89, lower < -1e-4 0, within 246
probe-off - before: |diff| > 1e-8 0, max |diff| 0
multistart message on 89 rows; of them gain over probe-off < 1e-2: 0
```

(The review's count on the first build: 48 short, 17 silent.) The 14 left:
Michaelis-Menten under binomial and bernoulli (7), cumulative `pow` 22,
gompertz bernoulli 33 and 22 (0.0125 and 0.0017, warned), bigx 11, 44 and
55 (2.3e-4 to 3.6e-3, warned), sig4 bernoulli 55 (1.3e-4, warned). The 4
silent: binomial mm 11 (4.32) and 44 (0.92), bernoulli mm 55 (1.84, with
the multistart message), cumulative pow 22 (0.60). The review found no
pattern that fixes the Michaelis-Menten rows either.

Time on / off on the fresh rows (same process, on first, not
interleaved), medians by body: bigx 1.67, five 2.90, gomp 1.87, logis3
1.63, mm 2.71, pow 1.67, sat_re 2.40, sat_sig 1.54, sat_z 2.34, sig4 1.21,
two_sat 3.49, two_sat_rev 3.46; total 2.92. The two-term bodies have four
curved coefficients and now run 16 patterns.

n = 20000, rotating off / on / off2, min of 3 (`dev/review-nlstart-cost.R`,
log `dev/nlstart-out/fix/cost-big.log`):

```
mm gaussian            on/off 1.55 control 1.00
logis3 poisson         on/off 1.86 control 0.95
sig4 gaussian          on/off 5.41 control 1.00   (four varied, 16 patterns)
five gaussian          on/off 46.88 control 0.98  ll off -20234.96 on -4329.16
exp RE gaussian        on/off 2.75 control 0.94
exp bernoulli          on/off 2.23 control 0.89
a + b*x (not probed)   on/off 1.10 control 1.07
ode nlf, default start on/off 2.27 control 0.88
```

The five-coefficient body takes 20.86 s against 0.45 s and gains 15906.
That is the price of 16 patterns there; the review measured 15.0 times
with the first build's 4.

### Tests in the fix round

`test-nl-start.R` gains 4 tests (order, wide covariate, prior unit, RE
restore) and asserts no `frmtmb_nl_restart` message on a win: 43 pass on
the final build. On rellib-r18 (`dev/nlstart-out/nsbefore-files/`):
`pass=21 fail=16 err=2`; the order test fails (-163.5 against -162.0), the
wide-covariate test fails (-404.4 against -385.4, no message); the prior
and RE-restore tests pass there, as base has no unit and no probe. On the
review's `norestore` mutant (`dev/nlstart-out/nsmut-files/`): `pass=33
fail=10`; the RE restore test fails 4 times (last.par.best 1.6442452117
against 1.6442448642, last.par, sdreport par.fixed and cov.fixed). That
mutant is the first lane build minus the restore, so the order,
wide-covariate and prior tests fail there too, on the first build's
behavior: the prior test reads code 1, a warning and -298.07061 against
-298.07016.

The review's 8 files plus the census, final build
(`dev/nlstart-out/fix1.log`): nl-start 43, nl 92, nl-family 323, nlf 77,
optstart 69, par-template 51, diagnostics-ux 130, se-flags 68,
conditions-census 11; 0 fail, 0 error, 0 skip, 0 warning; every `lib:`
line the lane's.

Core and frmtmb.ode, one file per process, final build
(`dev/nlstart-out/fix2.log`, 320 of 320, frmtmb from `wt-nlstart-lib` on
every `lib:` line): pass 23106, fail 0, err 0, skip 240 (ungated), warn
0. frmtmb.ode: test-ode-nlf.R 6, test-ode-tv.R 27, test-ode.R 72,
test-ode-events.R 142, test-lincmt.R 65, test-records.R 204, the rest
passing or skipped as on base.

R CMD check --as-cran of core on the fix build (`dev/nlstart-check.sh`,
log `dev/nlstart-check/check.log`, `dev/nlstart-out/check-00check-fix.log`):
`Status: 1 NOTE`, the V8 math-rendering note on the HTML manual; tests 31
min, OK. `.Rcheck` and the tarball deleted.

## Second fix round (review's short check, verdict FIX)

### What changed

- **R1, uncovered prior.** nl_prior_precision() returns `NA` when a prior
  with a curved density it has no formula for (gamma, lognormal, beta,
  inverse gamma, ...) or one on another scale reaches the coefficient; an
  exponential, uniform or flat prior counts 0 (no curvature). nl_par_units()
  then gives that coefficient no unit. Review's prior-class rows
  (`dev/review-nlstart-fix-prior2.R`, seeds 6201 and 6203, log
  `dev/nlstart-out/fix3/prior2-after.log`, compared by
  `dev/nlstart-out/fix3/prior2cmp.R` with the review's
  `fix-prior2-before.log`): `rows 66; below base by > 1e-6: 0; above by >
  1e-6: 4; code differs: 0`. The 4 above are `uniform(0, 1e21)` at s =
  1e-8 and 1e-20, 22.7 and 26.4 higher, as the review found. The 63 rows
  of `dev/review-nlstart-prior.R` (`fix3/prior-after.log`,
  `fix3/priorsum.R`): `rows 63; lane body lower than base by > 1e-6: 0,
  higher: 0`. NEWS's two sentences are corrected: the units bullet names
  the three covered prior families and says other curved priors keep the
  old steps; "No fit is lower than before" is scoped to the 785 fits.
- **R2, sizes read with the other coefficients at 1.** nl_probe_scale() now
  sizes in two stages, every coefficient at once (no order). Stage 1: a
  coefficient the body reads inside a curve gets 3 over its largest
  absolute derivative where it is 0 and the others are 1 (the curve's
  argument moves by 3 on the row it moves most). Then each coefficient the
  body is linear in (second difference 0 at the stage-1 sizes) gets 2
  response units over its derivative's RMS, read with the curve
  coefficients at their stage-1 starts. Not the suggested fixed-point pass:
  for `b1 * (1 - exp(-b2 * x))` the body reads `b1 b2` as one product near
  0, so taking each size from the other's trades the two sizes back and
  forth (reasoned from the derivatives, not run). Changing the curve size
  from `2 resp_u / rms` to `3 / max|d|` keeps the overflow fix (both give
  2^-10 = 9.8e-4 on x in 0 to 3000) and gives b2 about 1e4 on x about
  1e-4. Review's rows (`dev/review-nlstart-fix-small.R`, all four
  families, log `dev/nlstart-out/fix3/small-after.log`, compared by
  `dev/nlstart-smallcmp.R` with the review's reference and base logs):

  ```
  gaussian   exp  short > 1e-4: now 0 (silent 0), base 0; lower than base: 0
  gaussian   sat  short > 1e-4: now 0 (silent 0), base 5; lower than base: 0
  poisson    exp  short > 1e-4: now 0 (silent 0), base 0; lower than base: 0
  poisson    sat  short > 1e-4: now 0 (silent 0), base 0; lower than base: 0
  bernoulli  exp  short > 1e-4: now 0 (silent 0), base 2; lower than base: 0
  bernoulli  sat  short > 1e-4: now 0 (silent 0), base 0; lower than base: 0
  cumulative exp  short > 1e-4: now 0 (silent 0), base 0; lower than base: 0
  cumulative sat  short > 1e-4: now 0 (silent 0), base 3; lower than base: 0
  ```

  The 5 overflow seeds (`dev/review-nlstart-bigx.R`, log
  `dev/nlstart-out/fix3/bigx.log`) read as before: -375.6325, -385.4386,
  -381.3314, -397.5533, -392.9816, each with the multistart message.
- **Ties.** Ties in the ranking key go by `-abs(value)`, `-value` (the fitted
  values), then the coefficient's name, not by formula order. Review's
  `dev/review-nlstart-fix-order5.R` (36 rows: 6 and 5 curved coefficients,
  and a body whose three terms read x, 4 families, 3 seeds; log
  `dev/nlstart-out/fix3/order5.log`): 9 of 36 rows differ by more than
  1e-6 across the four spellings, the same 9 and values as the review's
  run: six cumulative 503 0.029, tie poisson 502 7.2e-4 and 503 4.7e-3,
  tie bernoulli 501 1.2e-5, 502 0.147 and 503 0.0149, tie cumulative 501
  2.3e-4, 502 0.12 and 503 0.0685. CORRECTION: an earlier version of this
  record said 0 of 36; its count read the wrong column of the log.
  Renaming changes no row; reversing the terms changes these 9, all
  warned fits whose probe-off first fits already differ by term order
  (review: six cumulative 503 -434.13 against -440.58). The tie-break
  did not change them. `dev/review-nlstart-order.R` (`fix3/order.log`):
  the six rows equal on all three spellings, 0 warnings, the same values
  as before.
- **Cost.** The 16-pattern default stays: the review measured no cheaper
  trigger that keeps the wins (a trigger on a warning loses 10 of 33
  four-coefficient gains of 0.29 to 15.3). A `frmtmb_control()` switch for
  the probe is filed for after the round (it was not added now because
  another lane edits frmtmb_control()'s signature this round).
- Two tests added to `test-nl-start.R`: the gamma row (code 0, no warning,
  no unit on `a`) and the gaussian sat body at x about 1e-4, seed 1 (the
  default reaches the start at the truth, with the multistart message).
  The gamma row failed on the first fix build (review's
  `fix-prior2-after.log`: code 1, 2 warnings, 992.09 below base); that
  build is overwritten, so its failure is the review's record, not rerun.

### A regression caught in this round, and its fix

The first build of this round (sizes `3 / max|d|` for a curve coefficient,
the 16-fold zone read after rounding) made the review's fresh rows worse
(`dev/nlstart-out/freshsum-fix5.txt`): short 19 (9 silent) against 14 (4)
on the previous fix build. The 5 new ones were all `logis3` (rate form
`Asym / (1 + exp(-k * (x - xmid)))`) under `Gamma("log")`, 9.7 to 23.4
short at code 0 with no warning and no probe win. Cause: `xmid`'s size
came to 3 / 0.25 = 12, which rounds to 16 and passed the zone, so `xmid`
started at +-16, off the data. Fix: read the zone before rounding (12 is
inside it). On that build (`dev/nlstart-out/fix4/logis3-gamma.log`) the 5
rows read -302.566394, -287.846894, -296.837260, -300.751527 and
-258.658399, each the reference, each with the multistart message. The
overflow rows (`fix4/bigx.log`) and the small-covariate rows
(`fix4/small-after.log`: 0 short in all 8 family-body cells) are
unchanged. The rest of this round's runs (450 rows, fresh rows, order5,
core and frmtmb.ode suites, R CMD check) were repeated on this build, in
the next section. The scan, suite and check of the build
before it (`scan-after5`, `fix3.log`, check `Status: 1 NOTE`) are
superseded.

### Second fix round: final-build numbers

All on the build with the zone read before rounding (install12).

- 450 rows (`scan-after6`, `dev/nlstart-out/scansum-after6.txt`): short
  > 1e-4 0 (silent 0); higher than 0.79.1 on 16, lower on 0, within 1e-4
  on 434.
- Review's 335 fresh rows (`freshsum-fix6.txt`, as before):
  `short > 1e-4: after 14 (silent 4), before 97 (silent 38)`;
  `short > 1: after 5 (silent 2)`; higher than base 89, lower 0; probe
  off equals base on every row (max |diff| 0); 89 multistart messages, 0
  for a gain below 0.01. The same 14 rows as the previous fix build.
  Time on / off, same process, on first: total 2.99; medians by body 1.41
  (sig4) to 3.50 (two_sat).
- Review's order5 (`fix4/order5.log`): 9 of 36 rows differ by more than
  1e-6, as above (corrected; an earlier version said 0).
- Small covariate and overflow rows: as in the regression note above.
- Prior rows: nl_par_units() is unchanged since `fix3/prior2-after.log`
  and `fix3/prior-after.log` (the zone fix is in nl_probe_scale() only,
  and those bodies are added terms, which the probe does not vary).
- `test-nl-start.R` on rellib-r18 (`nsbefore2-files/`): pass=24 fail=18
  err=2; the new small-covariate test fails there (-615.9 against -540.0,
  no message). The gamma-prior test passes on base, which had no unit; it
  failed on the first fix build, per the review's log.
- Core and frmtmb.ode, one file per process (`dev/nlstart-out/fix4.log`,
  320 of 320, every `lib:` line frmtmb from `wt-nlstart-lib`): pass 23111,
  fail 0, err 0, skip 240, warn 0. test-nl-start.R 48, test-nl.R 92,
  test-nl-family.R 323, test-ode-tv.R 27, test-ode-nlf.R 6.
- R CMD check --as-cran, core (`dev/nlstart-check/check.log`,
  `dev/nlstart-out/check-00check-fix4.log`): `Status: 1 NOTE` (V8); tests
  26 min, OK.
- Not rerun on this build: the 27 sessions, the extension suites other
  than frmtmb.ode, OpenBLAS.

## Third fix round (review's second short check, R3 and the order5 record)

### What changed

- **R3, curve coefficients stepped in absolute units.** nl_par_units()
  takes `curve = TRUE` for the refit of a probe's winner (fit_assembled()
  reads the attribute `frmtmb_nl_probe` that nl_sign_probe() puts on the
  winner's start) and then also gives a coefficient inside a curve
  (nl_probe_scale()'s `curve` row: the body is not linear in it) the
  probe's start size as its unit, where that is 16 or more times from 1
  and no prior reaches it. The probe's own runs use the same factor
  (`pu_run`). A first build of this round gave the unit to every fit, the
  first fit and user starts included: the scale sweep was the same as
  below, but the core and frmtmb.ode suites failed 6 times
  (`dev/nlstart-out/fix5.log`): frmtmb.ode test-ode-tv.R 4 (the late
  shift of `lk` again ran to a flat phase from its given start) and
  test-open-issues.R's nlme check 1 (`xmid` at a gradient of 59380), with
  test-nl-start.R 1 (its first fit no longer needed the probe), and R CMD
  check stopped on the same two. So the unit is now in the probe's runs and
  the refit only, which the review asked for; a first fit and a user's
  start step as before.
- **Silent refits.** nl_sign_probe() records whether the winner's refit
  settled: it warned (its own warnings are said), or it stopped at code 0
  no lower than the pattern's run on the tape. If not, nl_fit_probed()
  warns with class `frmtmb_nl_multistart_unsettled`. Tested by mocking
  nl_sign_probe() to report an unsettled refit (the class is raised beside
  the multistart message); I found no real case that triggers it.
- Test added: `b1 * (1 - exp(-b2 * xs))`, poisson, seed 122 (the review's
  seed 3), xs = x * 1e-6, equals the fit on x to 1e-6, no warning, b2
  rescaled to 1e-4.

### Scale sweep (`dev/review-nlstart-fix2-inv.R`)

On this build (`dev/nlstart-out/fix6/inv.log`, compared by
`dev/nlstart-invcmp.R` with the review's previous-build log):

```
rows 40; 1e-4 vs 1 differ > 1e-4: 0; 1e-6 vs 1 differ > 1e-4: 0
l1: below the best of every arm by > 1e-4: now 0, previous build 0; now lower than previous by > 1e-4: 0; warned now 0
l4: below the best of every arm by > 1e-4: now 0, previous build 0; now lower than previous by > 1e-4: 0; warned now 0
l6: below the best of every arm by > 1e-4: now 0, previous build 17; now lower than previous by > 1e-4: 0; warned now 0
```

Small-covariate rows at 1e-4 (`fix6/small-after.log`): 0 short in all 8
family-body cells. Overflow rows (`fix6/bigx.log`): all 5 now at code 0
with no warning, and higher than the review's best of 16 starts: -371.3312
(best -375.6322), -385.4149 (-385.4386), -381.2896 (-381.3345), -395.0881
(-397.5503), -392.4971 (-392.9815). Before this round they stopped at code
1 with 1 or 2 warnings.

### Lane bigx's nl_basis_descent() (read only, not run)

In `frmtmb-wt-bigx/R/fit.R`, nl_basis_descent() runs a second descent in
coordinates that whiten the body's Jacobian at the first descent's stop,
when that Jacobian is ill-conditioned, and keeps it if lower. At x about
1e-6 the Jacobian's columns differ by about 1e6, so on bigx alone it would
likely rescale b2 for the first fit. After the merge, three points from
the code: (1) it runs after the first descent, from a stop where b2 may
not have moved, so it is a repair, not a scale; (2) it sets
`obj$env$.frmtmb_opt_basis`, which bigx's optimize_obj() reads, so the
probe's runs on that tape would step in that basis too; (3) the review
found that opt_coords() replaces `par_units` inside a basis block, so for a
coefficient inside such a block the curve unit of this lane is replaced by
the basis. Whether it covers this case alone after the merge was not run.
This lane does not depend on it: the sweep above is on this lane's build.

### Third fix round: final-build numbers

Build install15 (the R3 design above); install17 adds only the help text
of par_template() and the NEWS wording, and is what the library holds.

- Tests, one file per process (`dev/nlstart-out/r4t.log`): test-nl-start.R
  53, test-nl.R 92, test-nl-family.R 323, test-open-issues.R 46,
  frmtmb.ode test-ode-tv.R 27 and test-ode-nlf.R 6; 0 fail, 0 error, 0
  skip, 0 warning.
- 450 rows (`scan-after8`, `scansum-after8.txt`): 0 short; higher than
  0.79.1 on 16, lower on 0.
- Review's 335 fresh rows (`freshsum-fix8.txt`): `short > 1e-4: after 11
  (silent 4), before 97 (silent 38)`; `short > 1: after 5 (silent 2)`;
  higher than base 89, lower 0; probe off equals base on every row (max
  |diff| 0); 89 multistart messages, 0 for a gain below 0.01. The 11: the
  same 8 Michaelis-Menten and power-law rows, gompertz bernoulli 33 and 22
  (0.0125 and 0.0017, warned), sig4 bernoulli 55 (1.3e-4, warned); the 3
  bigx rows are no longer short. Time on / off, same process, on first:
  3.20.
- Order5 (`fix6/order5.log`): 9 of 36 rows differ by more than 1e-6, the
  same 9 as before (warned fits; renaming changes none).
- Core and frmtmb.ode, one file per process (`dev/nlstart-out/fix6.log`,
  320 of 320, every `lib:` line frmtmb from `wt-nlstart-lib`): pass 23116,
  fail 0, err 0, skip 240, warn 0.
- R CMD check --as-cran, core, on the source with the final help text
  (`dev/nlstart-check/check.log`, `dev/nlstart-out/check-00check-fix6.log`):
  `Status: 1 NOTE` (V8); tests OK.
- Not rerun: the 27 sessions, the extension suites other than frmtmb.ode,
  OpenBLAS, the prior rows (added-term bodies, which neither the probe
  nor the curve unit touches).
